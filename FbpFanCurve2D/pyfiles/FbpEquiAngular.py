import ctypes as ct
import math
import numpy as np
import matplotlib.pyplot as plt
import yaml

from pyfiles.CreateHSP import CreateHSP
from pyfiles.LoadConfig import LoadConfig

# Init ctypes types
DOUBLE = ct.c_double
PtrDOUBLE = ct.POINTER(DOUBLE)
PtrPtrDOUBLE = ct.POINTER(PtrDOUBLE)


class TestStruct(ct.Structure):
    _fields_ = [
                ("ScanR", ct.c_double),         # Radius of the scanning trajectory
                ("DecFanAng", ct.c_double),     # Fan angle coverage of the detector
                ("startangle", ct.c_double),     # Fan angle coverage of the detector
                ("YL", ct.c_int),               # Detector cell number
                ("YL_Offset", ct.c_double),     # Detector ccenter offset in pixels (e.g. 0.25)
                ("AngleNumber", ct.c_int),      # View number along the the scanning trajectory
                ("Radius", ct.c_double),        # Radius of the reconstructed phantom/object 
                ("RecSizeX", ct.c_int),         # Reconstructed image size on horizontal axis
                ("RecSizeY", ct.c_int),         # Reconstructed image size on vertical axis
                ("centerX", ct.c_int),          # Rotational center coordiante in the reconstructed image along the horizontal axis (pixel
                ("centerY", ct.c_int),          # Rotational center coordinate in the Reconstructed image along the vertical axis (pixel)
                ("FOILength", ct.c_int),        # Physical length of the reconstructed image along the horizontal axis
                ("FOIWidth", ct.c_int),         # Physical length of the reconstructed image along the vertical axis
                ("GF", PtrPtrDOUBLE),           # Projection/Sinogram data in the Radon space
                ("RecIm", PtrPtrDOUBLE)         # Reconstruction image data
                ]


def Double2DArray2Pointer(array):
    # Converts a 2D numpy into a array ctypes 2D array.
    arr_dimx = DOUBLE * array.shape[1]
    arr_dimy = PtrDOUBLE * array.shape[0]
    arr_ptr = arr_dimy()

    for i, row in enumerate(array):
        arr_ptr[i] = arr_dimx()
        for j, val in enumerate(row):
            arr_ptr[i][j] = val

    return arr_ptr


def Double2DPointer2Array(ptr, n, m):
    # Converts ctypes 2D array into a 2D numpy array.
    arr = np.zeros(shape=(n, m))
    for i in range(n):
        for j in range(m):
            arr[i][j] = ptr[i][j]
    return arr


def FbpEquiAngular(yaml_path):
    # Load the config
    ProjData, DecFanAng, ScanR, Radius, Number, YL, input_path, dll_path = LoadConfig(yaml_path)

    DeltaY = DecFanAng/(YL-1)
    DeltaY2 = 2*DeltaY
    Dgy = np.array(ProjData, dtype=np.float32)

    ### pre-weighting for ramp-filter
    for Yindex in range(1, YL-1):
        Dgy[:,Yindex]=(ProjData[:,Yindex+1]-ProjData[:,Yindex-1])/DeltaY

    Dgy[:,0] = Dgy[:,1]
    Dgy[:,YL-1]= Dgy[:,YL-2]

    ### filtering
    WindowType=1
    nn=int(math.pow(2,(math.ceil(math.log2(abs(YL)))+1)))
    HS=CreateHSP(nn,WindowType)
    nn2= nn*2
    k = int(nn/2)
    TempF=np.zeros(nn2)
    TempF[0:k]=HS[k:nn]
    TempF[k+nn:nn2]=HS[0:k]
    HS=TempF*complex(0,1)
    FFT_F=np.fft.fft(HS)

    GF=Dgy
    for ProjIndex in range(0,Number):
        TempData=np.zeros(nn2)
        for k in range(0, YL):
            TempData[k]=Dgy[ProjIndex,k]
        FFT_S=np.fft.fft(TempData,nn2)
        TempData=np.fft.ifft(FFT_S*FFT_F).imag
        for k in range(0,YL):
            GF[ProjIndex,k]=-TempData[k]

    # init the struct
    t = TestStruct()

    t.ScanR = ScanR
    t.DecFanAng = DecFanAng
    t.startangle = 0
    t.YL = YL
    t.YL_Offset = 0
    t.AngleNumber = Number
    t.Radius = Radius

    # These are flexible parameters.
    t.RecSizeX = 256
    t.RecSizeY = 256
    t.centerX = 128
    t.centerY = 128
    t.FOILength = t.RecSizeX
    t.FOIWidth = t.RecSizeY
    # Generate a 2D ctypes array from numpy array
    GF = GF.T
    GF_ptr = Double2DArray2Pointer(GF)
    t.GF = GF_ptr

    RecIm = np.zeros(shape=(t.FOILength, t.FOIWidth))
    RecIm_ptr = Double2DArray2Pointer(RecIm)
    t.RecIm = RecIm_ptr

    # Load the compiled library
    recon = ct.CDLL(dll_path)
    # Define arguments of the C function
    recon.fbp.argtypes = [ct.POINTER(TestStruct)]
    # Define the return type of the C function
    recon.fbp.restype = None

    # interface with C function
    recon.fbp(ct.byref(t))

    # Convert ctypes 2D arrays to numpy arrays
    RecImg = Double2DPointer2Array(RecIm_ptr, *RecIm.shape)
    RecImg = RecImg.T

    return RecImg
    
