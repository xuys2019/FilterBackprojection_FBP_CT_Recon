import ctypes as ct
import numpy as np
import math

import matplotlib.pyplot as plt

from pyfiles.load_config import load_config
from pyfiles.CreateHSP import CreateHSP
# Init ctypes types
FLOAT = ct.c_float
PtrFLOAT = ct.POINTER(FLOAT)
PtrPtrFLOAT = ct.POINTER(PtrFLOAT)
PtrPtrPtrFLOAT = ct.POINTER(PtrPtrFLOAT)


class TestStruct(ct.Structure):
    _fields_ = [
                ("ScanR", ct.c_float),
                ("DistD", ct.c_float),
                ("YL", ct.c_int),
                ("ZL", ct.c_int),
                ("dectorYoffset", ct.c_float),  # Detector along the horizontal direction (pixel, e.g. quarter pixel)
                ("dectorZoffset", ct.c_float),  # Detector offset along the vertical direcion (pixel, e.g. quarter pixel)
                ("XOffSet", ct.c_float),  # recon offset along the x axis(mm)
                ("YOffSet", ct.c_float),  # recon offset along the y axis(mm)
                ("ZOffSet", ct.c_float),  # recon offset along the z axis(mm)
                ("phantomXOffSet", ct.c_float),  # phantom offset along the x axis(mm)
                ("phantomYOffSet", ct.c_float),  # phantom offset along the y axis(mm)
                ("phantomZOffSet", ct.c_float),  # phantom offset along the z axis(mm)
                ("DecFanAng", ct.c_float),
                ("DecHeight", ct.c_float),
                ("DecWidth", ct.c_float),
                ("h", ct.c_float),
                ("BetaS", ct.c_float),
                ("BetaE", ct.c_float),
                ("AngleNumber", ct.c_int),
                ("N_2pi", ct.c_int),
                ("Radius", ct.c_float),
                ("RecSize", ct.c_int),
                ("RecSizeZ", ct.c_int),
                ("delta", ct.c_float),
                ("HSCoef", ct.c_float),
                ("k1", ct.c_float),
                ("GF", PtrPtrPtrFLOAT),
                ("RecIm", PtrPtrPtrFLOAT)
                ]


def double3darray2pointer(arr):
    # Converts a 3D numpy to ctypes 3D array.
    arr_dimx = FLOAT * arr.shape[2]
    arr_dimy = PtrFLOAT * arr.shape[1]
    arr_dimz = PtrPtrFLOAT * arr.shape[0]

    arr_ptr = arr_dimz()

    for i, row in enumerate(arr):
        arr_ptr[i] = arr_dimy()
        for j, col in enumerate(row):
            arr_ptr[i][j] = arr_dimx()
            for k, val in enumerate(col):
                arr_ptr[i][j][k] = val
    return arr_ptr


def double3dpointer2array(ptr, n, m, o):
    # Converts ctypes 3D array into a 3D numpy array.
    arr = np.zeros(shape=(n, m, o))

    for i in range(n):
        for j in range(m):
            for k in range(o):
                arr[i, j, k] = ptr[i][j][k]

    return arr


def load_C_recon_lib():
    # add recon lib path to environment value "PATH" for depending DLLs
    # # # # recon_lib = my_path.find_dir("top", os.path.join("reconstruction", "lib"))
    # # # # my_path.add_dir_to_path(recon_lib)

    #  my_path.find_dir doesn't have the key "reconstruction", use the temp solution below:
    recon_lib = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib")

    # load C/C++ lib
    ll = ct.cdll.LoadLibrary
    if os.name == "nt":
        lib_file = "fdk_equiAngle.dll"
    else:
        lib_file = "fdk_equiAngle.so"
    clib = ll(os.path.join(recon_lib, lib_file))

    return clib


def helical_equiAngular(yaml_path):
    # Load the config
    Proj, SO, DD, YL, ZL, DecAngle, DecHeight, N_Turn, N_2pi, h, dectorYoffset, dectorZoffset, \
            imageSize, sliceCount, sliceThickness, ObjR, kernelType, centerOffset, \
            k1, delta, HSCoef, input_path, save_path, dll_path = load_config(yaml_path)

    PI =3.14159265358979
    YLC= (YL-1)*0.5+dectorYoffset  # Detector center along the horizontal direction of detector array
    ZLC= (ZL-1)*0.5+dectorZoffset  # Detector center along the vertical direction of detector array

    BetaS = -N_Turn*PI
    BetaE =  N_Turn*PI
    ViewN =  N_Turn*N_2pi+1
    DecWidth = math.tan(DecAngle*0.5)*(DD)*2
    dYL   =  DecWidth/YL
    dZL   =  DecHeight/ZL
    DeltaFai= (BetaE-BetaS)/(ViewN-1)
    DeltaTheta = DeltaFai
    DeltaT     = math.tan(DecAngle*0.5)*(DD)*2/YL
    dYA = DecAngle/YL
    PProj    = np.zeros((ViewN,YL,ZL))

    ## rebinning the projection
    print("* Rebinning the projection...")
    for i in range(ViewN):
        Theta=(i)*DeltaTheta                   # the view for the parallel projection
        for j in range(YL):
            t      = (j-YLC)*DeltaT     # the distance from origin to ray for parallel beam
            Beta   = math.asin(t/(SO))            # the fan_angle for cone_beam projection
            Fai    = Theta+Beta              # the view for cone_beam projecton
            a      = math.atan(t/math.sqrt(SO**2-t**2))  # the position of this ray on the flat detector
            FaiIndex        =  (Fai/DeltaFai)
            UIndex          =  (a/dYA)+YLC
            FI              =  math.ceil(FaiIndex)
            UI              =  math.ceil(UIndex)
            coeXB           =  FI-FaiIndex
            coeXU           =  1-coeXB
            coeYB           =  UI-UIndex
            coeYU           =  1-coeYB
            if (FI<=0):
                IndexXU = 0
                IndexXB = 0
            elif(FI > ViewN-1):
                IndexXU = ViewN-1
                IndexXB = ViewN-1
            else:
                IndexXU = FI
                IndexXB = FI-1

            if (UI<=0):
                IndexYU = 0
                IndexYB = 0
            elif(UI>YL-1):
                IndexYU = YL-1
                IndexYB = YL-1
            else:
                IndexYU=UI
                IndexYB=UI-1
            PProj[i,j,:]=coeXB*coeYB*Proj[IndexXB,IndexYB,:]+coeXU*coeYB*Proj[IndexXU,IndexYB,:]+coeXB*coeYU*Proj[IndexXB,IndexYU,:]+coeXU*coeYU*Proj[IndexXU,IndexYU,:]


    #Preweighting the conebeam projections
    print("* Pre-weighting the filter...")
    for j in range(YL):
        t=(j-YLC)*dYL
        for k in range(ZL):
              b=(k-ZLC)*dZL
              Proj[:,j,k] = PProj[:,j,k]*SO*SO/np.sqrt(SO**4+(SO*b)**2-(b*t)**2)

    Proj = Proj.transpose(1,2,0)

    #Perform Ramp filtering
    print("* Applying the filter...")
    Dg=Proj
    nn = int(math.pow(2, (math.ceil(math.log2(abs(YL))) + 1)))
    nn2 = nn*2
    FFT_F = CreateHSP(nn, kernelType)

    GF = Proj

    for ProjIndex in range(0, ViewN):
        for j in range(ZL):
            TempData = np.ones(YL)
            for k in range(YL):
                TempData[k] = Dg[k, j, ProjIndex]
            FFT_S = np.fft.fft(TempData, nn2)
            TempData = np.fft.ifft(FFT_S * FFT_F).imag
            for k in range(YL):
                GF[k, j, ProjIndex] = TempData[k]
    GF = GF/dYL

    #Backproject the filtered data into the 3D space
    # Load the compiled library
    recon = ct.CDLL(dll_path)
    # Define arguments of the C function
    recon.fbp.argtypes = [ct.POINTER(TestStruct)]
    # Define the return type of the C function
    recon.fbp.restype = None

    # init the struct
    t = TestStruct()

    t.ScanR = SO
    t.DistD = DD
    t.YL = YL
    t.ZL = ZL
    t.DecFanAng = DecAngle
    t.startangle = 0
    t.DecHeight = DecHeight
    t.DecWidth = DecWidth
    t.h = h
    t.BetaS = BetaS
    t.BetaE = BetaE
    t.dectorYoffset = 0
    t.dectorZoffset = 0
    t.AngleNumber = ViewN
    t.N_2pi = N_2pi

    t.Radius = ObjR
    t.RecSize = imageSize
    t.sliceThickness = sliceThickness
    t.RecSizeZ = sliceCount
    t.delta = delta
    t.HSCoef = HSCoef
    t.k1 = k1

    t.XOffSet = centerOffset[0]
    t.YOffSet = centerOffset[1]
    t.ZOffSet = centerOffset[2]
    t.phantomXOffSet = 0
    t.phantomYOffSet = 0
    t.phantomZOffSet = 0

    print("* Converting projection data from a numpy array to a C array...")
    GF_ptr = double3darray2pointer(GF)
    t.GF = GF_ptr

    # RecIm = np.zeros(shape=(t.RecSize, t.RecSize, t.RecSize))
    print("* Allocating a C array for the recon results...")
    RecIm = np.zeros(shape=(t.RecSize, t.RecSize, t.RecSizeZ))
    RecIm_ptr = double3darray2pointer(RecIm)
    t.RecIm = RecIm_ptr

    # interface with C function
    print("* In C...")
    recon.fbp(ct.byref(t))

    # Convert ctypes 2D arrays to numpy arrays
    print("* Converting the recon results from a C array to a numpy array...")
    RecImg = double3dpointer2array(RecIm_ptr, *RecIm.shape)

    return RecImg
