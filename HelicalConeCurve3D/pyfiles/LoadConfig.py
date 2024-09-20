import yaml
import scipy

def LoadConfig(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
## scanner geomery
    SO = config['geometry']['SO']
    DD = config['geometry']['DD']
    YL = config['geometry']['YL']
    ZL = config['geometry']['ZL']
    DecAngle = config['geometry']['DecAngle']
    DecHeight = config['geometry']['DecHeight']
    N_Turn = config['geometry']['N_Turn']
    N_2pi = config['geometry']['N_2pi']
    h = config['geometry']['h']
    dectorYoffset = config['geometry']['dectorYoffset']
    dectorZoffset = config['geometry']['dectorZoffset']
## recon params
    imageSize = config['recon']['imageSize']
    sliceCount = config['recon']['sliceCount']
    sliceThickness = config['recon']['sliceThickness']
    ObjR = config['recon']['ObjR']
    kernelType = config['recon']['kernelType']
    centerOffset = config['recon']['centerOffset']

    k1 = config['recon']['k1']
    delta = config['recon']['delta']
    HSCoef = config['recon']['HSCoef']
## data path
    input_path = config['data']['input_path']
    save_path = config['data']['save_path']
    data = scipy.io.loadmat(input_path)
    ProjData = data['Proj']

## dll path
    dll_path = config['dll']['dll_path']


    return  ProjData, SO, DD, YL, ZL, DecAngle, DecHeight, N_Turn, N_2pi, h, dectorYoffset, dectorZoffset, \
            imageSize, sliceCount, sliceThickness, ObjR, kernelType, centerOffset, \
            k1, delta, HSCoef, input_path, save_path, dll_path

        
