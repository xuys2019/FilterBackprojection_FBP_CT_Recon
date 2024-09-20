import yaml
import scipy

def load_config(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)

    DecFanAng = config['geometry']['DecFanAng']
    ScanR = config['geometry']['ScanR']
    Radius = config['geometry']['Radius']
    Number = config['geometry']['Number']
    YL = config['geometry']['YL']
    dll_path = config['dll']['dll_path']
    input_path = config['data']['input_path']
    input_data = scipy.io.loadmat(input_path)
    ProjData = input_data['proj']
    
    return ProjData, DecFanAng, ScanR, Radius, Number, YL, input_path, dll_path

        
