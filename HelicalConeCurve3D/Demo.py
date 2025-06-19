import matplotlib.pyplot as plt
import scipy
from pyfiles.HelicalEquiAngular import HelicalEquiAngular

## configuration path
yaml_path = './config/config.yaml'

## Reconstruction 
RecImg = HelicalEquiAngular(yaml_path)

# # save result
dataNew = './data/reconstruction/HelicalReconCurve.mat'
scipy.io.savemat(dataNew,
             {'RecImg': RecImg})

plt.figure()
plt.imshow(RecImg[:,:,32], cmap='gray')
plt.show()