import matplotlib.pyplot as plt
import scipy
from pyfiles.FbpEquiAngular import FbpEquiAngular
import numpy as np

yaml_path = './config/config.yaml'

## Reconstruction 
RecImg = FbpEquiAngular(yaml_path)

# # save result
dataNew = './data/reconstruction/ReconCurve.mat'
scipy.io.savemat(dataNew,
             {'RecImg': RecImg})
             
plt.figure()
plt.imshow(RecImg, cmap='gray')
plt.show()

