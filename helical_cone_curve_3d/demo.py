import matplotlib.pyplot as plt
import scipy
from pyfiles.helical_equiAngular import helical_equiAngular


yaml_path = './config/config.yaml'

## Reconstruction 
RecImg = helical_equiAngular(yaml_path)

# # save result
dataNew = './data/reconstruction/helical_recon_curve.mat'
scipy.io.savemat(dataNew,
             {'RecImg': RecImg})

plt.figure()
plt.imshow(RecImg[:,:,32], cmap='gray')
plt.show()