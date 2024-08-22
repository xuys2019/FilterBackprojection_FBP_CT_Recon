import matplotlib.pyplot as plt
import scipy
from pyfiles.fbp_equiAngular import fbp_equiAngular


yaml_path = './config/config.yaml'

## Reconstruction 
RecImg = fbp_equiAngular(yaml_path)

# # save result
dataNew = './data/recon_curve.mat'
scipy.io.savemat(dataNew,
             {'RecImg': RecImg})

plt.figure()
plt.imshow(RecImg, cmap='gray')
plt.show()