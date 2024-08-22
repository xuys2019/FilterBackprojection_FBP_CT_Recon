import matplotlib.pyplot as plt
from CreateHSP import CreateHSP
from recon_fbp_equiAngular import recon_fbp_equiAngular


yaml_path = './config.yaml'

## Reconstruction 
RecImg = recon_fbp_equiAngular(yaml_path)

# # save result
dataNew = './data/recon_curve.mat'
scio.savemat(dataNew,
             {'RecImg': RecImg})

plt.figure()
plt.imshow(RecImg, cmap='gray')
plt.show()