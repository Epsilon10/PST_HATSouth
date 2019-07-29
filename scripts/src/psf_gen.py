from dll_loader import PSFDllLoader
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.ticker import LinearLocator, FormatStrFormatter



fig = plt.figure()

ax = fig.gca(projection='3d')


c_src = PSFDllLoader("./libpsf.so")
c_src.execute(
    0.6, 11.0, 0.1, 0.1, 1, 1, int(5e3), int(5e3), 3,3
)


X, Y = np.meshgrid(c_src.kx_vals, c_src.ky_vals)

Z = c_src.psi_amps

ax.set_xlabel('x wavelengths')
ax.set_ylabel('y wavelengths')
ax.set_zlabel('PSF')
surf = ax.plot_surface(X[:240,:240], Y[:240,:240], Z.transpose()[:240,:240], cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('plots/psf')

