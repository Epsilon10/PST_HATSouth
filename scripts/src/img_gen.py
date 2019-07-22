from dll_loader import PSFDllLoader
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import math
import sys

wavelength = 0.6
focal_ratio = 11.0
s_x = s_y = 1.0
subpixel_res_x = 10
subpixel_res_y = 10
map_res_x = 6000
map_res_y = 6000
max_wavenumbner_x = 6
max_wavenumber_y = 6

fig = plt.figure()

ax = fig.gca(projection='3d')


c_src = PSFDllLoader("./libpsf.so")
c_src.execute(
    wavelength, focal_ratio, s_x, s_y, subpixel_res_x,subpixel_res_y, map_res_x, map_res_y, max_wavenumbner_x,max_wavenumber_y
)


X, Y = np.meshgrid(c_src.x_vals, c_src.y_vals)
Z = c_src.I_vals

surf = ax.plot_surface(X[:30, :30],Y[:30, :30],Z[:30, :30], cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('plots/I_map')

