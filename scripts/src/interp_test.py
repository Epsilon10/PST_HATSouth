from dll_loader import PSFDllLoader
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib import cm


fig = plt.figure()

ax = fig.gca(projection='3d')


c_src = PSFDllLoader("./libpsf.so")


c_src.execute(
    0.6, 11.0, 0.01, 0.01, 10,10, 2 * int(2e3),2 * int(2e3), 2 * 3,2 * 3
)


X1, Y1 = np.meshgrid(c_src.x_vals[:-1] , c_src.y_vals[:-1] )
V1 = c_src.I_vals

c_src.execute(
    0.6, 11.0, 0.01, 0.01, 10,10, int(2e3), int(2e3), 3,3
)

func = RegularGridInterpolator((c_src.x_vals,c_src.y_vals), c_src.I_vals)
res = func((X1, Y1))
diff = np.absolute(res - V1[:-1, :-1])
surf = ax.plot_surface(X1,Y1,diff, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
print('max diff:', np.amax(diff), 'avg diff:', np.mean(diff))
print('I max:', np.amax(V1), 'I min: ', np.amin(V1), 'Mean: ', np.mean(V1))
print('thing:', np.mean(diff / (np.absolute(res) + np.absolute(V1[:-1, :-1]))))

#fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('plots/interp_graph')
plt.cla()

min = 0
max = 150
Z = diff / (np.absolute(res) + np.absolute(V1[:-1, :-1]))
surf = ax.plot_surface(X1[min:max, min:max],Y1[min:max, min:max],Z[min:max, min:max], cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.savefig('plots/ratio_graph')
plt.cla()
