import numpy
from dll_loader import PSFDllLoader
from textwrap import dedent
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline

PIXEL_SIZE = 9.0 # microns

wavelength = 1.0 / PIXEL_SIZE
focal_ratio = 11.0
s_x = s_y = 1.0 / PIXEL_SIZE
subpixel_res_x = 1
subpixel_res_y = 1
map_res_x = 6000
map_res_y = 6000 
max_wavenumbner_x = 6 * PIXEL_SIZE 
max_wavenumber_y = 6 * PIXEL_SIZE

class PRF():
    def __init__(self, wavelength, focal_ratio, variance, subpixel_res, 
        map_res, max_wavenumber):
        c_src = PSFDllLoader("./libpsf.so")
        print(dedent(f"""\
        Building PSF Model with the following parameters:
        Focal Ratio:                {focal_ratio},
        Centroid Stddev (x,y):      {variance[0], variance[1]},
        Subpixel Resolution (x,y):  {subpixel_res[0], subpixel_res[1]},
        Map Resolution (x,y):       {map_res[0], map_res[1]},
        Max Wavenumbers (x,y):      {max_wavenumber[0], max_wavenumber[1]} 
        """))
        print(">>>> Loaded C library <<<<")
        c_src.execute(
            wavelength, focal_ratio, 
            variance[0], variance[1],
            subpixel_res[0], subpixel_res[1],
            map_res[0], map_res[1],
            max_wavenumber[0], max_wavenumber[1]
        )

        self.intensities = RegularGridInterpolator((c_src.x_vals, c_src.y_vals), c_src.I_vals)
        #self.intensities = RectBivariateSpline(c_src.x_vals, c_src.y_vals,c_src.I_vals)
    def __call__(self, x, y):
        return self.intensities((abs(x),abs(y)))

if __name__ == "__main__":
    prf = PRF(
        wavelength, focal_ratio, (s_x, s_y), (subpixel_res_x,subpixel_res_y), (map_res_x, map_res_y), (max_wavenumbner_x,max_wavenumber_y)
    )
    print(prf(numpy.array([[3.17918]]),numpy.array([[0.967098]])))