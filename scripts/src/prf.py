import numpy
from dll_loader import PSFDllLoader
from textwrap import dedent
from scipy.interpolate import RegularGridInterpolator

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

    def __call__(self, x, y):
        return self.intensities((abs(x),abs(y)))