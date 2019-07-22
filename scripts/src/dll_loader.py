from ctypes import cdll, c_double, c_uint
import numpy as np

class BaseDllLoader:
    def __init__(self, path):
        self.path = path
    def load_library(self):
        self.lib = cdll.LoadLibrary(self.path) # TODO: fix this

class PSFDllLoader(BaseDllLoader):
    def __init__(self, path):
        super().__init__(path)
        self.load_library()
        self.lib.fill_I_vals.argtypes = [
            c_double, c_double, c_double, c_double, c_uint, c_uint, c_uint, c_uint, 
            c_double, c_double, np.ctypeslib.ndpointer(dtype=c_double, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=c_double, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=c_double, ndim=1, flags='C_CONTIGUOUS')
        ]

        
    def execute(self, wavelength, focal_ratio, variance_x,variance_y, subpixel_res_x, subpixel_res_y, 
        map_res_x, map_res_y, max_x_wavenumber, max_y_wavenumber):
        self.I_vals = np.empty(shape=(map_res_x, map_res_y), dtype=c_double)
        self.kx_vals = np.empty(shape=map_res_x, dtype=c_double)
        self.ky_vals = np.empty(shape=map_res_y, dtype=c_double)
        self.psi_amps = np.empty(shape=(map_res_x, map_res_y), dtype=c_double)
        self.x_vals = np.empty(shape=map_res_x, dtype=c_double)
        self.y_vals = np.empty(shape=map_res_y, dtype=c_double)
        self.lib.fill_I_vals(
            wavelength, focal_ratio, variance_x,variance_y, subpixel_res_x, subpixel_res_y, 
            map_res_x, map_res_y, max_x_wavenumber, max_y_wavenumber, self.I_vals, self.kx_vals, self.ky_vals, self.psi_amps, self.x_vals, self.y_vals
        )

    