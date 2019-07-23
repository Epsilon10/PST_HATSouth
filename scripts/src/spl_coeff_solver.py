from ctypes import cdll, c_double, c_uint
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from dll_loader import PSFDllLoader



class SubpixelCoeffSolver:
    """
    Solves for the coeffecients of the subpixel sensitivity map.
    Args:
        - wavelength
        - focal_ratio 
        - variance_x: standard deviation of the x coordinates of the spot centers
        - variance_y: standard deviation of the y coordinates of the spot centers
        - num_subpxl_x: number of subpixels to split each pixel along the x axis
        - num_subpxl_y: number of subpixels to split each pixel along the y axis
        - subpxl_res_x: resolution of the ifft along the x axis
        - subpxl_res_y: resolution of the ifft along the y axis
        - x_0: x coordinate of the spot center
        - y_0: y coordinate of the spot center
        - x_min: lower x coordinate of pixel
        - x_max: upper x coordinate of pixel
        - y_min: lower y coordinate of pixel
        - y_max: upper y coordinate of pixel
    """
    def __init__(self):
        self.c_src = PSFDllLoader("./libpsf.so")
    
    """
    Allocates memory for matrices
    Args:
        - max_x_wavenumber
        - max_y_wavenumber
    """
    
    def _get_response(self, pixel_x, pixel_y, row):
        for nx in range(self.subpixel_res_x):
            for ny in range(self.subpixel_res_y):
                x = pixel_x + nx*self.dx + self.dx/2.0 - self.x_0
                y = pixel_y + ny*self.dy + self.dy/2.0 - self.y_0
                print(x,y)
                self.sums[row][ny + nx * self.subpixel_res_y] = self.interp((abs(x),abs(y)))
        return np.sum(self.sums[row])
    
    def solve(self, wavelength, focal_ratio, variance_x,variance_y, subpixel_res_x, subpixel_res_y, 
        map_res_x, map_res_y, max_x_wavenumber, max_y_wavenumber, x_0, y_0, x_min, x_max, y_min, y_max):
        self.c_src.execute(
            wavelength, focal_ratio, variance_x,variance_y, subpixel_res_x, subpixel_res_y, 
            map_res_x, map_res_y, max_x_wavenumber, max_y_wavenumber
        )
        num_x_pixels = x_max - x_min + 1
        num_y_pixels = y_max - y_min + 1

        self.subpixel_res_x = subpixel_res_x
        self.subpixel_res_y = subpixel_res_y

        self.dx = 1 / subpixel_res_x
        self.dy = 1/subpixel_res_y
        self.x_0 = x_0
        self.y_0 = y_0
        self.responses = np.empty(shape=num_x_pixels * num_y_pixels, dtype = c_double)
        self.sums = np.empty(shape = (num_x_pixels * num_y_pixels, subpixel_res_x * subpixel_res_y), dtype=c_double)
        self.interp = RegularGridInterpolator((self.c_src.x_vals,self.c_src.y_vals), self.c_src.I_vals)
        for pixel_x in range(x_min, x_max):
            for pixel_y in range(y_min, y_max):
                row = pixel_y + pixel_x * num_y_pixels
                self.responses[row] = self._get_response(pixel_x, pixel_y, row)
        return np.linalg.lstsq(self.sums, self.responses, rcond = None)

solver = SubpixelCoeffSolver()
x = solver.solve(0.6, 2.0, 0.01, 0.01, 10,10, int(6e3), int(6e3), 6,6, 3,3, -5, 5, -5 ,5)
