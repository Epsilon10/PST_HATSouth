from ctypes import cdll, c_double, c_uint
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from dll_loader import PSFDllLoader
from responses import get_responses_for_spot
from scipy.interpole import RegularGridInterpolator

class SubpixelCoeffSolver:
    """
    Solves for the coeffecients of the subpixel sensitivity map .
    Args:
        - wavelength
        - focal_ratio 
        - variance_x: standard deviation of the x coordinates of the spot centers
        - variance_y: standard deviation of the y coordinates of the spot centers
        - subpixel_res_x: number of subpixels to split each pixel along the x axis
        - subpixel_res_y: number of subpixels to split each pixel along the y axis
        - map_res_x: resolution of the ifft along the x axis
        - map_res_y: resolution of the ifft along the y axis
        - max_x_wavenumber: boundary of the fourier transform along x
        - max_y_wavenumber: boundary of the fourier transform along y
        - x_0: x coordinate of the spot center
        - y_0: y coordinate of the spot center
        - x_min: lower x coordinate of pixel
        - x_max: upper x coordinate of pixel
        - y_min: lower y coordinate of pixel
        - y_max: upper y coordinate of pixel
        - pixel_width: width of pixel
        - pixel_height: height of pixel
    """
    def __init__(self, wavelength, focal_ratio, variance, subpixel_res, 
        map_res, max_wavenumber,spots, all_pixels):
        c_src = PSFDllLoader("./libpsf.so")
        c_src.execute(
            wavelength, focal_ratio, 
            variance[0], variance[1],
            subpixel_res[0], subpixel_res[1],
            map_res[0], map_res[1],
            max_wavenumber[0], max_wavenumber[1]
        )
        self.all_pixels = all_pixels
        self.spots = spots
        self.intensities = RegularGridInterpolator((c_src.x_vals,c_src.y_vals), c_src.I_vals)
    
    def iterate(spot, solve_for_amplitude, amplitude=None, subpixel):
        responses, median_stddev, pixel_coords = get_responses_for_spot(spot, all_pixels)
        rotate_vector = numpy.ones_like(responses).T
        intensity_matrix = self.intensities((pixel_coords['x'], pixel_coords['y']))

        if solve_for_amplitude:
            u, s, vh = numpy.linalg.svd(rotate_vector, full_matrices = True)
            amplitudes = numpy.linalg.lstsq(intensity_matrix, 
            return amplitudes
        else:
            pass

    def solve_for_amplitude(spot, spl_map=None):
        responses, median_stddev, pixel_coords = get_responses_for_spot(spot, all_pixels)
        if spl_map is None:
            spl_map = numpy.ones_like(responses)
        
        stddev = 


