from ctypes import cdll, c_double, c_uint
import numpy
from scipy.interpolate import RegularGridInterpolator
from dll_loader import PSFDllLoader
from responses import get_responses_for_spot
from scipy.interpolate import RegularGridInterpolator
from utils import read_spot_sources, get_pixel_values
wavelength = 0.6
focal_ratio = 11.0
s_x = s_y = 1.0
subpixel_res_x = 10
subpixel_res_y = 10
map_res_x = 6000
map_res_y = 6000
max_wavenumbner_x = 6
max_wavenumber_y = 6

dx = 1.0 / subpixel_res_x
dy = 1.0 / subpixel_res_y

PIXEL_SIZE = 9.0 # microns

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
        map_res, max_wavenumber,spots, all_pixels, all_pixel_stddevs):
        c_src = PSFDllLoader("./libpsf.so")
        c_src.execute(
            wavelength, focal_ratio, 
            variance[0], variance[1],
            subpixel_res[0], subpixel_res[1],
            map_res[0], map_res[1],
            max_wavenumber[0], max_wavenumber[1]
        )
        self.all_pixels = all_pixels
        self.all_pixel_stddevs = all_pixel_stddevs
        self.spots = spots
        self.intensities = RegularGridInterpolator((c_src.x_vals / PIXEL_SIZE,c_src.y_vals / PIXEL_SIZE), c_src.I_vals)
        # print(self.intensities((1,1)))
    
    def iterate(self, spot, solve_for_amplitude, amplitude=None, spl_map=None):
        responses = get_responses_for_spot(spot, self.all_pixels, self.all_pixel_stddevs)
        # X, Y = numpy.meshgrid(responses.pixel_coords_x[:-1], responses.pixel_coords_y[:-1])
        intensity_matrix = numpy.empty((responses.data.shape[0], subpixel_res_x * subpixel_res_y))

        for idx, pixel in enumerate(responses.pixel_coords):
            intensity_matrix[idx,:] = self._get_pixel_intensities(pixel)

        stddev = numpy.sqrt(numpy.square(responses.pixel_stddevs) + responses.median_stddev**2)
        weighted_intensities = intensity_matrix / stddev[:,None]
        weighted_responses = responses.data / stddev

        if solve_for_amplitude:
            if spl_map is None:
                spl_map = numpy.full((subpixel_res_x * subpixel_res_y,), 1)

            amplitude = numpy.linalg.lstsq(weighted_intensities.dot(spl_map).reshape(responses.data.size,1), weighted_responses, rcond=None)
            return amplitude[0][0]
        else:
            assert amplitude is not None
            rotate_vector = numpy.full((1,intensity_matrix.shape[1]),1)
            u, s, vh = numpy.linalg.svd(rotate_vector, full_matrices = True)
            new_intensities = (weighted_intensities * amplitude).dot(vh[:,1:])

            rotated_spl_map = numpy.linalg.lstsq(new_intensities, weighted_responses - weighted_intensities.dot(vh[:,0]), rcond=None)
            spl_map = vh[:,1:].dot(rotated_spl_map[0])
            return spl_map
        
    def _get_pixel_intensities(self,pixel):
        half_dx = dx/2.0
        half_dy = dy/2.0
        
        pixel_x, pixel_y = pixel

        start_x, end_x = pixel_x + half_dx, 1.0 + pixel_x + half_dx
        start_y, end_y = pixel_y + half_dy, 1.0 + pixel_y + half_dy

        x_coords, y_coords = abs(numpy.arange(start_x, end_x, dx)), abs(numpy.arange(start_y, end_y, dy))
        x_coords = x_coords[:subpixel_res_x]
        y_coords = y_coords[:subpixel_res_y]
        #print(x_coords.size, y_coords.size)

        spl_x, spl_y = numpy.meshgrid(x_coords,y_coords)
        #print(spl_x, spl_y)

        return self.intensities((spl_x, spl_y)).flatten()



    def solve(self):
        test_amp = self.iterate(self.spots[0], True)
        test_spl = self.iterate(self.spots[0], False, test_amp)
        print("AMPLITUDE 1st ITER: ", test_amp)
        print("SPL MAP 1st ITER: ", test_spl)

spot_sources = read_spot_sources("/home/epsilon/Documents/photometry/images/projected/calibrated_nontilted_scan_projected_0068.fistar", 4,5)[:51]
pixel_values, pixel_stddevs = get_pixel_values("/home/epsilon/Documents/photometry/images/cal/calibrated_nontilted_scan_0068.fits.fz")

solver = SubpixelCoeffSolver(
    wavelength, focal_ratio, (s_x, s_y), (subpixel_res_x,subpixel_res_y), (map_res_x, map_res_y), (max_wavenumbner_x,max_wavenumber_y),
    spot_sources, pixel_values, pixel_stddevs
)

solver.solve()





