from astropy.io import fits
import numpy
import scipy
import scipy.spatial
from math import sqrt
from background_extractor import get_adjusted_responses 
from utils import read_spot_sources, get_pixel_values, clamp
from dataclasses import dataclass
import matplotlib.pyplot as plt

spot_coord_type=[('x', float), ('y', float)]

@dataclass(unsafe_hash=True)
class PixelResponses():
    data: numpy.ndarray
    median_stddev: float
    pixel_coords: numpy.ndarray
    pixel_stddevs: numpy.ndarray

def get_surrounding_pixels(source,pixels_shape, radius=5):
    x,y = source
    x_flr = int(x)
    y_flr = int(y)
    offset = radius + 1

    radius_squared = radius**2

    y_min = clamp(y_flr - offset, 0, pixels_shape[0] - 1)
    y_max = clamp(y_flr + offset, 0, pixels_shape[0] - 1)

    x_min = clamp(x_flr - offset, 0, pixels_shape[1] - 1)
    x_max = clamp(x_flr + offset, 0, pixels_shape[1] - 1)
    
    subgrid_shape = (y_max - y_min + 1, x_max - x_min + 1)

    indices = numpy.empty((81,), dtype=[('x', int), ('y', int)])

    idx = 0

    #print(f"S shape {subgrid_shape[0], subgrid_shape[1]}")

    for pixel_y in range(0, subgrid_shape[0]):
        for pixel_x in range(0, subgrid_shape[1]):
            img_y = y_min + pixel_y + 0.5
            img_x = x_min + pixel_x + 0.5
            dist_squared = (img_x - x)**2 + (img_y - y)**2
            if dist_squared <= radius_squared:
                indices[idx] = (img_x - 0.5, img_y -0.5)
                idx += 1
    return indices[:idx]

def get_center_row(source, radius):
    x,y = source
    row = round(y)
    indices = numpy.empty((radius * 2 + 1,),dtype=[('x', int), ('y', int)])
    indices['x'] = numpy.arange(int(x) - radius, int(x) + radius + 1)
    indices['y'] = row
    return indices

def get_responses_for_spot(spot, all_pixel_values, all_pixel_stddevs):
    matched_indices = get_center_row(spot, 4)

    raw_responses = all_pixel_values[matched_indices['y'], matched_indices['x']]
    pixel_stddevs = all_pixel_stddevs[matched_indices['y'], matched_indices['x']]

    pixel_coords_x = matched_indices['x'] - spot[0]
    pixel_coords_y = matched_indices['y'] - spot[1]

    response_data, median_stddev = get_adjusted_responses(spot, raw_responses, all_pixel_values, all_pixel_stddevs)

    pixel_coords = numpy.column_stack((pixel_coords_x, pixel_coords_y))

    return PixelResponses(
        data=response_data,
        median_stddev=median_stddev,
        pixel_coords=pixel_coords,
        pixel_stddevs=pixel_stddevs
    )


def test():
    spot_sources = read_spot_sources("/home/epsilon/Documents/photometry/images/projected/calibrated_nontilted_scan_projected_0068.fistar", 4,5)[:51]
    pixel_values, pixel_stddevs = get_pixel_values("/home/epsilon/Documents/photometry/images/cal/calibrated_nontilted_scan_0068.fits.fz")
    for spot in spot_sources:
        responses = get_responses_for_spot(spot, pixel_values, pixel_stddevs)
        #print("ADJ RESP:", responses)
        #print("PIXEL COORDS:", pixel_coords)
        print(responses.median_stddev)


if __name__ == "__main__":
    test()