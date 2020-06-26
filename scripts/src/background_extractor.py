import itertools
import numpy
from math import sqrt
from dataclasses import dataclass
from utils import clamp
MAX_ITER = 10

def get_pixels_in_range(source, bound, pixels_shape):
    x,y = source

    x_flr = int(source[0])
    y_flr = int(source[1])

    lo_bound, hi_bound = bound

    lo_squared, hi_squared = lo_bound**2, hi_bound**2

    offset = hi_bound + 1
    print(pixels_shape)
    y_min = clamp(y_flr - offset, 0, pixels_shape[0] - 1)
    y_max = clamp(y_flr + offset, 0, pixels_shape[0] - 1)

    x_min = clamp(x_flr - offset, 0, pixels_shape[1] - 1)
    x_max = clamp(x_flr + offset, 0, pixels_shape[1] - 1)
    
    subgrid_shape = (y_max - y_min + 1, x_max - x_min + 1)

    indices = numpy.empty((448,), dtype=[('x', int), ('y', int)])

    #print(pixel_values)
    idx = 0

    for pixel_y in range(0, subgrid_shape[0]):
        for pixel_x in range(0, subgrid_shape[1]):
            img_y = y_min + pixel_y
            img_x = x_min + pixel_x
            dist_squared = (img_x - x)**2 + (img_y - y)**2
            if dist_squared <= hi_squared and dist_squared >= lo_squared:
                indices[idx] = (img_x, img_y)
                idx += 1

    return indices[:idx]
    

def get_adjusted_responses(source, responses, pixel_values, pixel_stddevs):
    matched_indices = get_pixels_in_range(source=source, bound=(6,13), pixels_shape=pixel_values.shape)
    bg_pixels = pixel_values[matched_indices['y'], matched_indices['x']]
    bg_pixel_stddevs = pixel_stddevs[matched_indices['y'], matched_indices['x']]

    bg_pixels = numpy.sort(bg_pixels)
    num_bg_pixels = bg_pixels.shape[0]

    new_bg_pixels = None
    median_stddev = None
    prev_bg_pixels = bg_pixels

    for i in range(MAX_ITER):
        new_bg_pixels, median_stddev = _get_adjusted_responses(prev_bg_pixels)

        if prev_bg_pixels.shape[0]==new_bg_pixels.shape[0]:
            break
        
        prev_bg_pixels = new_bg_pixels

    median = numpy.median(new_bg_pixels)

    return responses - median, median_stddev

def _get_adjusted_responses(bg_pixels, N=5):
    num_bg_pixels = bg_pixels.shape[0]

    median = numpy.median(bg_pixels)
    # root(mean(square(pixels - median))/N)
    stddev = numpy.sqrt(numpy.mean(numpy.square(bg_pixels - median)))
    print("PIXEL STD", stddev)

    return bg_pixels[abs(bg_pixels - median) < stddev * N], stddev / sqrt(num_bg_pixels)