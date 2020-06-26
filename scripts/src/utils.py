import numpy
from astropy.io import fits

spot_coord_type=[('x', float), ('y', float)]

def read_spot_sources(filename, xcol, ycol) :
    """Reads the projected sources from a fistar file"""

    result=[]
    with open(filename, 'r') as srcfile :
        for line in srcfile :
            if line.strip()[0]=='#' : continue
            entries=line.split()
            result.append((float(entries[xcol]), float(entries[ycol])))
    return numpy.asarray(result, dtype=spot_coord_type)

def get_pixel_values(filename, compressed=True):
    """Returns the pixel values from a compressed fits file (fits.fz), data first, then stddev"""
    hdu_list = fits.open(filename)
    idx = 1 if compressed else 0
    return hdu_list[idx].data, hdu_list[idx + 1].data

def clamp(value, min_v, max_v):
    return max(min_v, min(value, max_v))