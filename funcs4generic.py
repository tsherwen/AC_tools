#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Generic functions for use with GEOS-Chem/Data Analysis.

Use help(<name of function>) to get details on a particular function.

NOTE(S):
 - This module is underdevelopment vestigial/inefficient code is being removed/updated.
 - Where external code is used credit is given.
"""

# ----------------------------- Section 0 -----------------------------------
# -------------- Required modules:

# --- compatibility with both python 2 and 3


# -- I/O / Low level
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pandas import DataFrame

# -- time
import time
import datetime as datetime

# -- math
from math import radians, sin, cos, asin, sqrt, pi, atan2

# --  This needs to be updated, imports should be specific and in individual functions
# import tms modules with shared functions
from . funcs4core import *
from . funcs_vars import *

# -------------------------- Section 7 -------------------------------
# -------------- Generic Processing
#

# -------------
# X.XX -  Split arrays into chunks (e.g. days)
# -------------


def chunks(l, n):
    """
    Split list in chunks - useful for controlling memory usage
    """
    if n < 1:
        n = 1
    return [l[i:i + n] for i in range(0, len(l), n)]

# -------------
# X.XX - Get len of file
# -------------


def file_len(fname):
    """
    Get length of file
    """
    return sum(1 for line in open(fname))

# -------------
# X.XX - My round - credit: Alok Singhal
# -------------


def myround(x, base=5, integer=True, round_up=False):
    """
    Round up values - mappable function for pandas processing
    NOTES:
     - credit: Alok Singhal
    """
#    round_up=True # Kludge - always set this to True

    # round to nearest base number
    rounded = base * round(float(x)/base)

    if round_up:
        # ensure rounded number is to nearest next base
        if rounded < x:
            rounded += base

    if integer:
        return int(rounded)
    else:
        return rounded

# --------------
# X.XX - Count number of files/the filenames that contain certain phrase/number conbinations (e.g. dates) - tms
# -------------


def counter_directory_contains_files(model_path, must_contain):
    """
    Count number of files in directory
    """
    model_ouput_file_counter = len(glob.glob1(model_path, must_contain))
    return model_ouput_file_counter


# ----
# X.XX - Rename vars/strs in files
# ----
def replace_strs_in_files(wd, input_str, output_str, debug=False):
    """
    replace text in files
    """
    print((wd, input_str, output_str))
    for f in os.listdir(wd):
        if not f.startswith('.'):
            print(f)
            os.rename(wd + f, wd + f.replace(input_str, output_str))
            print((f.replace(input_str, output_str)))

# ----
# X.XX - Get X and Y coordinates for a given grid - Credit: Eric Sofen
# ----


def get_xy(Lon, Lat, lon_edges, lat_edges, debug=False):
    """
    Takes lon,lat values for a point (e.g. station) and arrays of longitude
    and latitudes indicating edges of gridboxes.

    Returns the index of the gridbox corresponding to this point.
    NOTES:
     - Credit: Eric Sofen
     - Could be easily extended to return data values corresponding to points.
    """
    hasobs, lonrange, latrange = np.histogram2d(
        [Lon], [Lat], [lon_edges, lat_edges])
    gridindx, gridindy = np.where(hasobs >= 1)
    if not gridindx:
        if debug:
            print(('Lat, lon outside of x,y range.  Assigning -1 for', Lon, Lat))
        return -1, -1
    else:
        #print Lon, Lat, gridindx, gridindy
        return gridindx[0], gridindy[0]

# --------
# X.XX - Save as pdf.
# --------


def plot2pdf(title='new_plot', fig=None, rasterized=True, dpi=320,
             justHH=False, no_dstr=True, save2png=True, wd=None,
             save2eps=False, transparent=True, debug=False):
    """
    Save figures (e.g. matplotlib) to pdf file
    """
    # set save directory ( using default directory dictionary )
    if isinstance(wd, type(None)):
        wd = './'

    # Set pdf name
    if justHH:
        date_str = time.strftime("%y_%m_%d_%H")
    else:
        date_str = time.strftime("%y_%m_%d_%H_%M")
    if no_dstr:
        npdf = wd+title
    else:
        npdf = wd+date_str+'_'+title

    # setup pdf
    pdf = PdfPages(npdf + '.pdf')

    # Rasterise to save space?
    if rasterized:
        plt.gcf().set_rasterized(True)

    # save and close
    file_extension = 'PDF'
    pdf.savefig(dpi=dpi, transparent=transparent)
    pdf.close()

    # Also export to png and eps?
    if save2png:
        file_extension += '/PDF'
        plt.savefig(npdf+'.png', format='png', dpi=dpi,
                    transparent=transparent)
    if save2eps:
        file_extension += '/EPS'
        plt.savefig(npdf+'.eps', format='eps', dpi=dpi,
                    transparent=transparent)
    print((file_extension+' saved & Closed as/at: ', npdf))

# --------
# X.XX - Save as mulitple page pdf.
# --------


def plot2pdfmulti(pdf=None, title='new_plot', rasterized=True, wd=None,
                  dpi=320, open=False, close=False, justHH=False, no_dstr=True):
    """
    Save figures (e.g. matplotlib) to pdf file with multiple pages
    """
    # set save directory ( using default directory dictionary )
    if isinstance(wd, type(None)):
        wd = './'

    # Set pdf name
    if justHH and (not no_dstr):
        date_str = time.strftime("%y_%m_%d_%H")
    if not no_dstr:
        date_str = time.strftime("%y_%m_%d_%H_%M")
    if no_dstr:
        npdf = wd+title+'.pdf'
    else:
        npdf = wd+date_str+'_'+title+'.pdf'

    # If 1st call ( open ==True), setup pdf
    if open:
        pdf = PdfPages(npdf)
        print(('pdf opened @: {}'.format(npdf)))
        return pdf

    # Rasterise to save space?
    if rasterized:
        plt.gcf().set_rasterized(True)

    # save and close or keep open to allow additions of plots
    if close:
        pdf.close()
        print(('PDF saved & Closed as/at: ', npdf))
    else:
        pdf.savefig(dpi=dpi)
        print(('pdf is still open @: {}'.format(npdf)))


# --------------
# X.XX -  Return mesh and indics for which obs are within
# -------------
def obs2grid(glon=None, glat=None, galt=None, nest='high res global',
             sites=None, debug=False):
    """
    values that have a given lat, lon and alt
    """
    if isinstance(glon, type(None)):
        glon, glat, galt = get_latlonalt4res(nest=nest, centre=False,
                                             debug=debug)

    # Assume use of known CAST sites... unless others given.
    if isinstance(sites, type(None)):
        loc_dict = get_loc(rtn_dict=True)
        sites = list(loc_dict.keys())

    # pull out site location indicies
    indices_list = []
    for site in sites:
        lon, lat, alt = loc_dict[site]
        vars = get_xy(lon,  lat, glon, glat)
        indices_list += [vars]
    return indices_list

# --------
# X.XX - Sort (GAW) sites by Latitude and return
# --------


def sort_sites_by_lat(sites):
    """
    Order given list of GAW sties by latitudes
    """

    # Get info
    vars = [gaw_2_loc(s) for s in sites]  # lat, lon, alt, TZ

    # Sort by lat, index orginal sites list and return
    lats = [i[0] for i in vars]
    slats = sorted(lats)[::-1]
    return [sites[i] for i in [lats.index(ii) for ii in slats]]

# --------
# X.XX - Find nearest
# --------


def find_nearest(array, value):
    """
    Find nearest number in array to given value.

#    NOTES:
#    - Adapted from (credit:) HappyLeapSecond's Stackoverflow answer.
#    ( http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array )
#    """
    idx = (np.abs(array-value)).argmin()
    return idx

# --------
# X.XX - Get suffix for number
# --------


def get_suffix(n):
    """  Add the appropriate suffix (th/st/rd) to any number given  """
    def ordinal(n): return "%d%s" % (
        n, "tsnrhtdd"[(n/10 % 10 != 1)*(n % 10 < 4)*n % 10::4])
    return ordinal(n)

# --------
# X.XX - Get shortest distance on a sphere ( for finding closest point )
# --------


def get_shortest_in(needle, haystack, r_distance=False):
    """
    needle is a single (lat,long) tuple. haystack is a numpy array to find the point in
    that has the shortest distance to needle

    NOTES:
     - adapted from stackoverflow (Credit: jterrace):
    (http://stackoverflow.com/questions/6656475/python-speeding-up-geographic-comparison)
    """
    # set Earth's radius
    earth_radius_miles = 3956.0

    # convert to radians
    dlat = np.radians(haystack[:, 0]) - radians(needle[0])
    dlon = np.radians(haystack[:, 1]) - radians(needle[1])

    # get the vector
    a = np.square(np.sin(dlat/2.0)) + cos(radians(needle[0])) * np.cos(
        np.radians(haystack[:, 0])) * np.square(np.sin(dlon/2.0))

    # convert to distance
    great_circle_distance = 2 * np.arcsin(np.minimum(np.sqrt(a),
                                                     np.repeat(1, len(a))))
    d = earth_radius_miles * great_circle_distance

    if r_distance:
        return np.min(d)
    # return the index
    else:
        return list(d).index(np.min(d))

# --------
# X.XX - Get logarithmically spaced integers
# --------


def gen_log_space(limit, n):
    """
    Get logarithmically spaced integers

    NOTES:
    - credit: Avaris
    ( http://stackoverflow.com/questions/12418234/logarithmically-spaced-integers)
    """
    result = [1]
    if n > 1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result) < n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value
            # by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining
            # values will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    return np.array([round(x)-1 for x in result], dtype=np.uint64)

# --------
# X.XX - Get indices in array where change in value of x occurs
# --------


def get_arr_edge_indices(arr, res='4x5', extra_points_point_on_edge=None,
                         verbose=True, debug=False):
    """
    Find indices in a lon, lat (2D) grid, where value does not equal a given
    value ( e.g. the edge )
    """

    if verbose:
        print(('get_arr_edge_indices for arr of shape: ', arr.shape))

    # initialise variables
    lon_c, lat_c, NIU = get_latlonalt4res(res=res, centre=True)
    lon_e, lat_e, NIU = get_latlonalt4res(res=res, centre=False)
    lon_diff = lon_e[-5]-lon_e[-6]
    lat_diff = lat_e[-5]-lat_e[-6]
    nn, n, = 0, 0
    last_lat_box = arr[nn, n]
    coords = []
    last_lon_box = arr[nn, n]
    need_lon_outer_edge, need_lat_outer_edge = False, False
    if debug:
        print((lon_e, lat_e))

    # ---- Loop X dimension ( lon )
    for nn, lon_ in enumerate(lon_c):

        # Loop Y dimension ( lat ) and store edges
        for n, lat_ in enumerate(lat_c):

            if debug:
                print((arr[nn, n], last_lat_box, last_lon_box,
                       arr[nn, n] == last_lat_box, arr[nn, n] == last_lon_box))

            if arr[nn, n] != last_lat_box:

                # If 1st lat, selct bottom of box
                point_lon = lon_e[nn]+lon_diff/2
                if need_lat_outer_edge:
                    point_lat = lat_e[n+1]
                else:
                    point_lat = lat_e[n]
                    need_lat_outer_edge = True
                need_lat_outer_edge = False

                # Add mid point to cordinates list
                if isinstance(extra_points_point_on_edge, type(None)):
                    mid_point = [point_lon, point_lat]
                    coords += [mid_point]

                # Add given number of points along edge
                else:
                    coords += [[lon_e[nn]+(lon_diff*i), point_lat] for i in
                               np.linspace(0, 1, extra_points_point_on_edge,
                                           endpoint=True)]

            # temporally save the previous box's value
            last_lat_box = arr[nn, n]

    # ---- Loop Y dimension ( lat )
    for n, lat_ in enumerate(lat_c):

        if debug:
            print((arr[nn, n], last_lat_box, last_lon_box,
                   arr[nn, n] == last_lat_box, arr[nn, n] == last_lon_box))
        # Loop X dimension ( lon ) and store edges
        for nn, lon_ in enumerate(lon_c):

            # If change in value at to list
            if arr[nn, n] != last_lon_box:
                point_lat = lat_e[n]+lat_diff/2

                # Make sure we select the edge lon
                if need_lon_outer_edge:
                    point_lon = lon_e[nn+1]
                else:
                    point_lon = lon_e[nn]
                    need_lon_outer_edge = True
                need_lon_outer_edge = False

                # Add mid point to coordinates list
                if isinstance(extra_points_point_on_edge, type(None)):
                    mid_point = [point_lon, point_lat]
                    coords += [mid_point]

                # Add given number of points along edge
                else:
                    coords += [[point_lon, lat_e[n]+(lat_diff*i)] for i in
                               np.linspace(0, 1, extra_points_point_on_edge,
                                           endpoint=True)]

            # temporally save the previous box's value
            last_lon_box = arr[nn, n]

    return coords

# --------
# X.XX - Get data binned by uniques days in data
# --------


def split_data_by_days(data=None, dates=None, day_list=None,
                       verbose=False, debug=False):
    """
    Takes a list of datetimes and data and returns a list of data and
    the bins ( days )
    """
    if verbose:
        print('split_data_by_days called')

    # Create DataFrame of Data and dates
    df = DataFrame(data, index=dates, columns=['data'])
    # Add list of dates ( just year, month, day ) <= this is mappable, update?
    df['days'] = [datetime.datetime(*i.timetuple()[:3]) for i in dates]
    if debug:
        print(df)

    # Get list of unique days
    if isinstance(day_list, type(None)):
        day_list = sorted(set(df['days'].values))
    # Loop unique days and select data on these days
    data4days = []
    for day in day_list:
        print((day, df[df['days'] == day]))
        data4days += [df['data'][df['days'] == day]]
    # Just return the values ( i.e. not pandas array )
    data4days = [i.values.astype(float) for i in data4days]
    print([type(i) for i in data4days])
#    print data4days[0]
#    sys.exit()

    if debug:
        print(('returning data for {} days, with lengths: '.format(
            len(day_list)), [len(i) for i in data4days]))

    # Return as list of days (datetimes) + list of data for each day
    return data4days, day_list

# -------------------------- Section 2 -------------------------------
# -------------- Maskes for data analysis
#

# --------
# X.XX - Ocean mask
# --------


def ocean_unmasked(res='4x5', debug=False):
    """
    Get ocean mask from GEOS-Chem LWI ( at given resolution)

    NOTES:
     - this currently returns mask with all except ocean masked
     - NEEDS UPDATE
    """

    from .funcs4GEOSC import get_land_map
    if debug:
        print(('ocean_mask called for: ', res))

    # Create a mask from land/water/ice indices
    m = np.ma.masked_not_equal(get_land_map(res=res), 0)
    if debug:
        print((mask, mask.shape))
    return m.mask

# --------
# X.XX - Land mask
# --------


def land_unmasked(res='4x5', debug=False):
    """
    Get land mask from GEOS-Chem LWI ( at given resolution)
    """
    from .funcs4GEOSC import get_land_map  # Kludge, use GEOS-Chem LWI

    # Create a np.ma mask
    if debug:
        print(('land_mask called for: ', res))
    m = np.ma.masked_not_equal(get_land_map(res=res), 1)
    if debug:
        print((mask, mask.shape))
    return m.mask

# --------
# X.XX - Ice mask
# --------


def ice_unmasked(res='4x5', debug=False):
    """
    Get mask for all areas but ice surface. Ice mask is from GEOS-Chem
    LWI ( at given resolution)

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a np.ma mask
    m = np.logical_not((land_unmasked(res)*ocean_unmasked(res)))
    if debug:
        print((mask, mask.shape))
    return m

# --------
# X.XX - Surface mask
# --------


def surface_unmasked(res='4x5', trop_limit=False, mask2D=False,
                     debug=False):
    """
    Get mask for all areas but surface

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a np.ma mask
    m = np.ma.array(np.zeros(get_dims4res(res)), mask=False)
    m[..., 0] = 1
    m = np.ma.masked_not_equal(m, 1)

    if trop_limit:
        m = m[..., :38]
    if debug:
        print((mask, mask.shape))

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask
# --------
# X.XX  - Tropical Mask
# --------


def tropics_unmasked(res='4x5', saizlopez=False, pradosroman=False,
                     mask2D=False):
    """
    Get mask for all areas but tropics

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    if saizlopez or pradosroman:
        lats = np.arange(-20, 20, 1)
    else:
        lats = np.arange(-22, 22, 1)
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1

    # Create a np.ma mask and return solely mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask

# --------
# X.XX  - Mid Lats Mask
# --------


def mid_lats_unmasked(res='4x5', saizlopez=False, pradosroman=False,
                      mask2D=False):
    """
    Get mask for all areas but mid latitudes

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    if saizlopez:
        lats = np.concatenate((np.arange(-50, -20, 1), np.arange(20, 50, 1)))
    elif pradosroman:
        lats = np.concatenate((np.arange(-60, -20, 1), np.arange(20, 61, 1)))
    if ((not saizlopez) and (not saizlopez) and (res == '4x5')):
        lats = np.concatenate((np.arange(-50, -26, 1), np.arange(26, 51, 1)))
    if ((not saizlopez) and (not saizlopez) and (res == '2x2.5')):
        lats = np.concatenate((np.arange(-50, -24, 1), np.arange(24, 51, 1)))

    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1

    # Create a np.ma mask and return solely mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask
# --------
# X.XX - 40N to 40S Mask
# --------


def mask_lat40_2_40(res='4x5', mask2D=False):
    """
    Get mask for all areas but 40S-40N latitudes (40 to 40 deg latitude)

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    lats = np.arange(-42, 42, 1)  # use 42 due to 4x5 grid
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1

    # Create a np.ma mask and return solely mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask
# --------
# X.XX - Extra Tropics Mask
# --------


def extratropics_unmasked(res='4x5', mask2D=False):
    """ Get mask for all areas but extratropics
    NOTES:
        - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    lats = np.concatenate((np.arange(-89, -26, 1), np.arange(26, 90, 1)))
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1
    # Create a np.ma mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask
# --------
# X.XX - Create unmask array ( for ease of dataprocessing )
# --------


def all_unmasked(res='4x5', mask2D=False):
    """ Get un-masked mask of size GEOS-Chem dimensions
    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    m = np.ma.array(np.ones(get_dims4res(res)), mask=False)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask

# --------
# X.XX - Maskes Regions by Pressure
# --------


def mask_3D(hPa, sect, MBL=True, res='4x5', extra_mask=None,
            M_all=False, use_multiply_method=True, trop_limit=False,
            verbose=True, debug=False):
    """
    Creates Maskes by pressure array  (required shape: 72,46,47),
    with conditions (lower and upper bounds) set by given cases for
    MBL, UT, FT

    Parameters
    -------
    sect (Str): section of the atmosphere of interest (e.g. MBL, UT...)
    hPa (array): array for pressures ( in hPa)
    MBL (boolean): apply a mask for the marine boundary layer
    res (str): the resolution of required output/input arrays (e.g. '4x5' )
    use_multiply_method (boolean): Create arrays of ones and zeros
    trop_limit (boolean): limit 3D arrays to troposphere
    debug (boolean): legacy debug option, replaced by python logging
    verbose (boolean): legacy debug option, replaced by python logging
    extra_mask (str): name of additional region (e.g. ocean) to mask
    M_all (boolean): apply oceanic masking to all regions

    Returns
    -------
    (np.ma.mask)

    NOTES:
     - originally written to generate masks for mulitplication
    (i.e. use_multiply_method = True ), but can also be use to make
    more pythonic masks ( if use_multiply_method=False )
    """
    if verbose:
        print(('mask_3D called for sect={}, use_multiply_method={}'.format(
            sect, use_multiply_method) + ', M_all={}, and '.format(M_all) +
            'with debug: {}, verbose:{}'.format(sect, debug, verbose)))

    # Get atmospheric region as case defining lower and upper bounds
    cases = {
        'BL': [1200., 900.], 'MBL': [1200., 900.], 'FT': [900., 350.],
        'UT': [350., 75.], 'All': [1200., 75.]
    }
    l, h = cases[sect]

    # --- Mask between upper and lower values
    m = np.ones(get_dims4res(res))
    m[(hPa >= l)] = 0
    m[(hPa < h)] = 0
    logging.debug('Sect={}, l={}, h={}'.format(sect, l, h))
    logging.debug('{}'.format(
        *[[i.min(), i.max(), i.mean(), i.sum(), i.shape] for i in [m]]))

    # Mask off the 'sect' area that still equals 1
    m = np.ma.masked_equal(m, 1)

    # --- Remove above the "chemical tropopause" from GEOS-Chem (v9-2)
    if trop_limit:
        m = m[..., :38]

    if not isinstance(extra_mask, type(None)):
        # only consider MBL
        if (MBL and sect == 'BL') or (sect == 'MBL') or M_all:
            if use_multiply_method:
                return m.mask * extra_mask * land_unmasked(res)
            else:
                print("WARNING: needs 3D arrays, use 'mask_all_but' instead")
                sys.exit()
        else:
            if use_multiply_method:
                return m.mask * extra_mask
            else:
                print("WARNING: needs 3D arrays, use 'mask_all_but' instead")
                sys.exit()

    # --- Only consider MBL (or MFT/MFT for Saiz-Lopez 2014 comparison)
    if (MBL and sect == 'BL') or (sect == 'MBL') or M_all:
        if use_multiply_method:
            return m.mask * land_unmasked(res)
        else:
            land_unmasked_ = mask_all_but('Land', mask3D=True, res=res,
                                          use_multiply_method=False, trop_limit=trop_limit)

            # MBL unmasked
            m = np.logical_or(np.logical_not(m.mask),
                              np.logical_not(land_unmasked_))

            # Invert as function expected to return oposite
            m = np.logical_not(m)

            return m   # m is a mask here

    return m.mask


# --------
# X.XX  - Custom 2D (Lat) Mask
# --------
def lat2lat_2D_unmasked(lowerlat=None, higherlat=None, res='2x2.5',
                        debug=False):
    """
    Takes a lower and higher latitude value and then creates mask of area boxes outside
    given laitidue limits.
    """

    # Get vars
    lon_c, lat_c, NIC = get_latlonalt4res(res=res, centre=True)

    # mask between upper and lower values
    lats = [i for i in lat_c if ((i >= lowerlat) and (i < higherlat))]
    lats = [get_gc_lat(i, res=res) for i in lats]

    # fill all lat and lon True or False
    m = np.zeros(get_dims4res(res))[:, :, 0]
    print((m.shape, np.sum(m)))
    for i in lats:
        m[:, i] = 1
    m = np.ma.masked_not_equal(m, 1)
    return m.mask

# --------
# X.XX - South Pole mask
# --------


def southpole_unmasked(res='4x5', mask2D=False):
    """
    Return global np.ma with southpole masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    # adjust for resolution at grid start points at 62
    if res == '4x5':
        lats = np.arange(-89, -62, 1)  # define S pole as > 60S
    else:
        lats = np.arange(-89, -60, 1)  # define S pole as > 60S
#    lats = np.arange(-89, -80,1 ) # define S pole as > 80S
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1

    # Create a np.ma mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask

# --------
# X.XX - North Pole mask
# --------


def northpole_unmasked(res='4x5', mask2D=False):
    """
    Return global np.ma with northpole masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a dummy array of zeros
    m = np.zeros(get_dims4res(res))
    # adjust for resolution at grid start points at 62
    if res == '4x5':
        lats = np.arange(62, 90, 1)  # define N pole as > 60N
    else:
        lats = np.arange(60, 90, 1)  # define N pole as > 60N
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :] = 1

    # Create a np.ma mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask

# --------
# X.XX - North Hemisphere mask
# --------


def NH_unmasked(res='4x5', mask2D=False):
    """
    Return global np.ma with north hemisphere masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a dummy array of ones, with all locations masked
    m = np.ma.array(np.ones(get_dims4res(res)), mask=True)
    if res == '4x5':
        lats = np.arange(1, 91, 1)
    elif res == '2x2.5':
        lats = np.arange(0, 89, 1)
        print('CHECK (NH) mask for non 4x5 resolutions')
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :].mask = False
    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


# --------
# X.XX - South Hemisphere mask
# --------
def SH_unmasked(res='4x5', mask2D=False):
    """
    Return global np.ma with south hemisphere masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a dummy array of ones, with all locations masked
    m = np.ma.array(np.ones(get_dims4res(res)), mask=True)
    if res == '4x5':
        lats = np.arange(-89, 0, 1)
    if res == '2x2.5':
        lats = np.arange(-90, 0, 1)
        print('CHECK (SH) mask for non 4x5 resolutions')
    lats = [get_gc_lat(i, res=res) for i in lats]
    for i in lats:
        m[:, i, :].mask = False
    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


# --------
# X.XX - North Hemisphere mask
# --------
def location_unmasked(res='4x5', lat=None, lon=None, mask2D=False):
    """
    Return global np.ma.mask with all locations apart from that containing the
    location (lat, lon) masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a dummy array of ones, with all locations masked
    m = np.ma.array(np.ones(get_dims4res(res)), mask=True)
    # Get location index
    assert all([type(i) == float for i in (lat, lon)]
               ), 'lat & lon must be floats'
    lat_ind = get_gc_lat(lat=lat, res=res)
    lon_ind = get_gc_lon(lon=lon, res=res)
    # Unmask grid box for location
    m.mask[lon_ind, lat_ind] = False
    # Return 2D or 3D?
    if mask2D:
        return m[..., 0].mask
    else:
        return m.mask


# --------
# X.XX  - Get Analysis maskes
# --------
def get_analysis_masks(masks='basic',  hPa=None, M_all=False, res='4x5',
                       saizlopez=False, r_pstr=True, wd=None, trop_limit=True, mask4D=False,
                       use_multiply_method=True, debug=False):
    """
    Return list of mask arrays for analysis

    NOTES:
    - For comparisons with Saiz-Lopez et al. 2014, set M_all to True,
    this masks the UT and FT as well as the BL (therefore gives MBL, MFT, MUT)

    """
    # --- Check hPa has been provided as arg.
    if isinstance(hPa, type(None)):
        print('ERROR: Please provide array of hPa to get_analysis_masks')

    if masks == 'full':
        # ---- List of masks
        mtitles = [
            'All', 'Ocean', 'Land', 'Ice', 'All Sur.',  'Ocean Sur.', 'Land Sur.', 'Ice Sur.', 'NH', 'SH', 'Tropics', 'Ex. Tropics', 'Mid Lats',
            'Ocn. 50S-50N',  '50S-50N'
        ]
        tsects3D = ['MBL', 'BL', 'FT',  'UT']
        # get MBL, FT and UT maskes
        sects3D = [
            mask_3D(hPa, i, MBL=False, M_all=M_all, res=res,
                    use_multiply_method=use_multiply_method) for i in tsects3D
        ]
        # ---- Use non-pythonic mulitply method?
        if use_multiply_method:
            maskes = [mask_all_but(i, trop_limit=trop_limit, mask3D=True,
                                   use_multiply_method=True, res=res) for i in mtitles]
            # if comparison with saiz-lopez 2014,
            if M_all:
                ind = [n for n, i in enumerate(mtitles) if not ('MBL' in i)]
                for n in ind:
                    maskes[n] = maskes[n]*land_unmasked(res=res)
        # --- Use pythonic approach
        else:
            maskes = [mask_all_but(i, trop_limit=trop_limit, mask3D=True,
                                   use_multiply_method=False, res=res) for i in mtitles]
            # If not 'use_multiply_method', then invert hPa masks
            sects3D = [np.logical_not(i) for i in sects3D]
        # Add to mask and mask title lists
        maskes = maskes + sects3D
        mtitles = mtitles + tsects3D
        # Also create print strings...
        npstr = '{:<12}'*len(maskes)
        pstr = '{:<12,.2f}'*len(maskes)
    if masks == 'basic':
        tsects3D = ['All', 'MBL', 'BL', 'FT',  'UT']
        mtitles = [i+' (All)' for i in tsects3D] + \
            [i+' (Tropics)' for i in tsects3D] +   \
            [i+' (Mid Lats)' for i in tsects3D]
        # Standard maskes none, tropics, mid-lats (3)
        maskes = [
            np.logical_not(i) for i in (all_unmasked(res=res),
                                        tropics_unmasked(
                                            res, saizlopez=saizlopez),
                                        mid_lats_unmasked(res))]
        # Additional masks - tsects3D (4+1) * standard maskes (3)
        dmaskes = [
            [mask_3D(hPa, i, MBL=False, extra_mask=mask, M_all=M_all, res=res)
             for i in tsects3D] for mask in maskes]
        # Unpack and set as maskes (12) to list (3)
        dmaskes = [j for i in dmaskes for j in i]
        print([len(i) for i in (maskes, dmaskes, mtitles, tsects3D)])
        maskes = dmaskes
        print([len(i) for i in (maskes, dmaskes, mtitles, tsects3D)])
        # If comparison with saiz-lopez 2014 appli marine mask to all...
        if M_all:
            ind = [n for n, i in enumerate(mtitles) if not 'MBL' in i]
            for n in ind:
                maskes[n] = maskes[n]*land_unmasked(res=res)
        if debug:
            print([len(i) for i in (maskes, dmaskes, mtitles, tsects3D)])
        # Also create print strings...
        npstr = '{:<15}'*len(maskes)
        pstr = '{:<15,.2f}'*len(maskes)
    if masks == 'trop_regions':
        mtitles = ['BL', 'FT',  'UT']
        maskes = [mask_3D(hPa, i, M_all=M_all, MBL=False, res=res)[:, :, :38]
                  for i in mtitles]
        # Also create print strings...
        npstr = '{:<15}'*len(maskes)
        pstr = '{:<15,.2f}'*len(maskes)
    # Only consider the "chemical troposphere" - according v9-2
    if trop_limit:
        maskes = [i[:, :, :38] for i in maskes]
    # Create 4D array by concatenating through time dimension
    # ( assuming year long array of 1 months )
    if mask4D:
        for n, mask in enumerate(maskes):
            if any([(mask.shape[-1] == i) for i in [12]]):
                pass
            else:  # concatenate dimensions
                maskes[n] = np.concatenate([mask[..., None]]*12, axis=3)
    if r_pstr:
        return maskes, mtitles, npstr, pstr
    else:
        return maskes, mtitles

# --------
# X.XX -  Retrieve individual 4D mask of locations except region given
# --------


def mask_all_but(region='All', M_all=False, saizlopez=False,
                 res='4x5', trop_limit=True, mask2D=False, mask3D=False, mask4D=False,
                 use_multiply_method=True, lat=None, lon=None,
                 verbose=False, debug=False):
    """
    Mask selector for analysis. global mask provided for with given region
        unmasked

    Parameters
    -------
    res (str): the resolution if wd not given (e.g. '4x5' )
    M_all (boolean): maask all marine areas?
    saizlopez (boolean): use tropics definition from Saiz-Lopez er al 2014
    trop_limit (boolean): limit 4D arrays to troposphere
    mask2D/mask3D/mask4D(booolean): ensure mask returned is 2D/3D/4D
    use_multiply_method (boolean): return array of ones, that can be mulitpled
    through an array to set data to zero
    verbose (boolean): legacy debug option, replaced by python logging
    debug (boolean): legacy debug option, replaced by python logging
    loc (str): location
    lon, lat (float): lat/lon locations to leave nearest grid box unmasked

    Returns
    -------
    (np.ma.mask) or (np.array) (later if use_multiply_method==True)

    Notes
    -----
    "unmask_all" yeilds completely unmasked array
    function was oringialyl used to mulitple masks, however, this approch is
    unpythonic and therefore reccomended against.
    """
    logging.info('mask_all_but called for region {}'.format(region))
    # --- Setup cases...
    # ( except None, unmask_all and global to retrive no mask )
    case = {
        'Tropics': 0,
        'tropics': 0,
        'mid_lats': 1,
        'Mid Lats': 1,
        'Mid lats': 1,
        'south_pole': 2,
        'south pole': 2,
        'north_pole': 3,
        'north pole': 3,
        None: 4,
        'unmask_all': 4,
        'All': 4,
        'global': 4,
        # NEED TESTING ...
        'Extratropics': 5,
        'Ex. Tropics': 5,
        'Oceanic': 6,
        'Ocean': 6,
        'NH': 7,
        'SH': 8,
        'Ice': 10,
        'Land': 11,
        'lat40_2_40': 12,
        'Ocean Tropics': 13,
        'Oceanic Tropics': 13,
        'Ocn. Trop.': 13,
        'Land Tropics': 14,
        'All Sur.': 15,
        'surface': 15,
        'Ocean Sur.': 16,
        'Land Sur.': 17,
        'Ice Sur.': 18,
        'lat50_2_50': 19,
        '50S-50N': 19,
        #    'Oceanic lat50_2_50': 20,
        'Ocn. 50S-50N': 20,
        #     'South >60': 2,
        #      'North >60': 3
        'North Sea': 21,
        'Med. Sea': 22,
        'Mediterranean Sea': 22,
        'Black Sea': 23,
        'Irish Sea': 24,
        'Europe': 25,
        'EU': 25,
        #    'Surface BL': 26,
        'Land Tropics Sur.': 27,
        'Boreal Land': 28,
        'Alps':  29,
        'loc': 30,
        'location': 30,
        'France': 31,
    }[region]

    # --- This is a simple way of using masks ( as multiplers )
    # i.e. all (future) functions should have use_multiply_method=False
    # and not use the code below
    if use_multiply_method:  # Kludge
        print(('!'*50, 'WARNING: using mulitply method for masking. '))

        # For case, pull mask from case list
        if case == 0:
            mask = tropics_unmasked(res=res, saizlopez=saizlopez)
        elif case == 1:
            mask = mid_lats_unmasked(res=res)
        elif case == 2:
            mask = southpole_unmasked(res=res)
        elif case == 3:
            mask = northpole_unmasked(res=res)
        elif case == 4:
            #        mask = np.logical_not( all_unmasked( res=res ) )
            mask = all_unmasked(res=res)
        elif case == 5:
            mask = extratropics_unmasked(res=res)
        elif case == 6:
            mask = ocean_unmasked(res=res)
        elif case == 7:
            mask = NH_unmasked(res=res)
        elif case == 8:
            mask = SH_unmasked(res=res)
        elif case == 10:
            mask = ice_unmasked(res=res)
        elif case == 11:
            mask = land_unmasked(res=res)
        elif case == 12:
            mask = mask_lat40_2_40(res=res)
        elif case == 13:  # 'Oceanic Tropics'
            mask = np.ma.mask_or(ocean_unmasked(res=res),
                                 tropics_unmasked(res=res, saizlopez=saizlopez))
        elif case == 14:  # 'Land Tropics'
            mask = np.ma.mask_or(land_unmasked(res=res),
                                 tropics_unmasked(res=res, saizlopez=saizlopez))
        elif case == 15:  # 'All Sur.'
            mask = surface_unmasked(res=res)
        elif case == 16:  # 'Ocean Sur.'
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 ocean_unmasked(res=res))
        elif case == 17:  # 'Land Sur.':
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 land_unmasked(res=res))
        elif case == 18:  # 'Ice Sur.'
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 ice_unmasked(res=res))
        elif case == 19:  # '50S-50N'
            mask = lat2lat_2D_unmasked(lowerlat=-50, higherlat=50,
                                       res=res)
        elif case == 20:  # 'Ocn. 50S-50N'
            mask = np.ma.mask_or(lat2lat_2D_unmasked(lowerlat=-50,
                                                     higherlat=50, res=res), ocean_unmasked(res=res)[..., 0])
        elif case == 21:
            mask = get_north_sea_unmasked(res=res)
        elif case == 25:
            mask = get_EU_unmasked(res=res)
#        if case == 26:
#            mask = get_2D_BL_unmasked( res=res )
        elif case == 27:  # 'Land Tropics Sur.':
            tmp = np.ma.mask_or(surface_unmasked(res=res),
                                land_unmasked(res=res))
            mask = np.ma.mask_or(tmp, tropics_unmasked(res=res))
        else:
            print('WARNING - Mask not setup for case={}'.format(case))
            sys.exit()
        # Invert mask to leave exception unmasked if used to multiply
        mask = np.logical_not(mask)

    # --- This is a more pythonic way of using masks (Use as preference)
    else:
        # For case, pull mask from case list
        if case == 0:
            mask = tropics_unmasked(res=res, saizlopez=saizlopez)
        elif case == 1:
            mask = mid_lats_unmasked(res=res)
        elif case == 2:
            mask = southpole_unmasked(res=res)
        elif case == 3:
            mask = northpole_unmasked(res=res)
        elif case == 4:
            #        mask = np.logical_not( all_unmasked( res=res ) )
            mask = all_unmasked(res=res)
        elif case == 5:
            mask = extratropics_unmasked(res=res)
        elif case == 6:
            mask = ocean_unmasked(res=res)
        elif case == 7:
            mask = NH_unmasked(res=res)
        elif case == 8:
            mask = SH_unmasked(res=res)
        elif case == 10:
            mask = ice_unmasked(res=res)
        elif case == 11:
            mask = land_unmasked(res=res)
        elif case == 12:
            mask = mask_lat40_2_40(res=res)
        elif case == 13:
            mask = np.ma.mask_or(ocean_unmasked(res=res),
                                 tropics_unmasked(res=res, saizlopez=saizlopez))
        elif case == 14:
            mask = np.ma.mask_or(land_unmasked(res=res),
                                 tropics_unmasked(res=res, saizlopez=saizlopez))
        elif case == 15:  # 'All Sur.'
            mask = surface_unmasked(res=res)
        elif case == 16:  # 'Ocean Sur.'
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 ocean_unmasked(res=res))
        elif case == 17:  # 'Land Sur.':
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 land_unmasked(res=res))
        elif case == 18:  # 'Ice Sur.'
            mask = np.ma.mask_or(surface_unmasked(res=res),
                                 ice_unmasked(res=res))
        elif case == 19:
            mask = lat2lat_2D_unmasked(lowerlat=-50, higherlat=50,
                                       res=res)
        elif case == 20:
            mask = np.ma.mask_or(lat2lat_2D_unmasked(lowerlat=-50,
                                                     higherlat=50, res=res), ocean_unmasked(res=res)[..., 0])
        elif case == 21:
            mask = get_north_sea_unmasked(res=res)
        elif case == 22:
            mask = get_mediterranean_sea_unmasked(res=res)
        elif case == 23:
            mask = get_unmasked_black_sea(res=res)
        elif case == 24:
            mask = get_unmasked_irish_sea(res=res)
        elif case == 25:
            mask = get_EU_unmasked(res=res)
#        if case == 26:
#            mask = get_2D_BL_unmasked( res=res )
        elif case == 27:  # 'Land Tropics Sur.':
            tmp = np.ma.mask_or(surface_unmasked(res=res),
                                land_unmasked(res=res))
            mask = np.ma.mask_or(tmp, tropics_unmasked(res=res))
        elif case == 28:
            mask = np.ma.mask_or(lat2lat_2D_unmasked(lowerlat=50,
                                                     higherlat=80, res=res), land_unmasked(res=res)[..., 0])
        elif case == 29:  # Alps
            # Alps mask
            lowerlat = 43
            higherlat = 47
            lowerlon = 5
            higherlon = 15
            # Get a mask for lat and lon range, then combine
            mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                        higherlat=higherlat)
            mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                        higherlon=higherlon)
            mask = np.ma.mask_or(mask1, mask2)
        elif case == 30:  # Location ('loc' )
            # Alps
            mask = location_unmasked(lat=lat, lon=lon, res=res)
        elif case == 31:  # Rough(!) France map
            # mask
            mask = get_France_unmasked(res=res)

        else:
            print('WARNING - Mask not setup for case={}'.format(case))
            sys.exit()

    logging.debug('prior to setting dimensions: {}'.format(mask.shape))
    # Apply Saiz-Lopez Marine MFT/MUT? <= should this be before multiply op.?
    if M_all:
        if use_multiply_method:  # Kludge
            mask = mask*land_unmasked(res=res)
        else:
            # check this!!!
            mask = np.ma.mask_or(mask, land_unmasked(res=res))

    # Ensure returned arrays are 2D
    if mask2D:
        if len(mask.shape) == 2:
            pass
        elif len(mask.shape) == 3:
            mask = mask[..., 0]
        elif len(mask.shape) == 4:
            mask = mask[..., 0, 0]

    # Create 3D array by concatenating through altitude dimension
    if mask3D:
        if any([(mask.shape[-1] == i) for i in (38, 47)]):
            pass
        else:  # concatenate dimensions
            if len(mask.shape) == 3:
                mask = np.concatenate([mask]*47, axis=2)
            elif len(mask.shape) == 2:
                mask = np.concatenate([mask[..., None]]*47, axis=2)

    # Remove above the "chemical tropopause" from GEOS-Chem (v9-2)
    if trop_limit:
        if (len(mask.shape) == 2) or mask2D:
            pass
        else:
            mask = mask[..., :38]

    # Create 4D array by concatenating through time dimension
    # ( assuming year long array of 1 months )
    if mask4D:
        if any([(mask.shape[-1] == i) for i in [12]]):
            pass
        else:  # concatenate dimensions
            mask = np.concatenate([mask[..., None]]*12, axis=3)
    logging.debug('post to setting dimensions: {}'.format(mask.shape))
    logging.info("returning a 'mask' of type:{}".format(type(mask)))
    return mask


# --------
# X.XX - Custom 2D (Lon) Mask
# --------
def lon2lon_2D_unmasked(lowerlon, higherlon, res='2x2.5', debug=False):
    """
    Takes a lower and higher latitude value and then creates
    mask to given given limits.
    """

    # Get vars
    lon_c, lat_c, NIU = get_latlonalt4res(res=res, centre=True)
    if debug:
        print((lon_c, lowerlon, higherlon))

    # Mask between upper and lower values
    lons = [i for i in lon_c if ((i >= lowerlon) and (i < higherlon))]
    lons = [get_gc_lon(i, res=res) for i in lons]

    # Fill all lat and lon True or False
    m = np.zeros(get_dims4res(res))[:, :, 0]
    print((m.shape, np.sum(m)))
    for i in lons:
        m[i, :] = 1
    m = np.ma.masked_not_equal(m, 1)
    return m.mask

# --------
# X.XX - EU mask
# --------


def get_EU_unmasked(res='1x1'):
    """
    Mask 'EU' as defined by GEOS-Chem EU grid the grid of "'0.5x0.666" resolution is
    used by default, but any list of lat and lons could be provided and the extremities
    would be used as the mask edges
    """
    EU_resolutions = ['0.25x0.3125', '0.5x0.666']
    if res not in EU_resolutions:
        EU_res = EU_resolutions[0]  # default = '0.25x0.3125'
#        EU_res = EU_resolutions[1] # default = '0.5x0.666'
    else:
        EU_res = res
    # Get GEOS-Chem EU lat and lons
    lon, lat, NIU = get_latlonalt4res(res=EU_res)
    # mask lats
    m1 = lat2lat_2D_unmasked(lowerlat=lat.min(), higherlat=lat.max(), res=res)
    # Mask lons
    m2 = lon2lon_2D_unmasked(lowerlon=lon.min(), higherlon=lon.max(), res=res)
    # Combine maskes
    m = m1 + m2
    return m

# --------
# X.XX - Cruise track mask
# --------


def get_cruise_track_mask(max_lon=None, min_lon=None, max_lat=None,
                          min_lat=None, unmask_water=True, res='4x5', trop_limit=True):
    """
    Mask whole area of ship based research campaigns for bulk comparison
    """
    # only look at surface
    m = surface_unmasked(res=res, trop_limit=trop_limit)
    # apply ocean mask
    if unmask_water:
        m = m + ocean_unmasked(res=res)
    # Mask over given longitude range, if provided
    if not isinstance(max_lon, type(None)):
        m = m + lon2lon_2D_unmasked(lowerlon=min_lon, higherlon=max_lon,
                                    res=res)[:, :, None]
    # Mask over given latitude range, if provided
    if not isinstance(max_lat, type(None)):
        m = m + lat2lat_2D_unmasked(lowerlat=min_lat, higherlat=max_lat,
                                    res=res)[:, :, None]
    # Invert
    m = np.logical_not(m)
    return m


# --------
# X.XX - Get mask of mediterranean sea
# --------
def get_France_unmasked(res='4x5'):
    """
    A rough Mask of France for use with 2x2.5 / 4x5 model output
    """
    # France mask
    lowerlat = 42.5
    higherlat = 51
    lowerlon = -4.441
    higherlon = 7.7577
    # Get a mask for lat and lon range, then combine
    mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                higherlat=higherlat)
    mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                higherlon=higherlon)
    mask = np.ma.mask_or(mask1, mask2)
    # Only consider land grid boxes
    mask = np.ma.mask_or(mask, land_unmasked(res=res)[..., 0])
    return mask


# --------
# X.XX - Get mask of Mediterranean sea
# --------
def get_mediterranean_sea_unmasked(res='0.25x0.3125'):
    """
    A rough Mask of the Mediterranean Sea for use with ~0.5/~0.25 mdodel output.
    """
    # West South corner (south Jordan) =
    lowerlat = 34
    lowerlon = -6.5
    # East North corner (~Ukraine)
    higherlat = 47
    higherlon = 38
    # Get a mask for lat and lon range, then combine
    mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                higherlat=higherlat)
    mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                higherlon=higherlon)
    mask = np.ma.mask_or(mask1, mask2)
    # Add mask for water
    mask = np.ma.mask_or(mask, ocean_unmasked(res=res)[..., 0])
    # Also remove black sea ( by removing an inverted unmasked mask )
    mask3 = get_black_sea_unmasked(res=res, unmask_water=False)
    mask = np.ma.mask_or(mask, np.logical_not(mask3))
    # Also remove bay of biscay
    # Also remove black sea
    mask4 = get_bay_of_biscay_unmasked(res=res)
    mask = np.ma.mask_or(mask, np.logical_not(mask4))
    return mask


# --------
# X.XX - Get mask of black sea
# --------
def get_black_sea_unmasked(res='0.25x0.3125', unmask_water=True):
    """
    A rough Mask of the Mediterranean Sea for use with ~0.5/~0.25 mdodel output.
    """
    # West South corner (south Jordan) =
    lowerlat = 41
    lowerlon = 26.8
    # East North corner (~Ukraine)
    higherlat = 50
    higherlon = 43
    # Get a mask for lat and lon range, then combine
    mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                higherlat=higherlat)
    mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                higherlon=higherlon)
    mask = np.ma.mask_or(mask1, mask2)
    # Add mask for water
    if unmask_water:
        mask = np.ma.mask_or(mask, ocean_unmasked(res=res)[..., 0])
    return mask

# --------
# X.XX - Get mask of bay of biscay
# --------


def get_bay_of_biscay_unmasked(res='0.25x0.3125', unmask_water=True):
    """
    A rough Mask of the Mediterranean Sea for use with ~0.5/~0.25 mdodel output.
    """
    # West South corner (south Jordan) =
    lowerlat = 42.5
    lowerlon = -10
    # East North corner (~Ukraine)
    higherlat = 51
    higherlon = 0
    # Get a mask for lat and lon range, then combine
    mask1 = lat2lat_2D_unmasked(res=res, lowerlat=lowerlat,
                                higherlat=higherlat)
    mask2 = lon2lon_2D_unmasked(res=res, lowerlon=lowerlon,
                                higherlon=higherlon)
    mask = np.ma.mask_or(mask1, mask2)
    # Add mask for water
    if unmask_water:
        mask = np.ma.mask_or(mask, ocean_unmasked(res=res)[..., 0])
    # Also remove black sea
    return mask


# --------
# X.XX - Get mask of north sea
# --------
# def get_north_sea_unmasked( res='0.25x0.3125' ):
#     """
#     A rough Mask of the North Sea for use with ~0.5/~0.25 mdodel output.
#     (inc. English channel. )
#
#     """
#     mask latitudes of north sea
#     Drawing a box that includes all of North sea and English channel
#     Near Brest in France 48.3669927,-4.7560745
#     Lillehammer 61.1122408,10.4386779
#
#     mask lats
#     m1 = lat2lat_2D_unmasked( lowerlat=48.3669927, higherlat=61.1122408, res=res )
#
#     mask lons
#     m2 = lon2lon_2D_unmasked(lowerlon=-4.7560745, higherlon=10.4386779, res=res )
#
#      combine maskes
#     m = m1 + m2
#
#     remove all water.
#    m = np.ma.mask_or( m, ocean_unmasked( res=res)[...,0] ) # old approach
#     the below will work for muliple options.
#    m = m +
#    m = ocean_unmasked( res=res)[...,0]
#
#
#     remove irish sea
#
#     return m.mask
#
#
# --------
# X.XX - Get mask of black sea
# --------
# def get_unmasked_black_sea( res='0.25x0.3125'):
#     """
#     A rough Mask of the Black Sea for use with ~0.5/~0.25 mdodel output.
#     """
#
#     Drawing a box that include all of Black. Sea
#     East South corner (south Jordan) = 43.4623268,33.3809392
#     West North corner  = 44.17207,25.30604
#     pass
#
# --------
# X.XX - Get mask of irsh sea
# --------
# def get_unmasked_irish_sea( res='0.25x0.3125', unmasked_oceans=None):
#     """
#     A rough Mask of the Irish Sea for use with ~0.5/~0.25 mdodel output.
#     """
#
#     NE corner Glasgow - 55.8553803,-4.3725463
#     SW corner Cork - 51.8959842,-8.5332609
#     mask lats
#     m1 = lat2lat_2D_unmasked( lowerlat=51.8959842, higherlat=55.855380, res=res )
#
#     mask lons
#     m2 = lon2lon_2D_unmasked(lowerlon=8.5332609, higherlon=4.3725463, res=res )
#
#      combine maskes
#     m = m1 + m2
#
#     only consider oceans
#     m = np.ma.mask_or( ocean_unmasked( res=res)[...,0], m )
#
#
#     return mask

# --------
# X.XX - Get orthogonal distance regression (ODR) of two datasets
# --------
def get_linear_ODR(x=None, y=None, job=10, maxit=5000, beta0=(0, 1),
                   xvalues=None, return_model=True, debug=False, verbose=False):
    """
    Wrapper to run ODR for arrays of x and y

    NOTES
    -----
    adapted from example in manual
    (https://docs.scipy.org/doc/scipy/reference/odr.html)
    """
    import scipy.odr
    # Setup linear model to fit

    def f(B, x):
        '''Linear function y = m*x + b'''
        # B is a vector of the parameters.
        # x is an array of the current x values.
        # x is in the same format as the x passed to Data or RealData.
        #
        # Return an array in the same format as y passed to Data or RealData.
        return B[0]*x + B[1]
    # Create a model
    linear = scipy.odr.Model(f)
    # Create a Data or RealData instance.:
    mydata = scipy.odr.Data(x, y)
    # Instantiate ODR with your data, model and initial parameter estimate.:
#    myodr = scipy.odr.ODR(mydata, linear, beta0=[1., 2.])
    myodr = scipy.odr.ODR(mydata, linear, beta0,  maxit=maxit, job=job)
    # Run the fit.:
    myoutput = myodr.run()
    # Examine output.:
    if verbose:
        myoutput.pprint()
    if return_model:
        return myoutput
    else:
        if isinstance(xvalues, type(None)):
            xvalues = np.arange(min(x), max(x), (max(x)-min(x))/100.)
        return xvalues, f(myoutput.beta, xvalues)


# --------
# X.XX - Convert ug per m3 to 2 ppbv
# --------
def convert_ug_per_m3_2_ppbv(data=None,  spec='O3', rtn_units=False,
                             units='ug m$^{-3}$'):
    """
    Converts units of ugm^-3 to ppbv for a given species assuming SATP
    """
    # --- Get constants
    RMM_air = constants('RMM_air')  # g/mol
    # assume standard air density
    # At sea level and at 15 °C air has a density of approximately 1.225 kg/m3
    # (0.001225 g/cm3, 0.0023769 slug/ft3, 0.0765 lbm/ft3) according to
    # ISA (International Standard Atmosphere).
    AIRDEN = 0.001225  # g/cm3
    # moles per cm3
    #  (1/(g/mol)) = (mol/g) ; (mol/g) * (g/cm3) = mol/cm3
    MOLS = (1/RMM_air) * AIRDEN

    # --- Convert
    # convert ug/m3 to ppbv
    # get g per cm3, ( instead of ug/m3)
    data = data / 1E6 / 1E6
    # get mol/cm3 (mass/ RMM ) = ( ug/m3  /  g/mol )
    data = data/species_mass(spec)
    # convert to ppb
    data = data/MOLS * 1E9
    # update unit string
    units = 'ppbv'
    if rtn_units:
        return data, units
    else:
        return data


# --------
# X.XX - Convert mg per m3 to 2 ppbv
# --------
def convert_mg_per_m3_2_ppbv(data=None,  spec='O3', rtn_units=False,
                             units='mg m$^{-3}$'):
    """
    Converts units of ugm^-3 to ppbv for a given species assuming SATP
    """
    # --- Get constants
    RMM_air = constants('RMM_air')  # g/mol
    # assume standard air density
    # At sea level and at 15 °C air has a density of approximately 1.225 kg/m3
    # (0.001225 g/cm3, 0.0023769 slug/ft3, 0.0765 lbm/ft3) according to
    # ISA (International Standard Atmosphere).
    AIRDEN = 0.001225  # g/cm3
    # moles per cm3
    #  (1/(g/mol)) = (mol/g) ; (mol/g) * (g/cm3) = mol/cm3
    MOLS = (1/RMM_air) * AIRDEN

    # --- Convert
    # convert mg/m3 to ppbv
    # get g per cm3, ( instead of mg/m3)
    data = data / 1E6 / 1E3
    # get mol/cm3 (mass/ RMM ) =  (mg/m3  /  g/mol )
    data = data/species_mass(spec)
    # convert to ppb
    data = data/MOLS * 1E9
    # update unit string
    units = 'ppbv'
    if rtn_units:
        return data, units
    else:
        return data


# --------
# X.XX - Get 2D (lat, lon) mask of night time for a given datetime
# --------
def get_2D_nighttime_mask4date_pd(date=None, ncfile=None, res='4x5',
                                  mask_daytime=False, buffer_hours=0, debug=False):
    """
    Creates 2D (lon,lat) masked (1=Masked) for nighttime for a given list of
    dates

    Parameters
    -------
    date (datetime): date to use (UTC)
    mask_daytime (boolean): mask daytime instead of nightime
    ncfile (str): location to netCDF file - not implemented...
    res (str): resolution, if using resolutions listed in get_latlonalt4res
    buffer_hours (float/int): number of hours to buffer subrise/sunset with
     (This will act to increase the size of the mask - e.g. if masking
      nightime, then an extra hour of nightime would be added to sunrise, and
      removed from sunset. )

    Returns
    -------
    (np.array) with conditional values masked (1=masked)

    ncfile (NetCDF file): NetCDF file to extract lat and lon metadata from

    Notes
    -----
     - if ncfile provide programme will work for that grid.
    """
    # Astronomical math
    import ephem
    from ephem import AlwaysUpError, NeverUpError
    # And functions in other AC_tools modules
    from .funcs4time import add_days, add_hrs
    logging.info('get_2D_nighttime_mask4date_pd called for {}'.format(date))

    # Profile function...
    if debug:
        start_time = time.time()

    # --- Local variables?
    # reference data for ephem (number of days since noon on 1899 December 31)
    ref_date = datetime.datetime(1899, 12, 31, 12)

    # --- Get LON and LAT variables
    if isinstance(ncfile, type(None)):
        # extract from refence files
        lons, lats, alts = get_latlonalt4res(res=res)
    else:
        # TODO - allow any lat, lon grid to be used by taking input lats and
        # lons from ncfile file/arguments.
        print('Not implemented')
        sys.exit()
    if debug:
        print(("--- (start-1) %s seconds ---" % (time.time() - start_time)))

    # --- setup function to mask based on date, lat and lon
    def mask_nighttime(lon, lat, date=date, mask_daytime=mask_daytime,
                       ref_date=datetime.datetime(1899, 12, 31, 12),
                       buffer_hours=buffer_hours, debug=False):
        """
        sub-function to mask if nightime for a given date at a specific lat/lon
        """
        # --- get lat and lon values from columns
        if debug:
            print(("--- (s4-1) %s seconds ---" % (time.time() - start_time)))
        # --- get sunrise and sunset for location
        o = ephem.Observer()
        # set lat (decimal?), lon (decimal?), and date (UTC)
        o.lat = str(lat)
        o.long = str(lon)
        o.date = date
        # planetary body
        s = ephem.Sun()
        if debug:
            print(("--- (s4-2) %s seconds ---" % (time.time() - start_time)))

        # Compute sun vs observer
        s.compute()
        if debug:
            print(("--- (s4-3) %s seconds ---" % (time.time() - start_time)))

        # Work out if day or night based on sunrises and sunsets
        mask_value = 0
        try:

            # get sunrise time and date
            next_rising = o.next_rising(s)
            next_setting = o.next_setting(s)

            # convert to datetime.datetime
            next_rising = add_days(ref_date, next_rising)
            next_setting = add_days(ref_date, next_setting)

            # did the sun last rise or set? (inc. any adjustments)
            sun_last_rose = False
            if next_setting < next_rising:
                sun_last_rose = True

            # Add buffer to rising/setting if provided with buffer_hours
            if buffer_hours != 0:

                # Calculate last rise
                previous_rising = o.previous_rising(s)
                # convert to datetime.datetime
                previous_rising = add_days(ref_date, previous_rising)
                # Calculate last setting
                previous_setting = o.previous_setting(s)
                # convert to datetime.datetime
                previous_setting = add_days(ref_date, previous_setting)

                # Calculate absolute difference
                time_from_rise = (date-previous_rising).total_seconds()
                time_till_set = (date-next_setting).total_seconds()
                time_from_set = (date-previous_setting).total_seconds()
                time_till_rise = (date-next_rising).total_seconds()

                # If absolutely difference less than buffer
                if abs(time_from_rise)/60./60. <= buffer_hours:
                    mask_value = 1
                elif abs(time_till_set)/60./60. <= buffer_hours:
                    mask_value = 1
                elif abs(time_from_set)/60./60. < buffer_hours:
                    mask_value = 1
                elif abs(time_till_rise)/60./60. < buffer_hours:
                    mask_value = 1

            # --- Check if daytime or nighttime and mask if condition met.
            if sun_last_rose:
                if mask_daytime:
                    # ... and has not set yet, it must be daytime
                    if (date < next_setting):
                        mask_value = 1

            # if the sun last set... (mask nighttime is default)
            else:
                # if mask nighttime (aka not mask_daytime)
                if not mask_daytime:
                    # ... and has not risen yet, it must be nighttime
                    if (date < next_rising):
                        mask_value = 1

        # Add gotcha for locations where sun is always up.
        except AlwaysUpError:
            if mask_daytime:
                mask_value = 1

        # Add gotcha for locations where sun is always down.
        except NeverUpError:
            if not mask_daytime:
                mask_value = 1

        except:
            print('FAIL')
            sys.exit()

        # Mask value in array
        return mask_value

    # --- Setup an unstack(ed) pandas dataframe to contain masked values
    if debug:
        print(("--- (2) %s seconds ---" % (time.time() - start_time)))
    # Use list comprehension to setup list of indices for lat and lon
    # Better way of doing this? (e.g. pd.melt?)
    ind_lat_lons_list = [[lon_, lat_] for lat_ in lats for lon_ in lons]
    if debug:
        print(("--- (3) %s seconds ---" % (time.time() - start_time)))
    # Make this into a pd.DataFrame and label columns.
    df = pd.DataFrame(ind_lat_lons_list)
    df.columns = ['lons', 'lats']
    if debug:
        print(("--- (4) %s seconds ---" % (time.time() - start_time)))
    # Apply function to calculate mask value
#    df['mask'] = df.apply(mask_nighttime, axis=1)
    df['mask'] = df.apply(lambda x: mask_nighttime(
        x['lons'], x['lats']), axis=1)
    if debug:
        print(("--- (5) %s seconds ---" % (time.time() - start_time)))
    # Re-index by lat and lon
    df = pd.DataFrame(df['mask'].values, index=[df['lats'], df['lons']])
    if debug:
        print(("--- (6) %s seconds ---" % (time.time() - start_time)))
    # Unstack and return just as array
    df = df.unstack()
    marr = df.values
    if debug:
        print(("--- (end-7) %s seconds ---" % (time.time() - start_time)))

    return marr


# --------
# X.XX - Get 2D array of solar time.
# --------
def get_2D_solartime_array4_date(date=None, ncfile=None, res='4x5',
                                 lons=None, lats=None, varname='SolarTime',  debug=False):
    """
    Creates 2D (lon,lat) masked (1=Masked) for nighttime for a given list of
    dates

    Parameters
    -------
    date (datetime): date to use (UTC)
    mask_daytime (boolean): mask daytime instead of nightime
    ncfile (str): location to netCDF file - not implemented...
    res (str): resolution, if using resolutions listed in get_latlonalt4res
    lons (array): array of longditudes (optional)
    lats (array): array of lattiudes (optional)

    Returns
    -------
    (np.array)

    ncfile (NetCDF file): NetCDF file to extract lat and lon metadata from

    Notes
    -----
     - if ncfile provide programme will work for that grid.
    """
    # Astronomical math
    import ephem
    from ephem import AlwaysUpError, NeverUpError
    # And functions in other AC_tools modules
    from .funcs4time import add_days, add_hrs
    logging.info('get_2D_solartime_array4_dates called for {}'.format(date))

    # Profile function...
    if debug:
        start_time = time.time()

    # --- Local variables?
    # reference data for ephem (number of days since noon on 1899 December 31)
    ref_date = datetime.datetime(1899, 12, 31, 12)

    # --- Get LON and LAT variables (if lons/lats not provdided)
    if any([not isinstance(i, type(None)) for i in (lats, lons)]):
        pass
    else:
        if isinstance(ncfile, type(None)):
            # extract from refence files
            lons, lats, NIU = get_latlonalt4res(res=res)
        else:
            # TODO - allow any lat, lon grid to be used by taking input lats and
            # lons from ncfile file/arguments.
            print('Not implemented')
            sys.exit()

    # --- setup function to get solartime based on date, lat and lon
    def _solartime(lon, lat, date=date, rtn_as_epoch=True):
        """
        Get solartime for location (lat, lon) and date

        Returns
        ---
        (float) Epoch time
        """
        # --- get sunrise and sunset for location
        o = ephem.Observer()
        # set lat (decimal?), lon (decimal?), and date (UTC)
        o.lat = str(lat)
        o.long = str(lon)
        o.date = date
        # planetary body
        s = ephem.Sun()
        # Compute sun vs observer
        s.compute(o)
        # below code was adapted from stackoverflow (Credit: J.F. Sebastian)
        # http://stackoverflow.com/questions/13314626/local-solar-time-function-from-utc-and-longitude
        # sidereal time == ra (right ascension) is the highest point (noon)
        hour_angle = o.sidereal_time() - s.ra
        s_time = ephem.hours(
            hour_angle + ephem.hours('12:00')).norm  # norm for 24h
        # ephem.hours is a float number that represents an angle in radians
        # and converts to/from a string as "hh:mm:ss.ff".
        s_time = "%s" % s_time
        if len(s_time) != 11:
            s_time = '0'+s_time
        # return as datetime
        s_time = datetime.datetime.strptime(s_time[:8], '%H:%M:%S')
        if rtn_as_epoch:
            s_time = unix_time(s_time)
        return s_time

    # --- Setup an unstack(ed) pandas dataframe to contain masked values
    # Use list comprehension to setup list of indices for lat and lon
    # Better way of doing this? (e.g. pd.melt?)
    ind_lat_lons_list = [[lon_, lat_] for lat_ in lats for lon_ in lons]
    # Make this into a pd.DataFrame and label columns.
    df = pd.DataFrame(ind_lat_lons_list)
    df.columns = ['lons', 'lats']
    # Apply function to calculate mask value
    df[varname] = df.apply(lambda x: _solartime(x['lons'], x['lats']), axis=1)
    # Re-index by lat and lon
    df = pd.DataFrame(df[varname].values, index=[df['lats'], df['lons']])
    # Unstack and return just as array
    return df.unstack().values


# --------
# X.XX - Save 2D arrays (lat, lon) to 3D netCDF (3rd dim=time)
# --------
def save_2D_arrays_to_3DNetCDF(ars=None, dates=None, res='4x5', lons=None,
                               lats=None, varname='MASK', Description=None, Contact=None,
                               filename='misc_output', var_type='f8', debug=False):
    """
    makes a NetCDF from a list of dates and list of (lon, lat) arrays

    Parameters
    -------
    ars (list of np.array): list of (lon, lat) arrays
    dates (list of datetime.datetime)
    res (str): reoslution (opitional )
    lons (list): lons for use for NetCDF cooridinate variables
    lats (list): lats for use for NetCDF cooridinate variables
    filename (str): name for output netCDF file
    varname (str): name for variable in NetCDF
    var_type (str): variable type (e.g. 'f8' (64-bit floating point))

    Returns
    -------
    None

    Notes
    -----
     - needs updating to take non-standard reoslutions
     (e.g. those in funcs_vars)

    """
    logging.info('save_2D_arrays_to_3DNetCDF called')
    from .funcs4time import unix_time, dt64_2_dt
    print((locals()))

    # ---  Settings
    ncfilename = '{}_{}.nc'.format(filename, res)
    # Get lons and lats...
    if any([isinstance(i, type(None)) for i in (lats, lons)]):
        lons, lats, NIU = get_latlonalt4res(res=res)
        logging.debug('Using offline lons/lats, as either lats/lons==None')
#    else:
#        print 'WARNING: non-standard lats/lons not implemented!!! - TODO. '
#        sys.exit()

    # --- Setup new file
    # write file
    logging.debug('setting up new NetCDF file: {}'.format(ncfilename))
    ncfile = Dataset(ncfilename, 'w', format='NETCDF4')
    ncfile.createDimension('lat', len(lats))
    ncfile.createDimension('lon', len(lons))
    ncfile.createDimension('time', None)

    # Define the coordinate variables. They will hold the coordinate
    # information, that is, the latitudes and longitudes.
    time = ncfile.createVariable('time', 'f4', ('time',))
    lat = ncfile.createVariable('lat', 'f4', ('lat',))
    lon = ncfile.createVariable('lon', 'f4', ('lon',))

    # --- Add meta data
    # Assign units attributes to coordinate var data. This attaches a
    # text attribute to each of the coordinate variables, containing the
    # units.
    lat.units = 'degrees_north'
    lat.long_name = 'Latitude'
    lat.standard_name = 'Latitude'
    lat.axis = "Y"

    lon.units = 'degrees_east'
    lon.long_name = 'Longitude'
    lon.standard_name = 'Longitude'
    lon.axis = "X"

    time.units = 'seconds since 1970-01-01 00:00:00'
    time.calendar = "standard"
    time.standard_name = 'Time'
    time.axis = "T"

    # --- Set global variables
    if not isinstance(Description, type(None)):
        ncfile.Description = Description
    if not isinstance(Contact, type(None)):
        ncfile.Contact = Contact
    ncfile.Grid = 'lat: {}-{}, lon: {}-{}'.format(lats[0], lats[-1],
                                                  lons[0], lons[-1])

    # Write values to coordinate variables (lat, lon)
    lon[:] = lons
    lat[:] = lats

    # ---  Setup time dimension/variables (as epoch)
    # Convert to Epoch time
    def format(x): return unix_time(x)
    df = pd.DataFrame({'Datetime': dates})
    df['Epoch'] = df['Datetime'].map(format).astype('i8')
    # store a copy of dates at datetime.datetime, then remove from DataFrame
    dt_dates = dt64_2_dt(df['Datetime'].values.copy())
    del df['Datetime']
    dates = df['Epoch'].values
    # Assign to time variable
    time[:] = dates

    # --- Create new NetCDF variable (as f8) with common dimensions
    # (e.g. 'f8' = 64-bit floating point, 'i8'=(64-bit singed integer) )
    ncfile.createVariable(varname, var_type, ('time', 'lat', 'lon'), )
    # Close NetCDF
    ncfile.close()

    # --- Now open and add data in append mode
    for n, date in enumerate(dt_dates):

        if debug:
            fmtstr = "%Y%d%m %T"
            logging.debug('saving array date:{}'.format(date.strftime(fmtstr)))
            pcent = (float(n)+1)/float(len(dt_dates))*100
            logging.debug('array #: {} (% complete: {:.1f})'.format(n, pcent))

        # Open NetCDF in append mode
        ncfile = Dataset(ncfilename, 'a', format='NETCDF4')

        # Add data to array
        try:
            ncfile.variables[varname][n] = ars[n]
        except ValueError:
            err_msg = '>total size of new array must be unchanged<'
            err_msg += '(new arr shape {})'.format(str(ars[n].shape))
            print(err_msg)
            logging.info(err_msg)
            sys.exit()

    # Close NetCDF
    ncfile.close()
    if debug:
        logging.debug('saved NetCDF file:{}'.format(ncfilename))


# --------
# X.XX - Interpolate values from subset of 2D array
# --------
def interpolate_sparse_grid2value(Y_CORDS=None, X_CORDS=None,
                                  X=None, Y=None, XYarray=None, buffer_CORDS=5,
                                  verbose=True, debug=False):
    """
    Get an interpolated value for a location (X,Y) from surrounding array values

    Parameters
    -------
    X_CORDS (np.array): coordinate values for X axis of 2D array (e.g. lon)
    Y_CORDS (np.array): coordinate values for Y axis of 2D array (e.g. lat)
    X (float): value of interest (in same terms at X_CORDS, e.g. lon)
    Y (float): value of interest (in same terms at Y_CORDS, e.g. lat)
    XY (array): array of values with shape (X_CORDS, Y_CORDS)
    buffer_CORDS (int): number of units of X_CORDS and Y_CORDS to interpolate
        around (great this value, greater the cost.)
    verbose (boolean): print out extra infomation
    debug (boolean): print out debugging infomation

    Returns
    -------
    (float)
    """
    # X_CORDS=file_lonc; Y_CORDS=file_latc; X=lon; Y=lat; XYarray=file_data; buffer_CORDS=10
    import scipy.interpolate as interpolate
    # ---  Select a sub grid.
    # WARNING THIS APPRAOCH WILL NOT WORK NEAR ARRAY EDGES!!!
    # (TODO: add functionality for this.)
    # Get indices buffer sub-selection of grid.
#    assert X_CORDS[-1] > X_CORDS[0], 'This -180=>180 and
    if X_CORDS[0] < X_CORDS[-1]:
        s_low_X_ind = find_nearest_value(X_CORDS, X-buffer_CORDS)
        s_high_X_ind = find_nearest_value(X_CORDS, X+buffer_CORDS)
    else:  # Y_CORDS[0] > Y_CORDS[-1]
        s_low_X_ind = find_nearest_value(X_CORDS, X+buffer_CORDS)
        s_high_X_ind = find_nearest_value(X_CORDS, X-buffer_CORDS)
        XYarray[::-1, :]
    # WARNING: This approach assumes grid -90=>90 (but flips slice if not)
    if Y_CORDS[0] < Y_CORDS[-1]:
        s_low_Y_ind = find_nearest_value(Y_CORDS, Y-buffer_CORDS)
        s_high_Y_ind = find_nearest_value(Y_CORDS, Y+buffer_CORDS)
    else:  # Y_CORDS[0] > Y_CORDS[-1]
        s_low_Y_ind = find_nearest_value(Y_CORDS, Y+buffer_CORDS)
        s_high_Y_ind = find_nearest_value(Y_CORDS, Y-buffer_CORDS)
        XYarray[:, ::-1]
    # Print out values in use
    if verbose:
        prt_str = 'Y={}, subrange=({}(ind={}),{}(ind={}))'
        print(prt_str.format(Y, Y_CORDS[s_low_Y_ind], s_low_Y_ind,
                             Y_CORDS[s_high_Y_ind], s_high_Y_ind))
        prt_str = 'X={}, subrange=({}(ind={}),{}(ind={}))'
        print(prt_str.format(X, X_CORDS[s_low_X_ind], s_low_X_ind,
                             X_CORDS[s_high_X_ind], s_high_X_ind))
    # Select sub array and get coordinate axis
    subXY = XYarray[s_low_X_ind:s_high_X_ind, s_low_Y_ind:s_high_Y_ind]
    subX = X_CORDS[s_low_X_ind:s_high_X_ind]
    subY = Y_CORDS[s_low_Y_ind:s_high_Y_ind]
    # Debug (?) by showing 2D grid prior to interpolation
    if debug:
        print(([i.shape for i in (subX, subY, subXY)], XYarray.shape))
        plt.pcolor(subX, subY, subXY.T)
        plt.colorbar()
        plt.show()

    # ---  interpolate over subgrid.
    M = subXY
    rr, cc = np.meshgrid(subX, subY)
    # fill masked values with nans
    M = np.ma.filled(M, fill_value=np.nan)
    # only consider non nan values as values to interpolate with
    vals = ~np.isnan(M)
    if debug:
        print(vals)
    # interpolate
    f = interpolate.Rbf(rr[vals], cc[vals], M[vals], function='linear')
    # extract interpolation...
    interpolated = f(rr, cc)

    # ---  Select value of interest
    # Debug (?) by showing 2D grid post to interpolation
    if debug:
        print((interpolated.shape))
        plt.pcolor(subX, subY, interpolated.T)
        plt.colorbar()
        plt.show()

    # indix value from subgrid?
    Y_ind = find_nearest_value(subY, Y)
    X_ind = find_nearest_value(subX, X)

    return interpolated[X_ind, Y_ind]


# --------
# X.XX - Split a NetCDF by month
# --------
def split_NetCDF_by_month(folder=None, filename=None, ext_str='',
                          file_prefix='ts_ctm'):
    """
    Split a NetCDF file by month into new NetCDF files using xarray

    Parameters
    -------
    folder (str): the directory to search for files in
    filename (str): the NetCDF filename (e.g. ctm.nc)
    file_prefix (str): prefix to attach to new saved file
    ext_str (str): extra string for new filenames
    """
    import xarray as xr
    # --- Open data
    ds = xr.open_dataset(folder+filename)
    months = list(sorted(set(ds['time.month'].values)))

    # --- Loop months
    for month_ in months:

        # Custom mask
        def is_month(month):
            return (month == month_)
        # Make sure time is the dimension not module
        time = ds.time
        # Now select for month
        ds_tmp = ds.sel(time=is_month(ds['time.month']))

        # --- Save as NetCDF
        # Name of new NetCDF?
        year_ = list(set(ds_tmp['time.year'].values))[0]
        file2save = '{}_{}_{}_{:0>2}.nc'.format(file_prefix, ext_str, year_,
                                                str(month_))
        logging.info('saving month NetCDF as: {}'.format(file2save))
        # Save the file...
        ds_tmp.to_netcdf(folder+file2save)
        # Delete temporary dataset
        del ds_tmp


# --------
# X.XX - stack a 2D table (DataFrame) of lat/lon coords
# --------
def get_2D_df_of_lon_lats_and_time(res='4x5', df_lar_var='lat',
                                   df_lon_var='lon', df_time_var='month', add_all_months=False,
                                   lons=None, lats=None, verbose=True, month=9):
    """ stack a 2D table (DataFrame) of lat/lon coords """
    # Get lon and lat resolution (Add other ways to get lat and lon here...
    lons_not_provided = isinstance(lons, type(None))
    lats_not_provided = isinstance(lats, type(None))
    if (lons_not_provided) or (lats_not_provided):
        try:
            assert(type(res) == str), 'Resolution must be a string!'
            lons, lats, NIU = get_latlonalt4res(res=res)
        except:
            print('please provide lons/lats or a res in get_latlonalt4res')
    else:
        pass
    # Make Table like array
    b = np.zeros((len(lons), len(lats)))
    df = pd.DataFrame(b)
    # Use lats and lons as labels for index and columns
    df.index = lons
    df.columns = lats
    # Stack, then reset index to obtain table structure
    df = df.stack()
    df = df.reset_index()
    # Set column headers
    df.columns = df_lon_var, df_lar_var, df_time_var
    # Add time dims...
    if add_all_months:
        dfs = []
        for month in np.arange(1, 13):
            df_tmp = df.copy()
            df_tmp[df_time_var] = month
            dfs.append(df_tmp)
        df = pd.concat(dfs)
    else:  # Just select a single month (September is default )
        print('WARNING: Only September considered')
        df[df_time_var] = month
    return df


# --------
# X.XX - stack a 2D table (DataFrame) of lat/lon coords
# --------
def get_vars_from_line_printed_in_txt_file(filename=None, folder=None,
                                           prefix=' HOUR:', var_names=None, type4var_names=None):
    """
    Get variables from a line printed to non-binary file with a given prefix

    Parameters
    ----------
    filename (Str): name of non binary file (e.g. geos.log)
    folder (str): name of directory where file ("filename") is located
    prefix (str): the string that lines containing data begin with
    var_names (list): optional. names for variables in line (must be # in line)
    type4var_names (dict): dictionary of type to convert variables to

    Returns
    -------
    (pd.DataFrame)

    Notes
    ----------
     - Can be used for debugging or tabulating output. e.g. lines of
     grid index locations and values printed to GEOS-chem's geos.log file
    """
    import pandas as pd
    import sys
    # --- Local vars
    lines_with_var = []
    # --- Open file and loop vars
    with open(folder+filename) as f:

        for line in f:
            #			if prefix in line:
            #            print line
            if line.startswith(prefix):
                var_ = line.split()[1:]
                if not isinstance(var_names, type(None)):
                    try:
                        assert len(var_) == len(var_names)
                    except AssertionError:
                        print((len(var_), var_))
                        sys.exit()
                # Save data until later
                lines_with_var.append(var_)
    # If lines found, return as a DataFrame
    if len(lines_with_var) > 0:
        # Make DataFrame
        df = pd.DataFrame(lines_with_var)
        # Apply provided names to columns?
        if not isinstance(var_names, type(None)):
            df.columns = var_names
        # Covert data type of variable?
        if not isinstance(type4var_names, type(None)):
            for key_ in list(type4var_names.keys()):
                try:
                    df[key_] = df[key_].astype(type4var_names[key_])
                except KeyError:
                    err_str = 'key_ ({}) not in df'.format(key_)
                    logging.info(err_str)

        return df
    else:
        err_str = 'No lines with prefix ({})'.format(prefix)
        logging.info(err_str)

# --------
# X.XX - remove the spaces and extra vars from strings
# --------


def rm_spaces_and_chars_from_str(input_str, remove_slashes=True,
                                 replace_brackets=True, replace_quotes=True, replace_dots=True,
                                 remove_plus=True, swap_pcent=True, replace_braces=True):
    """ remove the spaces and extra vars from strings"""
    input_str = input_str.replace(' ', '_')
    if replace_brackets:
        input_str = input_str.replace('(', '_')
        input_str = input_str.replace(')', '_')
    if replace_braces:
        input_str = input_str.replace('{', '_')
        input_str = input_str.replace('}', '_')
    if replace_quotes:
        input_str = input_str.replace("'", '_')
    if replace_dots:
        input_str = input_str.replace(".", '_')
    if remove_slashes:
        input_str = input_str.replace("\\", '_')
        input_str = input_str.replace("/", '_')
    if remove_plus:
        input_str = input_str.replace("+", '_plus_')
    if swap_pcent:
        input_str = input_str.replace("%", 'pcent')
    return input_str
