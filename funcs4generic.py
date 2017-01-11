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
from funcs4core import *
from funcs_vars import *


# -------------------------- Section 7 -------------------------------
# -------------- Generic Processing
#

# -------------
# 1.01 -  Split arrays into chunks (e.g. days)
# -------------
def chunks(l, n):
    """ 
    Split list in chunks - useful for controlling memory usage 
    """
    if n < 1:
        n = 1
    return [l[i:i + n] for i in range(0, len(l), n)]
    
# -------------
# 1.02 - Get len of file
# -------------
def file_len(fname):
    """ 
    Get length of file 
    """
    return sum(1 for line in open(fname))

# -------------
# 1.03 - My round - credit: Alok Singhal
# -------------
def myround(x, base=5, integer=True, round_up=False):
    """ 
    Round up values - mappable function for pandas processing 
    NOTES:
     - credit: Alok Singhal
    """

    round_up=True # Kludge - always set this to True

    # round to nearest base number
    rounded = base * round(float(x)/base) 
        
    if round_up:
        # ensure rounded number is to nearest next base
        if rounded < x:
            rounded += base

    if integer:
        return int( rounded )
    else:        
        return rounded

# --------------
# 1.07 - Count number of files/the filenames that contain certain phrase/number conbinations (e.g. dates) - tms 
# -------------
def counter_directory_contains_files(model_path,must_contain):
    """ 
    Count number of files in directory 
    """
    model_ouput_file_counter = len(glob.glob1(model_path,must_contain))
    return model_ouput_file_counter
            

# ----
# 1.09 - Rename vars/strs in files
# ----
def replace_strs_in_files( wd, input_str, output_str, debug=False ):
    """ 
    replace text in files
    """
    print wd, input_str, output_str
    for f in os.listdir(wd) :
        if not f.startswith('.'):
            print f
            os.rename(wd+ f, wd+ f.replace(input_str, output_str)) 
            print f.replace(input_str, output_str)

# ----
# 1.10- Get X and Y coordinates for a given grid - Credit: Eric Sofen
# ----
def get_xy(Lon,Lat, lon_edges, lat_edges, debug=False):
    """
    Takes lon,lat values for a point (e.g. station) and arrays of longitude
    and latitudes indicating edges of gridboxes.

    Returns the index of the gridbox corresponding to this point.  
    NOTES:
     - Credit: Eric Sofen
     - Could be easily extended to return data values corresponding to points.
    """
    hasobs,lonrange,latrange=np.histogram2d([Lon],[Lat],[lon_edges, lat_edges])
    gridindx,gridindy=np.where(hasobs>=1)
    if not gridindx:
        if debug:
            print 'Lat, lon outside of x,y range.  Assigning -1 for', Lon, Lat
        return -1, -1
    else:
        #print Lon, Lat, gridindx, gridindy
        return gridindx[0],gridindy[0]
        
# --------   
# 1.12 - Save as pdf.
# --------
def plot2pdf(title='new_plot', fig=None, rasterized=True, dpi=160,\
        justHH=False, no_dstr=True, save2png=True, \
        save2eps=False, transparent=True, debug=False ):
    """ 
    Save figures (e.g. matplotlib) to pdf file 
    """

    # set save directory ( using default directory dictionary )
    from funcs_vars import get_dir 
    wd = get_dir('ppwd')

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
    pdf = PdfPages(npdf +'.pdf' )

    # Rasterise to save space?
    if rasterized:
        plt.gcf().set_rasterized(True)    

    # save and close
    type = 'PDF' 
    pdf.savefig( dpi=dpi, transparent=transparent )
    pdf.close()

    # Also export to png and eps?
    if save2png:
        type += '/PDF' 
        plt.savefig(npdf+'.png', format='png', dpi=dpi, \
            transparent=transparent )
    if save2eps:
        type += '/EPS' 
        plt.savefig(npdf+'.eps', format='eps', dpi=dpi, \
            transparent=transparent)
    print type+' saved & Closed as/at: ', npdf

# --------   
# 1.13- Save as mulitple page pdf.
# --------
def plot2pdfmulti(pdf=None, title='new_plot', rasterized=True, \
        dpi=160, open=False, close=False, justHH=False, no_dstr=False ):
    """ 
    Save figures (e.g. matplotlib) to pdf file with multiple pages 
    """

    # set save directory ( using default directory dictionary )
    from funcs4core import get_dir 
    wd = get_dir('ppwd')

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
        print 'pdf opened @: {}'.format( npdf )
        return pdf

    # Rasterise to save space?
    if rasterized:
        plt.gcf().set_rasterized(True)    

    # save and close or keep open to allow additions of plots
    if close:
        pdf.close()
        print 'PDF saved & Closed as/at: ', npdf
    else:
        pdf.savefig(dpi=dpi)
        print 'pdf is still open @: {}'.format( npdf )

        
# --------------
# 1.14 -  Return mesh and indics for which obs are within 
# -------------
def obs2grid(  glon=None, glat=None, galt=None, nest='high res global', \
        sites=None, debug=False ):
    """ 
    values that have a given lat, lon and alt
    """
    if isinstance(glon, type(None)):
        glon, glat, galt = get_latlonalt4res( nest=nest, centre=False, \
            debug=debug )

    # Assume use of known CAST sites... unless others given.
    if isinstance(sites, type(None)):
        loc_dict = get_loc( rtn_dict=True )
        sites = loc_dict.keys()            

    # pull out site location indicies
    indices_list=[] 
    for site in sites :
        lon, lat, alt = loc_dict[site]
        vars = get_xy(  lon,  lat , glon, glat ) 
        indices_list += [ vars ]
    return indices_list

# --------   
# 1.15 - Sort (GAW) sites by Latitude and return
# --------
def sort_sites_by_lat(sites): 
    """ 
    Order given list of GAW sties by latitudes 
    """ 

    # Get info
    vars = [ gaw_2_loc( s ) for s in sites ] # lat, lon, alt, TZ
    
    # Sort by lat, index orginal sites list and return
    lats = [ i[0] for i in vars ]
    slats = sorted( lats )[::-1]
    return [ sites[i] for i in [ lats.index(ii) for ii in slats ] ] 

# --------   
# 1.16 - Find nearest
# --------
def find_nearest(array,value):
    """
    Find nearest number in array to given value. 
    
#    NOTES:
#    - Adapted from (credit:) HappyLeapSecond's Stackoverflow answer. 
#    ( http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array )
#    """
    idx = (np.abs(array-value)).argmin()
    return idx

# --------   
# 1.17 - Get suffix for number
# --------
def get_suffix(n):
    """ 
    Add the appropriate suffix (th/st/rd) to any number given 
    """
    ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n/10%10!=1)*(n%10<4)*n%10::4])
    return ordinal(n)

# --------   
# 1.18 - Get shortest distance on a sphere ( for finding closest point )
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
    dlat = np.radians(haystack[:,0]) - radians(needle[0])
    dlon = np.radians(haystack[:,1]) - radians(needle[1])

    # get the vector
    a = np.square(np.sin(dlat/2.0)) + cos(radians(needle[0])) * np.cos(np.radians(haystack[:,0])) * np.square(np.sin(dlon/2.0))

    # convert to distance
    great_circle_distance = 2 * np.arcsin(np.minimum(np.sqrt(a), \
        np.repeat(1, len(a))))
    d = earth_radius_miles * great_circle_distance

    if r_distance:
        return np.min(d)
    # return the index
    else:
        return list(d).index( np.min(d) )

# --------   
# 1.19 - Get logarithmically spaced integers
# --------
def gen_log_space(limit, n):
    """
    Get logarithmically spaced integers
    
    NOTES:
    - credit: Avaris
    ( http://stackoverflow.com/questions/12418234/logarithmically-spaced-integers)
    """
    result = [1]
    if n>1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result)<n:
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
    return np.array(map(lambda x: round(x)-1, result), dtype=np.uint64)

# --------   
# 1.20 - Get indices in array where change in value of x occurs
# --------
def get_arr_edge_indices( arr, res='4x5', extra_points_point_on_edge=None, \
            verbose=True, debug=False):
    """ 
    Find indices in a lon, lat (2D) grid, where value does not equal a given 
    value ( e.g. the edge )    
    """
    
    if verbose:
        print 'get_arr_edge_indices for arr of shape: ', arr.shape

    # initialise variables
    lon_c, lat_c, NIU = get_latlonalt4res( res=res, centre=True )
    lon_e, lat_e, NIU = get_latlonalt4res( res=res, centre=False )
    lon_diff = lon_e[-5]-lon_e[-6]
    lat_diff = lat_e[-5]-lat_e[-6]
    nn, n, = 0,0
    last_lat_box = arr[nn, n]
    coords = []
    last_lon_box = arr[nn, n]
    need_lon_outer_edge, need_lat_outer_edge =  False, False
    if debug:
        print lon_e, lat_e

    # ---- Loop X dimension ( lon )
    for nn, lon_ in enumerate( lon_c ):
    
        # Loop Y dimension ( lat ) and store edges
        for n, lat_ in enumerate( lat_c ):

            if debug:
                print arr[nn, n], last_lat_box, last_lon_box,  \
                    arr[nn, n]==last_lat_box, arr[nn, n]==last_lon_box

            if arr[nn, n] != last_lat_box:

                # If 1st lat, selct bottom of box        
                point_lon = lon_e[nn]+lon_diff/2
                if need_lat_outer_edge:
                    point_lat = lat_e[n+1] 
                else:
                    point_lat = lat_e[n] 
                    need_lat_outer_edge = True
                need_lat_outer_edge=False

                # Add mid point to cordinates list
                if isinstance( extra_points_point_on_edge, type(None) ):
                    mid_point = [ point_lon, point_lat ] 
                    coords += [ mid_point ]

                # Add given number of points along edge
                else: 
                    coords += [ [lon_e[nn]+(lon_diff*i), point_lat ] for i in \
                        np.linspace(0, 1, extra_points_point_on_edge, \
                        endpoint=True) ]

            # temporally save the previous box's value
            last_lat_box = arr[nn, n]

    # ---- Loop Y dimension ( lat )
    for n, lat_ in enumerate( lat_c ):

        if debug:
            print arr[nn, n], last_lat_box, last_lon_box,  \
                arr[nn, n]==last_lat_box, arr[nn, n]==last_lon_box
        # Loop X dimension ( lon ) and store edges
        for nn, lon_ in enumerate( lon_c ):
        
            # If change in value at to list
            if arr[nn, n] != last_lon_box:
                point_lat = lat_e[n]+lat_diff/2

                # Make sure we select the edge lon
                if need_lon_outer_edge:
                    point_lon = lon_e[nn+1]
                else:
                    point_lon = lon_e[nn]
                    need_lon_outer_edge = True
                need_lon_outer_edge=False

                # Add mid point to coordinates list
                if isinstance( extra_points_point_on_edge, type(None) ):
                    mid_point = [ point_lon, point_lat ] 
                    coords += [ mid_point ]

                # Add given number of points along edge
                else:
                    coords += [ [point_lon, lat_e[n]+(lat_diff*i) ] for i in \
                        np.linspace(0, 1, extra_points_point_on_edge, \
                        endpoint=True) ]
                
            # temporally save the previous box's value
            last_lon_box = arr[nn, n]

    return coords

# --------
# 1.21 - Get data binned by uniques days in data
# --------
def split_data_by_days( data=None, dates=None, day_list=None, \
            verbose=False, debug=False ):
    """ 
    Takes a list of datetimes and data and returns a list of data and
    the bins ( days )  
    """
    if verbose:
        print 'split_data_by_days called'

    # Create DataFrame of Data and dates
    df = DataFrame( data, index=dates, columns=['data'] )
    # Add list of dates ( just year, month, day ) <= this is mappable, update?
    df['days'] = [datetime.datetime(*i.timetuple()[:3]) for i in dates ]
    if debug:
        print df

    # Get list of unique days
    if isinstance( day_list, type(None) ):
        day_list = sorted( set( df['days'].values ) )
    # Loop unique days and select data on these days 
    data4days = []
    for day in day_list:
        print day, df[df['days']==day] 
        data4days += [ df['data'][ df['days']==day ] ]
    # Just return the values ( i.e. not pandas array )
    data4days = [i.values.astype(float) for i in data4days ]
    print [ type(i) for i in data4days ]
#    print data4days[0]
#    sys.exit()

    if debug:    
        print 'returning data for {} days, with lengths: '.format( \
            len( day_list) ), [len(i) for i in data4days ]

    # Return as list of days (datetimes) + list of data for each day
    return data4days, day_list

# -------------------------- Section 2 -------------------------------
# -------------- Maskes for data analysis
#

# --------
# 2.01 - Ocean mask
# --------
def ocean_unmasked(res='4x5', debug=False):
    """ 
    Get ocean mask from GEOS-Chem LWI ( at given resolution) 
    
    NOTES: 
     - this currently returns mask with all except ocean masked 
     - NEEDS UPDATE 
    """

    from funcs4GEOSC import get_land_map 
    if debug:
        print 'ocean_mask called for: ', res

    # Create a mask from land/water/ice indices
    m = np.ma.masked_not_equal( get_land_map(res=res),0 )
    if debug:
        print mask, mask.shape, 
    return m.mask
 
# --------
# 2.02 - Land mask
# --------
def land_unmasked(res='4x5', debug=False):
    """ 
    Get land mask from GEOS-Chem LWI ( at given resolution) 
    """
    from funcs4GEOSC import get_land_map # Kludge, use GEOS-Chem LWI

    # Create a np.ma mask 
    if debug:
        print 'land_mask called for: ', res
    m = np.ma.masked_not_equal( get_land_map(res=res),1 )
    if debug:
        print mask, mask.shape, 
    return m.mask

# --------
# 2.03 - Ice mask
# --------
def ice_unmasked(res='4x5', debug=False):
    """ 
    Get mask for all areas but ice surface. Ice mask is from GEOS-Chem
    LWI ( at given resolution) 

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a np.ma mask 
    m= np.logical_not( (land_unmasked(res)*ocean_unmasked(res)) )
    if debug:
        print mask, mask.shape, 
    return m
 
# --------
# 2.04 - Surface mask
# --------
def surface_unmasked( res='4x5', trop_limit=False, mask2D=False, \
            debug=False ):
    """ 
    Get mask for all areas but surface

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """
    # Create a np.ma mask 
    m=np.ma.array( np.zeros( get_dims4res(res) ), mask=False)
    m[...,0] = 1
    m = np.ma.masked_not_equal(m, 1)

    if trop_limit:
        m = m[...,:38]
    if debug:
        print mask, mask.shape

    # Return 2D or 3D?
    if mask2D:
        return m[...,0].mask
    else:
        return m.mask 
 # --------
# 2.05 - Tropical Mask
# --------
def tropics_unmasked(res='4x5', saizlopez=False, pradosroman=False, \
        mask2D=False ):
    """ 
    Get mask for all areas but tropics

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m=np.zeros( get_dims4res( res ) )
    if saizlopez or pradosroman:
        lats = np.arange(-20, 20,1 )
    else:
        lats = np.arange(-22, 22,1 )
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:] =  1

    # Create a np.ma mask and return solely mask
    m = np.ma.masked_not_equal(m, 1)
    
    # Return 2D or 3D?
    if mask2D:
        return m[...,0].mask
    else:
        return m.mask

# --------
# 2.06 - Mid Lats Mask
# --------
def mid_lats_unmasked( res='4x5', saizlopez=False, pradosroman=False, \
            mask2D=False ):
    """ 
    Get mask for all areas but mid latitudes

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )    
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m=np.zeros(get_dims4res(res))
    if saizlopez:
        lats = np.concatenate((np.arange(-50, -20,1 ), np.arange(20, 50,1 )))
    elif pradosroman:
        lats = np.concatenate((np.arange(-60, -20,1 ), np.arange(20, 61,1 )))
    if ( (not saizlopez)  and (not saizlopez) and (res=='4x5') ):
        lats = np.concatenate((np.arange(-50, -26,1 ), np.arange(26, 51,1 )))
    if ( (not saizlopez)  and (not saizlopez) and (res=='2x2.5') ):
        lats = np.concatenate((np.arange(-50, -24,1 ), np.arange(24, 51,1 )))

    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:] =  1

    # Create a np.ma mask and return solely mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[...,0].mask
    else:
        return m.mask
# --------
# 2.07 - 40N to 40S Mask
# --------
def mask_lat40_2_40( res='4x5', mask2D=False ):
    """ 
    Get mask for all areas but 40S-40N latitudes (40 to 40 deg latitude) 

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )    
    """    

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m=np.zeros(get_dims4res(res))
    lats = np.arange(-42, 42,1 ) # use 42 due to 4x5 grid
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:] =  1

    # Create a np.ma mask and return solely mask
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[...,0].mask
    else:
        return m.mask    
# --------
# 2.08 - Extra Tropics Mask
# --------
def extratropics_unmasked( res='4x5', mask2D=False):
    """ Get mask for all areas but extratropics
    NOTES:
        - returns 3D array as default ( to return 2D set mask2D=True )    
    """    

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m=np.zeros(get_dims4res(res))
    lats = np.concatenate((np.arange(-89, -26,1 ), np.arange(26, 90,1 )))
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:] =  1
    # Create a np.ma mask 
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[...,0].mask
    else:
        return m.mask
# --------
# 2.09 - Create unmask array ( for ease of dataprocessing )
# --------
def all_unmasked( res='4x5', mask2D=False ):
    """ Get un-masked mask of size GEOS-Chem dimensions 
    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )    
    """

    m = np.ma.array( np.ones(get_dims4res(res) ), mask=False)

    # Return 2D or 3D?
    if mask2D:
        return m[...,0].mask
    else:
        return m.mask

# --------
# 2.10 - Maskes Regions by Pressure
# --------
def mask_3D( hPa, sect, MBL=True, res='4x5', extra_mask=None,    \
    M_all=False, use_multiply_method=True, trop_limit=False, \
    verbose=True, debug=False ):
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
        print 'mask_3D called for sect={}, use_multiply_method={}'.format( \
            sect, use_multiply_method ) + ', M_all={}, and '.format( M_all )+ \
            'with debug: {}, verbose:{}'.format( sect, debug, verbose )

    # Get atmospheric region as case defining lower and upper bounds
    cases = { 
    'BL': [1200., 900.], 'MBL': [1200., 900.], 'FT': [ 900., 350. ],  \
    'UT': [ 350., 75.], 'All' : [ 1200., 75.]
    }
    l, h = cases[ sect ] 

    # --- Mask between upper and lower values
    m=np.ones( get_dims4res(res) )
    m[ (hPa >=l) ] = 0
    m[ (hPa <h) ] = 0
    logging.debug( 'Sect={}, l={}, h={}'.format(sect, l, h) )
    logging.debug( '{}'.format( \
        *[ [i.min(), i.max(), i.mean(), i.sum(), i.shape] for i in [m] ]) )

    # Mask off the 'sect' area that still equals 1
    m = np.ma.masked_equal(m, 1 )

    # --- Remove above the "chemical tropopause" from GEOS-Chem (v9-2)
    if trop_limit:
        m = m[...,:38] 

    if not isinstance(extra_mask, type(None) ):
        # only consider MBL
        if ( MBL and sect == 'BL' ) or (sect == 'MBL') or M_all : 
            if use_multiply_method:
                return m.mask * extra_mask  * land_unmasked( res )    
            else:
                print "WARNING: needs 3D arrays, use 'mask_all_but' instead"
                sys.exit()
        else:    
            if use_multiply_method:
                return m.mask * extra_mask            
            else:
                print "WARNING: needs 3D arrays, use 'mask_all_but' instead"
                sys.exit()

    # --- Only consider MBL (or MFT/MFT for Saiz-Lopez 2014 comparison)
    if ( MBL and sect == 'BL' ) or (sect == 'MBL') or M_all: 
        if use_multiply_method:    
            return m.mask * land_unmasked( res )
        else:
            land_unmasked_ = mask_all_but( 'Land', mask3D=True, res=res, \
                use_multiply_method=False, trop_limit=trop_limit )

            # MBL unmasked 
            m = np.logical_or( np.logical_not( m.mask ), \
                np.logical_not( land_unmasked_ )  )

            # Invert as function expected to return oposite
            m = np.logical_not(  m )

            return  m   # m is a mask here

    return m.mask


# --------
# 2.11 - Custom 2D (Lat) Mask
# --------
def lat2lat_2D_unmasked( lowerlat=None, higherlat=None, res='2x2.5', \
            debug=False ):
    """ 
    Takes a lower and higher latitude value and then creates mask of area boxes outside 
    given laitidue limits.
    """

    # Get vars
    lon_c, lat_c, NIC = get_latlonalt4res( res=res, centre=True )

    # mask between upper and lower values
    lats =  [ i for i in lat_c if ( (i>=lowerlat) and (i<higherlat) )]
    lats = [ get_gc_lat(i, res=res) for i in  lats ]

    # fill all lat and lon True or False
    m=np.zeros(get_dims4res(res))[:,:,0]
    print m.shape, np.sum(m)
    for i in lats:
        m[:,i] = 1
    m = np.ma.masked_not_equal(m, 1)
    return m.mask

# --------
# 2.12 - South Pole mask
# --------
def southpole_unmasked(  res='4x5', mask2D=False ):
    """ 
    Return global np.ma with southpole masked

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )    
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    # adjust for resolution at grid start points at 62
    if res == '4x5': 
        lats  = np.arange(-89, -62,1 ) # define S pole as > 60S
    else:
        lats  = np.arange(-89, -60,1 ) # define S pole as > 60S
#    lats = np.arange(-89, -80,1 ) # define S pole as > 80S
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:] =  1

    # Create a np.ma mask 
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[...,0].mask
    else:
        return m.mask

# --------
# 2.13 - North Pole mask
# --------
def northpole_unmasked(  res='4x5', mask2D=False ):
    """ 
    Return global np.ma with northpole masked 

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )        
    """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m = np.zeros(get_dims4res(res))
    # adjust for resolution at grid start points at 62
    if res == '4x5': 
        lats = np.arange(62, 90,1 )  # define N pole as > 60N 
    else:
        lats = np.arange(60, 90,1 )  # define N pole as > 60N
#    lats = np.arange(80, 90,1 )  # define N pole as > 80N
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:] =  1

    # Create a np.ma mask 
    m = np.ma.masked_not_equal(m, 1)

    # Return 2D or 3D?
    if mask2D:
        return m[...,0].mask
    else:
        return m.mask

# --------
# 2.14 - North Hemisphere mask
# --------
def NH_unmasked(  res='4x5', mask2D=False ):
    """ 
    Return global np.ma with north hemisphere masked 

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )    
    """
    m=np.ma.array( np.ones( get_dims4res(res) ), mask=True)
    if res=='4x5':
        lats = np.arange(1, 91,1 )
    elif res == '2x2.5':
        lats = np.arange(0, 89,1 )
        print 'CHECK (NH) mask for non 4x5 resolutions'
#        sys.exit(0)
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:].mask =  False

    # Return 2D or 3D?
    if mask2D:
        return m[...,0].mask
    else:
        return m.mask

# --------
# 2.15 - South Hemisphere mask
# --------
def SH_unmasked(  res='4x5', mask2D=False ):
    """ 
    Return global np.ma with south hemisphere masked 

    NOTES:
     - returns 3D array as default ( to return 2D set mask2D=True )    
    """
    m=np.ma.array( np.ones( get_dims4res(res) ), mask=True)
    if res=='4x5':
        lats = np.arange(-89, 0,1 )        
    if res == '2x2.5':
        lats = np.arange(-90, 0,1 )
        print 'CHECK (SH) mask for non 4x5 resolutions'
#        sys.exit(0)
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:].mask =  False

    # Return 2D or 3D?
    if mask2D:
        return m[...,0].mask
    else:
        return m.mask
        
# --------
# 2.16 - Get Analysis maskes
# --------
def get_analysis_masks( masks='basic',  hPa=None, M_all=False, res='4x5',\
        saizlopez=False, r_pstr=True, wd=None, trop_limit=True, mask4D=False, \
        use_multiply_method=True, debug=False ):
    """
    Return list of mask arrays for analysis 

    NOTES:
    - For comparisons with Saiz-Lopez et al. 2014, set M_all to True, this maskes the 
    UT and FT as well as the BL  (therefore gives MBL, MFT, MUT) 
    
    """
    # --- Check hPa has been provided as arg.
    if isinstance( hPa, type(None) ):
        print 'ERROR: Please provide array of hPa to get_analysis_masks'

    if masks=='full':
        # ---- List of masks
        mtitles = [ 
        'All', 'Ocean', 'Land','Ice', 'All Sur.',  'Ocean Sur.', 'Land Sur.'
        , 'Ice Sur.', 'NH', 'SH', 'Tropics', 'Ex. Tropics', 'Mid Lats', 
        'Ocn. 50S-50N',  '50S-50N' 
        ]
        
        tsects3D= [  'MBL', 'BL','FT',  'UT']

        # get MBL, FT and UT maskes
        sects3D = [ 
            mask_3D(hPa, i, MBL=False, M_all=M_all, res=res, 
            use_multiply_method=use_multiply_method ) for i in tsects3D   
        ]
        
        # ---- Use non-pythonic mulitply method?
        if use_multiply_method:

            maskes = [ mask_all_but(i, trop_limit=trop_limit, mask3D=True, \
                use_multiply_method=True, res=res ) for i in mtitles ]    

            # if comparison with saiz-lopez 2014, 
            if M_all:
                ind = [n for n,i in enumerate( mtitles) if not ('MBL' in i) ]
                for n in ind:
                    maskes[n] = maskes[n]*land_unmasked(res=res)

        # --- Use pythonic approach
        else:
            maskes = [ mask_all_but(i, trop_limit=trop_limit, mask3D=True, \
                use_multiply_method=False, res=res ) for i in mtitles ]

            # If not 'use_multiply_method', then invert hPa masks
            sects3D = [ np.logical_not(i) for i in sects3D ]

        # Add to mask and mask title lists
        maskes = maskes + sects3D
        mtitles = mtitles + tsects3D

        # Also create print strings...    
        npstr ='{:<12}'*len(maskes)
        pstr ='{:<12,.2f}'*len(maskes)

    if masks=='basic':
        tsects3D= [ 'All', 'MBL', 'BL','FT',  'UT']
        mtitles = [ i+' (All)' for i in tsects3D ]    + \
        [ i+' (Tropics)' for i in tsects3D ]  +   \
        [i+' (Mid Lats)' for i in tsects3D ]     

        # standard maskes none, tropics, mid-lats (3)
        maskes = [ 
            np.logical_not( i)  for i in all_unmasked(res=res), \
            tropics_unmasked(res,saizlopez=saizlopez), \
            mid_lats_unmasked(res) ]

        # additional masks - tsects3D (4+1) * standard maskes (3)
        dmaskes = [ 
        [ mask_3D(hPa, i, MBL=False, extra_mask=mask, M_all=M_all, res=res ) 
        for i in tsects3D ]  for mask in  maskes ]

        # unpack and set as maskes (12) to list (3)
        dmaskes = [j for i in dmaskes for j in i ]
        print [ len(i) for i in maskes, dmaskes, mtitles, tsects3D ]
        maskes = dmaskes
        print [ len(i) for i in maskes, dmaskes, mtitles, tsects3D ]

        # if comparison with saiz-lopez 2014 appli marine mask to all...
        if M_all:
            ind = [n for n,i in enumerate( mtitles ) if not 'MBL' in i ]
            for n in ind:
                maskes[n] = maskes[n]*land_unmasked(res=res)

        print [ len(i) for i in maskes, dmaskes, mtitles, tsects3D ]
        # Also create print strings...    
        npstr ='{:<15}'*len(maskes)
        pstr ='{:<15,.2f}'*len(maskes)

    if masks=='trop_regions':
        mtitles = [ 'BL','FT',  'UT'] 
        maskes = [  mask_3D(hPa, i, M_all=M_all, MBL=False, res=res )[:,:,:38]  \
            for i in  mtitles ]
        # Also create print strings...    
        npstr ='{:<15}'*len(maskes)
        pstr ='{:<15,.2f}'*len(maskes)


    # Only consider the "chemical troposphere" - according v9-2
    if trop_limit:
        maskes = [i[:,:,:38] for i in maskes]

    # Create 4D array by concatenating through time dimension
    # ( assuming year long array of 1 months )
    if mask4D:
        for n, mask in enumerate( maskes ):
            if any( [ (mask.shape[-1] == i) for i in [12] ] ):
                pass
            else: # concatenate dimensions
                maskes[n] = np.concatenate( [ mask[...,None] ]*12, axis=3 )

    if r_pstr:
        return maskes, mtitles, npstr, pstr
    else:
        return maskes, mtitles

# --------
# 2.17 -  Retrieve individual 4D mask of locations except region given
# --------
def mask_all_but( region='All', M_all=False, saizlopez=False, \
        res='4x5', trop_limit=True, mask2D=False, mask3D=False, mask4D=False, \
        use_multiply_method=True, verbose=False, debug=False ):
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

    Returns
    -------
    (np.ma.mask) or (np.array) (later if use_multiply_method==True)

    Notes
    -----
    "unmask_all" yeilds completely unmasked array
    function was oringialyl used to mulitple masks, however, this approch is 
    unpythonic and therefore reccomended against. 
    """
    logging.info( 'mask_all_but called for region {}'.format(region) )
    # --- Setup cases...   
    # ( except None, unmask_all and global to retrive no mask )
    case = {
    'Tropics':0,
    'tropics' : 0,
    'mid_lats' : 1,
    'Mid Lats': 1, 
    'Mid lats':1, 
    'south_pole' : 2, 
    'south pole' : 2, 
    'north_pole' : 3,     
    'north pole' : 3,     
    None : 4, 
    'unmask_all': 4, 
    'All' : 4,
    'global': 4,
    # NEED TESTING ...
    'Extratropics':5, 
    'Ex. Tropics':5, 
    'Oceanic':6, 
    'Ocean':6,
    'NH': 7, 
    'SH':8,
    'Ice': 10,
    'Land' : 11, 
    'lat40_2_40':12, 
    'Ocean Tropics': 13, 
    'Oceanic Tropics': 13, 
    'Ocn. Trop.': 13, 
    'Land Tropics': 14,
    'All Sur.': 15,
    'surface': 15,
    'Ocean Sur.': 16, 
    'Land Sur.': 17 , 
     'Ice Sur.' : 18, 
    'lat50_2_50':19, 
    '50S-50N':19, 
#    'Oceanic lat50_2_50': 20, 
    'Ocn. 50S-50N' : 20, 
#     'South >60': 2,
#      'North >60': 3
    'North Sea' : 21, 
    'Med. Sea' : 22, 
    'Mediterranean. Sea' : 22, 
    'Black Sea' : 23, 
    'Irish Sea' : 24, 
    'Europe' : 25, 
    'EU' : 25, 
#    'Surface BL': 26, 
    }[region]


    # --- This is a simple way of using masks ( as multiplers )
    # i.e. all (future) functions should have use_multiply_method=False
    # and not use the code below 
    if use_multiply_method:  # Kludge
        print '!'*50, 'WARNING: using mulitply method for masking. '
    
        # For case, pull mask from case list 
        if case == 0:
            mask = tropics_unmasked( res=res, saizlopez=saizlopez )
        if case == 1:
            mask = mid_lats_unmasked( res=res )
        if case == 2:
            mask = southpole_unmasked(  res=res )
        if case == 3:
            mask = northpole_unmasked(  res=res )
        if case == 4:
#        mask = np.logical_not( all_unmasked( res=res ) )
            mask = all_unmasked( res=res ) 
        if case == 5:
            mask = extratropics_unmasked(res=res)
        if case == 6:
            mask = ocean_unmasked( res=res)
        if case == 7:
            mask = NH_unmasked(  res=res )
        if case == 8:
            mask = SH_unmasked(  res=res  )
        if case == 10:
            mask =ice_unmasked( res=res )
        if case == 11:
            mask =  land_unmasked( res=res )
        if case == 12:
            mask = mask_lat40_2_40( res=res )
        if case == 13:  #  'Oceanic Tropics'
             mask = np.ma.mask_or( ocean_unmasked( res=res),  \
                    tropics_unmasked(res=res, saizlopez=saizlopez ) )
        if case == 14: # 'Land Tropics'
            mask = np.ma.mask_or( land_unmasked( res=res ),  \
                    tropics_unmasked( res=res, saizlopez=saizlopez ) )
        if case == 15: # 'All Sur.'
            mask = surface_unmasked(res=res) 
        if case == 16: # 'Ocean Sur.'
            mask = np.ma.mask_or(  surface_unmasked(res=res) , \
                     ocean_unmasked( res=res )  )
        if case == 17: # 'Land Sur.':
            mask = np.ma.mask_or( surface_unmasked(res=res) ,  \
                     land_unmasked( res=res )  )
        if case == 18: # 'Ice Sur.' 
            mask = np.ma.mask_or( surface_unmasked(res=res) ,  \
                      ice_unmasked( res=res ) )
        if case == 19:  # '50S-50N' 
            mask = lat2lat_2D_unmasked( lowerlat=-50, higherlat=50, \
                res=res )
        if case == 20: # 'Ocn. 50S-50N' 
            mask = np.ma.mask_or( lat2lat_2D_unmasked( lowerlat=-50, 
                higherlat=50, res=res ), ocean_unmasked( res=res )[...,0]  )
        if case == 21:
            mask = get_north_sea_unmasked( res=res )

        if case == 25:
            mask = get_EU_unmasked( res=res )
#        if case == 26:
#            mask = get_2D_BL_unmasked( res=res )

        # Invert mask to leave exception unmasked if used to multiply
        mask = np.logical_not(mask)

    # --- This is a more pythonic way of using masks (Use as preference)
    else:
        # For case, pull mask from case list 
        if case == 0:
            mask = tropics_unmasked( res=res, saizlopez=saizlopez )
        if case == 1:
            mask = mid_lats_unmasked( res=res )
        if case == 2:
            mask = southpole_unmasked(  res=res )
        if case == 3:
            mask = northpole_unmasked(  res=res )
        if case == 4:
#        mask = np.logical_not( all_unmasked( res=res ) )
            mask = all_unmasked( res=res ) 
        if case == 5:
            mask = extratropics_unmasked(res=res)
        if case == 6:
            mask = ocean_unmasked( res=res )
        if case == 7:
            mask = NH_unmasked(  res=res )
        if case == 8:
            mask = SH_unmasked(  res=res  )
        if case == 10:
            mask = ice_unmasked( res=res )
        if case == 11:
            mask = land_unmasked( res=res )
        if case == 12:
            mask = mask_lat40_2_40( res=res )
        if case == 13:
            mask = np.ma.mask_or( ocean_unmasked( res=res),  \
                    tropics_unmasked(res=res, saizlopez=saizlopez ) )
        if case == 14:
            mask = np.ma.mask_or( land_unmasked( res=res ),  \
                    tropics_unmasked( res=res, saizlopez=saizlopez ) )
        if case == 15: # 'All Sur.'
            mask = surface_unmasked(res=res) 
        if case == 16: # 'Ocean Sur.'
            mask = np.ma.mask_or(  surface_unmasked(res=res) , \
                     ocean_unmasked( res=res )  )
        if case == 17: # 'Land Sur.':
            mask = np.ma.mask_or( surface_unmasked(res=res) ,  \
                     land_unmasked( res=res )  )
        if case == 18: # 'Ice Sur.' 
            mask = np.ma.mask_or( surface_unmasked(res=res) ,  \
                      ice_unmasked( res=res ) )
        if case == 19:
            mask = lat2lat_2D_unmasked( lowerlat=-50, higherlat=50, \
                res=res )
        if case == 20:
            mask = np.ma.mask_or( lat2lat_2D_unmasked( lowerlat=-50, 
                higherlat=50, res=res ), ocean_unmasked( res=res )[...,0]  )
        if case == 21:
            mask = get_north_sea_unmasked( res=res )
        if case == 22:
            mask = get_unmasked_mediterranean_sea( res=res )
        if case == 23:
            mask = get_unmasked_black_sea( res=res )
        if case == 24:
            mask = get_unmasked_irish_sea( res=res )
        if case == 25:
            mask = get_EU_unmasked( res=res )
#        if case == 26:
#            mask = get_2D_BL_unmasked( res=res )


    logging.debug( 'prior to setting dimensions: {}'.format(mask.shape) )

    # Apply Saiz-Lopez Marine MFT/MUT? <= should this be before multiply op.?
    if M_all:
        if use_multiply_method:  # Kludge
            mask = mask*land_unmasked(res=res)
        else:
            # check this!!!
            mask = np.ma.mask_or( mask, land_unmasked(res=res) )

    # Ensure returned arrays are 2D
    if mask2D:
        if len( mask.shape ) == 2:
            pass
        elif len( mask.shape ) == 3:
            mask = mask[..., 0]
        elif len( mask.shape ) == 4:
            mask = mask[..., 0, 0]

    # Create 3D array by concatenating through altitude dimension
    if mask3D:
        if any( [ (mask.shape[-1] == i) for i in 38, 47 ] ):
            pass
        else: # concatenate dimensions
            if len( mask.shape ) == 3:
                mask = np.concatenate( [ mask ]*47, axis=2 )
            elif len( mask.shape ) == 2:
                mask = np.concatenate( [ mask[...,None] ]*47, axis=2 )


    # Remove above the "chemical tropopause" from GEOS-Chem (v9-2)
    if trop_limit:
        if ( len( mask.shape ) == 2 ) or mask2D:
            pass
        else:
            mask = mask[...,:38] 

    # Create 4D array by concatenating through time dimension
    # ( assuming year long array of 1 months )
    if mask4D:
        if any( [ (mask.shape[-1] == i) for i in [12] ] ):
            pass
        else: # concatenate dimensions
            mask = np.concatenate( [ mask[...,None] ]*12, axis=3 )
    logging.debug('post to setting dimensions: {}'.format(mask.shape) )
    logging.info("returning a 'mask' of type:{}".format(type(mask)) )
    return mask
    
# --------
# 2.18 - Custom 2D (Lon) Mask
# --------
def lon2lon_2D_unmasked(lowerlon, higherlon, res='2x2.5', debug=False ):
    """
    Takes a lower and higher latitude value and then creates 
    mask to given given limits.
    """

    # Get vars
    lon_c, lat_c, NIU = get_latlonalt4res( res=res, centre=True )
    if debug:
        print lon_c, lowerlon, higherlon

    # Mask between upper and lower values
    lons =  [ i for i in lon_c if ( (i>=lowerlon) and (i<higherlon) )]
    lons = [ get_gc_lon(i, res=res) for i in  lons ]

    # Fill all lat and lon True or False
    m=np.zeros(get_dims4res(res))[:,:,0]
    print m.shape, np.sum(m)
    for i in lons:
        m[i,:] = 1
    m = np.ma.masked_not_equal(m, 1)
    return m.mask
    
# --------
# 2.19 - EU mask
# --------
def get_EU_unmasked( res='1x1'  ):
    """ 
    Mask 'EU' as defined by GEOS-Chem EU grid the grid of "'0.5x0.666" resolution is 
    used by default, but any list of lat and lons could be provided and the extremities 
    would be used as the mask edges  
    """

    EU_resolutions = [ '0.25x0.3125', '0.5x0.666' ]
    if res not in EU_resolutions:
        EU_res = EU_resolutions[0] # default = '0.25x0.3125'
#        EU_res = EU_resolutions[1] # default = '0.5x0.666'
    else:
        EU_res =res

    # Get GEOS-Chem EU lat and lons
    lon, lat, NIU = get_latlonalt4res( res=EU_res )

    # mask lats
    m1 = lat2lat_2D_unmasked( lowerlat=lat.min(), higherlat=lat.max(), res=res )
    
    # mask lons
    m2 = lon2lon_2D_unmasked(lowerlon=lon.min(), higherlon=lon.max(), res=res )

    #  combine maskes 
    m = m1 + m2 
    
    return m
    
# --------
# 2.20 - Cruise track mask
# --------
def get_cruise_track_mask(  max_lon=None, min_lon=None, max_lat=None, \
        min_lat=None, unmask_water=True, res='4x5', trop_limit=True ):
    """ 
    Mask whole area of ship based research campaigns for bulk comparison
    """
    
    # only look at surface
    m = surface_unmasked( res=res, trop_limit=trop_limit )
    
    # apply ocean mask
    if unmask_water:
        m = m + ocean_unmasked(res=res)
    
    # Mask over given longitude range, if provided
    if not isinstance( max_lon, type(None) ):
        m = m  + lon2lon_2D_unmasked(lowerlon=min_lon, higherlon=max_lon, \
                        res=res )[:,:,None]
    
    # Mask over given latitude range, if provided
    if not isinstance( max_lat, type(None) ):
        m = m + lat2lat_2D_unmasked( lowerlat=min_lat, higherlat=max_lat, \
                    res=res )[:,:,None]

    # Invert 
    m = np.logical_not( m )
    
    return m


def get_north_sea_unmasked( res='0.25x0.3125' ):
    """
    A rough Mask of the North Sea for use with ~0.5/~0.25 mdodel output. 
    (inc. English channel. )

    """
    # mask latitudes of mediterranean
    # Drawing a box that include all of North sea and English channel
    # Near Brest in France 48.3669927,-4.7560745
    # Lillehammer 61.1122408,10.4386779

    # mask lats
    m1 = lat2lat_2D_unmasked( lowerlat=48.3669927, higherlat=61.1122408, res=res )
    
    # mask lons
    m2 = lon2lon_2D_unmasked(lowerlon=-4.7560745, higherlon=10.4386779, res=res )

    #  combine maskes 
    m = m1 + m2     

    # remove all water.     
#    m = np.ma.mask_or( m, ocean_unmasked( res=res)[...,0] ) # old approach
    # the below will work for muliple options. 
#    m = m + 
#    m = ocean_unmasked( res=res)[...,0] 
    
    
    # remove irish sea

    return m.mask
    
    
def get_mediterranean_sea_unmasked( res='0.25x0.3125' ):
    """
    A rough Mask of the Mediterranean Sea for use with ~0.5/~0.25 mdodel output. 
    """
    
    # mask latitudes of mediterranean
    # Drawing a box that include all of Med. Sea
    # East South corner (south Jordan) = 43.4623268,33.3809392
    # West North corner (south Jordan) = 44.17207,25.30604
    
    # add mask for Black Sea
    
    # add mask for  
    pass

def get_unmasked_black_sea( res='0.25x0.3125'):
    """
    A rough Mask of the Black Sea for use with ~0.5/~0.25 mdodel output. 
    """
    
    # Drawing a box that include all of Black. Sea
    # East South corner (south Jordan) = 43.4623268,33.3809392
    # West North corner  = 44.17207,25.30604
    pass

def get_unmasked_irish_sea( res='0.25x0.3125', unmasked_oceans=None):
    """
    A rough Mask of the Irish Sea for use with ~0.5/~0.25 mdodel output. 
    """
    
    # NE corner Glasgow - 55.8553803,-4.3725463
    # SW corner Cork - 51.8959842,-8.5332609
    # mask lats
    m1 = lat2lat_2D_unmasked( lowerlat=51.8959842, higherlat=55.855380, res=res )
    
    # mask lons
    m2 = lon2lon_2D_unmasked(lowerlon=8.5332609, higherlon=4.3725463, res=res )

    #  combine maskes 
    m = m1 + m2     

    # only consider oceans     
    m = np.ma.mask_or( ocean_unmasked( res=res)[...,0], m )    
    
    
    return mask 

def get_ODR(x=None, y=None):
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

    # create a model 
    linear = scipy.odr.Model(f) 

    # Create a Data or RealData instance.:
    mydata = scipy.odr.Data(x, y)
    
    # Instantiate ODR with your data, model and initial parameter estimate.:
#    myodr = scipy.odr.ODR(mydata, linear, beta0=[1., 2.])
    myodr = scipy.odr.ODR(mydata, linear, [0., 1.],  maxit = 10000 )
    
    # Run the fit.:
    myoutput = myodr.run()

    # Examine output.:
    myoutput.pprint()
    
    return myoutput

       


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# ---------------- Section X -------------------------------------------
# -------------- Redundant Functions
# --------------------------------------------------------------------------
# 
# NOTE(s): 
# (1) These are retained even though they are redundant for back compatibility
# (2) It is not advised to use these. 

