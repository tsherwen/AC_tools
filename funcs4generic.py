""" Generic functions for use with GEOS-Chem/Data Analysis.
    
    These functions were written whilst learning python, vestigial 
    inefficient code is being removed/updated. 
    Where external code is used credit is given. """

# ------------------------- Section 0 -----------------------------------------
# -------------- Required modules:
#
#!/usr/bin/python
# -- I/O / Low level                                                                                
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# -- time                                                                                           
import time

# --- math
from math import radians, sin, cos, asin, sqrt, pi, atan2

# --- tms modules
from AC_tools.funcs4core import *
from AC_tools.funcs_vars import *

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 1 ----- Generic processing 
# 1.01 - Chunker (splits data into chunks )
# 1.02 - length of file
# 1.03 - round values to given base (stipulate if not integer )
# 1.04 - CSV opener - REDUDENT
# 1.05 - Data "binner" (splits x by y, with given conditions )
# 1.06 - Translate columns to rows  (in python 3)
# 1.07 - Count files in dir
# 1.08 - Bin range (redundent? )
# 1.09 - Renames files
# 1.10 - Get X and Y coordinates for a given grid
# 1.11 - hPa to Km 
# 1.12 - Save as pdf
# 1.13 - save as multi page pdf
# 1.14 - Obs finder for a given grid
# 1.15 - Sort (GAW) sites by latitude
# 1.16 - find nearest values in array
# 1.17 - get number suffix
# 1.18 - Get shortest distance on a sphere ( for finding closest point )
# 1.19 -  Get logritymcally spaced integers

# ------------- ------------- ------------- ------------- ------------- 
# ---- Section 2 ----- Masks, for data analysis involving maps
# 2.01 - Ocean Mask
# 2.02 - Mask land
# 2.03 - Mask Ice
# 2.04 - Mask Surface
# 2.05 - Tropical mask
# 2.06 - Mid Lats mask
# 2.07 - 40N to 40S Mask
# 2.08 - Extra Tropics Mask
# 2.09 - Mask 
# 2.10 - MBL, FT, UT - 3D mask by pressure
# 2.11 - Custom 2D mask by Lat (just give Lats) 
# 2.12 - South Pole mask
# 2.13 - North Pole mask
# 2.14 - North Hemisphere mask 
# 2.15 - South Hemisphere mask
# 2.16 - Get Analysis maskes (list of masks and titles)
# 2.17 - Get single mask ( single maks: e.g. tropical ) 
# 2.18 - Custom 2D mask by Lat (just give Lats) 
# 2.19 - EU mask
    
# -------------------------- Section 7 -------------------------------
# -------------- Generic Processing
#

# -------------
# 1.01 -  Split arrays into chunks (e.g. days)
# -------------
def chunks(l, n):
    """ Split list in chunks - useful for controlling memory usage """
    if n < 1:
        n = 1
    return [l[i:i + n] for i in range(0, len(l), n)]
    
# -------------
# 1.02 - Get len of file
# -------------
def file_len(fname):
    """ Get length of file """
    return sum(1 for line in open(fname))

# -------------
# 1.03 - My round - credit: Alok Singhal
# -------------
def myround(x, base=5, integer=True, round_up=False):
    """ Round up values - mappable function for pandas processing """

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
# 1.04 - csv, opener - upgrade to pandas makes this redundent
# -------------
#def csv_opener(fn_,big=None,delimiter=',',debug=False):
#    if debug:
#        print 'csv_opener'
#    reader = open(fn_, 'rb')
#    row_list=[]
#    for ii, row in enumerate(reader):
#        row = row.strip().split(delimiter)
#        if debug:
#            print row 
#        row_list.append(row)
#    return row_list

# -------------
# 1.05 - "Binner" - bins data by a given var for a given bin_size - tms
# -------------
def bin_data( data, bin_by, bin_size, set_min=None, set_max=None, debug=False ):
    """ Redundant program? - manual binner of data using list comprehension  """
    # set bin dimensions, round to get ranges and then extend if min and max are not captured.
    bmin, bmax = myround(np.ma.min(bin_by), bin_size, integer=False), myround( np.ma.max(bin_by), bin_size, integer=False)    
    if bmin > np.ma.min(bin_by):
        bmin =  bmin- bin_size
    if bmax < np.ma.max(bin_by):
        bmax =  bmax+ bin_size
    bins = np.arange( bmin, bmax+bin_size, bin_size ) 
    if all( [(not isinstance(i, type(None))) for i in set_min, set_max ]):
        bins = np.arange( set_min, set_max+bin_size, bin_size ) 
    if debug:
        print 'min: {}, max: {}, and # of bins: {}'.format(bmin, bmax, len(np.ma.arange( bmin, bmax, bin_size)) )

    # Loop and bin for given dimensions.
    for bin in bins:
        binned = [ d_ for n_, d_ in enumerate(data)  if ( (bin+bin_size) > bin_by[n_] >= ( bin ) )]
        if debug:
            print 'binning between {} and {}'.format( bin,  (bin+bin_size) ) 
        try :
            binned_list.append( binned )
        except:
            binned_list = [ binned ]
        if debug:
            print bin, [ len(i) for i in [ binned_list, data, binned ]]

    # remove any extra unfilled bins (if not using pre-defined min and max )
    if not all( [(not isinstance(i, type(None))) for i in set_min, set_max ]):
        empty_bins = [ n for n, i in enumerate(bins) if ( len(binned_list[n]) < 1) ]
        if len( empty_bins) > 0 :
            print 'removing empty bins: {}'.format( empty_bins )
            print [type(i) for i in  binned_list , bins  ]
            bins   = list(bins)
            [ [ i.pop(ii) for ii in empty_bins[::-1] ] for i in binned_list , bins ]
            bins   = np.ma.array(bins)
    return binned_list , bins

# --------------
# 1.06 - translsate columns to rows 
# ------------- 
# written for python 3! re-write
#def columns_to_rows(filename):
#        with open(filename) as f:
#            lis=[x.split() for x in f]
#
#        for x in zip(*lis):
#            for y in x:
#                print(y+'\t',end='')
#            print('\n')

# --------------
# 1.07 - Count number of files/the filenames that contain certain phrase/number conbinations (e.g. dates) - tms 
# -------------
def counter_directory_contains_files(model_path,must_contain):
    """ Count number of files in directory """
    model_ouput_file_counter = len(glob.glob1(model_path,must_contain))
    return model_ouput_file_counter
            
# ----
# 1.08 - bin by range
# ----
def bin_w_range(min, max, bin_size, data, press):
    """ REDUNDENT? ( see 1.05) - manual "binner" of data using list 
        comprehension """
    bins = []
    for bin in np.arange(min, max, bin_size):
        binned = [ d_ for n_, d_ in enumerate(data)  \
            if ( (bin+bin_size) > press[n_] >= ( bin ) )]
        try :  
            binned_list.append(binned)
        except:
            binned_list = [ binned ]
        bins.append( bin )
        print 'bin: {} to {}'.format(bin,bin+bin_size), \
             [ len(i) for i in [ binned_list, data, binned ]] 
    return binned_list#, np.array( bins )

# ----
# 1.09 - Rename vars/strs in files
# ----
def replace_strs_in_files( wd, input_str, output_str, debug=False ):
    """ replace text in files """
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
    """Takes lon,lat values for a point (e.g. station) and arrays of longitude
     and latitudes indicating edges of gridboxes.
    Returns the index of the gridbox corresponding to this point.  Could be
     easily extended to return data values corresponding
    to points."""
    hasobs,lonrange,latrange=np.histogram2d([Lon],[Lat],[lon_edges, lat_edges])
    gridindx,gridindy=np.where(hasobs>=1)
    if not gridindx:
        if debug:
            print 'Lat, lon outside of x,y range.  Assigning -1 for', Lon, Lat
        return -1, -1
    else:
        #print Lon, Lat, gridindx, gridindy
        return gridindx[0],gridindy[0]

# --------------
# 1.11 - Convert from hPa to km or vice versa.
# -------------
#=> mv's to funcs4core.py
        
# --------   
# 1.12 - Save as pdf.
# --------
def plot2pdf(title='new_plot', fig=None, rasterized=True, dpi=160,
                         justHH=False, no_dstr=True, save2png=True, 
                         save2eps=False, debug=False ):
    """ Save figures (matplotlib) to pdf file """

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
    pdf.savefig( dpi=dpi )
    pdf.close()

    # Also export to png and eps?
    if save2png:
        type += '/PDF' 
        plt.savefig(npdf+'.png', format='png', dpi=dpi)
    if save2eps:
        type += '/EPS' 
        plt.savefig(npdf+'.eps', format='eps', dpi=dpi)
    print type+' saved & Closed as/at: ', npdf

# --------   
# 1.13- Save as mulitple page pdf.
# --------
def plot2pdfmulti(pdf=None, title='new_plot', rasterized=True, 
                                    dpi=160, open=False, close=False, 
                                    justHH=False, no_dstr=False ):
    """ Save figures (matplotlib) to pdf file with multiple pages """

    # set save directory ( using default directory dictionary )
    from funcs_vars import get_dir # Kludge
    wd =get_dir('ppwd')

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
        """ values that have a given lat, lon and alt """
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
    """ Order given list of GAW sties by latitudes """ 

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
        Find nearest point. Adapted from HappyLeapSecond's Stackoverflow 
        answer. http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    idx = (np.abs(array-value)).argmin()
    return idx

# --------   
# 1.17 - Get suffix for number
# --------
def get_suffix(n):
    """ Add the appropriate suffix (th/st/rd) to any number given """
    ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n/10%10!=1)*(n%10<4)*n%10::4])
    return ordinal(n)

# --------   
# 1.18 - Get shortest distance on a sphere ( for finding closest point )
# --------
def get_shortest_in(needle, haystack, r_distance=False):
    """needle is a single (lat,long) tuple.
        haystack is a numpy array to find the point in
        that has the shortest distance to needle
    
        adapted from stackover (Credit: jterrace):
        http://stackoverflow.com/questions/6656475/python-speeding-up-geographic-comparison
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
# 1.19 - get logarithmically spaced integers
# --------
def gen_log_space(limit, n):
    """
    credit: Avaris
http://stackoverflow.com/questions/12418234/logarithmically-spaced-integers
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

# -------------------------- Section 2 -------------------------------
# -------------- Maskes for data analysis
#

# --------
# 2.01 - Ocean mask
# --------
def ocean_unmasked(res='4x5', debug=False):
    """ Get ocean mask from GEOS-Chem LWI ( at given resolution) 
        <= update 
        NOTE: this currently returns mask with all except ocean masked """

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
    """ Get land mask from GEOS-Chem LWI ( at given resolution) """
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
    """ Get ice mask from GEOS-Chem LWI ( at given resolution) """
    # Create a np.ma mask 
    m= np.logical_not( (land_unmasked(res)*ocean_unmasked(res)) )
    if debug:
        print mask, mask.shape, 
    return m
 
# --------
# 2.04 - Surface mask
# --------
def surface_unmasked( res='4x5', trop_limit=False, debug=False ):
    """ Get surface mask  """

    # Create a np.ma mask 
    m=np.ma.array( np.zeros( get_dims4res(res) ), mask=False)
    m[...,0] = 1
    m = np.ma.masked_not_equal(m, 1)

    if trop_limit:
        m = m[...,:38]
    if debug:
        print mask, mask.shape
    return m.mask
 
 # --------
# 2.05 - Tropical Mask
# --------
def tropics_unmasked(res='4x5', saizlopez=False, pradosroman=False):
    """ Get tropics mask  """

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
    return m.mask

# --------
# 2.06 - Mid Lats Mask
# --------
def mid_lats_unmasked( res='4x5', saizlopez=False, pradosroman=False ):
    """ Get mask of mid latitudes """

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
    return m.mask

# --------
# 2.07 - 40N to 40S Mask
# --------
def mask_lat40_2_40( res='4x5' ):
    """ Get mask of latitudes (40 to 40 deg latitude) """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m=np.zeros(get_dims4res(res))
    lats = np.arange(-42, 42,1 ) # use 42 due to 4x5 grid
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:] =  1

    # Create a np.ma mask and return solely mask
    m = np.ma.masked_not_equal(m, 1)
    return m.mask
    
# --------
# 2.08 - Extra Tropics Mask
# --------
def extratropics_unmasked( res='4x5'):
    """ Get mask of extratropics  """

    # Create a mask of 1s for chosen area and or 0s elsewhere
    m=np.zeros(get_dims4res(res))
    lats = np.concatenate((np.arange(-89, -26,1 ), np.arange(26, 90,1 )))
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:] =  1
    # Create a np.ma mask 
    m = np.ma.masked_not_equal(m, 1)
    return m.mask

# --------
# 2.09 - Create unmask array ( for ease of dataprocessing )
# --------
def all_unmasked( res='4x5' ):
    """ Get un-masked mask of size GEOS-Chem dimensions  """
    return np.ma.array( np.ones(get_dims4res(res) ), mask=False).mask

# --------
# 2.10 - Maskes Regions by Pressure
# --------
def mask_3D( hPa, sect, MBL=True, res='4x5', extra_mask=None,    \
    M_all=False, verbose=True, debug=False ):
    """
    Maskes by pressure array  (required shape: 72,46,47), with conditions,
    set by given cases for MBL, UT, FT 
    """
    if verbose:
        print 'mask_3D called for sect {}, with debug: {}, verbose:{}'.format(\
                    sect, debug, verbose )

    # Get case
    cases = { 
      'BL': [1200., 900.], 'MBL': [1200., 900.], 'FT': [ 900., 350. ]         \
     , 'UT': [ 350., 75.], 'All' : [ 1200., 75.]
    }
    l, h = cases[ sect ] 

    # mask between upper and lower values
    m=np.ones( get_dims4res(res) )
    m[ (hPa >=l) ] = 0
    m[ (hPa <h) ] = 0
    if debug:
        print sect, l, h, [ [np.ma.min(i), np.ma.max(i),   \
                    np.ma.mean(i), np.ma.sum(i), i.shape ] for i in [ m  ] ]
    m = np.ma.masked_equal(m, 1 )
    if not isinstance(extra_mask, type(None) ):
        # only consider MBL
        if ( MBL and sect == 'BL' ) or (sect == 'MBL') or M_all : 
            return m.mask * extra_mask  * land_unmasked( res )    
        else:    
            return m.mask * extra_mask

    # only consider MBL (or MFT/MFT for Saiz-Lopez 2014 comparison)
    if ( MBL and sect == 'BL' ) or (sect == 'MBL') or M_all: 
        return m.mask * land_unmasked( res )
    return m.mask

# --------
# 2.11 - Custom 2D (Lat) Mask
# --------
def lat2lat_2D_unmasked( lowerlat=None, higherlat=None, res='2x2.5', debug=False ):
    """
    Takes a lower and higher latitude value and then creates 
    mask to given given limits.
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
def southpole_unmasked(  res='4x5' ):
    """ Return global np.ma with southpole masked """

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
    return m.mask

# --------
# 2.13 - North Pole mask
# --------
def northpole_unmasked(  res='4x5' ):
    """ Return global np.ma with northpole masked """

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
    return m.mask

# --------
# 2.14 - North Hemisphere mask
# --------
def NH_unmasked(  res='4x5' ):
    """ Return global np.ma with north hemisphere masked """
    m=np.ma.array( np.ones( get_dims4res(res) ), mask=True)
    if res=='4x5':
        lats = np.arange(1, 91,1 )
    if res == '2x2.5':
        lats = np.arange(0, 89,1 )
        print 'CHECK mask for non 4x5 resolutions'
#        sys.exit(0)
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:].mask =  False
    return m.mask

# --------
# 2.15 - South Hemisphere mask
# --------
def SH_unmasked(  res='4x5' ):
    """ Return global np.ma with south hemisphere masked """
    m=np.ma.array( np.ones( get_dims4res(res) ), mask=True)
    if res=='4x5':
        lats = np.arange(-89, 0,1 )        
    if res == '2x2.5':
        lats = np.arange(-90, 0,1 )
        print 'CHECK mask for non 4x5 resolutions'
#        sys.exit(0)
    lats = [ get_gc_lat(i, res=res) for i in lats ]
    for i in lats:
        m[:,i,:].mask =  False
    return m.mask

# --------
# 2.16 - Get Analysis maskes
# --------
def get_analysis_masks( masks='basic',  hPa=None, M_all=False, res='4x5',\
        saizlopez=False, r_pstr=True, wd=None, trop_limit=True ):
    """
    Return list of mask arrays for analysis 

    For comparisons with Saiz-Lopez et al. 2014, set M_all to True,
    this maskes the UT and FT as well as the BL 
    (therefore gives MBL, MFT, MUT) 
    
    """
    # Check hPa has been provided as arg.
    if isinstance( hPa, type(None) ):
        print 'ERROR: Please provide array of hPa to get_analysis_masks'

    if masks=='full':
        mtitles = [ 
        'All', 'Ocean', 'Land','Ice', 'All Sur.',  'Ocean Sur.', 'Land Sur.'
        , 'Ice Sur.', 'NH', 'SH', 'Tropics', 'Ex. Tropics', 'Mid Lats'
        ]
        tsects3D= [  'MBL', 'BL','FT',  'UT']
        maskes = [ 
        np.logical_not( i) for i in  all_unmasked(res=res), ocean_unmasked(res=res) \
        , land_unmasked(res=res), ice_unmasked(res=res), surface_unmasked(res=res)] +  [                        
        np.logical_not( i)*np.logical_not( surface_unmasked(res=res) ) 
        for i in ocean_unmasked(res=res), land_unmasked(res=res)
        , ice_unmasked(res=res) ] + [  
        np.logical_not( i) for i in 
        NH_unmasked( res=res), SH_unmasked( res=res) ]    + [   
        np.logical_not( i)  for i in  tropics_unmasked(res, saizlopez=saizlopez)
        , extratropics_unmasked(res),mid_lats_unmasked(res, saizlopez=saizlopez) 
        ]

        # if comparison with saiz-lopez 2014, 
        if M_all:
            ind = [n for n,i in enumerate( mtitles) if not ('MBL' in i) ]
            for n in ind:
                maskes[n] = maskes[n]*land_unmasked(res=res)

        # get MBL, FT and UT maskes
        sects3D = [ 
        mask_3D(hPa, i, MBL=False, M_all=M_all, res=res ) for i in tsects3D   
        ]
        maskes = maskes + sects3D
        mtitles = mtitles + tsects3D

        # also create print strings...    
        npstr ='{:<12}'*len(maskes)
        pstr ='{:<12,.2f}'*len(maskes)

    if masks=='basic':
        tsects3D= [ 'All', 'MBL', 'BL','FT',  'UT']
        mtitles = [ i+' (All)' for i in tsects3D ]    + \
        [ i+' (Tropics)' for i in tsects3D ]  +   \
        [i+' (Mid Lats)' for i in tsects3D ]     

        # standard maskes none, tropics, mid-lats (3)
        maskes = [ 
            np.logical_not( i)  for i in all_unmasked(res=res), tropics_unmasked(res,\
            saizlopez=saizlopez), mid_lats_unmasked(res) ]

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
        npstr ='{:<15}'*len(maskes)
        pstr ='{:<15,.2f}'*len(maskes)

    if masks=='trop_regions':
        mtitles = [ 'BL','FT',  'UT'] 
        maskes = [ 
        mask_3D(hPa, i,M_all=M_all, res=res )[:,:,:38]  
        for i in  mtitles
                ]

    # Only consider the "chemical troposphere" - according v9-2
    if trop_limit:
        maskes = [i[:,:,:38] for i in maskes]

    if r_pstr:
        return maskes, mtitles, npstr, pstr
    else:
        return maskes, mtitles

# --------
# 2.17 -  Retrieve individual 4D mask of locations except region given
# --------
def mask_all_but( region='All', M_all=False, saizlopez=False, \
            res='4x5', mask3D=False, trop_limit=True, \
            use_multiply_method=True, verbose=False, debug=False ):
    """ Mask selector for analysis. global mask provided for with given region 
        unmasked 
        NOTE: 
            (1) "unmask_all" yeilds completely unmasked array
            (2) fucntion was oringialyl used to mulitple masks, however, 
                this approch is   """

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
    'unmask_all': 4, 
    'All' : 4,
    'global': 4,
    # NEED TESTING ...
    'Extratropics':5, 
    'Oceanic':6, 
    'Ocean':6,
    'NH': 7, 
    'SH':8,
    'Ice': 10,
    'Land' : 11, 
    'lat40_2_40':12, 
    'Oceanic Tropics': 13, 
    'Land Tropics': 14,
#     'South >60': 2,
#      'North >60': 3
    }[region]
    
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
        mask =ice_unmasked( res=res )
    if case == 11:
        mask =  land_unmasked( res=res )
    if case == 12:
        mask = mask_lat40_2_40( res=res )
    if case == 13:
         mask = ocean_unmasked( res=res)*tropics_unmasked(res=res, \
                            saizlopez=saizlopez )
    if case == 14:
        mask = land_unmasked( res=res )*tropics_unmasked( res=res, saizlopez=saizlopez )

    # Invert mask to leave exception unmasked if used to multiply
    # i.e. all (future) functions should have use_multiply_method=True 
    if use_multiply_method:  # Kludge
        mask = np.logical_not(mask)
    else:
        pass
    
    # Apply Saiz-Lopez Marine MFT/MUT? <= should this be before multiply op.?
    if M_all:
        mask = mask*land_unmasked(res=res)

    if mask3D:
        if any( [ (mask.shape[-1] == i) for i in 38, 47 ] ):
            pass
        else: # concatenate dimensions
            mask = np.concatenate( [mask]*47, axis=2 )

    # Remove above the "chemical tropopause" from GEOS-Chem (v9-2)
    if trop_limit:
        mask = mask[...,:38] 

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
def get_EU_unmasked( res='1x1',  ):
    """ Mask 'EU' as defined by GEOS-Chem EU grid
        the grid of "'0.5x0.666" resolution is used by default, but
        any list of lat and lons could be provided and the extremities 
        would be used as the mask edges  """

    # Get GEOS-Chem EU lat and lons
    lon, lat, NIU = get_latlonalt4res( res='0.5x0.666' )

    # mask lats
    m1 = lat2lat_2D_unmasked( lowerlat=lat.min(), higherlat=lat.max(), res=res )
    
    # mask lons
    m2 = lon2lon_2D_unmasked(lowerlon=lon.min(), higherlon=lon.max(), res=res )

    #  conbine maskes 
    m = m1 + m2 
    
    return m