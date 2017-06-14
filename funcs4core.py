#!/usr/bin/python
# -*- coding: utf-8 -*-
""" 
Core functions used to be used by all function levels of 
GEOS-Chem/Data analysis in AC_Tools 

NOTES: 
 - user specific functions need to be migrated to a seperate fuction
"""

# ----------------------------- Section 0 -----------------------------------
# -------------- Required modules:

import numpy as np
from netCDF4 import Dataset
from pandas import DataFrame
import platform
import sys
import logging
import os

from math import log10, floor
import math

# --------------                                                                                              
# X.XX - Store of dirs for earth0, atmosviz1, and tms mac                                                     
# -------------                                                                                               
def get_dir( input, loc='earth0' ):
    """
    Retrieves directories within structure on a given platform
        ( e.g. York computer (earth0), Mac, Old cluster (atmosviz1... etc)  )

    Notes
    -----
     - This function is not general enough to be transferrable.
     - Update to use $USER flag. 
    """
    import getpass
    import platform
    host = platform.node()
    user = getpass.getuser()
    tms_users = [ 'Tomas', 'tomassherwen', 'ts551' ]

    # Mac setup                                                                                               
    if (host == 'tomasmbp13.york.ac.uk') or ('Tomas-13-MBP') :
        home = '/Users/{}/'.format( user )
        if user in tms_users:
            d = { 
        'rwd'  : home+'PhD/Data/MUTD_iGEOS-Chem_output/',
        'dwd'  :  home+'PhD/Data/' ,
        'npwd' : home+'PhD/Data/np_arrs/' ,
        'tpwd' : home+'GITHub/PhD_progs/' ,
        'ppwd' : home+'Pictures/'  
            }
            d = d[input]
        else:
            d = home

    # Earth0 setup                                                                                            
    if host == 'earth0' :
        home =  '/work/home/{}/'.format( user )
        if user in tms_users:
            d = { 
        'rwd'  : home +'data/all_model_simulations/iodine_runs/',
        'dwd'  : home +'data/',
        'fwd'  : home +'labbook/Python_progs/d_fast-J_JX/data/',
        'lwd'  : home +'labbook/',
        'npwd' : home +'data/np_arrs/',
        'tpwd' : home +'labbook/Python_progs/' ,
        'ppwd' : home +'labbook/plots_images/'  
            }
            d = d[input]
        else:
            d = home

    # Atmosviz1 setup                                                                                         
    if host == 'atmosviz1' :
        home =  '/home/{}/'.format( user )
        if user in tms_users:
            d = { 
        'rwd'  : home +'data/model/',
        'dwd'  :  home +'data/',
        'lwd'  :  home +'labbook/',
        'npwd' : home +'data/np_arrs/',
        'tpwd' : home +'labbook/PhD_progs/'  
            }
            d = d[input]
        else:
            d = home

    return d
    
# ----                                                                                                                                                        
# X.XX -  Get Latitude as GC grid box number in dimension                                                                                                                  
# ----                                                                                                                                                        
def get_gc_lat(lat, res='4x5', wd=None, filename='ctm.nc', debug=False):
    """ 
    Get index of lat for given resolution 

    Parameters
    ----------
    lat (float): latitude to convert
    wd (str): the directory to search for file in
    filename (Str): name of NetCDF file (e.g. ctm.nc or ts_ctm.nc)
    res (str): the resolution if wd not given (e.g. '4x5' )
    debug (boolean): legacy debug option, replaced by python logging

    Returns
    -------
    (float)
    """
    NIU, lat_c, NIU = get_latlonalt4res( res=res, wd=wd, filename=filename )
    del NIU
    return find_nearest_value( lat_c, lat )

# ----                                                                                                                                                        
# X.XX -  Get Longitude as GC grid box number in dimension                                                                                                       
# ----                                                                                                                                                        
def get_gc_lon(lon, res='4x5', wd=None, filename='ctm.nc', debug=False):
    """ 
    Get index of lon for given resolution 

    Parameters
    ----------
    lon (float): longitude to convert
    wd (str): the directory to search for file in
    filename (Str): name of NetCDF file (e.g. ctm.nc or ts_ctm.nc)
    res (str): the resolution if wd not given (e.g. '4x5' )
    debug (boolean): legacy debug option, replaced by python logging

    Returns
    -------
    (float)
    """
    lon_c, NIU, NIU = get_latlonalt4res( res=res, wd=wd, filename=filename )
    del NIU
    return find_nearest_value( lon_c, lon )

# --------                                                                                          
# X.XX - Get model array dimension for a given resolution                                           
# --------                                                                                          
def get_dims4res(res=None, r_dims=False, invert=True, trop_limit=False, \
        just2D=False, full_vertical_grid=False, add_n_time_dims=None, 
        debug=False):
    """ 
    Get dimension of GEOS-Chem output for given resolution 

    Parameters
    ----------
    invert (boolean): invert dictionary keys and values
    trop_limit (boolean): limit 4D arrays to troposphere     
    r_dims (boolean): return dicionary of dimensions
    just2D (boolean): just return horizontal dimensions 
    debug (boolean): legacy debug option, replaced by python logging
    add_n_time_dims (int): number of time dimensions to add

    Returns
    -------
    (tuple)
    """
    # Dictionary of max dimensions of standard GEOS-Chem output
    dims = {
    '4x5'            : (72,46,47), 
    '2x2.5'          : (144,91,47), 
    '1x1'            : (360,181,47), 
    '0.5x0.5'        : (720,361,47), 
    '0.5x0.666'      : (121,81,47),
    '0.25x0.3125'    : (177,115,47), # tms - update to be '0.25x0.3125_EU' for consistancy?
    '0.25x0.3125_CH' : (225,161,47), 
    '0.25x0.3125_WA' : (145,89,47), 
    '0.5x0.625'      : (145,133,47),
    '0.083x0.083'    : (4320, 2160), # 9km resolution?
    }
    if debug:
        print(dims)

    # If full using full vertical 
    if full_vertical_grid:
        vals =[]
        for i in list(dims.values()):
            vals += [ ( i[0],i[1],72) ]
        dims = dict( list(zip( list(dims.keys()), vals )) )
        if debug:
            print(dims)

    # Only consider the GEOS-Chem chemical troposphere
    if trop_limit:
        vals =[]
        for i in list(dims.values()):
            vals += [ ( i[0],i[1],38) ]
        dims = dict( list(zip( list(dims.keys()), vals )) )
        if debug:
            print(dims)

    # Add a time dimension of length n (if add_n_time is an integer)
    if not isinstance( add_n_time_dims, type(None) ):
        vals =[]
        for i in list(dims.values()):
            vals += [ ( i[0],i[1],i[2], add_n_time_dims) ]
        dims = dict( list(zip( list(dims.keys()), vals )) ) 

    # Dictionary of lon, lat (e.g. for emissions and 2D datasets)
    if just2D:
        vals =[]
        for i in list(dims.values()):
            vals += [ ( i[0],i[1]) ]
        dims = dict( list(zip( list(dims.keys()), vals )) )
#        if debug:
#            print dims

#    if just2D:
#        _2Ddims = {}
#        for res in dims.keys():
#            _2Ddims[res] = dims[res][0,1] 
#        dims = _2Ddims
       


    if r_dims:
        if invert==True:
            return {v: k for k, v in list(dims.items())}
        else:
            return dims
    else:
        return dims[res ]

# ----                                                                                                                                                        
#  X.XX - Get grid values of lon, lat, and alt for a given resolution                                                                                                                         
# ----                                                                                                                                                        
def get_latlonalt4res( res=None, centre=True, hPa=False, nest=None, \
        dtype=None, wd=None, filename='ctm.nc', full_vertical_grid=False, \
        lat_bounds='latitude_bnds', lon_bounds='longitude_bnds',\
        lon_var='longitude', lat_var='latitude', \
#        lon_var=u'lon', lat_var=u'lat', 
        verbose=True, debug=False ):
    """ 
    Get lon, lat, and alt for a given model resolution. 

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    res (str): the resolution if wd not given (e.g. '4x5' )
    debug (boolean): legacy debug option, replaced by python logging
    lon_var, lat_var (str): variables names for lon and lat in the NetCDF 
    lon_bounds, lat_bounds (str): variables names for lon and lat bounds in the NetCDF 
    filename (str): name of NetCDF to use
    dtype (type): type for which data is return as, e.g. np.float64
    nest (str): manual override for retruned variables - vestigle?
    hPa (boolean): return altitudes in units of hPa, instead of km
    full_vertical_grid (boolean): use full vertical grid or reduced (47 vs. 72)

    Returns
    -------
    (list) variables for lon, lat, alt as arrays

    Notes
    -----
     - This function uses an updated version of gchem's variable dictionaries 
     - This function replaces most use dictionaries from "gchemgrid"
     - The update to using ctm.nc files has cause an bug linked to the lat 
     and lon variabe retrival. Just update to passing a wd with output at the 
     correct resolution to fix this.
    """ 
    logging.info("Calling get_latlonalt4res for res={}".format(res) )

    if isinstance( res, type(None) ):
        logging.warning("No resolution specified. Assuming 4x5!")
        res='4x5'

    if isinstance( wd, type(None) ):
        AC_tools_dir = os.path.dirname(__file__)
        dwd = os.path.join(AC_tools_dir, 'data/LM')
        dir_dict = {
        '4x5':'LANDMAP_LWI_ctm',  \
        '2x2.5': 'LANDMAP_LWI_ctm_2x25',  \
        '1x1' : 'work/data/GEOS/HEMCO/EMEP/v2015-03/',\
        # Kludge, use 1x1 for 0.5x0.5 <= remove this
        '0.5x0.5' :'work/data/GEOS/HEMCO/EMEP/v2015-03/',\
        '0.5x0.666' :'LANDMAP_LWI_ctm_05x0666',  \
        '0.25x0.3125' :'LANDMAP_LWI_ctm_025x03125',  \
        '0.25x0.3125_CH' :'LANDMAP_LWI_ctm_025x03125_CH',  \
        '0.25x0.3125_WA' :'LANDMAP_LWI_ctm_025x03125_WA',  \
        # Need to add a 0.5x0.625!
        # Temporary inclusion of local 0.083x0.083 file. 
        '0.083x0.083' : 'LANDMAP_LWI_ctm_0083x0083'
        }
        try:
            dir = dir_dict[res]
        except KeyError:
            logging.error("{res} not a recognised resolution!".format(res=res))
            raise KeyError
        wd = os.path.join(dwd, dir)
    
        if (res=='1x1') or (res=='0.5x0.5'):
            filename='EMEP.geos.1x1.nc'
            lat_var = 'lat'
            lon_var = 'lon'
            wd = '/work/data/GEOS/HEMCO/EMEP/v2015-03/'

    # Get the data file name
    data_fname = os.path.join(wd, filename)
    if not os.path.exists(data_fname):
        logging.error("Could not find {fn}".format(fn=data_fname))
        raise IOError("Could not find {fn}".format(fn=data_fname))

    if centre:
        try:
            # Extract lat and lon from model output data file
            with Dataset( data_fname, 'r' ) as d:
                lat = np.array( d[lat_var] )    
                lon = np.array( d[lon_var] )        
        except IOError:
            error = "Could not get {lat}, {lon} from {fn}"\
                    .format(fn=data_fname,lat=lat_bounds,
                            lon=lon_bounds)
            logging.error(error)
            if verbose:
                print("ERROR: are the refernces files in 'AC_tools/data/LM' ?")
                print("(To download just run AC_tools/Scripts/get_data_files.py)")
            raise IOError(error)        

    # Get edge values
    exception_res = ('1x1', '0.5x0.5', '0.083x0.083')
    if (not centre) and (res not in exception_res):
        # Extract lat and lon from model output data file
        try:
            with Dataset( data_fname, 'r' ) as d:
                lat = np.array( d[lat_bounds] )    
                lon = np.array( d[lon_bounds] )  

                # select lower edge of each bound, and final upper edge
                lat = [i[0] for i in lat ]+[ lat[-1][1] ]
                lon = [i[0] for i in lon ]+[ lon[-1][1] ]            
                lat, lon = [np.array(i) for i in (lat, lon) ]
        except:
            error = "Could not get {lat}, {lon} from {fn}"\
                    .format(fn=data_fname,lat=lat_bounds,
                            lon=lon_bounds)
            if verbose:
                print("ERROR: are the refernces files in 'AC_tools/data/LM' ?")
                print("(To download just run AC_tools/Scripts/get_data_files.py)")
            logging.error(error)
            raise IOError(error)

    # Manually set values for 0.25x0.3125_CH
#    if res=='0.25x0.3125_CH':
#        if centre:
#        lat = np.arange(15., 55.+.25, .25)
#        lon = np.arange(70, 140.+.3125, .3125)   
#        else:
#            lat = np.array( [-90]+list(np.arange(-89-0.5, 90+.5, 1))+[90] )
#            lon = np.arange(-180-0.5, 180+.5, 1)

    # Manually set values for 1.0x1.0
    if res=='1x1':
        step_size = 1.0
        if centre:
            lat = np.arange(-90, 90+step_size, step_size)
            lon = np.arange(-180, 180, step_size)     
        else:
            lat = np.array([-90]+list(np.arange(-89-0.5, 90+.5, step_size))+[90])
            lon = np.arange(-180-(step_size/2), 180+(step_size/2), step_size)

    # Manually set values for 0.5x0.5
    if res=='0.5x0.5':
        step_size = 0.5
        if centre:
            lat = np.array( [-90]+list(np.arange(-89, 90, step_size))+[90] )
            lon = np.arange(-180, 180, step_size)
        else:
            lat = np.array([-90]+list(np.arange(-89.75, 90+.25, step_size))+[90])
            lon = np.arange(-180-(step_size/2), 180+(step_size/2), step_size)
    # Manually set values for 0.1x0.1

    # Manually set values for 0.083x0.083
    if res=='0.083x0.083':
        step_size = 0.083333336

        if centre:
            lat = np.arange( -89.95833588, 89.95833588+step_size, step_size )
            lon = np.arange(-179.95832825, 179.95835876, step_size )        
            # adjust to center point
            lat = [i+step_size/2 for i in lat[:-1] ]
            lon = [i+step_size/2 for i in lon[:-1] ]
            
        else:
            lat = np.arange( -89.95833588, 89.95833588+step_size, step_size)
            lon = np.arange(-179.95832825, 179.95835876, step_size)        

    # Get dictionary variable name in Gerrit's GEOS-Chem dimensions list
    # ( now only doing this for alt, as alt values not in model output? )
    if hPa:
        alt = 'c_hPa_geos5'
    else:
        alt = 'c_km_geos5'
    # Use reduced vertical grid? (then add '_r')
    if not full_vertical_grid:
        alt+= '_r'
    d = gchemgrid(rtn_dict=True)
    alt = d[alt]

    # Also provide high resolution grid if requested from this function all
    if nest =='high res global':
        lon, lat = np.arange(-180, 180, 0.25), np.arange(-90, 90, 0.25)
        return lon, lat, alt

    if debug:
        print(lon, lat, alt)
    rtn_list = lon, lat, alt 
    if not isinstance( dtype, type( None ) ):
        return [ i.astype( dtype ) for i in rtn_list ]
    else:
        return rtn_list

# --------------
# X.XX - Convert from hPa to km or vice versa.
# -------------
def hPa_to_Km(input, reverse=False, debug=False):
    """ 
    hPa/km convertor

    Parameters
    ----------
    input (list): list of values (float) to convert
    reverse (boolean): Set "reverse" to True to convert Km to hPa 
    debug (boolean): legacy debug option, replaced by python logging

    Returns
    -------
    (list)
    """
    if reverse:
         return [ np.exp(  np.float(i) /-7.6)*1013. for i in input ]
    else:
        return [-7.6*np.log( float(i) / 1013.) for i in input ]

# --------   
# X.XX - Find nearest
# --------
def find_nearest_value( array, value ):
    """ 
    Find nearest point. 

    Parameters
    ----------
    arrary (np.array): 1D array in which to search for nearest value
    value (float): value to search array for closest point
    """
    # Adapted from (credit:) HappyLeapSecond's Stackoverflow answer. 
    # ( http://stackoverflow.com/questions/2566412/  )
    idx = (np.abs(array-value)).argmin()
    return idx

# -------------
# X.XX - Work out iGEOS-Chem version from working directory name
# ------------- 
def iGEOSChem_ver(wd, also_return_GC_version=False, verbose=True, debug=False):
    """ 
    Get iGEOS-Chem verson 

    NOTES:
     - These are not GEOSChem versions, but halogen branch versions
    (e.g. iGeosChem 1.1 or 1.2 from dir name ( wd )  )
    """
    # List iGEOSChem versions+ then DataFrame
    versions = [
    '1.1','1.2', '1.3', '1.4', '1.5', '1.6', '1.6.1', '1.6.2', \
     '1.6.3', '1.7', '2.0', '3.0', '4.0', '5.0', '6.0', 
    # Also hold a 'NOT FOUND' value to for ease of processing non-halogen code
     'NOT FOUND'  
    ]
    df = DataFrame( versions, columns=['Versions'] )
    if debug:
        print(wd, versions, df)

    # Which versions are in name?
    def element_in_str( element ): 
       return (element in wd)
    df['Run Version'] = df['Versions'].apply( element_in_str )

    # Select last value and return as string
    try:
        iGC_ver = df['Versions'][ df['Run Version'] ][-1:].values[0]
    except IndexError:
        err_msq='WARNING (i) GEOS-Chem ver. number not found in working dir. path'
        print('!'*15, err_msq, '!'*15)
        logging.debug(err_msq)
        iGC_ver='NOT FOUND'
#        sys.exit()

    if also_return_GC_version:
        # list GEOS-Chem versions (written with dashes and underscores)
        versions = [ 
        'v11-01', 'v11_01', 'v10-01', 'v10_01', 'v9-02', 'v9_02', 'v9-01-03', 
        'v9_01_03', 'v9-01-02', 'v9_01_02', 'v9-01-01', 'v9_01_01', 'v8-03-02',
        'v8_03_02', 'v8-03-01', 'v8_03_01', 'v8-02-04', 'v8_02_04', 'v8-02-03', 
        'v8_02_03', 'v8-02-02', 'v8_02_02', 'v8-02-01', 'v8_02_01', 'v8-01-04', 
        'v8_01_04', 'v8-01-03', 'v8_01_03', 'v8-01-02', 'v8_01_02', 'v8-01-01', 
        'v8_01_01', 'v7-04-13', 'v7_04_13', 'v7-04-12', 'v7_04_12'
        ]
        df = DataFrame( versions, columns=['Versions'] )
        if debug:
            print(wd, versions, df)
        #
        df['Run Version'] = df['Versions'].apply( element_in_str )
        # selection 
        try:
            GC_ver = df['Versions'][ df['Run Version'] ][-1:].values[0]
        except IndexError:    
            # map iGEOS-Chem versions to GEOS-Chem versions 
            dict_iGC_GC = {
        '1.4':'v9-2',
        '1.5':'v9-2',
        '1.6':'v9-2',
        '1.7':'v9-2',
        '1.1':'v9-2',
        '1.2':'v9-2',
        '1.3':'v9-2',
        '3.0':'v10-01',
        '6.0':'v11-01',
        '5.0':'v11-01',
        '4.0':'v10-01',
        '1.6.2':'v9-2',
        '1.6.3':'v9-2',
        '1.6.1':'v9-2',
        'NOT FOUND':'NOT FOUND', 
            }
            # 
            GC_ver = dict_iGC_GC[iGC_ver]
        return iGC_ver, GC_ver

    return iGC_ver


# --------------                                                                                 
# X.XX - Reference data (lon, lat, and alt) adapted from gchem - credit: GK (Gerrit Kuhlmann )             
# -------------                                                                                                                                      
def gchemgrid(input=None, rtn_dict=False, debug=False):
    """
    GeosChem grid lookup table. Can give the latitude, longitude or altitude 
    positions of geoschem in units to convert from model units. 

    Parameters
    ----------
    rtn_dict (boolean): return the lookup dictionary instead
    debug (boolean): legacy debug option, replaced by python logging

    Returns
    -------
    (np.array)

    Notes
    -----
    - Detail from original reposotiroy adapted from:
    Updated to dictionary from gchemgrid (credit: Gerrit Kuhlmann ) with 
    additional grid adds 

    This [function] contains (some) grid coordinates used within GEOS-Chem 
    as numpy arrays

    Distribution advice from gchem (credit: Gerrit Kuhlmann):

    Python Script Collection for GEOS-Chem Chemistry Transport Model (gchem)
    Copyright (C) 2012 Gerrit Kuhlmann

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.    
    """
    logging.info('gchemgrid called for:{} (rtn_dict={})'.format(input,rtn_dict))

    if isinstance(input, type(None)) and (rtn_dict==False):
        raise KeyError('gchemgrid requires an input or rtn_dict=True')
    d = {
    # 4x5                                                                                        
   'c_lon_4x5' : np.arange(-180, 175+5, 5) ,
   'e_lon_4x5' : np.arange(-182.5, 177.5+5, 5) ,
   'c_lat_4x5' : np.array( [-89]+ list(range(-86, 90, 4))+ [89] ),
   'e_lat_4x5' : np.array( [-90]+ list(range(-88, 92, 4))+ [90] ),

    # Altitude - Grid box level edges (eta coordinate):                                          
    'e_eta_geos5_r' : np.array([
            1.00179600e+00,   9.86769000e-01,   9.71665000e-01,
            9.56562000e-01,   9.41459000e-01,   9.26356000e-01,
            9.11253000e-01,   8.96152000e-01,   8.81051000e-01,
            8.65949000e-01,   8.50848000e-01,   8.35748000e-01,
            8.20648000e-01,   8.00515000e-01,   7.75350000e-01,
            7.50186000e-01,   7.25026000e-01,   6.99867000e-01,
            6.74708000e-01,   6.36974000e-01,   5.99251000e-01,
            5.61527000e-01,   5.23819000e-01,   4.86118000e-01,
            4.48431000e-01,   4.10759000e-01,   3.73114000e-01,
            3.35486000e-01,   2.85974000e-01,   2.42774000e-01,
            2.06167000e-01,   1.75170000e-01,   1.48896000e-01,
            1.26563000e-01,   1.07578000e-01,   9.14420000e-02,
            7.77260000e-02,   5.58200000e-02,   3.97680000e-02,
            2.80770000e-02,   1.95860000e-02,   9.19100000e-03,
            4.02600000e-03,   1.62500000e-03,   6.01000000e-04,
            1.99000000e-04,   5.50000000e-05,   0.00000000e+00]) ,

    # Grid box level edges [km]:                                                                      
    'e_km_geos5_r' : np.array([
            6.00000000e-03,   1.35000000e-01,   2.66000000e-01,
            3.99000000e-01,   5.33000000e-01,   6.69000000e-01,
            8.06000000e-01,   9.45000000e-01,   1.08600000e+00,
            1.22900000e+00,   1.37400000e+00,   1.52000000e+00,
            1.66900000e+00,   1.87100000e+00,   2.12800000e+00,
            2.39200000e+00,   2.66300000e+00,   2.94100000e+00,
            3.22800000e+00,   3.67300000e+00,   4.14000000e+00,
            4.63100000e+00,   5.14900000e+00,   5.69800000e+00,
            6.28300000e+00,   6.91000000e+00,   7.58700000e+00,
            8.32400000e+00,   9.41100000e+00,   1.05050000e+01,
            1.15780000e+01,   1.26330000e+01,   1.36740000e+01,
            1.47060000e+01,   1.57310000e+01,   1.67530000e+01,
            1.77730000e+01,   1.98550000e+01,   2.20040000e+01,
            2.42400000e+01,   2.65960000e+01,   3.17160000e+01,
            3.75740000e+01,   4.42860000e+01,   5.17880000e+01,
            5.99260000e+01,   6.83920000e+01,   8.05810000e+01]) ,

    # Grid box level edges [hPa]:                                                                     
    'e_hPa_geos5_r' : np.array([
            1.01181400e+03,   9.96636000e+02,   9.81382000e+02,
            9.66128000e+02,   9.50874000e+02,   9.35621000e+02,
            9.20367000e+02,   9.05114000e+02,   8.89862000e+02,
            8.74610000e+02,   8.59358000e+02,   8.44107000e+02,
            8.28856000e+02,   8.08522000e+02,   7.83106000e+02,
            7.57690000e+02,   7.32279000e+02,   7.06869000e+02,
            6.81458000e+02,   6.43348000e+02,   6.05247000e+02,
            5.67147000e+02,   5.29062000e+02,   4.90984000e+02,
            4.52921000e+02,   4.14873000e+02,   3.76851000e+02,
            3.38848000e+02,   2.88841000e+02,   2.45210000e+02,
            2.08236000e+02,   1.76930000e+02,   1.50393000e+02,
            1.27837000e+02,   1.08663000e+02,   9.23660000e+01,
            7.85120000e+01,   5.63880000e+01,   4.01750000e+01,
            2.83680000e+01,   1.97920000e+01,   9.29300000e+00,
            4.07700000e+00,   1.65100000e+00,   6.17000000e-01,
            2.11000000e-01,   6.60000000e-02,   1.00000000e-02]) ,

    # Grid box level centers (eta-coordinates)                                                        
    'c_eta_geos5_r' : np.array([
            9.94283000e-01,   9.79217000e-01,   9.64113000e-01,
            9.49010000e-01,   9.33908000e-01,   9.18805000e-01,
            9.03703000e-01,   8.88601000e-01,   8.73500000e-01,
            8.58399000e-01,   8.43298000e-01,   8.28198000e-01,
            8.10582000e-01,   7.87933000e-01,   7.62768000e-01,
            7.37606000e-01,   7.12447000e-01,   6.87287000e-01,
            6.55841000e-01,   6.18113000e-01,   5.80389000e-01,
            5.42673000e-01,   5.04968000e-01,   4.67274000e-01,
            4.29595000e-01,   3.91937000e-01,   3.54300000e-01,
            3.10730000e-01,   2.64374000e-01,   2.24471000e-01,
            1.90668000e-01,   1.62033000e-01,   1.37729000e-01,
            1.17070000e-01,   9.95100000e-02,   8.45840000e-02,
            6.67730000e-02,   4.77940000e-02,   3.39230000e-02,
            2.38320000e-02,   1.43890000e-02,   6.60900000e-03,
            2.82500000e-03,   1.11300000e-03,   4.00000000e-04,
            1.27000000e-04,   2.80000000e-05]) ,

    # Grid box level centers [km]                                                                     
    'c_km_geos5_r' : np.array([
            7.10000000e-02,   2.01000000e-01,   3.32000000e-01,
            4.66000000e-01,   6.01000000e-01,   7.37000000e-01,
            8.75000000e-01,   1.01600000e+00,   1.15700000e+00,
            1.30100000e+00,   1.44700000e+00,   1.59400000e+00,
            1.76900000e+00,   1.99900000e+00,   2.25900000e+00,
            2.52700000e+00,   2.80100000e+00,   3.08400000e+00,
            3.44800000e+00,   3.90400000e+00,   4.38200000e+00,
            4.88600000e+00,   5.41900000e+00,   5.98500000e+00,
            6.59100000e+00,   7.24100000e+00,   7.94700000e+00,
            8.84800000e+00,   9.93800000e+00,   1.10210000e+01,
            1.20860000e+01,   1.31340000e+01,   1.41700000e+01,
            1.51980000e+01,   1.62220000e+01,   1.72430000e+01,
            1.87270000e+01,   2.08360000e+01,   2.30200000e+01,
            2.53070000e+01,   2.86540000e+01,   3.40240000e+01,
            4.01660000e+01,   4.71350000e+01,   5.48340000e+01,
            6.30540000e+01,   7.21800000e+01]) ,

    # Grid box level centers [hPa]                                                                    
    'c_hPa_geos5_r' : np.array([
            1.00422500e+03,   9.89009000e+02,   9.73755000e+02,
            9.58501000e+02,   9.43247000e+02,   9.27994000e+02,
            9.12741000e+02,   8.97488000e+02,   8.82236000e+02,
            8.66984000e+02,   8.51732000e+02,   8.36481000e+02,
            8.18689000e+02,   7.95814000e+02,   7.70398000e+02,
            7.44984000e+02,   7.19574000e+02,   6.94163000e+02,
            6.62403000e+02,   6.24298000e+02,   5.86197000e+02,
            5.48105000e+02,   5.10023000e+02,   4.71952000e+02,
            4.33897000e+02,   3.95862000e+02,   3.57850000e+02,
            3.13844000e+02,   2.67025000e+02,   2.26723000e+02,
            1.92583000e+02,   1.63661000e+02,   1.39115000e+02,
            1.18250000e+02,   1.00514000e+02,   8.54390000e+01,
            6.74500000e+01,   4.82820000e+01,   3.42720000e+01,
            2.40800000e+01,   1.45420000e+01,   6.68500000e+00,
            2.86400000e+00,   1.13400000e+00,   4.14000000e-01,
            1.39000000e-01,   3.80000000e-02]) ,

    # GEOS-5 Native Vertical Grid (72 hybrid pressure-sigma levels)
    # http://acmg.seas.harvard.edu/geos/doc/archive/man.v9-01-02/appendix_3.html
    'c_hPa_geos5' : np.array([  
         1.00000000e-02,   1.50000000e-02,   2.00000000e-02,
         2.60000000e-02,   3.30000000e-02,   4.00000000e-02,
         4.80000000e-02,   5.70000000e-02,   6.60000000e-02,
         7.80000000e-02,   8.90000000e-02,   1.05000000e-01,
         1.20000000e-01,   1.40000000e-01,   1.59000000e-01,
         1.85000000e-01,   2.11000000e-01,   2.45000000e-01,
         2.79000000e-01,   3.22000000e-01,   3.65000000e-01,
         4.20000000e-01,   4.76000000e-01,   5.46000000e-01,
         6.17000000e-01,   7.06000000e-01,   7.95000000e-01,
         9.07000000e-01,   1.01900000e+00,   1.16000000e+00,
         1.30100000e+00,   1.47600000e+00,   1.65100000e+00,
         1.86800000e+00,   2.08500000e+00,   2.35300000e+00,
         2.62000000e+00,   2.94800000e+00,   3.27600000e+00,
         3.67700000e+00,   4.07700000e+00,   4.56200000e+00,
         5.04700000e+00,   5.63200000e+00,   6.21700000e+00,
         6.91800000e+00,   7.62000000e+00,   8.45600000e+00,
         9.29300000e+00,   1.02850000e+01,   1.12770000e+01,
         1.24600000e+01,   1.36430000e+01,   1.50500000e+01,
         1.64570000e+01,   1.81240000e+01,   1.97920000e+01,
         2.17610000e+01,   2.37300000e+01,   2.60490000e+01,
         2.83680000e+01,   3.10890000e+01,   3.38100000e+01,
         3.69930000e+01,   4.01750000e+01,   4.39100000e+01,
         4.76440000e+01,   5.20160000e+01,   5.63880000e+01,
         6.14960000e+01,   6.66030000e+01,   7.25580000e+01,
         7.85120000e+01,   8.54390000e+01,   9.23660000e+01,
         1.00514000e+02,   1.08663000e+02,   1.18250000e+02,
         1.27837000e+02,   1.39115000e+02,   1.50393000e+02,
         1.63661000e+02,   1.76930000e+02,   1.92587000e+02,
         2.08244000e+02,   2.26745000e+02,   2.45246000e+02,
         2.67087000e+02,   2.88927000e+02,   3.13966000e+02,
         3.39005000e+02,   3.58038000e+02,   3.77070000e+02,
         3.96112000e+02,   4.15155000e+02,   4.34212000e+02,
         4.53269000e+02,   4.72335000e+02,   4.91401000e+02,
         5.10475000e+02,   5.29550000e+02,   5.48628000e+02,
         5.67706000e+02,   5.86793000e+02,   6.05880000e+02,
         6.24967000e+02,   6.44054000e+02,   6.63146000e+02,
         6.82239000e+02,   6.94969000e+02,   7.07699000e+02,
         7.20429000e+02,   7.33160000e+02,   7.45890000e+02,
         7.58621000e+02,   7.71354000e+02,   7.84088000e+02,
         7.96822000e+02,   8.09556000e+02,   8.19743000e+02,
         8.29929000e+02,   8.37570000e+02,   8.45211000e+02,
         8.52852000e+02,   8.60493000e+02,   8.68135000e+02,
         8.75776000e+02,   8.83418000e+02,   8.91059000e+02,
         8.98701000e+02,   9.06342000e+02,   9.13984000e+02,
         9.21626000e+02,   9.29268000e+02,   9.36911000e+02,
         9.44553000e+02,   9.52195000e+02,   9.59837000e+02,
         9.67480000e+02,   9.75122000e+02,   9.82765000e+02,
         9.90408000e+02,   9.98051000e+02,   1.00565000e+03,
         1.01325000e+03]),
    # GEOS-5 Native Vertical Grid (72 hybrid pressure-sigma levels)
    'c_km_geos5' : np.array([  
        0.0905,   0.2215,   0.3535,   0.4875,   0.623 ,   0.7605,
         0.899 ,   1.0395,   1.182 ,   1.3265,   1.473 ,   1.6215,
         1.8095,   2.053 ,   2.3155,   2.5855,   2.862 ,   3.1465,
         3.552 ,   4.014 ,   4.499 ,   5.0105,   5.5525,   6.1285,
         6.745 ,   7.4095,   8.1315,   9.1275,  10.22  ,  11.2995,
        12.3595,  13.404 ,  14.438 ,  15.4645,  16.4875,  17.508 ,
        18.538 ,  19.582 ,  20.642 ,  21.721 ,  22.8195,  23.944 ,
        25.098 ,  26.2835,  27.502 ,  28.754 ,  30.0415,  31.3655,
        32.7365,  34.1605,  35.637 ,  37.1665,  38.7505,  40.388 ,
        42.078 ,  43.8205,  45.613 ,  47.454 ,  49.3395,  51.271 ,
        53.2445,  55.2555,  57.299 ,  59.37  ,  61.4615,  63.567 ,
        65.68  ,  67.8175,  70.0485,  72.496 ,  75.4755,  79.3635]),
    }
    # 
    Temp_arrays = ['c_km_geos5', 'c_hPa_geos5']
    if input in Temp_arrays:
        print('WARNING array ({}) is temporary! - Please check it!'.format(input))

    if rtn_dict:
        return d
    else:
        return d[input]


# --------------                                                                                 
# X.XX - Get significant figures 
# ------------- 
def get_sigfig( x, p=3 ):
    """
    Return a number with only the significant figures required.

    Parameters
    -------
    x (float): number to convert 
    p (int): The number of sig figs you want returned. Default=3

    Returns
    -------
    number with only significant figures.
    
    Notes
    -----    
     - Is this function outdated by get_scientific_number()?
    """

#####################
#    Found at https://github.com/randlet/to-precision/blob/master/to_precision.py
#
#    returns a string representation of x formatted with a precision of p
#    Based on the webkit javascript implementation taken from here:
#    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
##########

    import math
    x = float(x)

    if x == 0.:
        return float("0." + "0"*(p-1))

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1


    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return float("".join(out))

#    output round(x, -int(floor(log10(abs(x)))))
#    return output

# --------------                                                                                 
# X.XX - Get number in scientific form 
# ------------- 
def get_scientific_number( number, precision, string=False ): 
    """
    Gets a number in scientific notation with a given precision.
    Returns a rounded number by default, or can be returned as a string
    Recomended for plotting.
    Inputs:
    number (float) (number you want to change)
    precision (Integer) (How many significant figures you want)
    String = True (Boolian) Do you want the output returned as a string?
    Output: float(default) OR string(if string==True)
    """
    number = float(number)
    # Special case for 0 
    if number == 0.:
        if not string:
            return float("0." + "0"*(precision-1))
        else:
            return ("0." + "0"*(precision-1))

    # If negative prepend with - and use the absolute
    if number < 0:
        sign = "-"
#        precision = precision-1
        number = -number
    else:
        sign = ""

    # Get the exponent
    exponent = int(math.log10(number))
    mantisa = number / math.pow(10, exponent)

    # Fix for leading_zeroes
    if mantisa < 1:
        mantisa = mantisa * 10
        exponent = exponent - 1
 
     
#    if exponent < 0:                                 
#        precision = precision+1                      
    # Get only the sigfigs from the mantisa          
    mantisa = round(mantisa, precision)              
                                                     
    # Fix for leading 10                             
    if mantisa >= 10:                                
        print("hello from 10")                        
        mantisa = mantisa/10.0                       
        exponent = exponent+1                        
                                                     
    # Fix for mantisa=1                              
    if mantisa == 1.0:                               
        print("mantisa = 1.0")                        
        mantisa = "1." + "0"*(precision-1)           
                                                   

                                                                                
    # Create the string for the result.                                         
    out = sign + str(mantisa) 
    if not exponent==0:
        out = out + 'E' + str(exponent)                             
                                                                                
    # Return the result as a float unless asking for a string.                  
    if not string:                                                              
        return float(out)                                                       
    else:                      
	    return out

        
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
