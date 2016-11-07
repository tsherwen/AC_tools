#!/usr/bin/python
""" Core functions used to be used by all function levels of 
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

# --------------                                                                                              
# 1.01 - Store of dirs for earth0, atmosviz1, and tms mac                                                     
# -------------                                                                                               
def get_dir( input, loc='earth0' ):
    """
    Retrieves directories within structure on a given platform
        ( e.g. York computer (earth0), Mac, Old cluster (atmosviz1... etc)  )
    NOTES:
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
        'fwd'  : home +'labbook/PhD_progs/d_fast-J_JX/data/',
        'lwd'  : home +'labbook/',
        'npwd' : home +'data/np_arrs/',
        'tpwd' : home +'labbook/PhD_progs/' ,
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
# 1.02 -  Get Latitude as GC grid box number in dimension                                                                                                                  
# ----                                                                                                                                                        
def get_gc_lat(lat, res='4x5', wd=None, debug=False):
    """ 
    Get index of lat for given resolution 
    """
    NIU, lat_c, NIU = get_latlonalt4res( res=res, wd=wd )
    del NIU
    return find_nearest_value( lat_c, lat )

# ----                                                                                                                                                        
# 1.03 -  Get Longitude as GC grid box number in dimension                                                                                                       
# ----                                                                                                                                                        
def get_gc_lon(lon, res='4x5', wd=None, debug=False):
    """ 
    Get index of lon for given resolution 
    """
    lon_c, NIU, NIU = get_latlonalt4res( res=res, wd=wd )
    del NIU
    return find_nearest_value( lon_c, lon )

# --------                                                                                          
# 1.04 - Get model array dimension for a given resolution                                           
# --------                                                                                          
def get_dims4res(res=None, r_dims=False, invert=True, trop_limit=False, \
        just2D=False, debug=False):
    """ Get dimension of GEOS-Chem output for given resolution """

    # Dictionary of max dimensions of standard GEOS-Chem output
    dims = {
    '4x5' :  (72,46,47), 
    '2x2.5':(144,91,47) , 
    '1x1' :  (360,181,47), 
    '0.5x0.5' :  (720,361,47), 
    '0.5x0.666':(121,81,47) ,
    '0.25x0.3125':(177, 115, 47),
    '0.5x0.625':(145,133,47),
    }
    if debug:
        print dims

    # Only consider the GEOS-Chem chemical troposphere
    if trop_limit:
        vals =[]
        for i in dims.values():
            vals += [ ( i[0],i[1],38) ]
        dims = dict( zip( dims.keys(), vals ) )
        if debug:
            print dims

    # Dictionary of lon, lat (e.g. for emissions and 2D datasets)
#    if just2D:
#        vals =[]
#        for i in dims.values():
#            vals += [ ( i[0],i[1]) ]
#        dims = dict( zip( dims.keys(), vals ) )
#        if debug:
#            print dims

    if just2D:
        _2Ddims = {}
        for res in dims.keys():
            _2Ddims[ res] = dims[res][0:2] 
        dims = _2Ddims

    if r_dims:
        if invert==True:
            return {v: k for k, v in dims.items()}
        else:
            return dims
    else:
        return dims[res ]

# ----                                                                                                                                                        
#  1.05 - Get grid values of lon, lat, and alt for a given resolution                                                                                                                         
# ----                                                                                                                                                        
def get_latlonalt4res( res=None, centre=True, hPa=False, nest=None, \
            dtype=None, wd=None, filename='ctm.nc',\
            lat_bounds=u'latitude_bnds', lon_bounds=u'longitude_bnds',\
            lon_var=u'longitude', lat_var=u'latitude', debug=False ):
    """ Return lon, lat, and alt for a given resolution. 
        This function uses an updated version of gchem's variable 
        dictionaries 
    
    NOTE:
        - This function replaces most use dictionaries from "gchemgrid"
        - The update to using ctm.nc files has cause an bug linked to the lat 
        and lon variabe retrival. Just update to passing a wd with output at the 
        correct resolution to fix this.
    """ 
    logging.info("Calling get_latlonalt4res")
#    logging.debug( locals() )


    if res==None:
        logging.warning("No resolution specified. Assuming 4x5!")
        res='4x5'


# Kludge. Update function to pass "wd" 
# if model output directory ("wd") not provided use default directory
    if wd == None:
        AC_tools_dir = os.path.dirname(__file__)
        dwd = os.path.join(AC_tools_dir, 'data/LM')
        dir_dict = {
        '4x5':'LANDMAP_LWI_ctm',  \
        '2x2.5': 'LANDMAP_ctm_2x25',  \
        '1x1' : 'work/data/GEOS/HEMCO/EMEP/v2015-03/',\
        # Kludge, use 1x1 for 0.5x0.5 <= remove this
        '0.5x0.5' :'work/data/GEOS/HEMCO/EMEP/v2015-03/',\
        '0.5x0.666' :'LANDMAP_LWI_ctm_05x0666',  \
        '0.25x0.3125' :'LANDMAP_LWI_ctm_025x03125',  \
        # Need to add a 0.5x0.625!
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


    # Ge the data file name
    data_fname = os.path.join(wd, filename)
    if not os.path.exists(data_fname):
        logging.error("Could not find {fn}".format(fn=data_fname))
        raise IOError, "Could not find {fn}".format(fn=data_fname)

#    if debug:
#        print res, centre, hPa, nest, type( res)
#        print data_fname, res

    if centre:
        # Extract lat and lon from model output data file
        with Dataset( data_fname, 'r' ) as d:
            lat = np.array( d[lat_var] )    
            lon = np.array( d[lon_var] )        
            

    # Get edge values
    if (not centre) and ( not any([(res==i) for i in '1x1', '0.5x0.5' ]) ):
        # Extract lat and lon from model output data file
        try:
            with Dataset( data_fname, 'r' ) as d:
                lat = np.array( d[lat_bounds] )    
                lon = np.array( d[lon_bounds] )  

                # select lower edge of each bound, and final upper edge
                lat = [i[0] for i in lat ]+[ lat[-1][1] ]
                lon = [i[0] for i in lon ]+[ lon[-1][1] ]            
                lat, lon = [np.array(i) for i in lat, lon ]
        except:
            error = "Could not get {lat}, {lon} from {fn}"\
                    .format(fn=data_fname,lat=lat_bounds,
                            lon=lon_bounds)
            logging.error(error)
            raise IOError, error

    # Kludge - mannually give values for 1.0x1.0
    if res=='1x1':
        if centre:
            lat = np.arange(-90, 90+1, 1)
            lon = np.arange(-180, 180, 1)     
        else:
            lat = np.array( [-90]+list(np.arange(-89-0.5, 90+.5, 1))+[90] )
            lon = np.arange(-180-0.5, 180+.5, 1)

    # Kludge - mannually give values for 0.5x0.5
    if res=='0.5x0.5':
        if centre:
            lat = np.array( [-90]+list(np.arange(-89, 90, .5))+[90] )
            lon = np.arange(-180, 180, .5)             
        else:
            lat = np.array( [-90]+list(np.arange(-89.75, 90+.25, .5))+[90] )
            lon = np.arange(-180-0.25, 180+.25, .5)        


    # Get dictionary variable name in Gerrit's GEOS-Chem dimensions list
    # ( now only doing this for alt, as alt values not in model output? )
    if hPa:
        alt = 'c_hPa_geos5_r'
    else:
        alt='c_km_geos5_r'
    d = gchemgrid(rtn_dict=True)
    alt = d[alt]

    # Also provide high resolution grid if requested from this function all
    if nest =='high res global':
        lon, lat = np.arange(-180, 180, 0.25), np.arange(-90, 90, 0.25)
        return lon, lat, alt

    if debug:
        print lon, lat, alt
    rtn_list =  lon, lat, alt 
    if not isinstance( dtype, type( None ) ):
        return [ i.astype( dtype ) for i in rtn_list ]
    else:
        return  rtn_list

# --------------
# 1.06 - Convert from hPa to km or vice versa.
# -------------
def hPa_to_Km(input, reverse=False, debug=False):
    """ hPa/km convertor
        Set "reverse" to True to convert Km to hPa 
    """
    if reverse:
         return [ np.exp(  np.float(i) /-7.6)*1013. for i in input ]
    else:
        return [-7.6*np.log( float(i) / 1013.) for i in input ]

# --------   
# 1.07 - Find nearest
# --------
def find_nearest_value( array, value ):
    """ 
    Find nearest point. 
    
    NOTEs:
     - Adapted from (credit:) HappyLeapSecond's Stackoverflow answer.  
    ( http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array )
    """
    idx = (np.abs(array-value)).argmin()
    return idx

# -------------
# 1.08 - Work out iGEOS-Chem version 
# ------------- 
def iGEOSChem_ver(wd, verbose=True, debug=False):
    """ 
    Get iGEOS-Chem verson 

    NOTES:
     - These are not GEOSChem versions, but halogen branch versions
    (e.g. iGeosChem 1.1 or 1.2 from dir name ( wd )  )
    """
    # List iGEOSChem versions+ then DataFrame
    versions = [
    '1.1','1.2', '1.3', '1.4', '1.5', '1.6', '1.6.1', '1.6.2', \
     '1.6.3', '1.7', '2.0', '3.0', '4.0'  ]
    df= DataFrame( versions, columns=['Versions'] )
    if debug:
        print wd, versions, df

    # Which versions are in name?
    def element_in_str( element ): 
       return (element in wd)
    df['Run Version'] = df['Versions'].apply( element_in_str )

    # Select last value and return as string
    return df['Versions'][ df['Run Version'] ][-1:].values[0]


# --------------                                                                                 
# 1.99 - Reference data (lon, lat, and alt) adapted from gchem - credit: GK (Gerrit Kuhlmann )             
# -------------                                                                                                                                      
def gchemgrid(input=None, rtn_dict=False, debug=False):
    """
    GeosChem grid lookup table.
    
    Can give the latitude, longitude or altitude positions of geoschem in units to 
    convert from model units. Returns a numpy array.

    Example:
    gchemgrid('c_lon_4x5')
    [-180, -175, -170 ,..., 175, 180]

    ARGUEMENTS:
     - rtn_dict=False : return the lookup dictionary instead of only a numpy array.
     - debug=False : Prints extra infromaiton in the function for debugging.
    NOTES:
     - If the array begins with a 'c' then this denotes the center position of the gridbox.
    'e' denotes the edge.
     -  The following arrays ara available:
    c_lon_4x5, e_lon_4x5, c_lat_4x5, e_lat_4x5, 
    e_eta_geos5_r, e_km_geos5_r, e_hpa_geos5_r, 
    c_eta_geos5_r, c_km_geos5_r, c_hpa_geos5_r, 
    """

    if ((input==None) and (rtn_dict==False)):
        raise KeyError('gchemgrid requires an input or rtn_dict=True')

    """
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

    d = {
    # 4x5                                                                                        
   'c_lon_4x5' : np.arange(-180, 175+5, 5) ,
   'e_lon_4x5' : np.arange(-182.5, 177.5+5, 5) ,
   'c_lat_4x5' : np.array( [-89]+ range(-86, 90, 4)+ [89] ),
   'e_lat_4x5' : np.array( [-90]+ range(-88, 92, 4)+ [90] ),

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

    #Grid box level edges [km]:                                                                      
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

    #Grid box level edges [hPa]:                                                                     
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

    #Grid box level centers (eta-coordinates)                                                        
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

    #Grid box level centers [km]                                                                     
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

    #Grid box level centers [hPa]                                                                    
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
    }
    if debug:
        print 'gchemgrid called'
    logging.debug('gchemgrid called')

    if rtn_dict:
        return d
    else:
        return d[input]
        
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
