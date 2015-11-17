""" Core functions used to be used by all function levels of 
    GEOS-Chem/Data analysis in AC_Tools """

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 0 ----- Modules required

import numpy as np
from netCDF4 import Dataset
import platform

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 1 ----- Modules required
# 1.01 - Get file directory
# 1.02 - Get Latitude as GC grid box number 
# 1.03 - Get Longitude as GC grid box number in dimension
# 1.04 - Get model array dimension for a given resolution
# 1.05 - Get grid values of lon, lat, and alt for a given resolution 
# 1.06 - Convert from hPa to km or vice versa.
# 1.07 - find nearest
# 1.99 - Get Reference data (lon, lat, and alt)  for a given resolution/nest

# --------------                                                                                              
# 1.01 - Store of dirs for earth0, atmosviz1, and tms mac                                                     
# -------------                                                                                               
def get_dir( input, loc='earth0' ):
    """
        Retrieves directories within structure on a given platform
        ( e.g. York cluster, Mac, Old cluster  )
    """

    try:
        case = { 
        'Darwin-14.0.0-x86_64-i386-64bit' :1,
        'Darwin-14.1.0-x86_64-i386-64bit':1 ,
        'Darwin-14.3.0-x86_64-i386-64bit': 1, # updated ...
        'Darwin-14.5.0-x86_64-i386-64bit': 1, # updated 15/08
#        'Linux-3.0.101-0.46-default-x86_64-with-SuSE-13.1-x86_64':2,                                 
        'Linux-3.0.101-0.47.52-default-x86_64-with-SuSE-11-x86_64' :2, 
        # reverted on 15 03 10        
        'Linux-3.0.101-0.46-default-x86_64-with-SuSE-11-x86_64':2,  
        'Linux-3.2.0-56-generic-x86_64-with-debian-wheezy-sid':3
        }[platform.platform()]
    except:
        print '!'*50,' PLATFORM NOT IN LIST: >{}<'.format( \
            platform.platform()), '!'*50
        sys.exit(0)

    # Mac setup                                                                                               
    if case ==1 :
        home = '/Users/Tomas/'
        d = { 
        'rwd'  :home+'PhD/Data/MUTD_iGEOS-Chem_output/',
        'dwd'  :  home+'PhD/Data/' ,
        'npwd' : home+'PhD/Data/np_arrs/' ,
        'tpwd' : home+'GITHub/tools_progs/' ,
        'ppwd' : home+'Pictures/'  
        }

    # Earth0 setup                                                                                            
    if case ==2 :
        home =  '/work/home/ts551/'
        d = { 
        'rwd'  : home +'data/all_model_simulations/iodine_runs/',
        'dwd'  :  home +'data/',
        'fwd' : home+ 'labbook/tools_progs/d_fast-J_JX/data/',
        'lwd'  :  home +'labbook/',
        'npwd' : home +'data/np_arrs/',
        'tpwd' : home +'labbook/tools_progs/' ,
        'ppwd' : home +'labbook/plots_images/'  
        }

    # Atmosviz1 setup                                                                                         
    if case ==3 :
        home =  '/home/ts551/'
        d = { 
        'rwd'  : home +'data/model/',
        'dwd'  :  home +'data/',
        'lwd'  :  home +'labbook/',
        'npwd' : home +'data/np_arrs/',
        'tpwd' : home +'labbook/tools_progs/'  
        }
    return d[input]

# ----                                                                                                                                                        
# 1.02 -  Get Latitude as GC grid box number in dimension                                                                                                                  
# ----                                                                                                                                                        
def get_gc_lat(lat, res='4x5',debug=False):
    """ Get index of lat for given resolution """
    NIU, lat_c, NIU = get_latlonalt4res(res=res)
    del NIU
    return find_nearest( lat_c, lat )

# ----                                                                                                                                                        
# 1.03 -  Get Longitude as GC grid box number in dimension                                                                                                       
# ----                                                                                                                                                        
def get_gc_lon(lon, res='4x5',debug=False):
    """ Get index of lon for given resolution """
    lon_c, NIU, NIU = get_latlonalt4res( res=res )
    del NIU
    return find_nearest( lon_c, lon )

# --------                                                                                          
# 1.04 - Get model array dimension for a given resolution                                           
# --------                                                                                          
def get_dims4res(res=None, r_dims=False, invert=True, trop_limit=False, \
        just2D=False, debug=False):
    """ Get dimension of GEOS-Chem output for given resolution """

    dims = {
    '4x5' :  (72,46,47), 
    '2x2.5':(144,91,47) , 
    '1x1' :  (360,181,47), 
    '0.5x0,5' :  (720,361,47), 
    '0.5x0.666':(121,81,47) ,
    '0.25x0.3125':(177, 115, 47)
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
    if just2D:
        vals =[]
        for i in dims.values():
            vals += [ ( i[0],i[1]) ]
        dims = dict( zip( dims.keys(), vals ) )
        if debug:
            print dims

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
def get_latlonalt4res( res='4x5', centre=True, hPa=False, nest=None, \
            dtype=None, wd=None, \
            lat_bounds=u'latitude_bnds', lon_bounds=u'longitude_bnds',\
            lon_var=u'longitude', lat_var=u'latitude', debug=False ):
    """ Return lon, lat, and alt for a given resolution. 
        This function uses an updated version of gchem's variable 
        dictionaries """ 
    
    # Kludge. Update function to pass "wd" 
    # if model output directory ("wd") not provided use default directory
    if isinstance( wd, type(None) ):
        dwd = get_dir( 'dwd') + '/misc_ref/'
        dir = {
        '4x5':'/LANDMAP_LWI_ctm',  \
        '2x2.5': '/LANDMAP_ctm_2x25',  \
        '0.5x0.666' :'LANDMAP_LWI_ctm_05x0666',  \
        '0.25x0.3125' :'LANDMAP_LWI_ctm_025x03125',  \
        }[res]
        wd = dwd +dir

    if debug:
        print res, centre, hPa, nest, type( res)

    if centre:
        # Extract lat and lon from model output data file
        with Dataset( wd+'/ctm.nc', 'r' ) as d:
            lat = np.array( d[lat_var] )    
            lon = np.array( d[lon_var] )        

    # Get edge values
    else:
        # Extract lat and lon from model output data file
        with Dataset( wd+'/ctm.nc', 'r' ) as d:
            lat = np.array( d[lat_bounds] )    
            lon = np.array( d[lon_bounds] )  

            # select lower edge of each bound, and final upper edge
            lat = [i[0] for i in lat ]+[ lat[-1][1] ]
            lon = [i[0] for i in lon ]+[ lon[-1][1] ]            
            lat, lon = [np.array(i) for i in lat, lon ]

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
    """ hPa/km convertor"""
    if reverse:
         return [ np.exp(  np.float(i) /-7.6)*1013. for i in input ]
    else:
        return [-7.6*np.log( float(i) / 1013.) for i in input ]

# --------   
# 1.16 - Find nearest
# --------
def find_nearest(array,value):
    """
        Find nearest point. Adapted from HappyLeapSecond's 
        Stackoverflow answer.  http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    idx = (np.abs(array-value)).argmin()
    return idx

# --------------                                                                                 
# 1.99 - Reference data (lon, lat, and alt) adapted from gchem - credit: GK (Gerrit Kuhlmann )             
# -------------                                                                                  
"""
 Updated to dictionary from gchemgrid (credit: Gerrit Kuhlmann ) with 
     addition grid adds 

    This [function] contains (some) grid coordinates used within GEOS-Chem 
    as numpy arrays

# Distribution advice from gchem:

# Python Script Collection for GEOS-Chem Chemistry Transport Model (gchem)
# Copyright (C) 2012 Gerrit Kuhlmann
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""                                                                                        

def gchemgrid(input=None, rtn_dict=False, debug=False):
    d = {

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

    if rtn_dict:
        return d
    else:
        return d[input]

