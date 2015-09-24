
"""
    Core functions used to be used by all function levels of 
    GEOS-Chem/Data analysis
"""

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 0 ----- Modules required

import numpy as np
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
        'irwd' : home+'PhD/Data/MUTD_iGEOS-Chem_output/',
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
        'irwd' : home + 'data/all_model_simulations/iodine_runs/iodine_runs/',  #Error?               
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
    NIU, lat_c, NIU = get_latlonalt4res(res=res)
    del NIU
    return find_nearest( lat_c, lat )

# ----                                                                                                                                                        
# 1.03 -  Get Longitude as GC grid box number in dimension                                                                                                       
# ----                                                                                                                                                        
def get_gc_lon(lon, res='4x5',debug=False):
    lon_c, NIU, NIU = get_latlonalt4res( res=res )
    del NIU
    return find_nearest( lon_c, lon )

# --------                                                                                          
# 1.04 - Get model array dimension for a given resolution                                           
# --------                                                                                          
def get_dims4res(res=None, r_dims=False, invert=True, trop_limit=False, \
        just2D=False, debug=False):
    dims = {
    '4x5' :  (72,46,47), 
    '2x2.5':(144,91,47) , 
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
def get_latlonalt4res(res='4x5', centre=True, hPa=False, nest=None, \
            debug=False):
    if debug:
        print res, centre, hPa, nest

    # Get dictionary variable name in Gerrit's GEOS-Chem dimensions list
    if centre:
        lon, lat = { 
        '4x5': ['c_lon_4x5' , 'c_lat_4x5'  ],
        '2x2.5':  ['c_lon_2x25' , 'c_lat_2x25' ], 
        '0.5x0.666' :[ 'c_lon_05x0667_EU', 'c_lat_05x0667_EU' ], 
        '0.25x0.3125': ['c_lon_025x03125_EU', 'c_lat_025x03125_EU' ]
        }[res]
    else:
        lon, lat = { 
        '4x5': ['e_lon_4x5' , 'e_lat_4x5' ],
        '2x2.5': ['e_lon_2x25' , 'e_lat_2x25' ], 
        '0.5x0.666' : ['e_lon_05x0667_EU', 'e_lat_05x0667_EU'] , 
        '0.25x0.3125': ['e_lon_025x03125_EU', 'e_lat_025x03125_EU' ]
        }[res]
              
    # Override lat/lon variable name selection for nested grid (China)
    if nest == 'CH':
        if centre:
            lon, lat = 'c_lon_05x0667_CH', 'c_lat_05x0667_CH'
        else:
            lon, lat = 'e_lon_05x0667_CH', 'c_lat_05x0667_CH'
    if hPa:
        alt = 'c_hPa_geos5_r'
    else:
        alt='c_km_geos5_r'

    # Also provide high resolution grid if requested from this function all
    if nest =='high res global':
        lon, lat = np.arange(-180, 180, 0.25), np.arange(-90, 90, 0.25)
        return lon, lat, alt

    # Return values from gchemgrid
    d = gchemgrid(rtn_dict=True)
    if debug:
        print lon, lat, alt
    return  [ d[i] for i in lon, lat, alt ]

# --------------
# 1.06 - Convert from hPa to km or vice versa.
# -------------
def hPa_to_Km(input, reverse=False, debug=False):
    if reverse:
         return [ np.exp(  np.float(i) /-7.6)*1013. for i in input ]
    else:
        return [-7.6*np.log( float(i) / 1013.) for i in input ]

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

# --------------                                                                                 
# 1.99 - Reference data (lon, lat, and alt) adapted from gchem - credit: GK (Gerrit Kuhlmann )             
# -------------                                                                                  
"""
    Updated to dictionary from gchemgrid with addition grid adds
"""                                                                                        

def gchemgrid(input=None, rtn_dict=False, debug=False):
    d = {
    # 4x5                                                                                        
    'c_lon_4x5' : np.arange(-180, 175+5, 5) ,
    'e_lon_4x5' : np.arange(-182.5, 177.5+5, 5) ,
    'c_lat_4x5' : np.array( [-89]+ range(-86, 90, 4)+ [89] ),
    'e_lat_4x5' : np.array( [-90]+ range(-88, 92, 4)+ [90] ),

    # 2x2.5                                                                                      
    'e_lon_2x25' : np.arange(-181.25, 178.75+.5,2.5),
    'c_lon_2x25' : np.arange(-180, 177.5+2.5, 2.5),
    'e_lat_2x25' :np.array( [-90]+range(-89, 91, 2)+[90] ) ,
    'c_lat_2x25' : np.array( [-89.5]+range(-88, 90, 2)+[89.5] ),

    # China nested grid                                                                          
    'c_lon_05x0667_CH' : np.arange(70, 150+0.667, 0.667) ,
    'c_lat_05x0667_CH' : np.arange(-11, 55+.5, 0.5) ,
    'e_lon_05x0667_CH' : np.arange(70-(0.667/2), 150+0.667, 0.667) ,
    'e_lat_05x0667_CH' : np.arange(-11-(.5/2), 55+.5, 0.5) ,

    # Europe nested grid (0.5*0.5)                                                                    
    'c_lon_05x0667_EU'  : np.array([
    float(i) for i in np.arange(-30., 50+0.667, 0.667 ) 
    ]) ,
    'c_lat_05x0667_EU' : np.arange(30, 70+.5, 0.5),
    'e_lon_05x0667_EU' : np.arange( -30.333, 50.333+0.666, 0.667),
    'e_lat_05x0667_EU' : np.arange(29.75, 70.25+0.25, 0.5 ),

    # Europe nested grid (0.25*0.3125)                                                                    
    'c_lat_025x03125_EU' : np.arange(32.750, 61.25+.25, .25), 
    'c_lon_025x03125_EU' : np.arange(-15, 40+.3125, .3125), 
#    'e_lat_025x03125_EU' : np.arange(32.750-(.25/2), 61.25+(.25/2)+.25, .25), 
#    'e_lon_025x03125_EU' : np.arange(-15-(.3125/2), 40+(.3125/2)+.3125, .3125), 

    'e_lat_025x03125_EU': np.array([ 32.75,  33.  ,  33.25,  33.5 ,  33.75,  
        34.,  34.25,  34.5 ,
        34.75,  35.  ,  35.25,  35.5 ,  35.75,  36.  ,  36.25,  36.5 ,
        36.75,  37.  ,  37.25,  37.5 ,  37.75,  38.  ,  38.25,  38.5 ,
        38.75,  39.  ,  39.25,  39.5 ,  39.75,  40.  ,  40.25,  40.5 ,
        40.75,  41.  ,  41.25,  41.5 ,  41.75,  42.  ,  42.25,  42.5 ,
        42.75,  43.  ,  43.25,  43.5 ,  43.75,  44.  ,  44.25,  44.5 ,
        44.75,  45.  ,  45.25,  45.5 ,  45.75,  46.  ,  46.25,  46.5 ,
        46.75,  47.  ,  47.25,  47.5 ,  47.75,  48.  ,  48.25,  48.5 ,
        48.75,  49.  ,  49.25,  49.5 ,  49.75,  50.  ,  50.25,  50.5 ,
        50.75,  51.  ,  51.25,  51.5 ,  51.75,  52.  ,  52.25,  52.5 ,
        52.75,  53.  ,  53.25,  53.5 ,  53.75,  54.  ,  54.25,  54.5 ,
        54.75,  55.  ,  55.25,  55.5 ,  55.75,  56.  ,  56.25,  56.5 ,
        56.75,  57.  ,  57.25,  57.5 ,  57.75,  58.  ,  58.25,  58.5 ,
        58.75,  59.  ,  59.25,  59.5 ,  59.75,  60.  ,  60.25,  60.5 ,
        60.75,  61.  ,  61.25]) , 
    'e_lon_025x03125_EU' : np.array([-15.  , -14.69, -14.38, -14.06, -13.75, 
        -13.44, -13.12, -12.81,
       -12.5 , -12.19, -11.88, -11.56, -11.25, -10.94, -10.62, -10.31,
       -10.  ,  -9.69,  -9.38,  -9.06,  -8.75,  -8.44,  -8.12,  -7.81,
        -7.5 ,  -7.19,  -6.88,  -6.56,  -6.25,  -5.94,  -5.62,  -5.31,
        -5.  ,  -4.69,  -4.38,  -4.06,  -3.75,  -3.44,  -3.12,  -2.81,
        -2.5 ,  -2.19,  -1.88,  -1.56,  -1.25,  -0.94,  -0.62,  -0.31,
         0.  ,   0.31,   0.62,   0.94,   1.25,   1.56,   1.88,   2.19,
         2.5 ,   2.81,   3.12,   3.44,   3.75,   4.06,   4.38,   4.69,
         5.  ,   5.31,   5.62,   5.94,   6.25,   6.56,   6.88,   7.19,
         7.5 ,   7.81,   8.12,   8.44,   8.75,   9.06,   9.38,   9.69,
        10.  ,  10.31,  10.62,  10.94,  11.25,  11.56,  11.88,  12.19,
        12.5 ,  12.81,  13.12,  13.44,  13.75,  14.06,  14.38,  14.69,
        15.  ,  15.31,  15.62,  15.94,  16.25,  16.56,  16.88,  17.19,
        17.5 ,  17.81,  18.12,  18.44,  18.75,  19.06,  19.38,  19.69,
        20.  ,  20.31,  20.62,  20.94,  21.25,  21.56,  21.88,  22.19,
        22.5 ,  22.81,  23.12,  23.44,  23.75,  24.06,  24.38,  24.69,
        25.  ,  25.31,  25.62,  25.94,  26.25,  26.56,  26.88,  27.19,
        27.5 ,  27.81,  28.12,  28.44,  28.75,  29.06,  29.38,  29.69,
        30.  ,  30.31,  30.62,  30.94,  31.25,  31.56,  31.88,  32.19,
        32.5 ,  32.81,  33.12,  33.44,  33.75,  34.06,  34.38,  34.69,
        35.  ,  35.31,  35.62,  35.94,  36.25,  36.56,  36.88,  37.19,
        37.5 ,  37.81,  38.12,  38.44,  38.75,  39.06,  39.38,  39.69,  40.  ]), 

    # generic arrays                                                                             
    'e_lon_generic' : np.arange(-180.0,181.0),
    'c_lon_generic' : np.arange(-179.5,180.0),
    'e_lat_generic' : np.arange(-90.0,91.0),
    'c_lat_generic' : np.arange(-89.5,90.0),

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

