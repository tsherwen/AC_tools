# =========================================================================================================
# ------------------------------ tms - module of programs/code for re-use ----------------------------------
# -------------- 
# Section 1 - Plotters funcs
# Section 2 - Model Data Extractors
# Section 3 - Data Processing tools/drivers
# Section 4 - Data Analysis funcs
# Section 5 - Plotting Ancillaries 
# Section 6 - Time processing
# Section 7 - Generic processing 

# --------------- ------------- ------------- ------------- ------------- 
# --------------  Contents
# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 0 ----- Required modules

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 1 ----- Model Data Extractors
# 1.01 - open ctm.bpch 
# 1.02 - get np array (4D) of ctm.bpch ( lon, lat, alt, time)
# 1.03 - bulk ctm.bpch processor (Process mod output to ctms)
# 1.04 - Get surface area
# 1.05 - Get land map
# 1.06 - Get Species Gg
# 1.07 - Get air mass
# 1.08 - get OH mean (read from geos.log)
# 1.09 - Read X-Calc files (in need of update)
# 1.10 - BD file reader (from DB)
# 1.11 - SSA I2 file extractor
# 1.12 - Extract all pf data for a given site.
# 1.13 - Extract array from dir, saving the intermediate array to disc
# 1.14 - Extract OH and HO2 for a given ctm.bpch file
# 1.15 - Extract spec (surface) data from HDF for given time
# 1.16 - Open and re-structure data for a given 
# 1.17 - convert hdf (table to pandas panel) for spec
# 1.18 - Extract surface data for a given location from nested run HDF
# 1.19 - Extract spec ( EU surface) data from HDF for given time
# 1.20 - Extract spec ( EU surface) data from HDF for given time ( for a family)
# 1.21 - Process species for given arrays to (v/v) in respective scale + DU
# 1.22 - Get var data for model run ( PyGChem >0.3.0 version combatible )
# 1.23 - Get  surface data from HDF of surface plane flight data
# 1.24 - Get gc resolution from ctm.nc


# --------------- ------------- ------------- ------------- ------------- ------------- ------------- 
# ---- Section 2 ----- Data Processing tools/Extractors - GC...
# 2.01 - Retrieve model resolution 
# 2.02 - Get model array dimension for a given resolution - mv's to core
# 2.03 - Get GC ctm.bpch years 
# 2.04 - Get GC ctm.bpch months
# 2.05 - Get GEOS-Chem longitude index for a given latitude
# 2.06 - Get GEOS-Chem latitude index for a given latitude
# 2.07 - Get GC datetime
# 2.08 - Work out iGEOS-Chem version
# 2.09 - Get GC species burden of species (driver)
# 2.10 - Get GC species prod loss (driver)
# 2.11 - Get GC species emission (driver)
# 2.12 - Get CH4  lifetime (driver)
# 2.13 - Get GC Land/Water/Ice indices (LWI)
# 2.14 - Get GEOS-Chem altitude index for a given latitude
# 2.15 - Get Gg from v/v
# 2.16 - Get GC volume 
# 2.17 - Test is grid box for a given lat and lon is over water
# 2.18 - Deposition for a given species 
# 2.19 - Get Land map
# 2.20 - Convert 3D prod/loss molec/cm3/s2 to Gg/month (driver) 
# 2.21 - Convert 2D prod/loss tracer (molec/cm2/s2) to Gg/month (driver) 
# 2.22 - get Dry deposition for a given species
# 2.23 - Get GC run stats (driver) ***
# 2.24 - Get DU mean value (area weighted)
# 2.25 - Get O3 burden (redundent? - get_spec_burden does this already?)
# 2.26 - Get Ox prod(POX)/loss(LOX)
# 2.27 - Get Transport fluxes (driver)
# 2.28 - Get species wet deposition 
# 2.29 - Get species Boundary layer (BL) mixing
# 2.30 - Get cloud flux
# 2.31 - Volume weigh numbers
# 2.32 - Get GAW site info  (lat, lon, alt (press), timezone (UTC) ) 
# 2.33 - Return mesh and indics for which obs are within 
# 2.34 - takes indices or generates indices for a given set of sites, these are then substract from all arrays given
# 2.35 - Get Ox loss by route
# 2.36 - Split Tropospheric Ox by loss route families 
# 2.37 - Select only active wet dep tracers from given list
# 2.38 - split 4D ( lon, lat, alt, time) output by season

# --------------- 
# ---- Section 3 ----- Time processing
# => all generic funcs moved to time funcs
# 3.01 - Sort list of ctms into chronological order
# 3.02 - Takes monthly outputs and sort to chronological order

# ------------------------------------------- Section 0 -------------------------------------------
# -------------- Required modules:
#
#!/usr/bin/python
# -- I/O / Low level                                                                                
import os
import sys
import csv
import glob
from PIL import Image
import pygchem
if pygchem.__version__ == '0.2.0':
    import pygchem.diagnostics as gdiag
else:
    try:
        from pygchem import datasets
    except:
        import pygchem.datafields as datasets
import pandas as pd
import re
from bisect import bisect_left
from netCDF4 import Dataset
#import iris # Kludge as iris is not on Mac

# - Math/Analysis                                                                                   
import numpy as np
from time import mktime
import scipy.stats as stats
from pandas import HDFStore
from pandas import DataFrame
from pandas import Series
from pandas import Panel
import math

# -- Time                                                                                           
import time
import calendar
import datetime as datetime
from datetime import datetime as datetime_
#from iris.time import PartialDateTime # Kludge as iris is not on Mac

# Import tms modules with shared functions
from AC_tools.funcs4core import *
from AC_tools.funcs4generic import *
from AC_tools.funcs4time import *
#from funcs4obs import * #( need: get_CVO_DOAS_obs ... )
from AC_tools.funcs4pf import *
from AC_tools.funcs_vars import *

# --------------------------------- Section 2 ----------------------------------
# -------------- Model Data Extractors
#

# --------------
# 1.01 - open ctm.bpch using pygchem ( version '0.2.0' ) 
# -------------- 
def open_ctm_bpch(wd, fn='ctm.bpch', debug=False):
    """
        This is a vestigle programme, based on pychem version 0.2.0.
        Updates have made this incompatibile. 
        
        Update functions to use iris class. 
    """
    
    if ( debug ) :
        print 'fn: {}'.format(fn)
        print 'File to open: {}'. format(os.path.join(wd, fn))

    if pygchem.__version__ == '0.2.0':
        try:
            ctm_f = gdiag.CTMFile.fromfile(os.path.join(wd, fn))
        except:
            print 'Error @ open_ctm_bpch for {}'. format(os.path.join(wd, fn))  
            print gdiag.CTMFile.fromfile(os.path.join(wd, fn))
            sys.exit(0)
    else:
        print 'WARNING, using: {}'.format(pygchem.__version__)
        print 'Add a get data call, for specific data (spceies, category)' +\
                'using the get_GC_output function' 

    return ctm_f

# --------------
# 1.02 - get np array (4D) of ctm.bpch ( lon,lat , alt,time) 
# --------------
def get_gc_data_np(ctm, spec='O3', category="IJ-AVG-$", debug=False):
    """ REDUNDENT - This function extracts fields from diagnostics and 
        was written to work with - this is included for back compatibility"""

    if debug:
        print 'called get_np_gc_4D_diags'

    # Retrieve given diagnostics
    diags = ctm.filter(name=spec, category=category)
    if debug:
        print 'diagnostics', diags , spec, category

    # Extract diagnostics to np array    
    for diag in diags:
        ar = (diag.values[:,:,:])[...,None]
        if debug:
            print diag.name ,'len(ar)', len(ar), 'type(ar)', type(ar), \
                'diag.scale', diag.scale, 'ar.shape', ar.shape, 'diag.unit', \
                diag.unit
        try:
            arr = np.concatenate( (arr, ar), axis=3 )
        except NameError:
            arr = ar
        if debug:
            print 'arr' , type(arr), len(arr), arr.shape, 'ar', type(ar), \
                len(ar), ar.shape

    if len([d for d in diags]) <1 :
        print 'ERROR: No diags for {} and {} in {}'.format( spec, category, ctm)
        sys.exit(0)

    # if ctms len > 1, sort chronologically - this is to fix bug in pygchem
    if debug:
        print [ i[0] for i in get_gc_datetime(ctm ) ]
    if get_gc_datetime( ctm ) > 1:
        arr = np2chronological_fromctm([ ctm ], arr, debug=debug )
    if debug:
        print [ i[0] for i in get_gc_datetime(ctm ) ]

    return arr

# --------   
# 1.03 - Process mod output to ctms
# --------
def wd2ctms(wd, fn='ctm.bpch', cat_="IJ-AVG-$", spec='O3', 
                        start=None, end=None, debug=False):
    """ REDUNDENT - this function returns chronologically ordered 
        extracted ctm.bpch files - this is included for back compatibility"""

    # List Bpch files in dir and sort by name... 
    # (this assumes that months are labelled by YYYYMM etc...)
    if wd[-1] != '/':
        wd +=  '/'
    fns = list( sorted( glob.glob( wd + '*ctm*' ) ) )

    if pygchem.__version__ == '0.2.0':

        # Select for date  - form:  YYYYMM.ctm.bpch
        if ( (start != None) and  (end !=None) ):
            wds = [ i for i in fns ]# if ()
        else:
            wds = fns   
        if debug:
            print wds, [ i.split('/')[-1] for i in wds ], \
                ['/'.join(i.split('/')[:-1]) for i in wds ]

        # Open files
        ctms =  [ open_ctm_bpch('/'.join(i.split('/')[:-1]), i.split('/')[-1]) 
            for i in wds ]

        # What resolution
        res = mod_res(wd, fn=fn)
        if debug:
            print res, wd

        # Sort ctms to chronological order
        if debug:
            print ctms
        ctms = ctms2chronological( ctms )

        return ctms, res

    # Re-write will be need for programmes using this technique
    else:
        print 'WARNING: PyGChem {} '.format ( pygchem.__version__ ) + \
             'no longer combatible with this approach!'
        sys.exit(0 )
    
# ----
# 1.04 -Get surface area  ( m^2 )
# ----
def get_surface_area(res='4x5',time=None, debug=False):
    """ Get_surface_area of grid boxes for a given resolution
         - this function accessres previsouly run GEOS-Chem 
            1day files with just DXYP / DXYP diagnostic ouptuted
         - back compatibility with PyGChem 0.2.0 is retained  """

    if debug:
        print 'called get surface area'
    dwd = get_dir( 'dwd') + '/misc_ref/'
    dir = {
    '4x5':'/LANDMAP_LWI_ctm',  \
    '2x2.5': '/LANDMAP_ctm_2x25',  \
    '0.5x0.666' :'LANDMAP_LWI_ctm_05x0666',  \
    '0.25x0.3125' :'LANDMAP_LWI_ctm_025x03125',  \
    }[res]
    fd = dwd +dir
    if debug:
        print fd, res
    if pygchem.__version__ == '0.2.0' :
        ctm = gdiag.CTMFile.fromfile( fd +'/ctm.bpch' )
        diags = ctm.filter(name="DXYP",category="DXYP")#,time=date_list[0])     

        # extract diags
        first_time=True
        if debug:
            print diags
        while first_time:
            for diag in diags:
                if debug:
                    print diag.unit
                s_area = diag.values
            first_time=False
    else:
        s_area = get_GC_output( fd, vars=['DXYP__DXYP'] ) 

    return s_area

# ----
# 1.05 - Get Land map. 
# ----
def get_land_map(res='4x5', time=None, wd1=None,debug=False):
    """ Return land, water, ice indices (LWI ) from GEOS-Chem
        with integeres for Land (1) and Water (0). 
        Ice fraction is given as fractional values """

    if debug:
        print 'called get surface area, for ', res
    dwd = get_dir( 'dwd') + '/misc_ref/'
    dir = {
    '4x5':'/LANDMAP_LWI_ctm',  \
    '2x2.5': '/LANDMAP_LWI_ctm_2x25',  \
    '0.5x0.666' :'LANDMAP_LWI_ctm_05x0666',  \
    '0.25x0.3125' :'LANDMAP_LWI_ctm_025x03125',  \
    }[res]
    land_dir = dwd +dir
    if debug:
        print land_file

    # retain backwards compatibility
    if pygchem.__version__ == '0.2.0':
        ctm_f = gdiag.CTMFile.fromfile( land_dir + '/ctm.bpch' )
        diags = ctm_f.filter(name="LWI",category="LANDMAP")

        # Exctract from ctm.bpch
        first_time=True
        if (first_time):
            for diag in diags:
                scalar = diag.values[:,:,:]
                first_time=False
                if debug:
                    print diag.name ,'len(scalar)',len(scalar), 'type(scalar)',\
                        type(scalar), 'diag.scale', diag.scale, 'scalar.shape',\
                        scalar.shape,'diag.unit',diag.unit
                landmap=scalar

    # use new pygchem (>0.3.0) approach
    else:
        landmap = get_GC_output(  land_dir, species="LWI", \
            category="LANDMAP")

    return landmap

# -------------- 
# 1.06 - Get Gg of species 
# -------------
def get_Gg_species(ctm_f,species,air_mass,v_box, diagnostics=None,debug=False):
    diagnostics = ctm_f.filter(name=species, category="IJ-AVG-$")
    """ Get array of species in Gg
        Note: not compatible with PyGChem 0.3.0  """

    for i, diag in enumerate(diagnostics):
        scalar = diag.values[:,:,:]
        if debug:
            print diag.name ,'len(scalar)',len(scalar), 'type(scalar)',\
                type(scalar), 'diag.scale', diag.scale, 'scalar.shape', \
                scalar.shape, 'diag.unit', diag.unit

        #  number of moles in box
        moles = ( air_mass[0] * 1E3) / constants( 'RMM_air')     

        # g of I in species
        scalar = ( ( scalar * moles ) * 127 * spec_stoich(species) )[...,v_box] 

        try:
            np.append(scalar_list, scalar)
            scalar_total = np.array(scalar_total, scalar)
        except:
            scalar_list  = np.array(scalar)
            scalar_total = np.array(scalar)

    # average of values provided & convert to Gg
    scalar = ( scalar_total / len(diagnostics) )  /1E6 
    return scalar

# --------------
# 1.07 - get air mass (4D) np.array  in kg
# ------------- 
def get_air_mass_np( ctm_f=None, wd=None, times=None, trop_limit=True,\
            debug=False ):
    """ Get array of air mass in kg
        Note: not compatible with PyGChem 0.3.0  """

    if debug:
        print 'called get air mass'

    # retain PyGChem 0.2.0 version for back compatibility
    if pygchem.__version__ == '0.2.0':
        diagnostics = ctm_f.filter( name='AD', category="BXHGHT-$",time=times )
        if debug:
            print diagnostics
        for diag in diagnostics:
            # Extract data
            scalar = np.array( diag.values[:,:,:] )[...,None]              
            if debug:
                print diag.name ,'len(scalar)',len(scalar), 'type(scalar)',\
                    type(scalar), 'diag.scale', diag.scale, 'scalar.shape', \
                    scalar.shape,'diag.unit', diag.unit
            try:
                arr = np.concatenate( (np_scalar, scalar), axis=3 )
            except NameError:
                arr = scalar

    # Or use PyGChem 0.3.0 approach
    else:
        # Get air mass in kg
        arr = get_GC_output( wd=wd, vars=['BXHGHT_S__AD'], \
                            trop_limit=trop_limit, dtype=np.float64)

        if debug:
            print 'arr' , type(arr), len(arr), arr.shape
    return arr

# --------------
# 1.08 - Get OH mean value - [1e5 molec/cm3]
# -------------
def get_OH_mean( wd, debug=False ):
    """ Get mean OH concentration (1e5 molec/cm3) from geos.log files
         in given directory """

    if debug:
        print wd
    # find all geos log files...
    files = glob.glob( wd+'/*geos*log*' )

    # Loop and extract OH means, take an average if n>0
    z=[]
    for file in files:
        with open(file) as f:
            if debug:
                print file, f, z
            for line in f:
                if "Mean OH =    " in line:
                    z += [ float(line.split()[3] ) ]
    return np.mean(z)

# --------------
# 1.08 - Get CH4 mean value - v/v
# -------------
def get_CH4_mean( wd, rtn_global_mean=True, debug=False ):
    """ Get mean CH4 concentrtaion from geos.log files in given directory """

    if debug:
        print wd
    # find all geos log files...
    files = glob.glob( wd+'/*geos*log*' )

    # Loop and extract CH4 means, take an average if n>0
    regex = re.compile( 'CH4{}:'.format('.'*13) )
    z=[]
    for file in files:
        with open(file) as f:
            if debug:
                print file, f, z
            for line in f:
                if re.match(regex, line):
                    z += [ float(line.split()[-2])/1E9 ]
    if debug:
        print z, np.mean(z)
    if rtn_global_mean:
        rtn_value = np.mean(z)
    # return value for UK latitude
    else:
        rtn_value = np.array( z )[::4].mean()
    return rtn_value

# --------------
# 1.09 - Calculated X-section printer/
# -------------
def readfile_calc_Xsect(filename, species):
        """ Read X-sections from given file """

        for files in filename:
            #        print files
            reader=csv.reader(open(files,'rb'), delimiter=' ', \
                skipinitialspace = True)
            tmp = []
            count = 0
            for row in reader:
                #    count = 0
                if row[0] == species and count < 1:
                    #    print(species)
                    #  print(row[2:])
                    tmp = tmp + row[2:]
                    count = count + 1
                    for i in range(2): #print(row)
                        tmp = tmp + reader.next()
                #    print(reader.next())
                        GCbins = tmp[-7:]
                        GCbins = " ".join(GCbins)
                #   return GCbins
                        print("{0} {1}".format(species, GCbins)) 
                #   return GCbins
            #        tmp = str(tmp)
            #        tmp = ', '.join(tmp)
             #       print(tmp)
          #  else :
           #     pass
       # return GCbins  
#     print(tmp)
    #        return
      #  print(tmp)

# --------------
# 1.10 - DB's CV file reader
# -------------
def read_CV(filename,obs_data_name, obs_switch):
    """ Read binary of CV planeflight output
        credit: drb """
    data = np.genfromtxt(filename,skip_header = 428, names=True)
    names = data.dtype.names
    time = data['Day_of_yearstart']
    if obs_switch == 0:
        var1=data[obs_data_name]
    elif obs_switch == 1:
        var1=data[obs_data_name]+273.15

    return time, var1

# ----
# 1.11 - collate geos.log oupput for Aerosol release, and write to numpy array
# ----               
def get_SSA_I2( wd, s_total=0, debug=False ):
    """ get SSA I2 values from geos.log and write to np array
         - this is the binary writer for IGAC plots/analysis """

    np_wd = get_dir('np_wd')

    s_size =  48.
    arr, count, total = np.zeros((72,46,47,s_size)), 0., 0.

    for n, line in enumerate( open(wd+'/geos.log') ):
        if debug:
            print line.split()

        if 'PHYSPROC: Trop chemistry at' in line:
            count += 1.
            total += 1.

        if debug:
            print count, total
        if ( total > s_total ) :#or (count <144)):  # set start output point
            if 'tms check box SSA I2 EMISS :' in line:
                if debug:
                    print line.split()
                data = float( line.split()[7] )

            if 'tms check box - I, J, L : ' in line:
                if debug:
                    print line.split()
                I, J, L = [ int(line.split()[i]) for i in range(-3,0) ]
                if debug:
                    print [(i, type(i) ) for i in I, J, L, count, data ]
                arr[I-1, J-1, L-1,count ] = data
                if debug:
                    print I, J, L, count, total, data
            if (count == s_size):
                fp = np.memmap( np_wd+'SSA_I2_source{}_{:0>4}.npy'.format( \
                    wd.split('/')[-2].split('SSA_I2_source')[-1], total ), dtype='float32', \
                    mode='w+', shape=(72,46,47,s_size))

                fp[:] = arr[:]
                np.memmap.flush(fp)            
                count = 0.
                arr  = np.zeros((72,46,47,s_size))

                print np.float64(arr).shape
        if (count == s_size):
            count = 0.

# ----
# 1.12 - Extract all pf data for a given site. - REDUNDENT
# ----

# ----
# 1.13 - Extract array from dir, saving the intermediate array to disc
# ----
def wds2np( wds, spec='O3', res='4x5', category="IJ-AVG-$", trop_limit=False, 
                        months=12, dtype='float32', debug=False ):
    """ Read CTM output saved in Binary form
        REDUNDENT - now use automated NetCDF 2D/3D writer 
        https://github.com/tsherwen/MChem_tools/blob/master/bpch2netCDF.py
        """

    mod = []
    for wd in wds:        
        fstr = get_dir('npwd') +'GC_arr_store_{}_{}_{}_{}.npy'.format(spec, \
            wd.split('/')[-2], wd.split('/')[-1], res)
        if debug:
            print fstr, wd, res, category, trop_limit, months
        try:
            fp = np.memmap( fstr,  dtype=dtype, mode='r',\
                        shape=tuple(list(get_dims4res(res, \
                        trop_limit=trop_limit))+ [months] )  ) 
            print 'Read memory mapped np array', fp.shape,   \
                list(get_dims4res(res, trop_limit=trop_limit))+[months]
            mod += [  np.array( fp ) ]
            del fp
        except:
            print 'Writing memory mapped np array', fstr
#            ctms, res = wd2ctms( wd )#, cat_, spec, start, end )
#            arr = np.concatenate([get_gc_data_np( ctm, spec=spec, \
#                category=category, debug=debug) for ctm in ctms ], axis =3)
#            arr = np2chronological_fromctm(ctms, arr)
            # Use pygchem 0.3.0 to extract
            get_GC_output( wd,  species=spec, category=category )

            if trop_limit:
                arr=arr[:,:,:38,:]
            fp = np.memmap( fstr,  dtype=dtype, mode='w+', shape=arr.shape )
            fp[:] = arr[:]
            np.memmap.flush(fp)        
            mod += [ fp ] 
            del fp
    if debug:
        print np.ma.array(mod).shape, [np.ma.array(i).shape for i in mod]
    return mod, res

# ----
# 1.14 - Extract OH and HO2 for a given ctm.bpch file
# ---
def get_OH_HO2( ctm=None, t_p=None, a_m=None, vol=None, debug=False, \
            wd=None, HOx_weight=False, res='4x5', scale=1E5 ):
    """ Get OH/HO2 concentrations from ctm.bpch file """

    # Set/Get vars
    vars = ['OH','HO2']
    AVG= constants('AVG')
    RMM_air=constants('RMM_air')
    # load arrays if not provided...
    if not isinstance(t_p, np.ndarray):
        if pygchem.__version__ == '0.2.0':
            t_p = get_gc_data_np( ctm, spec='TIMETROP', category='TIME-TPS', \
                debug=debug)
        else:
            print 'WARNING!!! - provide time in trop diag. ' 
    if not isinstance(a_m, np.ndarray):
        if pygchem.__version__ == '0.2.0':
            a_m = get_air_mass_np( ctm, debug=debug )
        else:
            print 'WARNING!!! - provide air mass diag. ' 
    if not isinstance(vol, np.ndarray):
        if pygchem.__version__ == '0.2.0':
            vol = get_volume_np( ctm, res=res, debug=debug)
        else:
            print 'WARNING!!! - provide vol diag. ' 

    # Extract Data 
    if pygchem.__version__ == '0.2.0':
        OH, HO2 = [get_gc_data_np( ctm, \
            i, category='CHEM-L=$')*t_p[:,:,:38,:] /scale for i in vars ] 
    else:
        OH_HO2 = get_GC_output( wd, \
            vars=['CHEM_L_S__'+i for i in vars] )*t_p[:,:,:38,:]/scale
        OH, HO2 =  [OH_HO2[i,...] for i in range(len(vars )) ]

    # Convert data
    molecs = ( ( (a_m*1E3) / RMM_air ) * AVG )   # molecules
    moles =  a_m*1E3 / RMM_air # mols 

    # only consider troposphere
    molecs, a_m, vol, moles = [i[:,:,:38,:]*t_p[:,:,:38,:]  \
                for i in [molecs, a_m, vol, moles ]]  
    if debug:
        print [ (i.shape, np.mean(i) ) for i in [ OH, HO2, molecs, moles, a_m, \
            vol ]]

    # convert HO2 to molec/cm3
    # v/v * (mass total of air (kg)/ 1E3 (converted kg to g))  = moles of tracer => molecs.
    HO2 = HO2 * moles * AVG 
    HO2 =  HO2 / vol

    HO2, OH  = [ np.ma.masked_invalid( i) for i in HO2, OH   ]
    HO2, OH  = [ np.ma.masked_where( i<0, i) for i in HO2, OH   ]

    if debug:
        print [ ( np.ma.sum(i), np.ma.mean(i) ) for i in HO2, OH ]
        print 'volume weighted: ', [ np.ma.sum( i *vol) / np.ma.sum(vol) \
            for i in HO2, OH ] 

    if HOx_weight:
        HOx = HO2 + OH
        HO2, OH  = [ np.ma.sum( ( i* HOx )) /np.ma.sum(HOx) for i in HO2, OH   ]
        if debug:
            print [ i.shape for i in OH, HO2, moles, vol, HOx ]
            print  'HOx weighted: ',HO2, OH

    # volume weight
    else:
       HO2, OH = [ np.ma.sum( i *vol) / np.ma.sum(vol)  for i in HO2, OH ] 

    return OH, HO2

# ----
# 1.15 - Extract spec ( EU surface) data from HDF for given time
# ---

# ----
# 1.16 - Make 3D arrays from mass storage HDF of EU run
# ----
def mk_HDF_2D_from_table( spec='O3', nchunks=12, \
# Kludge for the partial year output. 
#            starttime=datetime.datetime(2005,2,1, 0, 0), 
#            endtime=datetime.datetime(2006,2, 1, 0, 0), \
#            res='2x2.5', region='ROW', \
#            run_names=[  'Just_Br' ],\
# updated full year dates  simulation
            starttime=datetime.datetime(2005,1,1,0,0), 
            endtime=datetime.datetime(2006,1,1,0,0), \
            res='0.5x0.666', region='EU', run_code='',  \
            run_names=[  'run' ], ver='1.7', \
            debug=False ):
    """ Make HDF file (table form ) to hold GEOS-Chem planeflight csv output 
        in binary form  
        REDUNDENT - Now use automated NetCDF writer
        https://github.com/tsherwen/MChem_tools/blob/master/bpch2netCDF.py""" 

    if debug:
        print 'mk_HDF_2D_from_table called'

#    run_names = [  'run', 'no_hal' ]    

    # --- Vars
    npwd=get_dir('npwd')
    
    if region == 'ROW':
#        run_name = ''
        nchunks = 11 # Kludge! - allow HDF writer to work for 11 month run

    if region == 'EU':
        nchunks = 12 # default, presume 1 year (12 months) of output

    # Kludge to allow analysis of EU EOH run / high res (0.25)
    nchunks=1

    # --- Loop Runs, creating a DataFrame per mont
#    for spec in slist:
    for run_name in run_names:
        big_dates, big_data = [], []
        for n in range(nchunks):         

            dates, data = hdf2panel( spec, starttime=add_months(starttime, n),                                                                                                                                  
                        endtime=add_months(starttime, n+1), run_name=run_name,
                        region=region, run_code=run_code, res=res, debug=debug  )
            big_dates += dates
            big_data +=  data
            if debug:
                print  [ i.shape for i in big_data ], len( big_dates  )
            print n, spec, run_name
        P = Panel(  dict(zip(big_dates, big_data))  )
        print P.shape

        # use standard naming 
        identifier1 = ver # ( e.g. '1.7', 'v10',  ... )
        identifier2 = region # ( e.g.  'EU', 'ROW' ... ) 
        identifier3 = res #( e.g. '4x5','2.5x2.5' ... )
        identifier4 = run_name #( e.g. 'run', 'no_hal' ... ) 
        identifier5 = run_code #( e.g.  'EOH', 'ICE' ... )
        identifier6 = spec #( e.g.  'EOH', 'ICE' ... )
        fname = npwd +'{}_surface_{}_{}_{}_{}_{}.h5'.format( identifier1, \
            identifier2, identifier3, identifier4, identifier5, identifier6 )
        print identifier1, identifier2, identifier3, identifier4, identifier5, identifier6
        print fname

        P.to_hdf( fname, 'd1', mode='w' )
                        
        del P 

# ----
# 1.17 - convert hdf (table to pandas panel) for spec
# ---
def hdf2panel( spec='O3', run_name='run',  region='EU',\
            res='0.5x0.666', run_code='', 
# Kludge for the partial year output. 
# region='ROW',\
#            starttime=datetime.datetime(2005,2,1, 0, 0),   \
#            endtime=datetime.datetime(2006,2,1, 0, 0),  \
# updated full year dates  simulation
            starttime=datetime.datetime(2005,1,1, 0, 0),   \
            endtime=datetime.datetime(2006,2,1, 0, 0),  \
            ver='1.7', debug=False ):
    """ stack tabulated planeflight output in HDF to pandas 3D panel """

#    debug=True

    # --- Vars
    # deal with run name issues... 
    # ACTION REQUIRED! - srename file to use EU instead of NP

#    if region == 'EU':
#        region ='NP'   # remove Kludge

    # this is no longer needed! 
#    if region == 'ROW':  
#        run_name =run_name

    # Set main vars
    npwd=get_dir('npwd')
    identifier1 = ver # ( e.g. '1.7', 'v10',  ... )
    identifier2 = region # ( e.g.  'EU', 'ROW' ... ) 
    identifier3 = res #( e.g. '4x5','2.5x2.5' ... )
    identifier4 = run_name #( e.g. 'run', 'no_hal' ... ) 
    identifier5 = run_code #( e.g.  'EOH', 'ICE' ... )

    hdf_filename = glob.glob(npwd +'/'+'pf*{}*{}*{}*{}*{}.h5'.format( identifier1, 
        identifier2, identifier3, identifier4, identifier5 )  )
    if debug:
        print npwd, region, run_name, hdf_filename
    hdf_filename = hdf_filename[0]
#    hdf_filename = npwd +'pf_iGEOSChem_1.7_v10_G5_EU_run.0.25x0.3125.2012.1week.II.h5'
#    hdf_filename = npwd +'pf_iGEOSChem_1.7_v10_G5_EU_run.0.25x0.3125.2012.1week.h5'
#    hdf_filename = npwd +'pf_iGEOSChem_1.7_v10_G5_EU_run.0.25x0.3125.2012.1week.I_I.h5'

    # transalte species to planeflight nomenclature
    GCFP_d2TRA_all = what_species_am_i( ver=ver, rtn_dict=True, invert=True)
    pspec = GCFP_d2TRA_all[ spec ]

    # Get datetimes between start and end
    dates = dt_hrs_a2b(a=starttime, b=endtime)
    if debug:
        print spec, pspec, dates[0], dates[-1]

    # Get Data in 3D arrays  (as pandas Dataframes )
    dates, dfs = read_HDF4spec( hdf_filename, dates, spec=pspec, 
        starttime=starttime, endtime=endtime )

    return dates, dfs 

def read_HDF4spec( hdf_filename, dates, spec='O3', \
        starttime=datetime.datetime(1985, 1, 1), 
        endtime=datetime.datetime(2100, 1, 1), debug=False ):
    """
        Reads HDF file for a given specifes and returns numpy array of 
        dates and list of data frames of spec's data contstructed into 3D 
        arrays 

        notes:
            - spec is required in pf format ( e.g. TRA_?? )
    """
    debug=True

    if debug:
        print hdf_filename , dates[0], dates[-1], spec
        print [type(i) for i in starttime, endtime ]

    # --- Open file
    df = pd.read_hdf( hdf_filename, 'd1', \
                where="Datetime>=starttime and Datetime<endtime", \
                columns=['LON','LAT', spec])        

    # Now select dates after querry.  - RESTORE preivsou method. 
#    df = df[starttime:endtime]
    if debug:
        print df.shape, df.index[0], df.index[1]

    # --- Split into list of DataFrames
    dfs = [  
    df[ (df.index>=i)  & (df.index<dates[n+1]) ]    \
    for n, i in enumerate( dates[:-1] ) 
    ] 
    if debug:
        print [ i.shape for i in dfs ]
    del df # is this wise?

    # --- Construct 3D arrays
    dfs = [ 
    DataFrame( i[spec].values, index=[i.LON, i.LAT] ).unstack()  for i in dfs 
    ]
    print [ len(i) for i in dates, set(dates), dfs ]

    return dates, dfs

# ----
# 1.18 - Extract surface data for a given location from nested run HDF
# ---
def get_surface_data( spec='O3', lname='Weyborne', time_span='year', \
            res='0.5x0.666', run_name='run', scale=1E9, region='EU', \
            ver='1.7', starttime=None, endtime=None, dates=None, debug=False):

    """ Extract surface data stored in HDF file 

        - set time_span to None to calculate for non standard 
            ( 2005 analysis periods - DOI: 10.5194/acpd-15-20957-2015 )
        
        REDUNDENT: now use NetCDF approach
        https://github.com/tsherwen/MChem_tools/blob/master/pf2NetCDF.py 
           """

    # Kludge, was is this need for func_vars extraction to work?
    from AC_tools.funcs_vars import gaw_2_name 

    # --- Vars
    lon, lat, alt = get_loc( lname )
    lat, lon = get_gc_lat(lat, res=res), get_gc_lon(lon, res=res)
    gaw_2_name = gaw_2_name()

    if isinstance( dates, type(None) ):
        dates = get_dt4run( time_span=time_span, a=starttime, b=endtime )  
    wd=get_dir('dwd')+'np_arrs/'

#    titles =  [ 'GEOS-Chem v9-2', 'Iodine Sim.']
#    run_names = ['run', 'no_hal']

    # Kludge 
#    starttime=datetime.datetime(2005,2,1, 0, 0) 
#    endtime=datetime.datetime(2006,1, 31, 23, 0)  # 1 year
#    dates = get_dt4run(a=starttime, b=endtime) 

    if debug:    
        print dates
        print  dates[0], dates[-1]

    # For both with an without halogens
    data = get_spec_pf_hdf_surface(spec=spec, run_name=run_name, \
                    starttime=dates[0], endtime=dates[-1], ver=ver, \
                    region=region, res=res, debug=debug )*scale
    data = data[:, lon, lat ]
    print [ i.shape for i in [data]]
    return dates, data

# ----
# 1.19 - Extract spec ( EU surface) data from HDF for given time
# ---
def get_spec_pf_hdf_surface( run_name='run', spec='O3', region='ROW',
            starttime=datetime.datetime(2005,1,1, 0, 0),  run_code='',  \
            endtime=datetime.datetime(2005,12, 31, 23, 0), \
            ver='1.7', res='2x2.5', just_return_values=True, debug=False ):

    """ Extracts species surface data from HDF store for a given date range in a 
    given region ( e.g. ROW or EU... ) """

    print 'get_spec_pf_hdf_surface called, for ', run_name
    npwd=get_dir('npwd')
    run_names = [run_name]

    # use standard naming 
    identifier1 = ver # ( e.g. '1.7', 'v10',  ... )
    identifier2 = region # ( e.g.  'EU', 'ROW' ... ) 
    identifier3 = res #( e.g. '4x5','2.5x2.5' ... )
    identifier4 = run_name #( e.g. 'run', 'no_hal' ... ) 
    identifier5 = run_code #( e.g.  'EOH', 'ICE' ... )
    identifier6 = spec #( e.g.  'EOH', 'ICE' ... )

    print identifier1, identifier2, identifier3, identifier4, identifier5, identifier6

    # is this already formed as a numpy?
    fn= npwd +'{}_surface_{}_{}_{}_{}_{}.h5'.format( identifier1, 
            identifier2, identifier3, identifier4, identifier5, identifier6 )

    try:
        print 'reading data, for ', fn
        P = pd.read_hdf( fn, 'd1', mode='r' )
        print P.shape#, [ (i.min(), i.max(), i.mean()) for i in [P] ]

    except:
        # extracts species surface data for a given date range  
        print 'Did not find/could not open: {}'.format(fn)
        print 'exacting and writing data '

        if debug:
            print npwd

        mk_HDF_2D_from_table( spec=spec, region=region, res=res,\
            starttime=starttime, endtime=endtime, run_names=run_names, 
            ver=ver, run_code=run_code, debug=debug ) 
        P = pd.read_hdf( fn, 'd1', mode='r' )
    
    # Select for given dats and return solely values as np array
    if debug:
        print P.shape, starttime, endtime
#    P= P.select(P.items>=starttime)#, axis=0 )  # This worked prior 15 05 26 
#    P = P.select( lambda x: x>= starttime, axis=0 )   # use Lambda instead.
#    P = P.select( lambda x: endtime>x>= starttime, axis=0 )   # control by start & end time
    P = P.select( lambda x: endtime>=x>= starttime, axis=0 )   # control by start & end time

    if debug:
        print P.shape

    # Just return values 
    if just_return_values:
        arr = P.values
        del P
    else:
        arr =P

    if debug:
        print [ [np.array(i).shape, np.mean(i), np.max(i), np.min(i)] for i in [arr ]]
    return arr 

# ----
# 1.20 - Extract spec ( EU surface) data from HDF for given time ( for a family)
# ---
def get_spec_pf_hdf_surface_procesor( run_name='run', spec='O3', \
            starttime=datetime.datetime(2005,1,1, 0, 0),   \
#            endtime=datetime.datetime(2005,11, 30, 23, 0), \
            endtime=datetime.datetime(2005,12, 31, 23, 0), \
            res='2x2.5', region='ROW', run_code='', just_return_values=True, \
            ver='1.7', debug=False ):

    """ Linker to allow for whole species families to be extracted 
        ( as well as single species ) """

    familes = 'NOy', 'Iy'
    # set stoich options to False
    N=False; I=False 

    if spec in familes:
        # get specs

        # NOy
        if spec == 'NOy':
            specs = GC_var('N_specs')
            # Kludge! - not all NOy was ouputted... 
            missing = 'R4N2', 'NH3',  'NH4',  'BrNO2', 'BrNO3', 'MPN', 'ISOPN',\
                'PROPNN', 'MMN'                
            [specs.remove(i) for i in missing ]
#            specs = specs[7:] # Kludge 
#           specs =  specs[3:] # Kludge 

            N = True
            
        # Iy 
        if spec == 'Iy':
            specs = GC_var( 'Iy' )
            I = True

        # loop for specs
        arr =[]
        for spec in specs:

            # GOTCHA - if scaled to pptv then divide by 1E3
#            units, scale = tra_unit( spec, scale=True )

            # Scale all to ppbv - N/A all in v/v  ( scale on bulk )
            units, scale = 'ppbv', 1E9

            # get N or I etc in spec - maintain stiochiometry
            stoich = spec_stoich( spec=spec, I=I, N=N )            

            # get data            
            arr += [get_spec_pf_hdf_surface(spec=spec, run_name=run_name, \
                res=res, starttime=starttime, endtime=endtime, run_code=run_code, \
                ver=ver, region=region, just_return_values=False, \
                debug=debug) *stoich* scale ]
        print type( arr )
        print [ type(i) for i in arr ]
        print np.array( [i.values for i in arr] ).shape
        arr = np.ma.array( [i.values for i in arr] ).sum(axis=0)
        print arr.shape  
            
    else:
        # scale species
        units, scale = tra_unit( spec, scale=True )

        # get data
        arr = get_spec_pf_hdf_surface(spec=spec, run_name=run_name, \
            res=res, starttime=starttime, endtime=endtime, \
            region=region, just_return_values=False, debug=debug) *scale
    
    return arr, units


# ----
# 1.21 - Process species for given arrays to (v/v) in respective scale + DU
# ---
def process_data4specs( specs=None, just_bcase_std=True, \
            just_bcase_no_hal=False, res='4x5', ver='1.6', diff=True, \
            pcent=True, tight_constraints=True, trop_limit=True, \
            NOy_family=False, Bry_family=False, Iy_family=False, \
            rtn_specs=False, debug=False ): 
    """ Return species values in v/v and DU. Also return datetimes for CTM 
        output and time in troposphere diagnostics  """

    # Get runs data and descriptors
    wds, titles = MUTD_runs( titles=True, just_bcase_std=just_bcase_std, \
        res=res, ver=ver, just_bcase_no_hal=just_bcase_no_hal )
    if debug:
        print wds, titles

    # --- Previous pygchem approach <= REDUNDENT
    if pygchem.__version__ == '0.2.0':
        runs  = [ wd2ctms( wd, debug=debug ) for wd in wds ]

        # Get species time Tropopause diagnostic
        t_ps  = [  np.concatenate(  [ get_gc_data_np( ctm, spec='TIMETROP', \
            category='TIME-TPS') for ctm in ctms ] , axis=3 ) \
            for ctms in [r[0] for r in runs] ]

        # get molecs in troposphere
        molecs  = [ arr*1E3/constants( 'RMM_air')*constants('AVG') \
            for arr in [np.concatenate([ get_air_mass_np(ctm) \
            for ctm in ctms ], axis=3 ) for ctms in [r[0] for r in runs] ] ]

        res= runs[0][1]
        s_area  = get_surface_area( res=res, debug=debug )[:,:,0]

        # Extract data for arrays once
        Vars = [ [ np.concatenate( [ get_gc_data_np( ctm, spec) \
            for ctm in ctms], axis=3) for ctms in [r[0] for r in runs] ] \
            for spec in specs ]
        del runs

        # convert to np arrays
        Vars, molecs, t_ps = [ np.ma.array(i) for i in Vars, molecs, t_ps ]

        # limit to GEOS-Chem's chemical troposphere
        Vars, molecs = [ i[...,:38,:] for i in Vars, molecs ]

        # Get GEOS-Chem datetime
        dlist = [ j for k in [get_gc_datetime(ctm ) 
            for ctm in  runs[0][0] ] for j in k ]
        dlist = sorted( set([  i[0] for i in dlist ])  )

    # ---- Updated pygchem approach, using iris cubes
    else:

        # Get GC dates as datetimes
        dlist = get_gc_datetime( wd=wds[-1] )

        # Get species time Tropopause diagnostic
        t_ps = [ get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], \
            trop_limit=trop_limit ) for wd in wds ]
    
        # Get air mass in kg
        a_m = [ get_GC_output( wd, vars=['BXHGHT_S__AD'], \
            trop_limit=trop_limit, dtype=np.float64) for wd in wds ]

        # Get molecs in troposphere
        molecs = [i*1E3/constants( 'RMM_air')*constants('AVG') for i in a_m ]

        # Get surface area
        s_area = get_surface_area( res=res, debug=debug )[:,:,0]

        # Get mixing ratio 
        Vars = [ get_GC_output( wd, vars=['IJ_AVG_S__'+i for i in specs ], \
            trop_limit=trop_limit ) for wd in wds ]
#        print [i[0].shape for i in Vars, molecs, t_ps ]
        print [i.shape for i in Vars + molecs + t_ps ]

        # convert to np arrays
        Vars, molecs, t_ps = [ np.ma.concatenate([ii[...,None] for ii in i], \
            axis=-1) for i in Vars, molecs, t_ps ]
        print [i.shape for i in Vars, molecs, t_ps ]

        # restore to shape used for previous analysis 
        molecs, t_ps= [ np.rollaxis( i, -1, 0 ) for i in molecs, t_ps ]
        print [i.shape for i in Vars, molecs, t_ps ]
        Vars=np.rollaxis( Vars, -1, -5 )
        print [i.shape for i in Vars, molecs, t_ps ]

    if debug:
        print [ ( i.shape, i.sum() ) for i in Vars, molecs, t_ps ]
    
    # apply tight constraints to troposphere? ( 
    if tight_constraints:
        t_ps = np.ma.masked_where( t_ps != 1, t_ps )

    if debug:
        print [ ( i.shape, i.sum() ) for i in Vars, molecs, t_ps ]
#    Vars, molecs = [ i*t_ps[None,...] for i in Vars, molecs ]
    # remove stratospheric values
    Vars = Vars*t_ps[None,...] 
    molecs  =  (molecs*t_ps)[None,...]

    if debug:
        print [ ( i.shape, i.sum() ) for i in Vars, molecs, t_ps ]
    
    # Get DU values for each array (2D not 3D )
    # molecules O3 in trop. summed
    DUarrs = Vars*molecs
    if debug:
        print [ ( i.shape, i.sum() ) for i in [DUarrs] ]
    DUarrs = DUarrs.sum(axis=-2)
    # Kludge back to 4D  ( this is due to the conversion to NetCDF for 0.3.0 )
    DUarrs=DUarrs[...,None,:]

    if debug:
        print [ ( i.shape, i.sum() ) for i in [DUarrs] ]

    # adjust to per unit area ( cm3 ) 
    tmp_s_area = s_area[None, None,...]
    tmp_s_area = tmp_s_area[...,None,None] 
    DUarrs =  DUarrs/ tmp_s_area

    if debug:
        print [ ( i.shape, i.sum() ) for i in DUarrs, tmp_s_area,s_area  ]

    # convert to DU   
    DUarrs = DUarrs/constants('mol2DU')

    if debug:
        print [ ( i.shape, i.sum() ) for i in DUarrs, tmp_s_area,s_area  ]

    if debug:
        for n in range( len(specs) ):
            for nn in range(2):
                print n,nn
                print '!'*30, (DUarrs[n,nn,...]*tmp_s_area[0,0,...]).sum()/  \
                    tmp_s_area[0,0,...].sum()

    # ---  Add familys to arrays by conbination of tracers
    # Shape of Vars, DUarrs = [(65, 2, 72, 46, 38, 12), (65, 2, 72, 46, 1, 12)]
    # (NOy_family, Bry_family, Iy_family ) 
    def add_family2arrays( fam=None, Vars=None, DUarrs=None, 
                I=False, N=False, Br=False, specs=None ):

        # get species in family
        d = { 'NOy':'N_specs', 'Iy':'Iy',  'Bry':'Bry' }
        fam_specs = GC_var( d[fam] )
        # Kludge - only consider the specs available in the NetCDF
        fam_specs = [i for i in fam_specs if (i in specs) ]

        print '!1.1'*200, [ i.shape for i in Vars, DUarrs], fam, d[fam], \
            fam_specs
    
        #  find indices and conctruct dictionary
        fam_specs = [(i,n) for n,i in enumerate( specs ) if ( i in fam_specs ) ]
        fam_specs = dict( fam_specs )

        # Add family name to specs 
        specs += [ fam ]

        # Extract values and adjust to stiochmetry  ( Vars )
        fam = [ Vars[fam_specs[i],...]*spec_stoich(i, I=I, N=N, Br=Br) \
                for i in fam_specs.keys() ]
        
        # sum species in family and add to returned array ( Vars )
        print np.array( [np.array(i) for i in fam] ).shape
        fam = np.ma.array( fam ).sum(axis=0)[None,...]
        print [i.shape for i in fam, Vars ]
        Vars = np.ma.concatenate( (Vars, fam) , axis=0 )

        print '!1.2'*200, [ i.shape for i in Vars, DUarrs]
        del fam

        # sum species in family and add to returned array ( DUarrs )
        fam = [ DUarrs[fam_specs[i],...]*spec_stoich(i, I=I, N=N, Br=Br) \
                for i in fam_specs.keys() ]

        # sum species in family and add to returned array ( DUarrs )
        fam = np.ma.array( fam ).sum(axis=0)[None,...]
        print [i.shape for i in fam, DUarrs ]
        DUarrs = np.ma.concatenate( (DUarrs, fam) , axis=0 )

        print '!1.3'*200, [ i.shape for i in Vars, DUarrs]
        del fam

        return Vars, DUarrs, specs

    print '!0'*200, [ i.shape for i in Vars, DUarrs]    

    if NOy_family:
        Vars, DUarrs, specs = add_family2arrays(  fam='NOy', N=True,\
                Vars=Vars, DUarrs=DUarrs, specs=specs )
    if Bry_family:
        Vars, DUarrs, specs = add_family2arrays(  fam='Bry', Br=True, \
                Vars=Vars, DUarrs=DUarrs, specs=specs )
    if Iy_family:
        Vars, DUarrs, specs = add_family2arrays(  fam='Iy', I=True, \
                    Vars=Vars, DUarrs=DUarrs, specs=specs )
    print '!2'*200, [ i.shape for i in Vars, DUarrs]

    # consider difference in v/v and DU
    if diff:
        if pcent:
            Vars, DUarrs = [ np.ma.array(i) for i in Vars, DUarrs ]
            Vars = (Vars[:,1,...] - Vars[:,0,...] )/ Vars[:,0,...] *100
            DUarrs = ( DUarrs[:,1,...] - DUarrs[:,0,...] ) / \
                DUarrs[:,0,...] *100
            Vars, DUarrs = [ np.ma.masked_invalid(i) for i in Vars, DUarrs ]

        else:
            Vars = Vars[:,1,...] - Vars[:,0,...]      
            DUarrs = DUarrs[:,1,...] - DUarrs[:,0,...] 
    else:
        Vars = Vars[:,1,...] 
        DUarrs = DUarrs[:,1,...]
    print [i.shape for i in DUarrs, Vars ]

    # Close/remove memory finished with
    del molecs, tmp_s_area, wds, titles

    # return requested variables
    rtn_vars =[  Vars, DUarrs, dlist, t_ps.mean(axis=0) ]
    if rtn_specs:
        rtn_vars += [specs ]
    return rtn_vars

# ----
# 1.22 - Get var data for model run
# ---
def get_GC_output( wd, vars=None, species=None, category=None, \
            r_cubes=False, r_res=False, restore_zero_scaling=True, \
            r_list=False, trop_limit=False, dtype=np.float32, use_NetCDF=True, \
            debug=False):
    """
        Data extractor for GEOS-Chem using pygchem (>= 0.3.0 )
        ( Credit: Ben Bovy -  https://github.com/benbovy/PyGChem )

        ctm.bpch files are extracted for a given directory, with only 
        the specific category and species returned. This can either be as
        a Iris Cube (retaining metadata) or as a numpy array.

            - To return Iris Cubes, set r_cubes to True
            - To get resolution of model run, set r_res to True

        For simplicity use variable names (species, category) as used in 
        the iris cube ( e.g. IJ_AVG_S__CO). if variable is unknown, just 
        print full dataset extracted to screen to see active diagnostics. 

        Species and category variables are maintained ( and translated ) 
        to allow for backwards compatibility with functions written for 
        pygchem version 0.2.0
        
        This replaces the now defunct functions: open_ctm_bpch 
        and get_gc_data_np.

        Examples and use of pygchem is discussed on Ben Bovy's GITHib 
        https://github.com/benbovy/PyGChem_examples/blob/master/Tutorials/Datasets.ipynb        
    """
    
    if debug:
        if not isinstance( vars, type(None) ):
            print 'Opening >{}<, for var: >{}<'.format( wd, ','.join(vars) ) + \
                '(additional) gamap variables provided: >{}< + >{}<'.format( \
                category, species )
                
    # Temporary fix for back compatibility: 
    # Convert gamap names ( species  + category ) to iris cube names
    # Just for use whilst updating functions written to use pygchem 0.2.0
    if not any( [ isinstance(i, type(None) ) for i in species, category ] ):
        # convert to Iris Cube name
        if (category == None) and ( vars ==  None ):
            category = "IJ-AVG-$"
        if (species == None) and ( vars ==  None ):
            species =  'O3'
    
        # remove scaling for 'IJ-AVG-$' - KLUDGE - needed for all ?
        if category == 'IJ-AVG-$':
            category = diagnosticname_gamap2iris( category )

        # setup var for Iris Cube, then remove known rm. chars.
        var = category+'__'+species
        vars=  [ var.replace('-', '_').replace('$', 'S') ]

    else:
        # Get default settings for reader
        if isinstance(vars, type(None)):
            vars = [ 'IJ_AVG_S__O3'  ]

    # ensure wd has a leading '/'
    if wd[-1] != '/':
        wd +=  '/'

    # Work with NetCDF. Convert ctm.bpch to NetCDF if not already done.
    if use_NetCDF:

            # Check for compiled NetCDF file
            # If not found, create NetCDF file from ctm.bpch files                
            fname = wd+ '/ctm.nc'
            import os.path
            if not os.path.isfile(fname):
                from bpch2netCDF  import convert_to_netCDF
                convert_to_netCDF( wd )

            print fname
            # "open" NetCDF + extract requested variables as numpy arr.
            with Dataset( fname, 'r' ) as rootgrp:
                arr = [ np.array(rootgrp[i]) for i in vars ]  

                # files are stored in NetCDF at GC scaling. 
                # ( This is different to ctm.bpch, rm for back compatibility. )
                if restore_zero_scaling:
                    try:
                        arr =[ arr[n]/get_unit_scaling( rootgrp[i].ctm_units ) \
                            for n, i in enumerate( vars ) ]
                    except:
                        print 'WARNING: SCALING NOT ADJUSTED TO' + \
                            ' PREVIOUS APPROACH'

    # Use Iris cubes via PyGChem to extract ctm.bpch files 
    else:
    
        # Get files in dir ( more than one? )
#        fns = sorted( glob.glob( wd+ '/*ctm*' ) )
#        if len(fns) == 0:
#            fns = glob.glob(wd + '*trac_avg*')
#            print 'using *trac_avg* bpch name convention: ', fns

#        if debug:
#            print fns

        # Load files into Iris Cube
#        print 'WARNING - with ctm.bpch files, all variables are loaded.'+\
#            'Convert to NetCDF or HDF for speed up or use use_NetCDF=True)'
#        cubes = datasets.load( fns, vars )

        # If no data extracted, print our variables
#        try:
#            [ i[0].data for i in cubes ]
#        except:
#            print datasets.load( fns[0] )
#            print 'WARNING: no vars found for >{}<'.format( ','.join(vars) )
#            sys.exit( 0 )
        
        # Temporary fix for back compatibility: 
        # Just extract as numpy 

#        if not r_cubes:
    
            # Extract data
#            try:
#                arr = [ cubes[i].data for i in range( len(vars) ) ]
#            except:
#                print 'WARNING: length of extracted data array < vars'
#                print 'vars: >{}<'.format( ','.join(vars) )
#                sys.exit( 0 )

            #  restore to zero scaling (v/v instead of pptv etc ) - only for GC tracers
#            if restore_zero_scaling: 
#                if (category == 'IJ_AVG_S') or ('IJ_AVG_S' in vars[0] ):
#                    print [ cubes[i].attributes['ctm_units'] for i in range( len(vars) ) ]
#                    arr = [ arr[n]/get_unit_scaling( cubes[n].attributes['ctm_units'])
#                             for n, var in enumerate( vars ) ]
#            del cubes
        print 'WARNING this approach has been removed due to time cost'
        sys.exit( 0 )

    # Process extracted data to gamap GC format and return as numpy 
    if not r_cubes:

        # Limit to GEOS-Chem "chemical troposphere'
        if trop_limit:
            arr = [ i[...,:38] for i in arr]

        # convert to GC standard 4D fmt. - lon, lat, alt, time 
        if len((arr[0].shape)) == 4:
            if debug:
                print 'prior to roll axis: ', [i.shape for i in arr]
            arr = [np.rollaxis(i,0, 4) for i in arr]
            if debug:
                print 'post to roll axis: ', [i.shape for i in arr]

        # convert to GC standard 3D fmt. - lon, lat, time   
        # two reasons why 3D (  missing time dim  or  missing alt dim ) 
        # <= update, there must be a better gotcha than this... 
        if ( (len((arr[0].shape)) == 3 ) \
            # '4x5'
            and  ( (72, 46, 47) != arr[0].shape ) \
            and ( (72, 46, 38) != arr[0].shape ) \
            # '2x2.5'
            and  ( (144, 91, 47) != arr[0].shape ) \
            and ( (144, 91, 38) != arr[0].shape ) \
            # '0.25x0.3125':
            and ( (177, 115, 47) != arr[0].shape ) \
            and ( (177, 115, 38) != arr[0].shape )  \
            # '0.5x0.666'
            and  ( (121, 81, 47) != arr[0].shape ) \
            and ( (121, 81, 38) != arr[0].shape ) ):

            if debug:
                print 'prior to roll axis: ', [i.shape for i in arr]
            arr = [np.rollaxis(i,0, 3) for i in arr]
            if debug:
                print 'post to roll axis: ', [i.shape for i in arr]

        # For multiple vars, concatenate to var, lon, lat, lat, time 
        if len(vars) >1:
            arr = np.concatenate( [ i[None,...] for i in arr ], axis=0 )
        else:
            arr =arr[0]
            
        # Add altitude dimension to 2D (lon, lat)
        # a bug might occur for emission etc ( where lon, lat, time are dims )
        if len((arr.shape)) == 2:
            arr =arr[...,None]

        if dtype != np.float32:
            arr = arr.astype( dtype )

    # Get res by comparing 1st 2 dims. against dict of GC dims.
    if r_res:
        res=get_dims4res( r_dims=True, trop_limit=trop_limit, \
                    just2D=True )[arr.shape[:2]]

    if r_list:
        arr = [ arr[i,...] for i in range( len(vars) ) ]

    # Sort output - return Cubes?
    if r_cubes:
        output =  cubes
    else:
        output = arr
#        return cubes.data

    # Return model resolution? 
    if r_res:
        return output, res
    else:
        return output


# ----
# 1.23 - Get  surface data from HDF of surface plane flight data
# ---
def get_surface_data_HDF( spec='O3', res='2x2.5', region='ROW', \
        actual_values=True, change=False, pcent=False, just_no_hal=False, \
        ref_name='Just_Br',main_name='run', ver='1.7', min_change=0.5,     \
#        starttime=datetime.datetime(2005,2,1, 0, 0),      # start Feb
        starttime=datetime.datetime(2005,1,1, 0, 0),      # start jan
        endtime=datetime.datetime(2005,12, 31, 23, 0),   # 1 year   
        run_code='', debug=False): 

    """ Get surface data for given region and given species. 
        
        This is a linker to allow for selection of data based on run name.
        
        Defunct. Now use 3D NetCDF made by:
        https://github.com/tsherwen/MChem_tools/blob/master/pf2NetCDF.py
    """
    
#    if debug:
#        print starttime, endtime
    
#    debug=True
    
    # Get data from iodine containing run. 
    if (actual_values or change) and (not just_no_hal):
        arr1, units = get_spec_pf_hdf_surface_procesor(spec=spec, \
            run_name=main_name, res=res, ver=ver, starttime=starttime, \
            endtime=endtime, run_code=run_code, debug=debug)

        if debug:
            print [ [ i.shape, i.min(), i.max(), i.mean()] for i in [arr1]  ]

        # plot up actual values
        arr, arr1  = [ np.ma.array(i) for i in arr1, arr1 ] 
        title = '{} Surface Concentration '.format( latex_spec_name(spec) )

    if change:
        # get ref name ( e.g no halogen, or Just Br ) values
        arr2, units = get_spec_pf_hdf_surface_procesor(spec=spec, \
            run_name=ref_name, res=res, ver=ver, starttime=starttime, \
            endtime=endtime, debug=debug)
        arr2 = np.ma.array( arr2 )

        if debug:
            print [ [ i.shape, i.min(), i.max(), i.mean()] for i in [arr2] ]
        
    # Get precentage change
    if pcent:
        arr = np.ma.subtract(arr1, arr2)
        if debug:
            print 0,'!'*5, 'detail on output: ', [ [ i.min(), i.max(), \
                i.mean(), type(i), i.shape ] for i in arr1, arr2 ]
        
        # Mask for tiny absolute changes
        arr = np.ma.masked_where( np.fabs(arr)<min_change, arr) 
        if debug:
            print 1.20,'!'*5, 'detail on output: ', [ [ i.min(), i.max(), \
                i.mean(), type(i), i.shape ] for i in arr1, arr2 ]

        arr = np.ma.divide( arr, arr2 )*100
        if debug:
            print 1.21,'!'*5, 'detail on output: ', [ [ i.min(), i.max(), i.mean() ] \
                for i in arr1, arr2 ]

    # Get actual change
    if (change) and (not pcent):
        if debug:
            print 1.22,'!'*5, 'detail on output: ', [ [ i.min(), i.max(), i.mean() ] \
                for i in arr1, arr2 ]
        arr = np.ma.subtract(arr1, arr2)

    # Just return no halogen run values
    if just_no_hal:
        arr = arr2

        if debug:
            print 2, 'detail on output: ', [ [ np.ma.min(i), np.ma.max(i), \
                np.ma.mean(i) ] for i in [ arr2 ] ] 

    title = '$\Delta$ {} Surface Concentration / {}'.format( \
            latex_spec_name(spec), units )

    return arr, title


# ----
# 1.24 - Get gc resolution from ctm.nc
# ---
def get_gc_res( wd ) :

    # create NetCDf if not created.
    fname = wd+ '/ctm.nc'
    if not os.path.isfile(fname):
        from bpch2netCDF  import convert_to_netCDF
        convert_to_netCDF( wd )

    # "open" NetCDF + extract time
    with Dataset( fname, 'r' ) as rootgrp:
        lon = rootgrp['longitude']
        lat = rootgrp['latitude']
#            lvls = rootgrp['model_level_number']
        lat, lon = [ np.array(i) for i in lat, lon ]

    # compare with dictionary to get resoslution    
    dims = (len(lon),len(lat))
    res = get_dims4res( r_dims=True, just2D=True )[dims]

    return res 

# ------------------------------------------- Section 3 -------------------------------------------
# -------------- Data Processing tools/drivers
#

# --------   
# 2.01 - Retrieve model resolution
# --------
def mod_res(wd, spec='O3', fn='ctm.bpch', debug=False):
    """ Extract model resolution
        note: this is not compatible with PyGChem 0.3.0 """

    if debug:
        print '>'*10, wd,  glob.glob(wd + '*ctm*'), glob.glob(wd + '*trac_avg*')
    # assume v9-2... ( ctm.bpch output ... )
    try:
        fn = glob.glob(wd + '*ctm*')[0].split('/')[-1]
    except IndexError:
        try:
            fn = glob.glob(wd + '*trac_avg*')[0].split('/')[-1]
        except:
            print 'ERROR: for wd: {}'.format( wd )
            sys.exit( 0 )
            
    ar = get_gc_data_np( open_ctm_bpch(wd, fn), spec, debug=debug )
    if debug:
        print ar.shape
    if (len(ar[:,0,0,0]) == 72 ):
        res = '4x5'
    elif (len(ar[:,0,0,0]) == 144 ):
        res='2x2.5'
    elif (len(ar[:,0,0,0]) == 121 ):
        res='0.5x0.666'        
    return res

# ----
# 2.03 - get GC years
# -----
def get_gc_years(ctm_f=None, wd=None, set_=True, debug=False):
    """ Return list of years in CTM output """

    if pygchem.__version__ == '0.2.0':
        diagnostics = ctm_f.filter(name='O3', category="IJ-AVG-$" )

        if (set_):
            return list(set( [diag.times[0].strftime('%Y')  \
                for diag in diagnostics ] ))
        else:
            return [ diag.times[0].strftime('%Y')  for diag in diagnostics ]

    # Use updated PyGChem ( >0.3.0 ) approach
    else:
        dates = get_gc_datetime( wd=wd )
        return [ i.year for i in dates ]         

# ----
# 2.04 - get GC months
# -----
def get_gc_months(ctm_f=None, wd=None, debug=False):
    """ Return list of months in CTM output """

    if pygchem.__version__ == '0.2.0':
        diagnostics = ctm_f.filter(name='O3', category="IJ-AVG-$" )
        return [diag.times[0].strftime('%m')  for diag in diagnostics ]

    # Use updated PyGChem ( >0.3.0 ) approach
    else:
        dates = get_gc_datetime( wd=wd )
        return [ i.month for i in dates ] 

# ----
# 2.07 - Get gc datetime
# -----
def get_gc_datetime(ctm_f=None, wd=None, spec='O3', cat='IJ-AVG-$', \
            debug=False):
    """ Return list of dates in datetime output from CTM output """

    # REDUNDENT - retain for backwards compatibility
    if pygchem.__version__ == '0.2.0':
        d  = ctm_f.filter(name=spec, category=cat)
        if debug:
            print '>'*30,d
            print '>.'*15, sorted(d)
        return [ i.times for i in d ]

    # Extract datetime from cube    
    else:
        if debug:
            print  '>'*5, wd, glob.glob( wd+ '/*ctm*' )

        # create NetCDf if not created.
        fname = wd+ '/ctm.nc'
        if not os.path.isfile(fname):
            from bpch2netCDF  import convert_to_netCDF
            convert_to_netCDF( wd )

        # "open" NetCDF + extract time
        with Dataset( fname, 'r' ) as rootgrp:
            dates = rootgrp['time']
            print dates, dates.units
            # Pull out units from cube, default is 'hours since 1985-01-01 00:00:00'
            starttime = time.strptime( str(dates.units),\
                    'hours since %Y-%m-%d %H:%M:%S' )
            starttime = time2datetime( [starttime] )[0]
            dates = np.array(dates) 

        if debug:
            print starttime

        # allow for single date output <= is there a better gotcha than this?
        if len(dates.shape)== 0 :
            dates = [float(dates)] 

        # convert to date time
        dates = [ add_hrs( starttime, i ) for i in dates ]
        if debug:
            print dates
        
        # return datetime objects
        return dates
    
# -------------
# 2.08 - Work out iGEOS-Chem version- e.g. iGeosChem 1.1 or 1.2 from dir name ( wd ) 
# ------------- 
def iGEOSChem_ver(wd, debug=False):
    """ get GEOS-Chem verson - iGEOS-Chem """

    vers = [
    '1.1','1.2', '1.3', '1.4', '1.5', '1.6', '1.7', '2.0'
    ]
    v = [ (i in wd) for i in vers ]
    if debug:
        print vers, v
    return [vers[n] for n, i in enumerate(v) if i==True ][0]


# --------------
# 2.09 - Get tropospheric Burden - 
# -------------
def get_trop_burden( ctm=None, spec='O3', wd=None, a_m=None, t_p=None, \
        Iodine=False, all_data=True, total_atmos=False , res='4x5',  \
        trop_limit=True, debug=False ):
    """ Get Trosposheric burden. If requested (total_atmos=False)
        limit arrays to "chemical troposphere"  (level 38) and boxes that are in 
        the troposphere removed by multiplication of "time in troposphere"
        diagnostic
        REDUNDENT - use get_gc_burden instead """

    if not isinstance(a_m, np.ndarray):
        a_m = get_air_mass_np( ctm_f=ctm, wd=wd, debug=debug, \
            dtype=np.float64 )[...,:38,:]
    if not isinstance(t_p, np.ndarray):
        # retain PyGChem 0.2.0 approach for back compatibility
        if pygchem.__version__ == '0.2.0':
            t_p = get_gc_data_np( ctm=ctm, spec='TIMETROP', \
                category='TIME-TPS', debug=debug)
        else:
            t_p = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], \
                            trop_limit=trop_limit  ) 

    # retain PyGChem 0.2.0 approach for back compatibility
    if pygchem.__version__ == '0.2.0':            
        ar = get_gc_data_np( ctm, spec, debug=debug )[:,:,:38,:]
    else:
        ar = get_GC_output( wd, vars=['IJ_AVG_S__'+ spec], \
                            trop_limit=trop_limit  ) 
    if debug:
        print [i.shape for i in ar, t_p, a_m ]

    # v/v * (mass total of air (kg)/ 1E3 (converted kg to g))  = moles of tracer
    ar = ar* ( a_m[...,:38,:]*1E3 / constants( 'RMM_air')) 
    if (Iodine):
        # convert moles to mass (* RMM) , then to Gg 
        ar = ar * float( species_mass('I') ) * spec_stoich(spec) /1E9 
    else:
        # convert moles to mass (* RMM) , then to Gg 
        ar = ar * float( species_mass(spec) ) /1E9 

    # Cut off at the "chemical troposphere" ( defined by GC integrator as 38th)
    if (not total_atmos):        
        ar = ar[...,:38,:] * t_p[...,:38,:]
    else:
        print '!'*200, 'TOTAL ATMOS'

    if debug:
        print 'Got buden for {} from {}'.format( spec, ctm )
    if (all_data):
        return ar
    else:
        return np.mean( ar, axis=3 )

# --------------
# 2.10 - Get prod/Loss (change) ( [molec/cm3/s] =>  Gg I /s => Gg/  per month )
# -------------
def get_pl_in_Gg(ctm_f, specs, years=None, months=None, monthly=False,\
            Iodine=True, IO=False, I=False, res='4x5', vol=None, spec='O3', \
            ver='1.6', debug=False  ):
    """ Return prod/loss diagnostic for griven p/l species/tags
        Note: this approach is not comparible with PyGChem >0.3.0"""

    if debug:
        print specs
    # -- Get Vars
    if not isinstance(vol, np.ndarray):
        vol = get_volume_np( ctm_f, res=res, debug=debug)
    if (  not isinstance(years, list) or  not isinstance(months, list)):
        years, months  = get_gc_years( ctm_f, set_=False ), get_gc_months(ctm_f)

    # --- Process data  # [molec/cm3/s] =>  Gg I /s 
    ars = [ molec_cm3_s_2_Gg_Ox_np( arr, specs[i], vol, ctm_f, Iodine=Iodine,\
             IO=IO, I=I, spec=spec, year_eq=False, debug=debug)    \
            for i, arr in enumerate( [ get_gc_data_np(ctm_f, PLO3_to_PD( rxn,  \
            ver=ver, fp=True ), category=GC_var('p_l')[0], debug=debug ) \
            for rxn  in specs ])  ]

    # adjust to per month.
    day_adjust = d_adjust( months, years)
    ars = [  ar_* day_adjust  for ar_ in [ i for i in ars ] ]  

    if (monthly):
        return np.concatenate(   [ i[...,None] 
                                  for ii ,i in enumerate(ars) ], axis=4 ) # concat. specs
    else:
        return np.concatenate(   [ i.sum(axis=3)[...,None] 
                                  for ii ,i in enumerate(ars) ], axis=3 ) # yearly sum. dep*, concat. 

# --------------
# 2.11 - Get Emission of species in Gg
# -------------
def get_emiss( ctm_f=None, spec=None, wd=None, years=None, \
            molec_cm2_s=False, nmonl_m2_d=False, kg_m2_s=False, \
            monthly=False, months=None, s_area=None, res='4x5', debug=False ):
    """ Extract given species emissions from BIOGSRCE category diagnostic
        Back compatibility maintained with PyGChem 0.2.0 """
            
    if debug:
        print 'get_emiss called for >{}<'.format(spec)
    if not isinstance(years, list):
        years = get_gc_years( ctm_f=ctm_f, set_=False, wd=wd )
    if not isinstance(months, list):
        months = get_gc_months( ctm_f=ctm_f, wd=wd)
                                             
    # Adjust to "I Gg / monthly" from "Kg/m2/ s"  
    m_adjust = d_adjust(months, years)

    #  get emissions in  Kg/m2/ s
    # retain compatibility with version 0.2.0
    if pygchem.__version__ == '0.2.0':
        arr = get_gc_data_np( ctm_f, spec, category="BIOGSRCE")[:,:,0,:] 
    else:
        arr = get_GC_output( wd=wd, species=spec, category="BIOGSRCE") 
        res=get_dims4res( r_dims=True, just2D=True )[arr.shape[:2]]

    if not isinstance(s_area, np.ndarray):        
        s_area = get_surface_area(res)  # m2 land map   

    if kg_m2_s:
        arr_ = arr
    else:
        # Kg/m2/ s => Kg/ s
        arr_ = arr*s_area

        # Kg/ s => "kg / monthly" 
        arr_ = arr_ * m_adjust
        # (Gg? - no ) I / month 
        arr_ = arr_*1E3/ species_mass(spec)* float(species_mass('I')) * \
                    float(spec_stoich(spec)) 

    if nmonl_m2_d:
        # convert to / m2
        arr_  =  arr_ / s_area 

        # convert to / day
        arr_  =  arr_ / (365./12.) 

        # convert to nmol 
        arr_  =  ( arr_ / float(species_mass('I') ) ) / float(spec_stoich(spec)) 
        arr_ = arr_*1E9

        if debug:
            print 'get_emiss - 2', arr_.shape

    if molec_cm2_s:
        # from  "I Gg/month" to "I Gg/month/cm/2" #convert to /m2 => cm/2
        arr_  =  arr_ / (s_area *10000.) 

        # convert to / day => hour => hour => minute => sec 
        arr_  =  arr_ / (365./12.) / 24. / 60. / 60. 

        # convert to molecules 
        arr_ = ( arr_ / float(species_mass('I') ) ) /float(spec_stoich(spec)) *\
                        constants('AVG')  

        if debug:
            print 'get_emiss - 3', arr_.shape

    return arr_

# --------
# 2.12 - Get CH4 Lifetime - 2.45E-12, -1775
# --------
def get_CH4_lifetime( ctm_f, wd, years=None, months=None, K=None, \
            vol=None, a_m=None, t_p=None, monthly=False, debug=False ):
    """ Get CH4 lifetime using reaction rate in globchem.dat and OH/CH4
         mean concentrations from geos.log """
    if not isinstance(vol, np.ndarray):
        vol = get_volume_np( ctm_f ) # cm^3 
    if not isinstance(a_m, np.ndarray):
        a_m = get_air_mass_np( ctm_f, dtype=np.float64 )  # Kg

    #  number of moles in box
    moles     = ( a_m *1E3 ) / constants( 'RMM_air') * constants( 'AVG')  

    # Get CH4 from geos.log
    arr = get_CH4_mean( wd )
    #*  kg  per grid box
    arr = arr * moles *  constants( 'AVG') / vol 

    # get loss rate - in kg/s
#    LCH4 = get_gc_data_np(ctm_f, spec='CH4Loss', category='CH4-LOSS') # only works for CH4 sim.

    # calc from reaction loss rate with OH
    K    = get_gc_data_np(ctm_f, spec='TMPU', category='DAO-3D-$')  # K
    OH   = get_OH_mean( wd ) * 1E5   # [molec/cm3] 
    LCH4 = 2.45E-12 * np.exp( (-1775. / K)  )   # cm^3  molec.^-1  s^-1    # rate per grid box 

    return ( 1 / ( np.mean( LCH4 * OH ) )) /60/60/24/365

# ----
# 2.13 - get Land / Water /Ice fraction
# ---- 
def get_LWI(lon, lat, res='4x5',debug=False):
    """ Return LWI for a given lon and lat """
    # lat
    lat=get_gc_lat(lat, res=res)
    lon=get_gc_lon(lon, res=res)
    LWI=get_land_map(res, res=res)
    return LWI[lon,lat,0]
    
# -------------- 
# 2.14 - get_gc_alt  ( KM => box num )   
# -------------
def get_gc_alt(alt):
    """ Return index of nearest altitude (km ) of GEOS-Chem box """

    alt_c = gchemgrid('c_km_geos5_r')
    return find_nearest( alt_c, alt )

# --------------
# 2.15 - species arrary (4D) in v/v to Gg  I/ O3/Species
# ------------- 
def species_v_v_to_Gg(arr, spec, a_m=None, Iodine=True, All =False, \
        Ox=False, ctm_f=None, debug=False):
    """ Convert array of species in v/v to Gg """

    if not isinstance(a_m, np.ndarray):
        a_m     =  get_air_mass_np( ctm_f, dtype=np.float64, debug=debug )  # kg
     #  number of moles in box ( g / g mol-1 (air) )
    moles     = ( a_m *1E3 ) / constants( 'RMM_air')   
    if debug:
        print 'moles.shape', moles.shape, 'arr.shape', arr.shape

    # I mass
    if (( Iodine ) and ( ( not Ox) and (not All) ) ):  
        arr = ( ( arr * moles )  * 127. ) /1E9 * spec_stoich(spec)

    # O3 mass 
    if (( Iodine ) and ( Ox) ):      
        arr = ( ( arr * moles )  * (16.*3.) ) /1E9 

    # O3 mass 
    if (( not Iodine ) and ( All) ):      
        arr = ( ( arr * moles )  * (species_mass( spec  )) ) /1E9 
    return arr

# --------------
# 2.16 - retrive volume for geos chem run ( in cm3 )
# -------------
def get_volume_np(ctm_f=None, box_height=None, s_area=None, res='4x5', \
            wd=None, trop_limit=False, debug=False):
    """ Get grid box volumes for CTM output in cm3 """
    if debug:
        print 'get_volume_np called'

    if not isinstance(box_height, np.ndarray):
        if pygchem.__version__ == '0.2.0':
            box_height =  get_gc_data_np(ctm_f, 'BXHEIGHT',\
                category="BXHGHT-$", debug=debug)  # ( m )
        else:
            print '!'*100
            box_height = get_GC_output( wd, ['BXHGHT_S__BXHEIGHT'], \
                debug=debug )

            # Gotcha for single (None) time dimension output:
            none_time_dim_shapes = (72, 46, 47), (72, 46, 38) 
            if any( [ box_height.shape == i for i in none_time_dim_shapes ]) :
                box_height= box_height[...,None]

    if not isinstance(s_area, np.ndarray):
        # m2 land map                                                
        s_area = get_surface_area(res=res)[...,None] 
        if debug:
            print '**'*100, s_area.shape

    # Gotcha for 2D s_area array 
    if len( s_area.shape) == 2:
        s_area = s_area[..., None, None]

    # m3 ( *1E6 ) => cm3 
    volume = ( box_height * s_area ) *1E6  

    if trop_limit:
        volume =  volume[...,:38,:]
    
    if debug:
        print  'box_height' , box_height.shape , 's_area' , s_area.shape,\
             'volume', volume.shape
    return volume

# --------------
# 2.17 - Test is grid box for a given lat and lon is over water
# ------------- 
def loc_is_water_grid_box( lat, lon, res='4x5' ):
    """ Return Boolean for if grid box is water or not """

    # Load masked array, where all ocean is non-masked
    o_mask = np.ma.equal( get_land_map(res=res)[...,0], 0) # Just look over oceans (use masked array)

    # Loop Lat and Lon
    for n, l in enumerate( lat ):
        glat, glon = get_gc_lat( l ), get_gc_lon( lon[n] )
       
        # check aagainst mask and store as boolean
        if o_mask[glon, glat] :
            bool = 'F'
        else:
            bool = 'T'
        try:
            booll.append( bool )
        except:
            booll = [ bool ]
    return booll

# --------------  
# 2.18 - dry dep for given spec
# -------------   
def spec_dep(ctm_f=None, wd=None, spec='O3', s_area=None, months=None, \
                years=None, res='4x5', vol=None, year_eq=True, debug=False, \
                trop_limit=True, Iodine=False):
    """ Get array of dry deposition values for a given species """

    if debug:
        print 'spec dep called'

    # Get surface area if not provided
    if not isinstance(s_area, np.ndarray):
        s_area =  get_surface_area(res)  # m2 land map                                                

    # Extract dry dep flux in  [molec/cm2/s]
    if pygchem.__version__ == '0.02':
        df = get_gc_data_np( ctm_f, spec=spec+'df', category='DRYD-FLX', \
            debug=debug) 
    else:
        df = get_GC_output( wd, category='DRYD-FLX', species=spec+'df' )

    if debug:
        print '*'*10, [ ( i.shape, np.sum(i), np.mean(i)) for i in [df] ], \
            len(df)

    #  convert to Gg "Ox" (Gg I /s)
    df = molec_cm2_s_2_Gg_Ox_np( df, spec, s_area=s_area, Iodine=Iodine, \
                res=res, debug=debug ) 

    if debug:
        print '0'*20, [ ( i.shape, np.sum(i), np.mean(i)) for i in [df] ]
    if ( ( not months) or (not years) ):
        years = get_gc_years( ctm_f, wd=wd, set_=False )
        months = get_gc_months( ctm_f, wd=wd )
        if debug:
            print years, months

    # adjust time dimension
    day_adjust = d_adjust( months, years)
    df = np.multiply(df, day_adjust)

    return df

# --------------  
# 2.19 - Land map - REDUNDENT
# -------------   

# --------------
# 2.20 - Gg Ox yr^-1 from individual PORL-L$ rxns ([molec/cm3/s] )
# -------------
def molec_cm3_s_2_Gg_Ox_np(arr, rxn=None, vol=None, ctm_f=None, \
            Iodine=False, I=False, IO=False, months=None, years=None, wd=None, \
            year_eq=True, month_eq=False, spec=None, res='4x5',debug=False ):
    """ Convert species/tag prod/loss from molec/cm3/s] to Gg (Ox) yr^-1.
        This function was originally used to process diagnostic outputs
        from PORL-L$  """ 

    if debug:
        print ' molec_cm3_s_2_Gg_Ox_np  called'
    if not isinstance(vol, np.ndarray):
        vol = get_volume_np( ctm_f=ctm_f, wd=wd )
        print 'WARNING: extracting volume online - inefficent'

    # --- Process conversion
    if (Iodine):
        if debug:
            print rxn , arr.shape, vol.shape 

        #  [molec/cm3/s] => * vol / moles * RAM *I in species / 1E9 
        arr = ( arr * vol[:,:,:38,:] ) /  constants( 'AVG')* (127.) /1E9 * \
            spec_stoich( rxn, I=I, IO=IO ) 

    else: # ELSE return Gg Ox yr^-s
        if ( rxn ==  'CH4' ):
            RMM =  species_mass(rxn) 
        else:
            RMM =  species_mass('O3')             
        if not isinstance( spec, type(None) ):
            RMM = species_mass( spec  )
        arr = ( arr * vol[:,:,:38,:] ) / constants( 'AVG') * RMM /1E9

    # --- Process time
    if year_eq: # convert /second to / year
        arr = arr * 60*60*24*365
    if month_eq:
        day_adjust = d_adjust( months, years)
        arr = arr * day_adjust

    if debug:
        print 'arr', arr.shape , 'vol', vol.shape
    return arr

# --------------
# 2.21 - Gg Ox yr^-1 from induvial spec (2D array - dry dep)  ( ([molec/cm2/s] ) to Gg X (or Gg X / s if not year_eq. )
# -------------
def molec_cm2_s_2_Gg_Ox_np( arr, spec='O3', s_area=None, ctm_f=None, \
            Iodine=False, res='4x5', year_eq=False, fix_Kludge=True, \
            debug=False):
    """ Convert 2D depositional array from [molec/cm2/s] to Gg Ox yr^-1 """
    
    if debug:
        print ' molec_cm2_s_2_Gg_Ox_np  called'

    # Kludge and anti-kludge (  REMOVE THIS once 
    s_area = s_area[...,None]  # tmp bug fix (loads all output surface areas, only requires a single time dimensions as does not chnage... )
    if fix_Kludge:
        s_area = s_area[...,0]
        
    if debug:
        print '_'*10, arr.shape, s_area.shape
    if ( Iodine ):
        # [molec/cm2/s] * cm2 (m*1E4) /  moles * RMM I /1E9   ( convert Gg I /s )
        arr  = ( arr * (s_area*1E4) ) / constants( 'AVG') * (127.) * \
            spec_stoich(spec) /1E9   
    elif spec == 'O3':
        arr  = ( arr * (s_area*1E4) )  / constants( 'AVG')* \
            (species_mass(spec)) /1E9  * Ox_in_species( spec)
    else:
        arr  = ( arr * (s_area*1E4) )  / constants( 'AVG') * \
            (species_mass(spec)) /1E9  
    if (year_eq):
        print "giving 'year_eq' "
        arr = arr* 60*60*24*365  

    if debug:
        print 'arr', arr.shape
    return arr

# --------------
# 2.22 - Get Dry dep
# -------------
def get_dry_dep( ctm_f, months=None, years=None, s_area=None, \
        res='4x5', debug=False):
    """ Get dry deposition from a given ctm.bpch file
        REDUNDENT - not compatible with PyGChem >0.3.0 """
    if not isinstance(months, list):
        months = get_gc_months( ctm_f )
    if not isinstance(years, list):
        years  = get_gc_years( ctm_f, set_=False )
    if not isinstance(s_area, np.ndarray):
        s_area      =  get_surface_area(res)  # m2 land map                                                
    d_dep = GC_var('d_dep')
#    Ox    = GC_var('Ox')
    Ox = GC_var('Ox_spec')
    Oxdf    = [i +'df' for i in Ox ]
    ind =[ n for n,i in enumerate ( Oxdf  )if any([(i==ii) for ii  in ['NO3df', 'HNO4df', 'BrOdf', 'BrNO2df', 'MPNdf', 'IOdf', 'OIOdf'] ] ) ]
    if debug:
        print ind, Ox, Oxdf
    [ l.pop(i) for i in ind[::-1]  for l in Ox, Oxdf ]

    print [ (len(i) , i) for i in [Oxdf ]]
    
    # Ox wet dep - get specs for wet dep within category, then load their data and conatenate    
    Oxdf = Ox_vars_in_cat_4_run(ctm_f, Oxdf, d_dep)[0]
    print [ (len(i) , i) for i in [Oxdf ]]

    
    ars_Ox_dep =  [ molec_cm2_s_2_Gg_Ox_np(array, Ox[i], s_area, ctm_f)                        # Ox dry dep  [molec/cm2/s]  => Gg Ox
                    for i, array in enumerate([ get_gc_data_np(ctm_f, rxn, \
                        category=d_dep[0] ) for rxn in Oxdf ]) ]

    day_adjust = d_adjust( months, years)
    ars_Ox_dep = [day_adjust * arr for arr in ars_Ox_dep ]

    return  np.sum( np.concatenate( [ (i[:,:,:,:]* \
            Ox_in_species(Ox[ii]))[:,:,:,:,np.newaxis] \
            for ii ,i in enumerate(ars_Ox_dep) ], axis=4 ), axis=4 )   

# --------------  
# 2.23 - Get run data for a given GC run - 'Mean surface (Oceanic) IO', 'Chem Ox loss', 'Chem Ox prod', 'O$_{3}$ Burden','DU Column', 'OH Mean Conc'
# -------------   
def get_GC_run_stats( wd, Iodine=True, HOx_weight=False, trop_limit=True, debug=False ):
    """ This is version of get_GC_run_stats that works with monthly outputted 
        data """

    # Open file(s)
#    ctms, res = wd2ctms( wd )
                                 
    # --- Get vars need for calcs
    # Get species time Tropopause diagnostic
    t_ps, res = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], \
        trop_limit=trop_limit, r_res=True ) 
  
    # Get surface area
    s_area = get_surface_area(res=res)[...,0] # m2 land map           
  
    # Get air mass
    a_m = get_GC_output( wd, vars=['BXHGHT_S__AD'], trop_limit=trop_limit, \
                    dtype=np.float64) 

    # Get volume
    vol  = get_volume_np( wd=wd, s_area=s_area[:,:,None, None], res=res )

    # Get O3 v/v array
    O3_arr  = get_GC_output( wd, species='O3' )[:,:,:38,:]

    # Get molecs in troposphere
#    molecs = a_m*1E3/constants( 'RMM_air')*constants('AVG') 

    # Get Months and years in model output
    months = get_gc_months( wd=wd ) 
    years = get_gc_years( wd=wd )

    # get model version    
    ver = iGEOSChem_ver( wd )

#    print months, years, ver
#    print [i.shape for i in  t_ps, a_m, vol, O3_arr, molecs  ]

    if any( [(ver ==i) for i in  '1.3' ,'1.4' ,'1.5', '1.6', '1.7' ]):
        vars =  [  'p_l', 'd_dep_specs', 'w_dep_specs', 'Iy','IOx', 'LIOx' ]
    else:
        print 'FAIL'
        sys.exit(0)

    p_l, d_dep_specs, w_dep_specs, Iy, IOx, LIOx = [ GC_var(i)  for i in vars ]

    # Get surface IO conc
    IO_arr =  get_GC_output( wd, species='IO' )
    sur_IO = np.ma.mean(  np.logical_not( ocean_mask(res)[:,:,0] ) *  \
        np.mean( IO_arr[:,:,0,:], axis=2) ) *1E12

    # get POx and LOx
    POx, LOx = get_POxLOx( wd=wd, vol=vol, t_p=t_ps, ver=ver, debug=debug)

    # get O3 burden
    O3_bud = get_O3_burden( a_m=a_m[:,:,:38,:], t_p=t_ps, O3_arr=O3_arr, \
        debug=debug ) /1E3

    # Get O3 dep
    O3_dep = spec_dep( spec='O3', wd=wd, s_area=s_area, months=months, \
        years=years, vol=vol, Iodine=False )/1E3 

    # get DU column
    DU_O3  = get_DU_mean( s_area=s_area, a_m=a_m[:,:,:38,:], t_p=t_ps, \
        O3_arr=O3_arr, debug=debug)

    # OH mean Conc
#    OH_mean = get_OH_mean( wd, debug=debug )

    # Get OH and HO2
    Hl = get_OH_HO2( wd=wd, t_p=t_ps, a_m=a_m, HOx_weight=HOx_weight, \
        vol=vol, res=res ) 

    # Get Loss routes for Iy and cal Iy lifetime Gg / Gg/year => days  (inc. wet & dry dep)
    # Kludge - tmp. rm for v10
#    d_dep = np.concatenate( [ [spec_dep(ctm, spec=i,months=months, years=years,                        \
#                                     s_area=s_area[...,None], debug=debug)  for i in                                                \
#                                    [ i.split('df')[0] for i in d_dep_specs[:-1] ] ]  or ctm in ctms ],  axis=-1 )
    d_dep =  np.zeros( get_dims4res(res=res) ) 
    # Kludge
    w_dep =  np.zeros( get_dims4res(res=res) )     
#    w_dep = np.concatenate( [ get_wet_dep( ctm, months=months, years=years, vol=vol,                \
#                                             Iodine=Iodine, debug=debug) for ctm in ctms ], axis=3)
#    print [ np.array(i).shape for i in d_dep, w_dep ] 
    
    try:
        Iy_loss = [ np.concatenate( [ get_pl_in_Gg(ctm, i, Iodine=Iodine, \
            res=res, ver=ver, debug=debug  ) for ctm in ctms ], axis=3)  \
            for i in [ ['L_Iy'] ]]
        Iy_burdens  = np.mean( [ np.concatenate( [get_trop_burden( ctm, i, \
            Iodine=Iodine, res=res,  debug=debug ) \
            for ctm in ctms ], axis=3) for i in Iy], axis =4 )
        ars =  [np.sum(i) for i in [ Iy_burdens, Iy_loss, d_dep, w_dep ] ]
        Iy_lifetime = ( ars[0] /( np.sum(ars[1:]) ) ) *365
    except:
        Iy_lifetime, Iy_burdens = 0, 0

    # Get Loss routes for IOx and cal IOx lifetime Gg / Gg/year => *365*24 => ms 
    try:
        IOx_loss = [  np.concatenate( [ get_pl_in_Gg(ctm, i, Iodine=Iodine, \
            res=res, ver=ver, debug=debug ) for ctm in ctms ], axis=3)  \
            for i in [ LIOx  ]] #
        burdens = np.mean( [ np.concatenate(  [ get_trop_burden( ctm, spec, \
             a_m=a_m, t_p=t_p, Iodine=Iodine, res=res, debug=debug ) \
             for n, ctm in enumerate( ctms ) ], axis=3) \
            for spec in IOx ], axis=4)
        ars = [ np.sum(i) for i in [ burdens, IOx_loss ]]
        IOx_lifetime = ( ars[0] /(np.sum(ars[1:])) ) *365*24*60

    except:
        IOx_lifetime = 0

    # Get CH4 Lifetime
#    CH4_lifetime =  np.mean( [ get_CH4_lifetime( ctm, wd, vol=vol[:,:,:38,:], \
#        a_m=a_m[:,:,:38,:], t_p=t_p ) for ctm in ctms ] )
    CH4_lifetime=0

    # remove files from memory
#    del ctms

    #'Mean surface (Oceanic) IO', 'Chem Ox loss', 'Chem Ox prod', 'O$_{3}$ Burden','DU Column', 'OH Mean Conc', 
    #                                'HO2 Mean Conc', 'HO2/OH Mean Conc', 'CH4 lifetime'
    return [ sur_IO, POx, LOx, POx -LOx, np.sum( O3_bud ), \
        np.sum( O3_dep ), np.mean( DU_O3 ), Hl[0],Hl[1] , Hl[1]/Hl[0],   \
        Iy_lifetime, np.mean(IOx_lifetime), np.sum(Iy_burdens), CH4_lifetime ]


# --------------
# 2.24 - Get DU mean value
# -------------
def get_DU_mean(s_area=None, a_m=None, t_p=None, O3_arr=None, \
            ctm_f=False, area_weight=True, res='4x5', debug=False ):
    """ Get mean DU value weighed by grid box array """

    if not isinstance(s_area, np.ndarray):
        print "Extracting s_area in 'get_DU_mean'" 
        s_area  =  get_surface_area( res )[...,0]  # m2 land map
    if not isinstance(a_m, np.ndarray):
        print "Extracting a_m in 'get_DU_mean'" 
        a_m = get_air_mass_np( ctm_f, debug=debug )
    if not isinstance(t_p, np.ndarray):
        print "Extracting t_p in 'get_DU_mean'" 
        t_p = get_gc_data_np( ctm_f, spec='TIMETROP', category='TIME-TPS', \
                        debug=debug)
    if not isinstance(O3_arr, np.ndarray):
        print "Extracting arr in 'get_DU_mean'"     
        arr  = get_gc_data_np( ctm_f, 'O3', debug=debug )
    else:
        arr = O3_arr
    if debug:
        print [ (i.shape, i.mean()) for i in s_area, a_m, t_p, arr ]

    # Get molec. of air to get moles of spec
    # molecs air
    molecs_air  =  (a_m*1E3) / constants( 'RMM_air') * constants('AVG')   
    # molecules O3
    arr = np.ma.array(arr) * molecs_air     

    # average of time within the troposphere
     # cut columns to tropopause
    arr = np.sum( arr[:,:,:38,:]* t_p, axis=2 )
    # average over time
    arr = arr.mean(axis=2) 

    # convert to DU
    arr = arr / s_area # adjust to per area
    arr  = arr / constants('mol2DU') # convert to DU       
    if area_weight:
        arr = np.sum(arr * s_area)/np.sum(s_area)    # weight by area
    return arr

# --------------
# 2.25 - Get O3 Burden
# -------------
def get_O3_burden(wd=None, spec='O3', a_m=None, t_p=None, O3_arr=None, \
                ctm_f=False, trop_limit=True,all_data=False,  debug=False ):
    """ Get O3 burden in CTM output 
        REDUNDENT? - this function is replicated in the get_burden function """

    if not isinstance(a_m, np.ndarray):
        print 'Getting a_m in get_O3_burden'
        if pygchem.__version__ == '0.2.0':
            a_m = get_air_mass_np( ctm_f, debug=debug )[:,:,:38,:]

        # use update pygchem
        else:
            a_m = get_GC_output( wd, vars=['BXHGHT_S__AD'], 
                trop_limit=trop_limit) 

    if not isinstance(t_p, np.ndarray):
        print 'Getting t_p in get_O3_burden'
        if pygchem.__version__ == '0.2.0':
            t_p = get_gc_data_np( ctm_f, spec='TIMETROP', \
                category='TIME-TPS', debug=debug)[:,:,:38,:]

        # use update pygchem
        else:
            t_p = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], \
                trop_limit=trop_limit ) 
            
    if not isinstance(O3_arr, np.ndarray):
        print 'Getting O3_arr in get_O3_burden'
        if pygchem.__version__ == '0.2.0':
            ar = get_gc_data_np( ctm_f, 'O3', debug=debug )[:,:,:38,:]
        else:
            ar  = get_GC_output( wd, species='O3' )[:,:,:38,:]
            
    else:
        ar  =  O3_arr[:,:,:38,:]

    # v/v * (mass total of air (kg)/ 1E3 (converted kg to g))  = moles of tracer
    ar = ar * ( a_m[:,:,:38,:]*1E3 / constants( 'RMM_air') ) 

    # convert moles to mass (* RMM) , then to Gg 
    ar = ar * species_mass(spec) /1E9             
    ar = ar[:,:,:38,:] * t_p[:,:,:38,:]
    if (all_data):
        return ar
    else:
        return ar.mean(axis=3 )

# --------------
# 2.26 - Get Prod / loss for O3
# -------------
def get_POxLOx( ctms=None, vol=None, all_data=False, t_p=None, ver='1.6', \
    wd=None, debug=False):
    """ Get production and loss terms for O3 from prod/loss diagnostic """

    specs = GC_var('POxLOx')

    if debug:
        print ver 

    # get prod/loss in [molec/cm3/s]
    if pygchem.__version__ ==  '0.2.0':
        arrs = [ np.concatenate( [get_gc_data_np( ctm, spec=PLO3_to_PD(spec, \
            ver=ver, fp=True),category="PORL-L=$", debug=debug) \
            for ctm in ctms ],  axis=3 )  for spec in specs  ]
    else:
        arrs = get_GC_output( wd, vars=['PORL_L_S__'+i for i in specs] )
        arrs = [ arrs[i,...] for i in range( len(specs) ) ]

    if (all_data):
        # [molec/cm3/s] => Gg Ox / month 
        months = [ get_gc_months( ctm) for ctm in ctms ]
        years =  [ get_gc_years( ctm, set_=False ) for ctm in ctms ]
        months = [j for i in months for j in i ]
        years = [j for i in years for j in i ]        
        arrs = [ molec_cm3_s_2_Gg_Ox_np(arr, specs[i], vol=vol, months=months, \
                    years=years, debug=debug, month_eq=True, year_eq=False) \
                    for i, arr in enumerate(arrs) ] 
        return [ arrs[i] for i in range(len(specs )) ]  # Gg

    else:
        # [molec/cm3/s] => Gg Ox / yr
        arrs = [ molec_cm3_s_2_Gg_Ox_np(arr, specs[i], vol=vol, \
            debug=debug) for i, arr in enumerate(arrs) ] 
        arrs = [ (arr*t_p).mean(axis=3) for arr in arrs ] # get yearly mean

        return [ int(np.sum( np.ma.masked_invalid( arrs[i] ) )/1E3)  \
            for i in range(len(specs )) ] # Tg

# --------------
# 2.27 - Get Transport fluxes
# -------------
def get_tran_flux( ctm_f, spec='O3', years=None, months=None, \
            s_area=None, all_data=False, debug=False):
    """ Get transport fluxes from UPFLX... diagnostics """
    
    if not isinstance(months, list):
        months = get_gc_months( ctm_f )
    if not isinstance(years, list):
        years  = get_gc_years( ctm_f, set_=False )
    f_var = GC_var('f_var')
    Ox    = GC_var('Ox_spec')    

    ars_f_var  =  [ [ get_gc_data_np(ctm_f, spec, category=var )*1E3  \
            for spec in Ox ] \
            for var in f_var ]   

    # list Ox specs, concat. and sum of rxns
    ars_f_var  = [ np.sum( np.concatenate( [ i[...,None]  \
                for ii ,i in enumerate(arrs_l) ], axis=4 ), axis=4 )  \
                for arrs_l in ars_f_var ]

    # Subtract values from adjoint box to get net flux
    if debug:
        print [(np.sum(i), i.shape) for i in ars_f_var ]
    ars_f_var[0]  =  [ ars_ - np.roll(ars_, 1, axis=0) \
        for ars_ in  [ ars_f_var[0]]  ][0] 
    ars_f_var[1]  =  [ ars_ - np.roll(ars_, 1, axis=1) \
        for ars_ in  [ ars_f_var[1]]  ][0] 
    ars_f_var[2]  =  [ ars_ - np.roll(ars_, 1, axis=2) \
        for ars_ in  [ ars_f_var[2]]  ][0] 
    if debug:
        print [(np.sum(i), i.shape) for i in ars_f_var ]

    day_adjust = d_adjust( months, years)
    return [day_adjust * arr/1E9 for arr in ars_f_var ]

# --------------
# 2.28 - Get wet dep
# -------------
def get_wet_dep( ctm_f=None, months=None, years=None, vol=None, \
            scale=1E9, s_area=None, res='4x5', wd=None, specs=None, \
            Iodine=False, all_wdep=False, sep_rxn=False, debug=False):
    """ Extract wet deposition for given species  """

    if not isinstance(months, list):
        months = get_gc_months( ctm_f, wd=wd )
    if not isinstance(years, list):
        years  = get_gc_years( ctm_f=ctm_f, wd=wd, set_=False )
    if not isinstance(vol, np.ndarray):
        vol = get_volume_np( ctm_f, s_area=s_area, wd=wd, res=res,\
             debug=debug )
    w_dep = GC_var('w_dep')
    
    
    if Iodine and isinstance( specs, type(None) ):
        if (all_wdep):
            specs = GC_var('w_dep_specs')
        else:
            specs = GC_var('w_dep_specs')[:-1] # skip AERI - Kludge

        # No longer nessesary as all wet is constant for the last few versions. 
#        specs      = [ [i for i in specs if  i in vars_in_w_dep ] for vars_in_w_dep in
#                       [ [ diag.name for diag in diagnostics ] for diagnostics in [ ctm_f.filter(category=i) for i in w_dep ] ]    ] 

#    elif not isinstance(species, type(None) ):
#        specs = [ species, species ]

    # Ox wet dep - get specs for wet dep within category, then load their data and conatenate    
    # REDUNDENT - this is fo v9-01-03 ( now using v10 as default + v9-2 )
#    else:  
#        Ox    = GC_var('Ox')
#        specs = Ox_vars_in_cat_4_run(ctm_f, Ox, w_dep)
    if debug:
        print specs , len(specs)

     # [kg/s] => g / s    of I equiv.
    if Iodine:

        # retain back compatibility
        if pygchem.__version__ == '0.2.0':
            dep_w  = [ [ get_gc_data_np(ctm_f, spec ,category=var, \
                debug=debug )*1E3 /species_mass(spec)*  \
                species_mass('I')*spec_stoich(spec)    \
                for spec in specs[ii] ]  \
                for ii,var in enumerate(w_dep)]  

        else:
            # Frontal rain?
            WETDLS_S__ = get_GC_output( wd=wd, vars=['WETDLS_S__'+i \
                for i in specs ] )

            # get convective scavegning 
            WETDCV_S__ = get_GC_output( wd=wd, vars=['WETDCV_S__'+i \
                for i in specs ] )
#            if debug:
            print [ i.shape for i in WETDLS_S__, WETDCV_S__  ]


            # convert to g/ s I equiv.  + conbine two lists
            WETDLS_S__ = [ i*1E3/species_mass( specs[n] )* species_mass('I')*spec_stoich( specs[n] )  
                    for n, i in enumerate( [WETDLS_S__[ii,...] for ii in range(len(specs)) ]  ) ]
            WETDCV_S__  = [ i*1E3/species_mass( specs[n] )* species_mass('I')*spec_stoich( specs[n] )  
                    for n, i in enumerate( [WETDCV_S__ [ii,...] for ii in range(len(specs)) ]  ) ]

            dep_w =  WETDLS_S__+ WETDCV_S__ 
            print [ i.shape for i in dep_w  ]
#            sys.exit()
                            
    # [kg/s] => g / s    # 
    else: 
        # retain back compatibility
        if pygchem.__version__ == '0.2.0':
            dep_w  = [ [ get_gc_data_np(ctm_f, spec ,category=var, debug=debug )*1E3  
                for spec in specs[ii] ]  for ii,var in enumerate(w_dep)]  

        else:
            # Frontal rain?
            WETDLS_S__ = get_GC_output( wd=wd, vars=['WETDLS_S__'+i \
                for i in specs ] )*1E3

            # get convective scavegning 
            WETDCV_S__ = get_GC_output( wd=wd, vars=['WETDCV_S__'+i \
                for i in specs ] )*1E3

            dep_w = [ WETDLS_S__, WETDCV_S__ ]

        # concatenate list of lists
        dep_w = [j for i in dep_w for j in i]  

    # adjust to monthly values...
    m_adjust = d_adjust( months, years) # => Gg / month 
    dep_w = [m_adjust * arr for arr in dep_w ]

    # list wet dep rxns, concat. and sum of rxns
    if (sep_rxn): 
        return np.concatenate( [ i[np.newaxis,:,:,:38,:] /scale \
                        for ii ,i in enumerate(dep_w) ], axis=0 )

    # list wet dep rxns and concat. 
    else:   
        return np.sum( np.concatenate( [ i[:,:,:38,:,None] /scale \
                        for ii ,i in enumerate(dep_w) ], axis=4 ), axis=4 )

# --------------
# 2.29 - BL mixing
# -------------
def get_BL_mix( ctm_f, months=None, years=None, spec='O3', debug=False):
    """ Get Boundary layer mixing """

    if not isinstance(months, list):
        months = get_gc_months( ctm_f )
    if not isinstance(years, list):
        years  = get_gc_years( ctm_f, set_=False )
    BL_m    = GC_var('BL_m')[0]
    Ox      = GC_var('Ox_spec')

    #   kg/s] =>  g / s    # BL mixing upwards of tracer
    ars_BL_m     = [ get_gc_data_np(ctm_f, spec=i, category=BL_m )*1E3 \
        for i in Ox ] 

    # list Ox specs, concat. and sum of rxns                                            
    ars_BL_m = np.sum( np.concatenate( [ i[:,:,:38,:,None]  \
                for ii ,i in enumerate(ars_BL_m) ], axis=4 ), axis=4 )

    ars_BL_m = [ ars_ - np.roll(ars_, 1, axis=2) \
        for ars_ in  [ ars_BL_m ]   ][0]

    day_adjust   = d_adjust( months, years)
    return day_adjust * ars_BL_m /1E9

# --------------
# 2.30 - Cloud Flux
# -------------
def get_cld_f( ctm_f, months=None, years=None, spec='O3',  debug=False ):
    """ Extract cloud flux"""

    if not isinstance(months, list):
        months = get_gc_months( ctm_f )
    if not isinstance(years, list):
        years  = get_gc_years( ctm_f, set_=False )
    Cld_flx = GC_var('Cld_flx')[0]
    Ox      = GC_var('Ox_spec')

    # kg/s =>  g / s    # UPWARD MASS FLUX DUE TO WET CONVECTION
    ars_cld_f    = [ get_gc_data_np(ctm_f, spec=i, category=Cld_flx )*1E3 \
        for i in Ox ]  

    # list Ox specs, concat. and sum of rxns                                     
    ars_cld_f    = np.sum( np.concatenate( [ i[:,:,:38,:,None]  \
                        for ii ,i in enumerate(ars_cld_f) ], axis=4 ), axis=4 )

    ars_cld_f = [ ars_ - np.roll(ars_, 1, axis=2) for ars_ in [ars_cld_f] ][0]
    day_adjust   = d_adjust( months, years)
    return day_adjust *  ars_cld_f /1E9 

# --------
# 2.31 - Vol weight array, sum, mean
# --------
def vol_weight( arr, vol=None, ctm_f=None ):
    """ Volume weight array """
    if not isinstance(vol, np.ndarray):
        vol = get_volume_np( ctm_f ) # cm^3 
    if (sum):
        return np.sum( arr/vol) *np.sum(vol) 

# ----
#  3.32 - returns  (lat, lon, alt (press), timezone (UTC) ) for a given site
# ----
# moved to vars mod.

# --------------
# 2.34 -  takes indices or generates indices for a given set of sites, these are then substract from all arrays given
# -------------
def rm_data_at_locs( data, lat, lon, alt, below_km=2, indiceslocs=None, \
            debug=False):
    """ This is a slow method quickly cobbled together, but is only used once, 
            if called more than once if needs to re-written for better speed """

    # Get indices of obs. on grid. for obs to fit to ...
    glon, glat, galt = get_latlonalt4res( nest='high res global', centre=False,\
        debug=debug )

    # get grids and obs,. indicies
    indices_list = obs2grid( glon=glon, glat=glat, galt=galt )
    if debug:
        print indices_list

    # check for these locations with given spatial values
    # for BL layer remove values near atolls of CHUUK, PILAU, GUAM
    coastal = []
    for n, lat_ in  enumerate( lat ):
        if alt[n] < below_km:
            loc  = get_xy(  lon[n], lat_ , glon, glat )
            if any( [ loc == i for i in indices_list ] ):
                if debug:
                    print '-'*20
                    print any( [ loc == i for i in indices_list ] ), n, lat_, \
                         lon[n], alt[n], loc, indices_list
                coastal += [n]
                if debug:
                    print '-'*20

    print [len(i) for i in coastal, data, lat, lon, alt]    
    arrs = data, lat, lon, alt 
    arrs = [ np.delete(i, coastal) for i in arrs ]
    data, lat, lon, alt  = arrs
    print [ i.shape for i in arrs ]
    print [ i.shape for i in data, lat, lon, alt ] 

    # remove values at given sties and return data
    return data, lat, lon, alt

# --------------
# 2.35 - Get Tropospheric Ox loss routes as list of arrays, + model resolution (res)
# -------------
def get_trop_Ox_loss( wd, pl_dict=None,  spec_l=None, ver='1.6' ,   \
                trop_limit=True, units='Gg Ox/yr', debug=False): 
    """Get Ox loss individually for all routes and return in "Gg Ox / yr".
        Can also return in  """

    # Create Variable lists/dictionaries is not defined.
    if isinstance( pl_dict, type(None)):
        pl_dict = get_pl_dict( wd, spec='LOX', ver=ver, rmx2=True )
    if isinstance( spec_l, type(None)):
        spec_l =  pl_dict.keys()

    # ---  create arrays from ctm files
    if pygchem.__version__ == '0.2.0':

        #  Find and open CTM.BPCH files 
        ctms, res = wd2ctms( wd )

        # Get time in troposphere and volume
        t_p =  np.concatenate( [ get_gc_data_np( ctm, spec='TIMETROP', \
            category='TIME-TPS') 
                        for ctm in ctms ] , axis=3 )
        vol = np.concatenate( [ np.around( get_volume_np( ctm, res=res ) )    \
                        for ctm in ctms ], axis=3 )[:,:,:38,:] 

        # Get prod loss in  [molec/cm3/s]
        ars = [ np.concatenate( [ \
        get_gc_data_np(ctm, PLO3_to_PD(s, fp=True,wd=wd, ver=ver ), \
                    "PORL-L=$") for ctm in ctms ], axis=3 )   for s in spec_l ]   
    else:
        # Get time in troposphere and volume
        t_p, res = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], \
            trop_limit=trop_limit, r_res=True ) 
        vol = get_volume_np( wd=wd, trop_limit=trop_limit )

        # Get prod loss in  [molec/cm3/s]
        ars = get_GC_output( wd, vars=[ 'PORL_L_S__'+ \
            PLO3_to_PD(i, fp=True,wd=wd, ver=ver )  for i in spec_l ] )

    # convert species arrays [molec/cm3/s] into Gg Ox / yr
    if units == 'Gg Ox/yr':
        ars = [ molec_cm3_s_2_Gg_Ox_np(ar, spec_l[i], vol=vol, debug=debug) \
                        for  i, ar in enumerate(ars) ]
    else:
        print 'units = {}'.format(units)
                                        
    # adjust for Ox difference within prod loss tracers
    ars = [ ar*pl_dict[spec_l[i] ][-1]  for i, ar in enumerate(ars) ]

    # Only consider troposphere
    ars = [i*t_p for i in ars]

    return ars, res

# --------------
# 2.36 - Split Tropospheric Ox loss routes 
# -------------
def split_Ox_loss_by_fam( wd, arr, r_t=None, pl_dict=None, \
            NOy_as_HOx=True, as_numpy=True, ver='1.6', debug=False ):
    """ Takes a n dimension array ( typically: 2D/4D) array, then splits this
    into a list of single arrays for each Ox family  """

    # Create Variable lists/dictionaries is not defined.
    if isinstance( pl_dict, type(None)):
        pl_dict = get_pl_dict( wd, spec='LOX', ver=ver, rmx2=True )
    if isinstance( r_t, type(None)):
        if NOy_as_HOx:
            r_t  = [ 'Photolysis','HOx ','Bromine', 'Iodine' ] #+ ['Total']
        else:
            r_t  = [ 'Photolysis','HOx ','NOy' ,'Bromine', 'Iodine' ] #+ ['Total']

    # generate indicies map to split Ox loss by route.   
    spec_l =  pl_dict.keys()
    if debug:
        print len( spec_l)
    r_, fam, spec_l = get_indicies_4_fam( spec_l, fam=True, IO_BrOx2=True,\
        rtnspecs=True )
    if debug:
        print len( spec_l)

    # split list of arrays and combine into Ox loss families
    ars = [ np.array( [ arr[r,...]  for r in r_[f] ] ).sum(axis=0) \
        for f in range(len(r_t)) ]                           

    # Kludge - Place NOy at front of list to "correct" order
    if not NOy_as_HOx:
        if debug:
            print 1, [ i.shape for i in ars ], type(ars), r_t.index( 'NOy' ),\
                            np.ma.array( ars).sum()
        ars = [ ars[ r_t.index( 'NOy' ) ]  ] + ars
        if debug:
            print 2, [ i.shape for i in ars ], type(ars), r_t.index( 'NOy' ),\
                            np.ma.array( ars).sum()
        ars.pop( r_t.index( 'NOy' ) +1 )
    if as_numpy:
        ars = np.array( ars )
        if debug:
            print 3, [ i.shape for i in [ ars ] ], np.ma.array( ars).sum()

    # manually rename plot labels
    r_t = GC_var('r_tn' )

    return ars, r_t

# --------------
# 2.37 - Return only dry/wet dep variaibles aloocated within model output arrays
# -------------       
def var_in_dep(ctm_f, specs, d_=False, w_=False, d_dep=False, \
        w_dep=False, debug=False):
    """ Check if deposition variable is in ctm.bpch  - REDUNDENT. """

    if ( ( not d_) or ( not  w_ ) ):
        d_dep, w_dep  = [ GC_var(i)  for i in [ 'd_dep','w_dep' ] ] # which species
    if ( w_ ):
        specs      = [ [i for i in specs if  i in vars_ ] for vars_ in
                         [ [ diag.name for diag in diagnostics ] for diagnostics in [ ctm_f.filter(category=i) for i in w_dep ] ]    ]
    if ( d_ ):
        specs      = [ [i for i in specs if  i in vars_ ] for vars_ in
                         [ [ diag.name for diag in diagnostics ] for diagnostics in [ ctm_f.filter(category=i) for i in d_dep ] ]    ][0]
    if debug:
        print specs
    return specs

# --------------
# 2.38 - split 4D ( lon, lat, alt, time) output by season
# -------------       
def split_4D_array_into_seasons( arr, annual_plus_seasons=True, \
            debug=False ):
    """ split 4D ( lon, lat, alt, time) output by season"""
    
    if debug:
        print arr.shape 
    
    if annual_plus_seasons:
        seasons  = ['Annual', 'DJF', 'MAM', 'JJA', 'SON']
        # assume calender month order 
        # ( this can be automated using get GC datetime )
        indices = [range(0,12), [11, 0, 1], [2, 3, 4], [5, 6, 7], [8, 9, 10] ]
             
    else:
        seasons  = ['DJF', 'MAM', 'JJA', 'SON']
        # assume calender month order 
        # ( this can be automated using get GC datetime )
        indices = [[11, 0, 1], [2, 3, 4], [5, 6, 7], [8, 9, 10]]

    ars = []
    # extract data by month
    for n, s in enumerate( seasons ):
        print s, n , indices[n] 
        ars += [ [ arr[...,i] for i in indices[n] ] ]
    
    # average by season
#    print [ np.array(i).shape for i in ars ], np.array(ars).mean(), seasons
    print [ np.array(i).shape for i in ars ],  seasons
    ars =  [ np.array(i).mean(axis=0) for i in ars]
    print [ i.shape for i in ars ], np.array(ars).mean(), seasons

    # return list array averaged by season 
    return ars, seasons

# --------------
# 2.39 - Convert v/v to ng/m^3
# -------------     
def convert_v_v2ngm3(  arr, wd=None, spec='AERI', trop_limit=True,\
            s_area=None, res='4x5', aerosol_mass=False, debug=False ):
    """ """

    # Get volume (m^3, adjusted (x1E6) from cm^3) 
    vol = get_volume_np( wd=wd, trop_limit=trop_limit, s_area=s_area, \
                    res=res ) /1E6  

    # Get mass
    a_m = get_GC_output( wd, vars=['BXHGHT_S__AD'], trop_limit=trop_limit, 
                    dtype=np.float64)

    # Get moles  ( converting airmass from kg 1st)
    mols = a_m*1E3/constants( 'RMM_air')

    # adjust to mols, then mass ( assume AERI = 
#    species_mass = {'SO4':96.0, 'AERI': 338.0/2, 'O3':48 }
    # Setup dictionary of species masses... 
    # ( this is separate from species_mass dictionary as assumptions differ)(
    if aerosol_mass:
        I2Ox_mass = (286.0+302.0+318.0/3)/2
        species_mass = {
        'SO4':96.0, 'SO4s':96.0, 'AERI': I2Ox_mass,'ISALA':I2Ox_mass, 'ISALC' : I2Ox_mass,  'O3':48 
        }
    else:
        species_mass = {'SO4':96.0, 'SO4s':96.0, 'AERI': 127.0, 'O3':48.0 }
    arr = arr*mols*species_mass[ spec ]

    print species_mass[ spec ], np.sum( arr )

    # convert to molecules of spec  / m^3   ( then to ng )
#    arr = arr/ vol #*1E9
    # convert to (nano)g/m3
    arr = arr*1E9/vol 
    
    return arr

# ------------------ Section 6 -------------------------------------------
# -------------- Time Processing
#

# --------------
# 3.01 - Takes monthly outputs and sort to chronological order
# -------------
def ctms2chronological( ctms, debug=False ):
    """ ensure list of ctms is chronological  """
    
    # get datetimes for month
    dts = [ get_gc_datetime(ctm ) for ctm in ctms ]
    if debug:
        print dts

    # if multiple months within ctm file
    if debug:
        print ctms
        print dts
        print len(dts), 
    if len(dts) > 1:
        dts = [ i[0] for i in dts ]
    else:
        dts = [ i[0] for i in dts[0] ]

#    print dts
    # Get indicies for sorted list and return chronological list
    ind = [dts.index(i) for i in sorted(dts)]
    if not ( len(dts) > 1 ):
        debug=True        
    if debug:
        print 'before: ', dts, len(dts), ind
        print ctms, ind

    if not ( len(dts) > 1 ):
        if debug:
            print [ctms[i] for i in ind ] 
        ctms =  [ctms[i] for i in ind ]
        if debug:
            print 'after: ', dts, len(dts)

    # deal with single files containing multiple months
    else:
        if debug:
            print ctms, dts

    return ctms

# --------------
# 3.02 - Takes monthly outputs and sort to chronological order
# -------------
def np2chronological_fromctm( ctms, arr, debug=False ):
    """ Ensure np array is in chronological order """

    # get datetimes for month and sort chron.
    dts = [ get_gc_datetime(ctm ) for ctm in ctms ]

    if debug:
        print 'before', dts, len(dts )

    # if multiple months within ctm file
    if len(dts) > 1:
        dts = [ i[0] for i in dts ]
    else:
        dts = [ i[0] for i in dts[0] ]

    sdts = sorted(dts)
    if sdts != sorted(dts):
        print 'sorted in chronological order of np in np2chronological_fromctm'
        if debug:
            print 'after',dts, sdts, len(sdts)

    # Get indicies for sorted list and return chronological list
    ind = [dts.index(i) for i in sdts ]
    if debug:
        print 2, ind

    return np.concatenate( [arr[...,i][...,None] for i in ind ], axis =3)
