#!/usr/bin/python
"""
Functions for use with the GEOS-Chem chemical transport model (CTM).

Use help(<name of function>) to get details on a particular function. 

NOTE(S):    
 - This module is underdevelopment vestigial/inefficient code is being removed/updated. 
 - Where external code is used credit is given. 
"""

# ----------------------------- Section 0 -----------------------------------
# -------------- Required modules:
#
# -- I/O / Low level                                                                                
import os
import sys
import csv
import glob
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
import iris 
import logging

# -- Math/Analysis                                                                                   
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
from iris.time import PartialDateTime 

# --  This needs to be updated, imports should be specific and in individual functions
# import tms modules with shared functions
from AC_tools.funcs4core import *
from AC_tools.funcs4generic import *
from AC_tools.funcs4time import *
from AC_tools.funcs4pf import *
from AC_tools.funcs_vars import *

from AC_tools.Scripts.bpch2netCDF import convert_to_netCDF



# --------------------------------- Section 2 ----------------------------------
# -------------- Model Data Extractors
#


    
# ----
# 1.04 -Get surface area  ( m^2 )
# ----
def get_surface_area(res=None,time=None, debug=False, wd=None):
    """ 
    Get_surface_area of grid boxes for a given resolution

    INPUTS:
    wd=None (Specify the wd to get the results from a run.)
    res='4x5' (Specify the resolution if wd not given.)
    time=None (Not used atall = Probably legacy? bjn)
    debug=False (legacy debug, replaced by logging)
    OUTPUTS:
    s_area (2d numpy array of surface area per gridbox)    

    NOTE(s):
	 - this function accesses previsouly run GEOS-Chem 
        1day files with just DXYP / DXYP diagnostic ouptuted
     - back compatibility with PyGChem 0.2.0 is retained  
    """
    logging.info( "Getting the surface area" ) 

    if res==None and wd==None:
        res = ('4x5')
        logging.warning("No res or wd specified. Assuming 4x5.")

    print locals()

    # if the wd has not been specified then use the previous runs
    if wd==None:

        # What is dwd? 
        # All of this might make sense to replace with example data?
        dwd = os.path.join( get_dir('dwd'),  'misc_ref/' )
        logging.debug("dwd = " + str( dwd) )

    #    dwd = get_dir( 'dwd') + '/misc_ref/'
        dir = {
        '4x5':'/LANDMAP_LWI_ctm',  \
        '2x2.5': '/LANDMAP_ctm_2x25',  \
        '0.5x0.666' :'LANDMAP_LWI_ctm_05x0666',  \
        '0.25x0.3125' :'LANDMAP_LWI_ctm_025x03125',  \
        }
        fd = os.path.join( dwd , dir[res])

        logging.debug( "resolution = {res}, lookup directory = {fd}".format( \
            res=res, fd=fd) )
    #    if debug:
    #        print fd, res
        wd = fd



    logging.debug("Trying to get surface area from {wd}".format(wd=wd))
    try:
        s_area = get_GC_output( wd, vars=['DXYP__DXYP'] ) 
    except ValueError:
        logging.debug("Failed getting surface area from wd")
        logging.debug("Trying to get surface area grom bpch file pygchem v0.2")
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
    except:
        logging.error("Could not get the surface area!")
        raise ValueError("Could not find the surface area")

    return s_area

def list_variables(wd=None):
    """
    Show a list of variables in a wd that we can get at.

    INPUTS:
    wd=None (Specify the working directory)

    OUTPUTS:
    prints a list of variables.

    #Note
    Only currently prints variables in the ctm.nc
    Should be expanded to include planeflight and HEMCO
    """

    logging.info("Listing variables in {wd}".format(wd=wd))
    if wd==None:
        raise ValueError("Please specify a working dir")

    # Try listing all bpch variables in ctm.nc
    try:
        logging.info("Listing all variables in netCDF file")
        ctm_nc = os.path.join(wd, 'ctm.nc')
        if not os.path.isfile( ctm_nc ):
            convert_to_netCDF( wd )

        
        for var in Dataset( ctm_nc ).variables:
            print var
    except:
        logging.info("No ctm.nc (bpch) data found")

    # Try listing all hemco variables
    try:
        logging.info("Listing all variables in HEMCO files")
        hemco_nc = os.path.join(wd, 'hemco.nc' )
        if not os.path.isfile( hemco_nc ):
            convert_to_netCDF( wd )
        for var in Dataset( hemco_nc ).variables:
            print var
    except:
        logging.info("No HEMCO data found")

    # Try listing all variables form planeflight
#    try:
#        logging.info("Listing all variables in planeflight")

    return
        
    

# ----
# 1.05 - Get Land map. 
# ----
def get_land_map(res='4x5', time=None, wd=None,debug=False):
    """ 
    Return land, water, and ice indices (LWI ) from GEOS-Chem with integers for Land (1) 
    and Water (0). Ice fraction is given as fractional values. 

	NOTES:
	 - This approach is inefficent and requires large files. Could this be improved with
	 on-line extract on inclusion of generic output for various resoltions as txt files? 
    """

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
# 1.07 - get air mass (4D) np.array in kg
# ------------- 
def get_air_mass_np( ctm_f=None, wd=None, times=None, trop_limit=True,\
            debug=False ):
    """ 
    Get array of air mass in kg 
    """
    logging.info( 'called get air mass')

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

        if trop_limit:
            arr = arr[...,:38,...]

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
    """ 
    Get mean OH concentration (1e5 molec/cm3) from geos.log files
    in given directory 
    """

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
# 1.0X - Get CH4 mean value - v/v
# -------------
def get_CH4_mean( wd, rtn_global_mean=True, rtn_avg_3D_concs=False, \
        res='4x5', debug=False ):
    """ 
    Get mean CH4 concentrtaion from geos.log files in given directory 
    """

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
    if rtn_avg_3D_concs:
        rtn_value = get_CH4_3D_concenctrations( z )
    # return value for UK latitude
    if (not rtn_global_mean) and (not rtn_avg_3D_concs):
        rtn_value = np.array( z )[::4].mean()
    return rtn_value


# ----
# 1.14 - Extract OH and HO2 for a given ctm.bpch file
# ---
def get_OH_HO2( ctm=None, t_p=None, a_m=None, vol=None, \
            wd=None, HOx_weight=False, res='4x5', scale=1E5, trop_limit=True, \
             molec_weight=True, time_averaged=True, debug=False ):
    """ 
    Get OH/HO2 concentrations from ctm.bpch/NetCDF of GEOS-Chem output 
    """

    if debug:
        print 'get_OH_HO2 called for ', wd

    # --- Set specs, get constants, and extract output variables
    specs = ['OH','HO2']
    AVG= constants('AVG')
    RMM_air=constants('RMM_air')
    # Load arrays if not provided...
    if not isinstance(t_p, np.ndarray):
        if pygchem.__version__ == '0.2.0':
            t_p = get_gc_data_np( ctm, spec='TIMETROP', category='TIME-TPS', \
                debug=debug)
        else:
            print 'WARNING!!! - provide time in trop diag. ' 
            t_p = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], \
                trop_limit=trop_limit ) 

    if not isinstance(a_m, np.ndarray):
        if pygchem.__version__ == '0.2.0':
            a_m = get_air_mass_np( ctm, debug=debug )
        else:
            print 'WARNING!!! - provide air mass diag. ' 
            a_m = get_GC_output( wd, vars=['BXHGHT_S__AD'],  \
                trop_limit=trop_limit, dtype=np.float64)
    if not isinstance(vol, np.ndarray):
        if pygchem.__version__ == '0.2.0':
            vol = get_volume_np( ctm, res=res, debug=debug)
        else:
            print 'WARNING!!! - provide vol diag. ' 
            vol = get_volume_np( wd=wd, res=res, trop_limit=trop_limit, \
                debug=debug)

    # --- Extract OH and HO2 Data ( OH in [molec/cm3], HO2 in [v/v])
    if pygchem.__version__ == '0.2.0':
        OH, HO2 = [get_gc_data_np( ctm, \
            i, category='CHEM-L=$') for i in specs ] 
    else:
        OH, HO2 = get_GC_output( wd, trop_limit=trop_limit, \
            vars=['CHEM_L_S__'+i for i in specs], r_list=True )

    # Mask for troposphere. 
    OH, HO2  = mask4troposphere( [OH, HO2 ], t_ps=t_p,  \
        use_time_in_trop=True, multiply_method=True)    
    
    # --- Process data
    molecs = ( ( (a_m*1E3) / RMM_air ) * AVG )   # molecules
    moles =  a_m*1E3 / RMM_air # mols 

    # Only consider troposphere
    print [ i.shape for  i in molecs, a_m, vol, moles  ]
    molecs, a_m, vol, moles = mask4troposphere( \
        [molecs, a_m, vol, moles ], t_ps=t_p,  \
        use_time_in_trop=True,  multiply_method=True)    

    if debug:
        print [ (i.shape, np.mean(i) ) for i in [ OH, HO2, molecs, moles, a_m, \
            vol ]]

    # convert HO2 [v/v] to [molec/cm3]
    HO2 = convert_v_v_2_molec_cm3( HO2, a_m=a_m, vol=vol, mols=moles )

    # Remove invalid values
    HO2, OH  = [ np.ma.masked_invalid( i) for i in HO2, OH   ]
    HO2, OH  = [ np.ma.masked_where( i<0, i) for i in HO2, OH   ]
    
    # Average over provided timesteps ( 4th dimension )
    if time_averaged:
        HO2, OH, molecs, vol = [ i.mean(axis=-1) for i in HO2, OH, molecs, vol ]

    if debug:
        print [ ( np.ma.sum(i), np.ma.mean(i) ) for i in HO2, OH ]
        print 'volume weighted: ', [ np.ma.sum( i *vol) / np.ma.sum(vol) \
            for i in HO2, OH ] 

    if HOx_weight: # weigh by [HO]+[HO2] molecules
        HOx = HO2 + OH
        HO2, OH  = [ np.ma.sum( ( i* HOx )) /np.ma.sum(HOx) for i in HO2, OH   ]
        if debug:
            print [ i.shape for i in OH, HO2, moles, vol, HOx ]
            print  'HOx weighted: ',HO2, OH

    elif molec_weight: # weight by # of molecules
        HO2, OH  = [ ( i*molecs).sum() /molecs.sum() for i in HO2, OH   ]
        
    else: # Volume weight
#       HO2, OH = [ np.ma.sum( i *vol) / np.ma.sum(vol)  for i in HO2, OH ] 
        print 'Please specify weighting in get_OH_HO2() '
        sys.exit()

    # Scale to set value ( e.g. 1E6 )
    HO2, OH = [i/scale for i in HO2, OH ]

    return OH, HO2

# ----
# 1.21 - Process species for given arrays to (v/v) in respective scale + DU
# ---
def process_data4specs( specs=None, just_bcase_std=True, preindustrial=False, \
            just_bcase_no_hal=False, res='4x5', ver='1.6', diff=True, \
            pcent=True, tight_constraints=True, trop_limit=True, \
            NOy_family=False, Bry_family=False, Iy_family=False, \
            rtn_specs=False, debug=False ): 
    """ 
    Return species values in v/v and DU. Also return datetimes for CTM output and time in 
    troposphere diagnostics  

	NOTE(s):
	 - This function is fairly inefficient, but is compbatible with pygchem 0.2.0 and 0.3.0
    """

    # Get runs data and descriptors
    wds, titles = MUTD_runs( titles=True, just_bcase_std=just_bcase_std, \
        res=res, ver=ver, just_bcase_no_hal=just_bcase_no_hal, \
        preindustrial=preindustrial )
    if debug:
        print wds, titles

    if preindustrial:
        # Force 'Cl+Br+I (PI)' to be base case and 'Cl+Br+I' to be second entry
        wds = [wds[3],  wds[1] ]
        titles = [ titles[3], titles[1] ]
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
            verbose=False, debug=False):
    """ 
    Data extractor for GEOS-Chem using PyGChem (>= 0.3.0 ). ctm.bpch files are extracted 
    for a given directory, with only the specific category and species returned. This can 
    either be as a Iris Cube (retaining metadata) or as a numpy array.
        
    ARGUMENTS:
     - vars(list): variables to extract (in NetCDF form, e.g. IJ_AVG_S__CO)
     ( alterately provide a single species and category through named input variables )
     - r_cubes (Boolean): To return Iris Cubes, set to True
     - r_res (Boolean): To return resolution of model NetCDF set to True
     - restore_zero_scaling(Boolean): restores scale to ctm.bpch standard (e.g. v/v not pptv)
     - trop_limit(boolean): limit to "chemical troposphere" (level 38 of model)
     - dtype: type of variable to be returned 
     - use_NetCDF(boolean)
     
     NOTES:
      - Credit for PyGChem: Ben Bovy - https://github.com/benbovy/PyGChem
      - Examples and use of pygchem is discussed on Ben Bovy's GITHib 
     ( https://github.com/benbovy/PyGChem_examples/blob/master/Tutorials/Datasets.ipynb )
      - This replaces the now defunct AC_tools functions: open_ctm_bpch and get_gc_data_np
      - For simplicity use variable names (species, category) as used in 
      the iris cube ( e.g. IJ_AVG_S__CO). if variable is unknown, just 
      print full dataset extracted to screen to see active diagnostics. 
      - Species and category variables are maintained ( and translated ) to allow for 
      backwards compatibility with functions written for pygchem version 0.2.0
    """

# bjn
# This function is not completly clear to me, and could do with a re-write
# The try command would probably be useful here for large parts.
# logging would also be good for replacing debug.
    logging.info("Called get_GC_output")
    logging.debug("get_GC_output inputs:")
    logging.debug(locals())
    
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

#    # ensure wd has a leading '/'
#    if wd[-1] != '/':
#        wd +=  '/'

    # Work with NetCDF. Convert ctm.bpch to NetCDF if not already done.
    if use_NetCDF:

        # Check for compiled NetCDF file
        # If not found, create NetCDF file from ctm.bpch files                
        import os.path
        fname = os.path.join(wd, 'ctm.nc')
        if not os.path.isfile(fname):
            from bpch2netCDF  import convert_to_netCDF
            convert_to_netCDF( wd )

        logging.debug("Opening netCDF file {fname}".format(fname=fname))
        # "open" NetCDF + extract requested variables as numpy arr.

        netCDF_data = Dataset( fname, 'r' )
        arr = []
        for var in vars:
            try:
                logging.debug("opening variabel {var}".format(var=var))
                var_data =  netCDF_data.variables[var] 
            except:
                logging.warning("Variable {var} not found in netCDF")\
                    .format(var=var)
                logging.warning("Will attempt renaming")
                abrv_var = get_ctm_nc_var( var )
                try:
                    var_data =  netCDF_data.varialbes[var] 
                except:
                    logging.error("Renamed variable {var} not found in netCDF")\
                        .format(var=abrv_var)




####################################################################################                
#####--- bjn - re-wrote (above) to make more understandable ---###
#
#            with Dataset( fname, 'r' ) as rootgrp:
#                try:
#                    arr = [ np.array(rootgrp[i]) for i in vars ]  
#                except IndexError:
#                    if verbose:
#                        print 'WARNING: VAR NOT IN NETCDF '
#                        print 'IndexError was found either due to lack of ' + \
#                            'variable in NetCDF file or error in name' + \
#                            ' species import from *.dat. Will attempt renaming'
#                    arr =[]
#                    for var_ in vars:
#                        if verbose:
#                            print 'Atempting indiviudual extraction of: ', var_
#                        try:
#                            arr += [ np.array( rootgrp[var_] ) ]
#                            if verbose:
#                                print 'successfull indiv. extraction of: ',var_
#                        except IndexError:
#                            logging.error('failed to find {var}'.format(var=var_))
##                            raise ValueError('failed to find {var} in netCDF file'\
##                                    .format(var=var_))
#
#                
#                            # If not found try an abreviation
#                            abrv_var_ = get_ctm_nc_var(var_)
#                            arr += [ np.array( rootgrp[ abrv_var_ ] )]
#                            if verbose:
#                                print 'using {} instead of {}'.format( \
#                                     abrv_var_, var_ )
#
##################################################################################


#                # files are stored in NetCDF at GC scaling. 
#                # ( This is different to ctm.bpch, rm for back compatibility. )

############################################################################
#    # This is not in a working state currently - needs work
            if restore_zero_scaling:
                try:
                    var_data = np.divide(var_data,get_unit_scaling(var_data.ctm_units))
                except:
                    logging.warning("Scaling not adjusted to previous approach")

            arr.append(var_data[:])

####--- The above re-write does not work so still using old version ---###

# Temp fix for out of place code if true:
#                arr.append(var_data[:]) # temp fix
#        if True:# temp fix
#                rootgrp = netCDF_data #temp fix
#
#                if restore_zero_scaling:
#                    try:
#                        arr =[ arr[n]/get_unit_scaling( rootgrp[i].ctm_units ) \
#                            for n, i in enumerate( vars ) ]
#                    except:
#                        print 'WARNING: SCALING NOT ADJUSTED TO' + \
#                            ' PREVIOUS APPROACH'
#############################################################################


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
        logging.error( 'WARNING this approach has been removed due to time cost')
        sys.exit( 0 )

    # Process extracted data to gamap GC format and return as numpy 
    if not r_cubes:

        # Limit to GEOS-Chem "chemical troposphere'
        if trop_limit:
            arr = [ i[...,:38] for i in arr]

        # Convert to GC standard 4D fmt. - lon, lat, alt, time 
        if len((arr[0].shape)) == 4:
            if debug:
                print 'prior to roll axis: ', [i.shape for i in arr]
            arr = [np.rollaxis(i,0, 4) for i in arr]
            if debug:
                print 'post roll axis: ', [i.shape for i in arr]

        # Convert to GC standard 3D fmt. - lon, lat, time   
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

        # Convert type if dtype not float32 
        # ( needed for some arrays e.g. air mass )
        if dtype != np.float32:
            arr = arr.astype( dtype )

    # Get res by comparing 1st 2 dims. against dict of GC dims.
    if r_res:
        res=get_dims4res( r_dims=True, trop_limit=trop_limit, \
                    just2D=True )[arr.shape[:2]]

    if r_list:
        # Make sure returned type is list of arrays
        if len( vars ) > 1:
            arr = [ arr[i,...] for i in range( len(vars) ) ]
        else:
            arr = [ arr ]

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
# 1.24 - Get gc resolution from ctm.nc
# ---
def get_gc_res( wd ):
    """
    Extract resolution of GEOS-Chem NetCDF
    """
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

# ----
# 1.25 - Get surface area ( m^2 ) for any global resolution
# ----
def calc_surface_area_in_grid( res='1x1', debug=False ):
    """ 
    Use GEOS-Chem appraoch for GEOS-Chem grid surface areas.

	NOTE(s):
     - update this to take any values for res... 
	 - Does this function need updating?

        Credit: Bob Yantosca
        Original docs from ( grid_mod ):
            !======================================================================
    ! Compute grid box surface areas (algorithm from old "input.f")
    !
    ! The surface area of a grid box is derived as follows:
    ! 
    !    Area = dx * dy
    !
    ! Where:
    !
    !    dx is the arc length of the box in longitude
    !    dy is the arc length of the box in latitude
    !  
    ! Which are computed as:
    !  
    !    dx = r * delta-longitude
    !       = ( Re * cos[ YMID[J] ] ) * ( 2 * PI / IIIPAR )
    !
    !    dy = r * delta-latitude
    !       = Re * ( YEDGE[J+1] - YEDGE[J] )
    !  
    ! Where:
    !    
    !    Re         is the radius of the earth
    !    YMID[J]    is the latitude at the center of box J
    !    YEDGE[J+1] is the latitude at the N. Edge of box J
    !    YEDGE[J]   is the latitude at the S. Edge of box J
    !
    ! So, the surface area is thus:
    ! 
    !    Area = ( Re * cos( YMID[J] ) * ( 2 * PI / IIIPAR ) *
    !             Re * ( YEDGE[J+1] - YEDGE[J] )
    !
    !    2*PI*Re^2    {                                            }      
    ! = ----------- * { cos( YMID[J] ) * ( YEDGE[J+1] - YEDGE[J] ) }
    !     IIIPAR      {                                            }
    !
    ! And, by using the trigonometric identity:
    !
    !    d sin(x) = cos x * dx
    !
    ! The following term:
    !
    !    cos( YMID[J] ) * ( YEDGE[J+1] - YEDGE[J] ) 
    !
    ! May also be written as a difference of sines:
    !
    !    sin( YEDGE[J+1] ) - sin( YEDGE[J] ) 
    ! 
    ! So the final formula for surface area of a grid box is:
    ! 
    !            2*PI*Re^2    {                                     }
    !    Area = ----------- * { sin( YEDGE[J+1] ) - sin( YEDGE[J] ) }
    !              IIIPAR     {                                     }
    !
    !
    ! NOTES:
    ! (1) The formula with sines is more numerically stable, and will 
    !      yield identical global total surface areas for all grids.
    ! (2) The units are determined by the radius of the earth Re.
    !      if you use Re [m], then surface area will be in [m2], or
    !      if you use Re [cm], then surface area will be in [cm2], etc.
    ! (3) The grid box surface areas only depend on latitude, as they
    !      are symmetric in longitude.  To compute the global surface
    !      area, multiply the surface area arrays below by the number
    !      of longitudes (e.g. IIIPAR).
    ! (4) At present, assumes that GEOS-Chem will work on a
    !      Cartesian grid.
    !
    ! (bmy, 4/20/06, 2/24/12)
    !======================================================================
        

    """

    logging.info('called calc surface area in grid')

    # Get latitudes and longitudes in grid    
    lon_e, lat_e, NIU = get_latlonalt4res( res=res, centre=False, debug=debug )    
    lon_c, lat_c, NIU = get_latlonalt4res( res=res, centre=True, debug=debug )    

    # Set variables values
    PI_180 = pi/ 180.0
    Re = np.float64( 6.375E6 ) # Radius of Earth [m] 
    lon_dim = get_dims4res(res=res)[0]
    lon_degrees = float( lon_dim )

    # Loop lats and calculate area    
    A1x1 = []
    for n, lat_ in enumerate( lat_e[:-1] ):

        # Lat at S and N edges of 1x1 box [radians]
        S = PI_180 * lat_e[n]
        N = PI_180 * lat_e[n+1]        

        # S to N extent of grid box [unitless]
        RLAT    = np.sin( N ) - np.sin( S )

        # 1x1 surface area [m2] (see [GEOS-Chem] "grid_mod.f" for algorithm)
        A1x1 += [ 2.0 * pi * Re * Re / lon_degrees  * RLAT ]

    A1x1 = np.array( A1x1 )
    if debug:
        print A1x1
        print [ (i.shape, i.min(),i.max() ) for i in [ A1x1 ] ]

    # convert to 2D array / apply to all longitudes 
    A1x1 = np.array( [ list( A1x1 ) ] * int(lon_dim) )

    return A1x1

# ----
# 1.26 - Process species for given family arrays to (v/v)
# ---
def get_chem_fam_v_v_X( wd=None, fam='Iy', res='4x5', ver='3.0' , specs=None, \
    trop_limit=True, N=False, I=False, Cl=False, Br=False, t_ps=None, a_m=None,\
    vol=None, verbose=True, rm_strat=False, debug=False ):
    """  
    Return array of family in mols of X ( e.g. Cl, Br, I, N ) equiv. in mol/mol.  

	NOTES:
	 - Is this function just a double up of fam_data_extractor? 
	 (which is more up to date)
    """

    # Get species time Tropopause diagnostic
    if isinstance( t_ps, type(None) ) and rm_strat:
        t_ps = get_GC_output( wd=wd, vars=['TIME_TPS__TIMETROP'], \
            trop_limit=trop_limit )

    # Get (fam) specs if not given
    if isinstance( specs, type(None) ):
        # get species in family
        d = { 
        'NOy':'N_specs', 'Iy':'Iy',  'Bry':'Bry', 'Cly':'Cly', 'HOx':'HOx',  \
        'SOx':'SOx', 'NOx':'NOx'
                }
        specs = GC_var( d[fam] )

    # Use correct stiochmetry 
    # ( This is no longer required as ref_spec is passed, however it is retained
    # for back compatibility )    if fam =='Bry':
        Br=True
    if fam =='Iy':
        I=True
    if fam =='Cly': 
        Cl=True        
    if any( [( fam ==i) for i in  ' NOy', 'NOx' ] ):
        N=True

    # Get mixing ratio 
    if fam == 'HOx':
        # OH ( in molec/cm3 )
        arr = get_GC_output( wd=wd, vars=['CHEM_L_S__'+'OH'], \
                trop_limit=trop_limit )    
        # Convert to  v/v 
        arr = convert_v_v_2_molec_cm3( [arr], a_m=a_m, vol=vol, wd=wd )[0] 
        # HO2 ( in v/v )
        arr = [ arr, get_GC_output( wd=wd, vars=['CHEM_L_S__'+'HO2'], \
                trop_limit=trop_limit ) ]

    else: # Just extract v/v 
        arr = get_GC_output( wd=wd, vars=['IJ_AVG_S__'+i for i in specs ],
            trop_limit=trop_limit, r_list=True ) 
    print [ i.shape for i in arr], len( arr), np.sum( arr ), specs

    # Adjust to stiochmetry  ( Vars )
    arr = [ arr[n]*spec_stoich(i, ref_spec=fam) \
        for n,i in enumerate( specs) ]
    print [ i.shape for i in arr], len( arr), np.sum( arr ), specs

    # Sum over stiochmertically adjusted list of specs
    arr = np.array( arr ).sum(axis=0)

    # Remove stratosphere by multiplication of time in trop. diag. 
    if rm_strat:
        arr = mask4troposphere( [arr], t_ps=t_ps )[0]
        
    return arr, specs

# ----
# 1.27 - Convert v/v array to DU array 
# ---
def convert_v_v_2_DU( arr, wd=None, \
     a_m=None, trop_limit=True, s_area=None, molecs=None, \
    verbose=True, debug=False):
    """ 
    Convert a 4D array of v/v for species (or family) to  DU
    """

    # Get DU values for each array (2D not 3D )
    if isinstance( molecs, type(None) ):
        # If 'a_m' not given, get air mass ('a_m') in kg
        if isinstance( a_m, type(None) ):
            a_m = get_GC_output( wd=wd, vars=['BXHGHT_S__AD'], 
                trop_limit=trop_limit, dtype=np.float64) 
        # Get molecs in troposphere
        molecs = a_m*1E3/constants( 'RMM_air')*constants('AVG')

    # Get surface area
    if isinstance( s_area, type(None) ):
        s_area = get_surface_area( res=res, debug=debug )

    # Molecules O3 in trop. summed
    DUarrs = arr*molecs
    if debug:
        print [ ( i.shape, i.sum() ) for i in [DUarrs, s_area] ]
    # sum over altitude in 4D array ( lon, lat, alt, time )
    DUarrs = DUarrs.sum(axis=-2)
    if debug:
        print [ ( i.shape, i.sum() ) for i in [DUarrs, s_area] ]

    # adjust to per unit area ( cm3 ) 
    DUarrs =  DUarrs/ s_area

    if debug:
        print [ ( i.shape, i.sum() ) for i in DUarrs, tmp_s_area,s_area  ]

    # convert to DU   
    DUarrs = DUarrs/constants('mol2DU')

    return DUarrs 
    
# ----
# 1.28 - Get common GC diagnostic arrays 
# ---   
def get_common_GC_vars( wd=None, trop_limit=True, res='4x5', \
        verbose=True, debug=False):
    """  
    Returns t_ps, a_m, molecs, s_area
    NOTES:
     - This appraoch is taken to avoid circular calls + reduces code repitions 
    """

    # Get species time Tropopause diagnostic
    t_ps = get_GC_output( wd=wd, vars=['TIME_TPS__TIMETROP'], \
            trop_limit=trop_limit )
    # Get air mass in kg
    a_m = get_GC_output( wd=wd, vars=['BXHGHT_S__AD'], 
        trop_limit=trop_limit, dtype=np.float64) 

    # Get molecs in troposphere
    molecs = a_m*1E3/constants( 'RMM_air')*constants('AVG')

    # Get surface area
    s_area = get_surface_area( res=res, debug=debug )

    return t_ps, a_m, molecs, s_area

# ----
# 1.29 - Get 3D CH4 concentrations
# ---   
def get_CH4_3D_concenctrations( z, res='4x5', trop_limit=True, debug=False ):
    """ 
    Takes monthly ( or any equllay spaced output) CH4 concentrations from geos.log 
    files. 4 value are given per monthly file, split by laititude, as list "z"
    
    Which apears in GC output ( geos.log) as: '
    CH4 (90N - 30N) :  X.X [ppbv]
    CH4 (30N - 00 ) :  X.X [ppbv]
    CH4 (00  - 30S) :  X.X [ppbv]
    CH4 (30S - 90S) :  X.X [ppbv]   '

    NOTE:   
     - this programme is typically called with z in units of v/v, not ppbv 
    """
    
    # latitudes between which CH4 conc is given
    lats = [90, 30, 0, -30, -90][::-1] 
    
    # make a array of zeros, including time dimension
    arr_shape  = get_dims4res( res=res, trop_limit=trop_limit ) + (len( z )/4, )
    arr  = np.zeros( arr_shape )

    # Get an average of each latitude band ( 4 outputs per geos.log file )
    z = np.array( [ z[i::4] for i in range( 4 ) ] )
    # Also reverse list to allow for increasing indice (for latitude )
    z = z[::-1]

    # Loop time stamps ( 4th dimension )
    for t in range( len( z[0,:] ) ):

        # Apply these values to the GC gird at given resolution
        for n, lat in enumerate( lats[:-1] ):
            lat_0 =  get_gc_lat( lats[n], res=res)
            lat_1 =  get_gc_lat( lats[n+1], res=res)
            arr[:,lat_0:lat_1,:,t] = z[n,t]

        # Check all values have been filled
        if debug:
            print z.shape, arr.shape, arr[ arr<0 ]        

    return arr

# ----
# 1.30 - Get Strat-Trop exchange (from geos.log files )
# ---  
def get_STRAT_TROP_exchange_from_geos_log( fn=None, ver='3.0', \
        rtn_date=False, rtn_Tg_per_yr=True, verbose=False, debug=False ):
    """ 
    Extract all tracer Stat-trop exchange values for a given geos.log file. 
    These are then returned as a Pandas Dataframe
    
    NOTEs:
     - Works by extracting all lines between start ("Strat-Trop Exchange")
         and end ("================") of section. 
     - file name (fn) is assumed to include directory as well as name
    """
    # --- Open the file
    file_ =  open( fn, 'rb' )    

    # Read in just the TROP-STRAT exchange section 
    start_line =  'Strat-Trop Exchange'
    end_line = '================'

    # --- Loop and extract lines of file with data on exchange
    readline = False
    lines = []
    for row in file_:

        # once at prod/loss section, start added to list
        if start_line in row:
            readline=True

        # if not at end of prod/loss section, add to list
        if end_line in row:
            readline=False
        if readline:
            try:
                lines.append( row )
            except:
                lines = [ row ]  

    # --- Process extracted lines
    # remove starting lines
    headers = [ i.strip() for i in lines[5].split('    ') ]
    # remove end line
    lines.pop()
    # What data range is this?
    date_range = lines[3]
    # split into columns
    lines = [ i.split() for i in lines[6:] ]
    # remove colon
    TRAs = [ i[0].split(':')[0] for i in lines  ]
    # loop columns and extract data
    vars = []
    if debug:
        print headers, date_range, lines
    # Extract data as floats by header 
    vars = [ [ float( i[n+1]) for i in lines ] for n in range( len(headers)-1) ]
    # make a DataFrame of the output
    if debug:
        print [ len(i) for i in vars ]
    d = dict( zip(headers[1:], vars) )
    df = pd.DataFrame(  d, index=TRAs )
    if rtn_Tg_per_yr:
        df = df['= [Tg a-1]' ]

    if rtn_date:
        return df, date_range
    else:
        return df

# ----
# 1.30 - Get Wind direction from U and V vectors 
# ---  
def get_mod_WIND_dir(  sdate=datetime.datetime(2012, 8, 1, 0 ), \
            edate = datetime.datetime(2012, 8, 8, 0 ), loc='KEN', \
            scale=1, adjustby=0, period = 'summer', \
			vars = ('GMAO_UWND', 'GMAO_VWND'), \
            verbose=False, debug=False):
    """ 
    Extract synoptic wind direction 

    NOTES:
	 - this function was written to work with GEOS-Chem planeflight output
	(thus U/V vector variables are set to pf varaibles), but alternates are accepted
	as arguements
    """
	
    # Extract U10, W10 ( use 10 m wind? )
    datal = []
    
    for spec in vars:
        # Extract model data
#        data, dates, units = get_pf_mod_data( sdate=sdate, edate=edate, \
#                loc=loc, spec=spec, scale=scale,adjustby=adjustby, \
#                period=period, units=None, debug=debug)
        data, dates, units = get_NetCDF_mod_data( sdate=sdate, \
                edate=edate, loc=loc, spec=spec, scale=scale,\
                adjustby=adjustby, period=period, EOHemisions=True, \
                units=None, debug=debug)

        datal += [data ]

    # Make dataframe to allow for function mapping 
    df = DataFrame( data=np.array(datal).T, columns=vars )

    # Calculate wind dir  " (270-atan2(V,U)*180/pi)%360  "
    def f(x):    
        return (270-atan2(x['GMAO_VWND'],x['GMAO_UWND'])*180/pi)%360
    df= df.apply( f, axis=1)

    # Return 
    return [ np.array(i) for i in df, dates ]

# ----------------------- Section 3 -------------------------------------------
# -------------- Data Processing tools/drivers
#

# --------   
# 2.01 - Retrieve model resolution
# --------
def mod_res(wd, spec='O3', fn='ctm.bpch', debug=False):
    """ 
    Extract model resolution
    NOTES: 
     - this is not compatible with PyGChem 0.3.0 
    """

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
    """ 
    Return list of years in CTM output 
    """

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
def get_gc_months(ctm_f=None, wd=None, verbose=False, debug=False):
    """ 
    Return list of months in CTM output 
    """
    if pygchem.__version__ == '0.2.0':
        diagnostics = ctm_f.filter(name='O3', category="IJ-AVG-$" )
        return [diag.times[0].strftime('%m')  for diag in diagnostics ]

    # Use updated PyGChem ( >0.3.0 ) approach
    else:
        dates = get_gc_datetime( wd=wd, debug=debug, verbose=verbose )
        return [ i.month for i in dates ] 

# ----
# 2.07 - Get gc datetime
# -----
def get_gc_datetime(ctm_f=None, wd=None, spec='O3', cat='IJ-AVG-$', \
            verbose=False, debug=False):
    """ 
    Return list of dates in datetime output from CTM output 
    """

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
            if verbose:
                print dates, dates.units
            # Get units from cube, default is 'hours since 1985-01-01 00:00:00'
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
# 2.08 - Work out iGEOS-Chem version 
# ------------- 
# NOTE: moved to funcs4core ... 

# --------------
# 2.09 - Get tropospheric Burden - 
# -------------
# Redundent - mv'd to bottom of this module



# --------------
# 2.11 - Get Emission of species in Gg
# -------------
def get_emiss( ctm_f=None, spec=None, wd=None, years=None, \
            molec_cm2_s=False, nmonl_m2_d=False, kg_m2_s=False, \
            monthly=False, months=None, s_area=None, res='4x5', \
            ref_spec='I', debug=False ):
    """ 
    Extract given species emissions from BIOGSRCE category diagnostic

    NOTES:
     - Back compatibility maintained with PyGChem 0.2.0 
     - Asumption on iodine emissions 
     ( set ref_spec to mass unit equivelnces wanted ( e.g. Br )  )
    """
            
    if debug:
        print 'get_emiss called for >{}<'.format(spec)
    if not isinstance(years, list):
        years = get_gc_years( ctm_f=ctm_f, set_=False, wd=wd )
    if not isinstance(months, list):
        months = get_gc_months( ctm_f=ctm_f, wd=wd)
                                             
    # Adjust to " Gg (X) / monthly" from "Kg/m2/ s"  
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
        # g ( e.g. I  ) / month 
        arr_ = arr_*1E3/ species_mass(spec)*species_mass(ref_spec) * \
            spec_stoich(spec)

    if nmonl_m2_d:
        # Convert to (g) / m2
        arr_  =  arr_ / s_area 

        # Convert to (g) / m2 / day
        arr_  =  arr_ / (365./12.) 

        # Convert to nmol ( /m2/day ) ( reverse normalisation to X mass equiv. )
        arr_ = arr_ / species_mass(ref_spec)  /  spec_stoich(spec) 
        arr_ = arr_*1E9

        if debug:
            print 'get_emiss - 2', arr_.shape

    if molec_cm2_s:
        # From  "I Gg/month" to "I Gg/month/cm/2" #convert to /m2 => cm/2
        arr_  =  arr_ / (s_area *10000.) 

        # Convert to / day => hour => hour => minute => sec 
        arr_  =  arr_ / (365./12.) / 24. / 60. / 60. 

        # Convert to molecules ( reverse normalisation to X mass equiv. )
        arr_ = ( arr_ / species_mass(ref_spec)  ) / spec_stoich(spec) *\
                        constants('AVG')  

        if debug:
            print 'get_emiss - 3', arr_.shape

    return arr_

# --------
# 2.12 - Get CH4 Lifetime in years
# --------
def get_CH4_lifetime( ctm_f=None, wd=None, res='4x5', \
            vol=None, a_m=None, t_ps=None, K=None, t_lvl=None, n_air=None, \
            years=None, months=None, monthly=False, trop_limit=True, \
            use_OH_from_geos_log=False, average_value=True, LCH4_Cl=None,  \
            use_time_in_trop=True, masktrop=False, include_Cl_ox=False, \
            verbose=True, debug=False ):
    """ 
    Get amtospheric methane (CH4) lifetime by calculating loss rate against OH.
    
    NOTES:
     - Multiple options for calculation of this value. 
     - uses diagnostic arrays for [OH] instead of goes.log values as default.
     - Calculates the total atmospheric lifetime of CH4 due to oxidation by 
     tropospheric OH as default ( aka does not mask  stratosphere. )
     - 1st approach ( set use_time_in_trop=True) uses reaction rate in 
     globchem.dat (K(o)= 2.45E-12, Ea/R= -1775 ) and OH/CH4 mean
     concentrations from geos.log. The tropospheric mask uses is defined by 
     the time in the troposphere diagnostic. the value is weighted by **air 
      mass** if  ( average_value=True  )
     - Other approach: Effectively the same approach, but uses the tropopause 
     level to define tropopause ( set use_time_in_trop=False to use this), 
     and weights output by **molecules** if ( average_value=True  ) 
    """
    if debug:
        print 'get_CH4_lifetime called ( using time in trop diag?={})'.format(
            use_time_in_trop ) 

    # --- Get shared variables that are not provided
    if not isinstance(vol, np.ndarray): # cm^3 
        vol = get_volume_np( ctm_f =ctm_f, wd=wd, trop_limit=trop_limit ) 
    if not isinstance(a_m, np.ndarray): # Kg
        a_m = get_air_mass_np( ctm_f=ctm_f, wd=wd, trop_limit=trop_limit )  
    if not isinstance(K, np.ndarray): # K 
        K = get_GC_output( wd=wd, vars=[u'DAO_3D_S__TMPU'], \
            trop_limit=trop_limit  ) 
    if not isinstance(t_ps, np.ndarray): 
        t_ps = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], \
            trop_limit=trop_limit ) 
    if not use_time_in_trop:
        if not isinstance(t_lvl, np.ndarray): 
            t_lvl = get_GC_output( wd, vars=['TR_PAUSE__TP_LEVEL'], \
                trop_limit=False ) 
        if not isinstance(n_air, np.ndarray): 
            n_air = get_GC_output( wd, vars=['BXHGHT_S__N(AIR)'], \
                trop_limit=trop_limit, dtype=np.float64 ) 
    
    # Get OH conc [molec/cm3] 
    if use_OH_from_geos_log:
        OH   = get_OH_mean( wd ) * 1E5   
    else: # Extract from GEOS-Chem fields
        OH = get_GC_output( wd=wd, vars=['CHEM_L_S__'+'OH'], \
            trop_limit=trop_limit )
#        #  What are the units for 'CHEM_L_S__OH' ? need to convert?
          # (NetCDF says molec/m3, and wiki says  [molec/cm3]  )
#        OH = OH *1E6
    if include_Cl_ox:
        # Get mixing ratio [v/v]
        Cl = get_GC_output( wd, vars=['IJ_AVG_S__'+'Cl'], \
            trop_limit=trop_limit) 
        # Convert units to [molec/cm3] 
        Cl = convert_v_v_2_molec_cm3( Cl, vol=vol, a_m=a_m )

        # Get global CH4 conc  ( v/v )
        CH4 = get_CH4_mean( wd, rtn_avg_3D_concs=True  )
        # Convert to [molec/cm3] from v/v
        CH4 = convert_v_v_2_molec_cm3( CH4, vol=vol, a_m=a_m )

        if use_time_in_trop and masktrop:
            # Mask non tropospheric boxes
            OH = mask4troposphere( [OH], t_ps=t_ps)[0]
            if include_Cl_ox:
                Cl = mask4troposphere( [Cl], t_ps=t_ps)[0]
                CH4 = mask4troposphere( [CH4], t_ps=t_ps)[0]

    if use_time_in_trop and masktrop:
        # Mask non tropospheric boxes
        K, vol, a_m = mask4troposphere( [K, vol, a_m], t_ps=t_ps )

        # get loss rate - in kg/s - only works for CH4 sim.
    #    LCH4 = get_gc_data_np(ctm_f, spec='CH4Loss', \
    #            category='CH4-LOSS') 

    # --- Shared processed variables
    # CH4 loss rate (OH) per grid box  ( cm^3  molec.^-1  s^-1  )
    KCH4 = 2.45E-12 * np.ma.exp( (-1775. / K)  )   
    # Fixing PD from smvgear tag (LR63) as PD354
    if include_Cl_ox:
        Lrate_CH4_Cl = get_GC_output( wd, vars=['PORL_L_S__'+'PD354'], \
                trop_limit=trop_limit) 
        
    # --- Now Calculate CH4 lifetimes
    if use_time_in_trop:

        # --- Calc from (CH4 +OH) reaction (loss) rate 
        # CH4 loss rate via reaction
    #    arr = LCH4 * CH4

        # Get CH4 lifetime with respect to OH 
        # (( cm^3  molec.^-1  s^-1  )* [molec/cm3]  =  [s^-1] )
        LCH4 = KCH4 * OH 

        # Get CH4 lifetime with respect to Cl        
        if include_Cl_ox:
            # ( 1/ ( [molec/cm3]/[molec/cm3/s] ) =  [s^-1]  )
#            LCH4_Cl =1/ np.ma.divide( Cl, Lrate_CH4_Cl  )
            LCH4_Cl =1/ np.ma.divide( CH4, Lrate_CH4_Cl  )

        if debug:
            print ( 1 / ( (LCH4*a_m).sum()/a_m.sum() )   ) /60/60/24/365
            print [ (i*a_m).sum()/a_m.sum() for i in OH, LCH4, KCH4 ]
            print get_OH_mean( wd ) * 1E5   
            print [ ( i.shape, i.sum(), i.mean(), i.min(), i.max()) \
                for i in  LCH4, LCH4_Cl, Lrate_CH4_Cl, Cl ]

        if average_value:
            # get mass weighed average
            LCH4 =  (LCH4*a_m).sum()/a_m.sum() 
    #        arr = np.ma.mean(  arr  )
            if include_Cl_ox:
                LCH4_Cl =  (LCH4_Cl*a_m).sum()/a_m.sum() 
                print [ ( i.shape, i.sum(), i.mean(), i.min(), i.max()) \
                    for i in  LCH4_Cl, Lrate_CH4_Cl, Cl ]

        # return CH4 lifetime 
        CH4_tau_s =  LCH4

        if include_Cl_ox:
            CH4_tau_s = np.ma.sum( [LCH4_Cl, LCH4] ) 
        print 'Using 1st apporach: ', LCH4, LCH4_Cl, CH4_tau_s

        return 1/CH4_tau_s /60/60/24/365 

    else: # Second Approach ( as per Schmidt et al 2015)

        # CH3CCl loss rate per grid box  ( cm^3  molec.^-1  s^-1  )
#        LCH3CCl = 1.64E-12 * np.ma.exp( (-1520. / K)  )   

        # Mask non tropospheric boxes ( use: tropopause level )
        print [i.shape for i in K, vol, a_m, n_air, KCH4, OH  ]
        if masktrop:
            K, vol, a_m, n_air, KCH4, OH = mask4troposphere( \
                [K, vol, a_m, n_air, KCH4, OH ], t_lvl=t_lvl,  \
                use_time_in_trop=use_time_in_trop )                
            if include_Cl_ox:
                Lrate_CH4_Cl = mask4troposphere( [Lrate_CH4_Cl], t_lvl=t_lvl,  \
                    use_time_in_trop=use_time_in_trop )[0]     

        # get # molecules per grid box
        molecs = n_air *vol

        # Get loss with respect to OH
        # (( cm^3  molec.^-1  s^-1  )* [molec/cm3]  =  [s^-1] )
        LCH4 = KCH4 * OH
        if include_Cl_ox: # Get loss with respect to Cl ? 
            # ( 1/ ( [molec/cm3]/[molec/cm3/s] ) =  [s^-1]  )
#            LCH4_Cl =  1/np.ma.divide( Cl, Lrate_CH4_Cl  )
            # ( 1/ ( [molec/cm3]/[molec/cm3/s] ) =  [s^-1]  )
            LCH4_Cl =  1/np.ma.divide( CH4, Lrate_CH4_Cl  )

        if debug:
            print [ ( i.shape, i.sum(), i.min(), i.max(), i.mean() )  \
                for i in [ LCH4, KCH4, molecs ] ]
            print [ ( i.mean(), i.min(), i.max()) \
                for i in  LCH4, LCH4_Cl, Lrate_CH4_Cl, Cl ]

        if average_value:

            # Get weighted values per time step
            LCH4_l, LCH4_Cl_l = [], []
            for n in range( LCH4.shape[-1] ):
    
                # Weight by molecules
                LCH4_l += [ ( LCH4[...,n] * molecs[...,n] ).sum()  \
                    / molecs[...,n].sum() ]
                if include_Cl_ox: # Get loss with respect to Cl ? 
                    LCH4_Cl_l += [ ( LCH4_Cl[...,n] * molecs[...,n] ).sum()  \
                        / molecs[...,n].sum() ]

            if debug:
                print LCH4_l

            # take mean
            LCH4 = np.ma.mean( LCH4_l )
            if include_Cl_ox: 
                LCH4_Cl = np.ma.mean( LCH4_Cl_l )

        if debug:
            print [ ( i.shape, i.sum(), i.min(), i.max(), i.mean() )  \
                for i in [ LCH4, KCH4, molecs ] ]
            
        # Convert units from /s to yrs
        CH4_tau_s = LCH4 

        if include_Cl_ox:
            CH4_tau_s = np.ma.sum( [LCH4_Cl, LCH4] ) 
        print 'Using 2nd apporach: ', LCH4, LCH4_Cl, CH4_tau_s

        return 1/ CH4_tau_s / (3600*24*365)

# ----
# 2.13 - get Land / Water /Ice fraction
# ---- 
def get_LWI(lon, lat, res='4x5',debug=False):
    """ Return LWI for a given lon and lat """
    lat=get_gc_lat(lat, res=res)
    lon=get_gc_lon(lon, res=res)
    LWI=get_land_map(res, res=res)
    return LWI[lon,lat,0]
    
# -------------- 
# 2.14 - get_gc_alt  ( KM => box num )   
# -------------
def get_gc_alt(alt):
    """ 
    Return index of nearest altitude (km ) of GEOS-Chem box 
    """
    alt_c = gchemgrid('c_km_geos5_r')
    return find_nearest( alt_c, alt )

# --------------
# 2.15 - species arrary (4D) in v/v to Gg  I/ O3/Species
# ------------- 
def species_v_v_to_Gg(arr, spec, a_m=None, Iodine=True, All =False, \
        Ox=False, ctm_f=None, wd=None, debug=False):
    """ 
    Convert array of species in v/v to Gg (of spec)
    
    NOTE(s):
     - The processing is not clear in this function, NEEDS UPDATE
     - To output in terms of provide spec (e.g. Bry ) set All=True and 
      Iodine=False. This should be default bevahiour, but is not for reasons 
      of back compatibility.
    """
    var_list = All, Iodine, Ox, spec
    print 'WARNING: Check settings in species_v_v_to_Gg: ',  var_list

    if not isinstance(a_m, np.ndarray):
        a_m = get_air_mass_np( ctm_f, wd=wd, debug=debug )  # kg
    #  number of moles in box ( g / g mol-1 (air) )
    moles = ( a_m *1E3 ) / constants( 'RMM_air')   
    if debug:
        print [ i.shape for i in moles, arr ]

    # I mass ( moles => mass (g) => Gg (in units of I ) )
    if ( (Iodine) and ( (not Ox) and (not All) ) ):  
        arr = ( ( arr * moles )  * 127. ) /1E9 * spec_stoich( spec )

    # O3 mass ( moles => mass (g) => Gg (in units of O3 ) )
    if ( (Iodine) and (Ox) ):      
        arr = ( ( arr * moles )  * (16.*3.) ) /1E9 

    #  In "species" mass terms ( moles => mass (g) => Gg (in units of spec ) )
    if ( ( not Iodine ) and ( All) ):      
        arr = ( ( arr * moles )  * (species_mass( spec  )) ) /1E9 
    return arr

# --------------
# 2.16 - retrive volume for geos chem run ( in cm3 )
# -------------
def get_volume_np(ctm_f=None, box_height=None, s_area=None, res='4x5', \
            wd=None, trop_limit=False, debug=False):
    """ 
    Get grid box volumes for CTM output in cm3 
    """
    if debug:
        print 'get_volume_np called'

    if not isinstance(box_height, np.ndarray):
        if pygchem.__version__ == '0.2.0':
            box_height =  get_gc_data_np(ctm_f, 'BXHEIGHT',\
                category="BXHGHT-$", debug=debug)  # ( m )
        else:
            box_height = get_GC_output( wd=wd, vars=['BXHGHT_S__BXHEIGHT'], \
                debug=debug )

            # Gotcha for single (None) time dimension output:
            # <= improve this - this only works for 4x5 resolution
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
    """ 
    Return Boolean for if grid box is water or not 
    """

    # Load masked array, where all ocean is non-masked
    # ( Just look over oceans (use masked array) )
    o_mask = np.ma.equal( get_land_map(res=res)[...,0], 0) 

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
# 2.18 - Get dry dep for given spec
# -------------   
def spec_dep(ctm_f=None, wd=None, spec='O3', s_area=None, months=None, \
                years=None, res='4x5', vol=None, debug=False, \
                trop_limit=True, Iodine=False):
    """ 
    Get array of dry deposition values for a given species

    NOTES:
     - Values are returned as a spatial arry with a time dimension 
    (e.g. shape= (72, 46, 12) )
    """
    if debug:
        print 'spec dep called for: ', spec

    # Get surface area if not provided
    if not isinstance(s_area, np.ndarray):
        s_area =  get_surface_area( res )  # m2 land map                                                

    # Extract dry dep flux in  [molec/cm2/s]
    if pygchem.__version__ == '0.02':
        df = get_gc_data_np( ctm_f, spec=spec+'df', category='DRYD-FLX', \
            debug=debug) 
    else:
        df = get_GC_output( wd, category='DRYD-FLX', species=spec+'df' )

    if debug:
        print '*'*10,[( i.shape, np.sum(i), np.mean(i)) for i in [df] ], len(df)

    # Convert to Gg "Ox" (Gg X /s)
    df = molec_cm2_s_2_Gg_Ox_np( df, spec, s_area=s_area, \
                Iodine=Iodine, res=res, debug=debug ) 

    if debug:
        print '0'*20, [( i.shape, np.sum(i), np.mean(i)) for i in [df]], len(df)

    if isinstance( months, type(None) ):
        months = get_gc_months( ctm_f=ctm_f, wd=wd )
    if isinstance( years, type(None) ):
        years = get_gc_years( ctm_f=ctm_f, wd=wd, set_=False )

    # Adjust time dimension (to per month)
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
    """ 
    Convert species/tag prod/loss from "molec/cm3/s" to "Gg (Ox) yr^-1".

	!!! This is redundent. use convert_molec_cm3_s_2_g_X_s instead. !!!

    NOTES:
     - This function was originally used to process diagnostic outputs from PORL-L$  
     - Alterate Output term apart from Ox are possible ( e.g. I, mass ... )
        Just stipulate input species  
     - Output can also be requested on a monthly basis 
        ( set  month_eq=True,  and year_eq=False )

    """ 

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

    else: # ELSE return Gg X yr^-s ( ASSUMING UNITY OF TAG! )
        if ( rxn ==  'CH4' ):
            RMM =  species_mass(rxn) 
        else:
            RMM =  species_mass('O3')             
        # Return mass in terms of species provide ( if given )
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
            Iodine=False, res='4x5', year_eq=False, #fix_Kludge=True, \
            debug=False):
    """ 
    Convert 2D depositional array from [molec/cm2/s] to Gg Ox yr^-1

	NOTE(s):
	 -  NEEDS UPDATE. The code is not clear. 
    """
    
    if debug:
        print ' molec_cm2_s_2_Gg_Ox_np  called'

    # Kludge (loads all output surface areas, 
    # only requires a single time dimensions as does not chnage... )
    s_area = s_area[...,None]  # tmp bug fix 
    #  anti-kludge 
#    if fix_Kludge:
#    s_area = s_area[...,0]
        
    if debug:
        print '_'*10, arr.shape, s_area.shape
    if ( Iodine ):
        # [molec/cm2/s] * cm2 (m*1E4) /  moles * RMM I /1E9 ( convert Gg I /s )
        arr  = ( arr * (s_area*1E4) ) / constants( 'AVG') * (127.) * \
            spec_stoich(spec) /1E9   
    elif spec == 'O3':
        arr  = ( arr * (s_area*1E4) )  / constants( 'AVG')* \
            species_mass(spec) /1E9  * Ox_in_species( spec)
    else:
        arr  = ( arr * (s_area*1E4) )  / constants( 'AVG') * \
            species_mass(spec) /1E9  
    if year_eq:
        print "giving 'year_eq' "
        arr = arr* 60*60*24*365  

    if debug:
        print 'arr', arr.shape
    return arr

# --------------
# 2.22 - Get Dry dep - rm as redundent.
# -------------
 
# --------------  
# 2.23 - Get run data for a given GC run - 
# -------------   
# Not generic enough for module - mv'd to bottom 

# --------------
# 2.24 - Get DU mean value
# -------------
def get_DU_mean(s_area=None, a_m=None, t_p=None, O3_arr=None, \
            ctm_f=False, area_weight=True, res='4x5', debug=False ):
    """ 
    Get mean DU value weighed by grid box area 
    
    ARGUMENTS:
     - t_p: time in the troposphere diganostic ( float values 0 to 1 )
     - res: resolution of the model input (e.g. 4x5, 2x2.5 ) 
     - ctm_f: ctm.bpch file object (vestigle from PyGChem <3.0)
     - s_area: Surface array (array of values in metres)
     - area_weight (boolean): weight by area of grid box?
    """

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
# REDUNDENT - mv'd to bottom of module


# --------------
# 2.26 - Get Prod / loss for O3
# -------------
def get_POxLOx( ctms=None, vol=None, all_data=False, t_p=None, ver='1.6', \
    wd=None, debug=False):
    """ 
    Get production and loss terms for O3 from prod/loss diagnostic 

    ARGUMENTS:
     - t_p: time in the troposphere diganostic ( float values 0 to 1 )
     - res: resolution of the model input (e.g. 4x5, 2x2.5 ) 
     - ver: GEOSChem version with halogens (default = 3.0), ignore if not using halogen code
     - all_data(boolean): return the full data (instead of two single numbers)
     - ctms: list of ctm.bpch file objects (vestigle from PyGChem <3.0)
    """

    specs = GC_var('POxLOx')

    if debug:
        print ver 

    # Get prod/loss in [molec/cm3/s]
    if pygchem.__version__ ==  '0.2.0':
        arrs = [ np.concatenate( [get_gc_data_np( ctm, spec=PLO3_to_PD(spec, \
            ver=ver, fp=True), category="PORL-L=$", debug=debug) \
            for ctm in ctms ], axis=3 ) for spec in specs ]
    else:
        arrs = get_GC_output( wd=wd, vars=['PORL_L_S__'+i for i in specs] )
        arrs = [ arrs[i,...] for i in range( len(specs) ) ]

    if all_data:
        # [molec/cm3/s] => Gg Ox / month 
        months = [ get_gc_months( ctm) for ctm in ctms ]
        years =  [ get_gc_years( ctm, set_=False ) for ctm in ctms ]
        months = [j for i in months for j in i ]
        years = [j for i in years for j in i ]        
        arrs = [ molec_cm3_s_2_Gg_Ox_np(arr, specs[i], vol=vol, months=months, \
                    years=years, debug=debug, wd=wd, \
                    month_eq=True, year_eq=False) \
                    for i, arr in enumerate(arrs) ] 
        return [ arrs[i] for i in range(len(specs )) ]  # Gg

    else:
        # [molec/cm3/s] => Gg Ox / yr
        arrs = [ molec_cm3_s_2_Gg_Ox_np(arr, specs[i], vol=vol, wd=wd, debug=debug)\
        	 for i, arr in enumerate(arrs) ]
        # get yearly mean + remove stratosphere
        # NEEDS UPDATE - update troposphere removal method?
        arrs = [ (arr*t_p).mean(axis=3) for arr in arrs ] 

        return [ int( np.ma.masked_invalid( arrs[i] ).sum()/1E3)  \
            for i in range(len(specs )) ] # Tg


# --------------
# 2.28 - Get wet dep
# -------------
def get_wet_dep( ctm_f=None, months=None, years=None, vol=None, \
        scale=1E9, s_area=None, res='4x5', wd=None, specs=None, \
        Iodine=False, all_wdep=False, sep_rxn=False, ref_spec=None, \
        ver='1.6', trop_limit=True, debug=False):
    """ 
    Extract wet deposition for given species in terms of X Gg of ref_spec. 

    ARGUMENTS:
     - ctm_f: ctm.bpch file object (vestigle from PyGChem <3.0)
     - ver: GEOSChem version with halogens (default = 3.0), ignore if not using halogen code
     - s_area: Surface array (array of values in metres)
     - vol: volumne of grid boxes
     - trop_limit: limit output to "chemical troposphere" (level 38 )
     - return in terms of mass of iodine
     - scale: scaleing to use (e.g. 1E9 = Gg )
     - sep_rxn: return list of arrays by "reaction"/process
     - res: resolution of the model input (e.g. 4x5, 2x2.5 )  

    NOTES: 
     - If Iodine=True, re_spec == 'I'
     - this fucntion expects list of species
    ( if a species list is not provided, then Iodine deposition is returned)
     - A reference species is expected to define the unit terms of the 
    returned values. 
    """

    logging.info( 'get_wet_dep called for {}, with ref_spec: {} (Iodine:{})'.format( 
            specs, ref_spec, Iodine ) )

    if not isinstance(months, list):
        months = get_gc_months( ctm_f, wd=wd )
    if not isinstance(years, list):
        years  = get_gc_years( ctm_f=ctm_f, wd=wd, set_=False )
    if not isinstance(vol, np.ndarray):
        vol = get_volume_np( ctm_f=ctm_f, s_area=s_area, wd=wd, res=res,\
             debug=debug )
    w_dep = GC_var('w_dep')
    
    if Iodine and isinstance( specs, type(None) ):
        if (all_wdep):
            specs = GC_var('w_dep_specs')
            if ver =='3.0':
                specs += ['ISALA', 'ISALC' ]
        else:
            specs = GC_var('w_dep_specs')[:-1] # skip AERI - Kludge
            if ver =='3.0':
                specs = GC_var('w_dep_specs')[:-2] # skip ISALA/ISALC
        ref_spec = 'I'

    # Return 
    if isinstance( ref_spec, type(None) ):
        print 'PLEASE PROVIDE REF SPEC'
        print 'CHECK implementation of ref_spec  - UPDATED'
        sys.exit()

    if debug:
        print specs , len(specs)

     # --- Extract and convert [kg/s] => g / s of ref_spec ( e.g. I equiv. )
    # Retain back compatibility
    if pygchem.__version__ == '0.2.0':
        dep_w = [ [ get_gc_data_np(ctm_f, spec ,category=var, \
            debug=debug )*1E3 /species_mass(spec)*  \
            species_mass(ref_spec)*spec_stoich(spec) \
            for spec in specs[ii] ]  \
            for ii,var in enumerate(w_dep)]  

    else:
        # Frontal rain? [kg/s]
        WETDLS_S__ = get_GC_output( wd=wd, vars=['WETDLS_S__'+i \
            for i in specs ], r_list=True, trop_limit=trop_limit )

        # Get convective scavenging  [kg/s]
        WETDCV_S__ = get_GC_output( wd=wd, vars=['WETDCV_S__'+i \
            for i in specs ], r_list=True, trop_limit=trop_limit )
        if debug:
            print [ [ i.shape for i in l ] for l in WETDLS_S__, WETDCV_S__  ]

        # Convert to g/ s X(ref_spec) equiv.  + conbine two lists
        dep_w = []
        for n, spec in enumerate( specs ):
            if debug:
                print n, spec, ref_spec, spec_stoich( spec, ref_spec=ref_spec) 

            # Combine convective scavenging and Frontal rain
            dep_w += [ WETDLS_S__[n] + WETDCV_S__[n] ] 

            # Convert [kg/s] to [g], then moles of species 
            dep_w[n] = dep_w[n]*1E3/species_mass( spec )
                    
            # Convert to unit terms of ref_spec
            dep_w[n] = dep_w[n]* species_mass(ref_spec)

            # Adjust to stoichiometry of ref_spec
            dep_w[n] = dep_w[n]*spec_stoich( spec, ref_spec=ref_spec) 

    # Adjust to monthly values...
    m_adjust = d_adjust( months, years) # => Gg / month 
    dep_w = [m_adjust *i for i in dep_w ]

    # List wet dep rxns, concat. and sum of rxns ( and convert to Gg )
    if sep_rxn: 
        return np.concatenate( [ i[None,...] /scale \
                            for i in dep_w ], axis=0 )

    # List wet dep rxns and concat. 
    else:   
        return np.concatenate( [ i[...,None] /scale \
                        for i in dep_w ], axis=-1 ).sum(axis=-1 )


# --------
# 2.31 - Vol weighted array average value
# --------
# REDUNDENT 

# ----
# 2.32 - molecule weighted array average value
# ----
def molec_weighted_avg( arr, wd=None, ctm_f=None, vol=None, t_p=None, n_air=None, \
        trop_limit=True, multiply_method=False, rm_strat=True, molecs=None,\
        weight_lon=False, weight_lat=False, LON_axis=0, LAT_axis=1, \
        annual_mean=True, debug=False ):
    """ Takes array and retuns the average (molecular weighted) value 

    ARGUMENTS:
     - weight over longitudinal (weight_lon=True) or latitudinal(weight_lat=True)
     - annual_mean (boolean): average the time axis?
     - n_air (array): number desnity of air
     - molecs (array): number of molecules in air
     - vol: volumne of grid boxes
     - trop_limit: limit output to "chemical troposphere" (level 38 )
     - res: resolution of the model input (e.g. 4x5, 2x2.5 )  
     - t_p: time in the troposphere diganostic ( float values 0 to 1 )
     - trop_limit: limit output to "chemical troposphere" (level 38 )
     ( (boolean) options for this include using 
     - definition of troposphere to use? 
        - use_time_in_trop (boolean): time a given box is in the troposphere
        ( if use_time_in_trop=False, the level of the troposphere is used )
     - multiply_method (boolean): use a multiplication method, rather than masking for
      the arrays (aka set stratosphere to have zero values)
     - Index of longitudinal (LON_axis) and latitudinal (LAT_axis) axis

    NOTES:
        - Uses mask of given array if given array is a numpy masked array
        - Assume axis dimensions are  (LON, LAT, ALT ), aka the array does not 
        have a TIME dimension. Any shape can be handled as long as given molecs 
        and arr have the same shape
        - Molecs is the same as n_air ( [molec air/m3] * [m3]  vs.   
       air mass [kg] / RMM [kg mol^-1] * Avogadros [mol^-1] )
    """

    if isinstance( molecs, type(None) ): 
        if not isinstance( n_air, np.ndarray): 
            n_air = get_GC_output( wd, vars=['BXHGHT_S__N(AIR)'], \
                trop_limit=trop_limit, dtype=np.float64 )  # [molec air/m3]
        if not isinstance(vol, np.ndarray):
            vol = get_volume_np( ctm_f=ctm_f, wd=wd, trop_limit=trop_limit ) 
            vol = vol /1E6 # [cm^3 ]
        # Calculate molecules per grid box
        molecs =  n_air * vol # [molec air]
        if annual_mean:
            molecs = molecs.mean(axis=-1)

    # Limit for troposphere? 
    if trop_limit:

        # Get species time Tropopause diagnostic
        if isinstance( t_p, type(None) ) and rm_strat:
            t_p = get_GC_output( wd=wd, vars=['TIME_TPS__TIMETROP'], \
                trop_limit=trop_limit ) 

            # mask for troposphere
            arr, molecs = mask4troposphere( [arr, molecs], t_ps=t_p,  \
                use_time_in_trop=True,  multiply_method=multiply_method)      

    # --- If masked array provided, applied same mask to molecules
    if isinstance( arr, np.ma.core.MaskedArray ):
        try:
            molecs  = np.ma.array( molecs, mask=arr.mask )
        except:# MaskError:
            print "MaskError for array shapes in 'molec_weighted_avg': ", \
                [ i.shape for i in molecs, arr ]

    if weight_lon and (not weight_lat):  # 1st axis  
        return (arr *molecs).sum(axis=LON_axis)/molecs.sum(axis=LON_axis)

    elif weight_lat and (not weight_lon):  # 2nd axis (LON, LAT, ALT, TIME )
        return (arr *molecs).sum(axis=LAT_axis)/molecs.sum(axis=LAT_axis)
    
    elif weight_lat and weight_lon:  # 1st+2nd axis (LON, LAT, ALT, TIME )
        return (arr *molecs).sum(axis=LAT_axis).sum(axis=LON_axis)/  \
            molecs.sum(axis=LAT_axis).sum(axis=LON_axis)
  
    else: # weight whole array to give single number
            return (arr *molecs).sum()/molecs.sum()



# --------------
# 2.35 - Get Tropospheric Ox loss routes as list of arrays, + model resolution (res)
# -------------
# Not generic enough for module - mv'd to bottom 

# --------------
# 2.36 - Split Tropospheric Ox loss routes 
# -------------
# Not generic enough for module - mv'd to bottom 

# --------------
# 2.37 - Return only dry/wet dep variaibles aloocated within model output arrays
# -------------       
# Redundant


# --------------
# 2.38 - split 4D ( lon, lat, alt, time) output by season
# -------------       
def split_4D_array_into_seasons( arr, annual_plus_seasons=True, \
            debug=False ):
    """ 
    Split 4D ( lon, lat, alt, time) output by season, then take 
    average of the seasons 

    NOTE(s): 
     - currently seasons index is mannual set assuming Jan-Dec 
    """
    
    if debug:
        print arr.shape 
    
    if annual_plus_seasons:
        seasons  = ['Annual', 'DJF', 'MAM', 'JJA', 'SON']
        # assume calender month order 
        # ( this can be automated using get_GC_datetime )
        indices = [range(0,12), [11, 0, 1], [2, 3, 4], [5, 6, 7], [8, 9, 10] ]
             
    else:
        seasons  = ['DJF', 'MAM', 'JJA', 'SON']
        # assume calender month order 
        # ( this can be automated using get GC datetime )
        indices = [[11, 0, 1], [2, 3, 4], [5, 6, 7], [8, 9, 10]]

    ars = []
    # extract data by month
    for n, s in enumerate( seasons ):
        if debug:
            print s, n , indices[n] 
        ars += [ [ arr[...,i] for i in indices[n] ] ]
    
    # average by season
    if debug:
    #    print [ np.array(i).shape for i in ars ], np.array(ars).mean(), seasons
        print [ np.ma.array(i).shape for i in ars ], seasons
    ars =  [ np.ma.array(i).mean(axis=0) for i in ars]
    if debug:
        print [ i.shape for i in ars ], np.ma.array(ars).mean(), seasons

    # return list array averaged by season 
    return ars, seasons

# --------------
# 2.39 - Convert v/v to ng/m^3
# -------------     
def convert_v_v2ngm3( arr, wd=None, spec='AERI', trop_limit=True, \
            s_area=None, vol=None, a_m=None, res='4x5', debug=False ):
    """ 
    Take v/v array for a species, and conver this to mass loading
    units used as standard are ng/m3

    ARGUMENTS:
     - spec: species to convert to (uses mass of this species)
     - a_m (array): array of air mass
     - vol: volumne of grid boxes
     - trop_limit: limit output to "chemical troposphere" (level 38 )
     - res: resolution of the model input (e.g. 4x5, 2x2.5 )  
    """

    # Get volume (m^3, adjusted (x1E6) from cm^3) 
    if not isinstance(vol, np.ndarray):
        vol = get_volume_np( wd=wd, trop_limit=trop_limit, s_area=s_area, \
                    res=res ) /1E6  

    # Get air mass ( kg )
    if not isinstance(a_m, np.ndarray):
        a_m = get_GC_output( wd, vars=['BXHGHT_S__AD'], trop_limit=trop_limit, 
                    dtype=np.float64)

    # Get moles  ( converting airmass from kg 1st)
    mols = a_m*1E3/constants( 'RMM_air')

    # Adjust to mols, then mass 
    arr = arr*mols*species_mass( spec )
    if debug:
        print species_mass( spec ), np.sum( arr )

    # Convert to (nano, x1E9)g/m3
    arr = arr*1E9/vol 
    
    return arr
    
# --------------
# 2.40 - Print weighted array 
# -------------   
def prt_seaonal_values( arr=None, res='4x5', area_weight=True, zonal=False, \
            region='All', monthly=False, mask3D=True, trop_limit=True, \
            prt_by_3D_region=False, hPa=None, wd=None, \
            verbose=True, debug=False ):
    """ 
    Provide zonal/surface area weighted values  
    """

    if verbose:   
        print 'function prt_seaonal_values called for region: ', region, \
            'debug: {},verbose: {}'.format( debug, verbose )

    print [ (i.min(), i.max(), i.mean() ) for i in  [arr] ]
    print [ (i.min(), i.max(), i.mean() ) for i in [ arr.mean(axis=-1)[...,0] ]]

    # Get surface area
    s_area =  get_surface_area( res=res ) # m2 land map
    s_area = s_area[...,0]

    months = range(12)

    # --- If region provided, mask elsewhere - else
    if ( 'asked' not in str( type( arr ) ) ):
        print 'WARNING: converting array to masked array'
        arr = np.ma.array( arr ) 
    if debug:
        print [ ( i.min(), i.max(),i.mean(), len(i.mask[i.mask==True].ravel()) ) 
                        for i in [arr] ]
    m = mask_all_but( region=region, mask3D=True, debug=debug,\
                use_multiply_method=False, \
                trop_limit=trop_limit )[...,:38]
    print [i.shape for i in m, arr.mask ]
    m = np.ma.concatenate( [m[...,None]] *len( months), axis=-1 )
    # Mask array individually 
    print [i.shape for i in m, arr.mask, arr ]
    arr = np.ma.array( arr, mask=np.ma.mask_or(m, arr.mask) ) 
    #    ars = [ np.ma.array( i, mask=m ) for i in ars ]
    if debug:
        print m, region, [ (type(i), i.shape) for i in [ m ] +[arr] ], \
                    arr.mask, [ ( i.min(), i.max(),i.mean() ) for i in [arr] ]

    if debug:
        print [ ( i.min(), i.max(),i.mean(), len(i.mask[i.mask==True].ravel()) ) 
                        for i in [arr] ]

    # Also mask surface to also for area weighting of data
    s_area = np.ma.array( s_area, mask=m[...,0,0] )

    # --- Split array by seasons ( on months if monthly==True)
    if monthly:
        # this is assuming monthly output 
        ars = [ arr[...,i] for i in range( 12 ) ]
        seasons = num2month( rtn_dict=True).values() 
        # Also plot annual value
        ars += [ arr.mean( axis=-1 ) ]
        seasons += ['Annual']
    else:
        ars, seasons = split_4D_array_into_seasons( arr, \
                                        annual_plus_seasons=True ) 

    # --- Print values by 3D region
    if prt_by_3D_region:

        # Get pressure levels
        if isinstance( hPa, type(None) ):
            hPa = get_GC_output( wd=wd, vars=[u'PEDGE_S__PSURF'] )
            hPa = hPa.mean(axis=-1)[...,:38]

        # Setup makes for BL, FT, UT
        regions = [ 'MBL','BL','FT',  'UT']
        sects3D = [ mask_3D( hPa, i, res=res )[:,:,:38] for i in regions]
        print [i.shape for i in sects3D ]
    
        # --- Print values by 3D region 
        pstr  = '{:<20}'*( len(regions) +1 )
        npstr  = '{:<20}'+ '{:<20,.3f}'*len(regions)
        print pstr.format( 'spec' , *regions )
        
        print seasons, ars[0].mean() 
        # Mean of 3D region
        for n, s in enumerate( seasons ):
            vars = [ float( (i*ars[n]).mean() ) for i in sects3D] 
            print npstr.format( s, *vars )

        # Min of 3D region
        for n, s in enumerate( seasons ):
            vars = [float( (i*ars[n]).min() ) for i in sects3D] 
            print npstr.format( s, *vars )

        # Max of 3D region
        for n, s in enumerate( seasons ):
            vars = [ float( (i*ars[n]).max() ) for i in sects3D] 
            print npstr.format( s, *vars )
        
    # ---- Get surface or Zonal array ( prt by 2D region )
    else:
        if zonal:
            ars = [i.mean(axis=0) for i in ars ]
        # If not zonal return surface values
        else:    
#            ars = [ i[...,0] for i in ars ]
            ars = [ np.ma.array( i[...,0], mask=m[...,0,0] ) for i in ars ]
        if debug:
            print [i.shape for i in ars ], s_area.shape, \
                    [type(i) for i in ars, ars[0],  s_area ]


        # --- Print out 
        header = [
         'Season./Month', 'Min', 'Max','5th','95th', 'Mean' , 'Median', \
         "wtg'd Mean"
        ]
        pstr =  '{:<15}'*(len( header))
        npstr =  '{:<15}'+'{:<15,.4f}'*(len( header)-1)
        print pstr.format( *header )
        for n, s in enumerate( seasons ):

            # Get vars for printing 
            vars = [ 
            ( float( i.min() ), float( i.max() ),  \
            float( np.percentile( i.compressed(), 5)),  \
            float( np.percentile( i.compressed(), 95) ), \
            float( i.mean() ), \
            float( np.median( i.compressed() ) ),   \
            float( (i*s_area).sum()/s_area.sum() )   ) \
            for i in [ ars[n] ] 
            ][0]
            # Print vars
            print npstr.format( s,  *vars )


# --------------
# 2.41 - Extract data by family for a given wd
# -------------   
def fam_data_extractor( wd=None, fam=None, trop_limit=True, ver='3.0', \
        annual_mean=True, t_ps=None, a_m=None, vol=None, \
        title=None, rtn_list=False, use_time_in_trop=True, multiply_method=True, \
        rtn_specs=False, verbose=False, debug=False ):
    """ 
    Driver to extract data requested ( as families have differing diagnostic units)
    
    ARGUMENTS:
     - fam: "family" to extract ( e.g. NOy, NOx, POx, CH4 loss rate, ... )
     - a_m (array): array of air mass
     - vol: volumne of grid boxes
     - trop_limit: limit output to "chemical troposphere" (level 38 )
     - res: resolution of the model input (e.g. 4x5, 2x2.5 )  
     - rtn_specs (boolean): return list of species extracted?
     - t_ps: time in the troposphere diganostic ( float values 0 to 1 )
     - trop_limit: limit output to "chemical troposphere" (level 38 )
     ( (boolean) options for this include using 
     - definition of troposphere to use? 
        - use_time_in_trop (boolean): time a given box is in the troposphere
        ( if use_time_in_trop=False, the level of the troposphere is used )
     - use a multiplication method, rather than masking for the arrays
     (aka set stratosphere to have zero values)
     - title: run title, ignore if not using halogen code
     - ver: GEOSChem version with halogens (default = 3.0), ignore if not using halogen code

    NOTES:
     - to return species as list ( not fully implimented ) set rtn_list=True
     - to return species extract, set  rtn_species=True
     - this function should be used in preference to other bulk output extractors 
      in this module.
    """
    if verbose:
        print 'fam_data_extractor called for ', fam, wd, title

    # --- Nitrogen Oxides NOx ( NO + NO2)
    if fam == 'NOx' :
        # Select species in family
        specs = [ 'NO2', 'NO' ] 
        scale = 1E12
#        units, scale = tra_unit(specs[0], IUPAC_unit=True, scale=True)
        # Extract data
        arr = get_GC_output( wd=wd, vars=['IJ_AVG_S__'+i for i in specs ], \
                    trop_limit=trop_limit )
        if debug:
            print [ ( i.shape, i.min(), i.max(), i.mean() ) for i in [arr ] ]
        if not rtn_list:
            arr = np.ma.concatenate( [ i[...,None] for i in arr ], axis=-1 )
            arr = arr.sum( axis=-1) * scale

    # --- OH ( in molec/cm3 )
    if fam == 'OH' :
        spec = 'OH'
        arr = get_GC_output( wd=wd, vars=['CHEM_L_S__'+spec], \
            trop_limit=trop_limit )

    # --- HOX 
    if fam == 'HOx' :
        # OH ( in molec/cm3 )
        arr = get_GC_output( wd=wd, vars=['CHEM_L_S__'+'OH'], \
            trop_limit=trop_limit )
        # Convert OH to  v/v 
        arr = convert_v_v_2_molec_cm3( arr, a_m=a_m, vol=vol )
        # Get HO2 
        arr2 = get_GC_output( wd=wd, vars=['CHEM_L_S__'+'HO2'], \
            trop_limit=trop_limit )
        # HOx ( HO + HO2)
        arr = arr + arr2

    # --- NOy
    if fam == 'NOy' :
        # Select species in family
        specs = GC_var('N_specs' )
        if ver == '3.0':
            if not any( [ (title != i) for i in 'NOHAL', 'BROMINE' ]):
                specs  += ['ClNO2', 'ClNO3'] 
            
#        units, scale = tra_unit(specs[0], IUPAC_unit=True, scale=True)
        scale =1E12
        # Extract data
        arr = get_GC_output( wd=wd, vars=['IJ_AVG_S__'+i for i in specs ], \
                    trop_limit=trop_limit, r_list=True  )
        # Adjust to stoichiometry
        arr = [ arr[n]*spec_stoich(i, ref_spec='N') \
             for n,i in enumerate( specs ) ]
                
        if debug:
            print [ ( i.shape, i.min(), i.max(), i.mean() ) for i in arr  ]
        if not rtn_list:
            arr = np.ma.concatenate( [ i[...,None] for i in arr ], axis=-1 )
            arr = arr.sum( axis=-1) * scale

    # --- Ozone (O3)
    if fam == 'O3' :
        # Select species in family
        units, scale = tra_unit(fam, IUPAC_unit=True, scale=True)
        # Extract data
        arr = get_GC_output( wd=wd, vars=['IJ_AVG_S__'+fam ], \
                    trop_limit=trop_limit )
        if debug:
            print [ ( i.shape, i.min(), i.max(), i.mean() ) for i in [arr ] ]
        arr = arr *scale

    # Get Ox prod (POx/POX) ([molec/cm3/s])
    if fam == 'POX' :
        # Select species in family
        arr = get_GC_output( wd=wd, vars=['PORL_L_S__'+fam], 
                    trop_limit=trop_limit )

    # Get Ox prod (POx/POX) ([molec/cm3/s])
    if fam == 'CH4 loss' :
        # Select species in family
        arr = get_CH4_lifetime( wd=wd, use_OH_from_geos_log=False,\
            average_value=False )
        # Want in units of yr^-1
        arr = 1/arr

    # ---  Inorganic bromine ( Bry )
    if fam == 'Bry' :
        # Select species in family
        specs = GC_var('Bry' )
        # Also consider Br- on SS
#        specs += [ 'BrSALA', 'BrSALC', ]

        # Extract data
        arr = get_GC_output( wd=wd, vars=['IJ_AVG_S__'+i for i in specs ], \
                    trop_limit=trop_limit, r_list=True  )
        # Adjust to stoichiometry
        arr = [ arr[n]*spec_stoich(i, ref_spec='Br') \
             for n,i in enumerate( specs ) ]
                
        if debug:
            print [ ( i.shape, i.min(), i.max(), i.mean() ) for i in arr  ]
        if not rtn_list:
            arr = np.ma.concatenate( [ i[...,None] for i in arr ], axis=-1 )
            arr = arr.sum( axis=-1) 

    # ---  Inorganic iodine ( Iy )
    if fam == 'Iy' :
        # Select species in family
        specs = GC_var('Iy' )

        # Extract data
        arr = get_GC_output( wd=wd, vars=['IJ_AVG_S__'+i for i in specs ], \
                    trop_limit=trop_limit, r_list=True  )
        # Adjust to stoichiometry
        arr = [ arr[n]*spec_stoich(i, ref_spec='I') \
             for n,i in enumerate( specs ) ]
                
        if debug:
            print [ ( i.shape, i.min(), i.max(), i.mean() ) for i in arr  ]
        if not rtn_list:
            arr = np.ma.concatenate( [ i[...,None] for i in arr ], axis=-1 )
            arr = arr.sum( axis=-1) 

    # ---  Inorganic iodine ( Cly )
    if fam == 'Cly' :
        # Select species in family
        specs = GC_var('Cly' )

        # Extract data
        arr = get_GC_output( wd=wd, vars=['IJ_AVG_S__'+i for i in specs ], \
                    trop_limit=trop_limit, r_list=True  )
        # Adjust to stoichiometry
        arr = [ arr[n]*spec_stoich(i, ref_spec='Cl') \
             for n,i in enumerate( specs ) ]
                
        if debug:
            print [ ( i.shape, i.min(), i.max(), i.mean() ) for i in arr  ]
        if not rtn_list:
            arr = np.ma.concatenate( [ i[...,None] for i in arr ], axis=-1 )
            arr = arr.sum( axis=-1) 

    # ---  Reactive bromine ( BrOx )
    if fam == 'ClOx' :
        # Select species in family
        specs = ['Cl', 'ClO', 'Cl2O2', 'ClOO', ]

        # Extract data
        arr = get_GC_output( wd=wd, vars=['IJ_AVG_S__'+i for i in specs ], \
                    trop_limit=trop_limit, r_list=True  )
        # Adjust to stoichiometry
        arr = [ arr[n]*spec_stoich(i, ref_spec='Cl') \
             for n,i in enumerate( specs ) ]
                
        if debug:
            print [ ( i.shape, i.min(), i.max(), i.mean() ) for i in arr  ]
        if not rtn_list:
            arr = np.ma.concatenate( [ i[...,None] for i in arr ], axis=-1 )
            arr = arr.sum( axis=-1) 


    if debug and (not rtn_list):
        print [ ( i.shape, i.min(), i.max(), i.mean() ) for i in [arr ] ]    
    # --- Mask for troposphere if t_ps provided (& trop_limit=True)
    if not isinstance( t_ps, type(None) ) and trop_limit:
        if rtn_list:
            arr = mask4troposphere( arr, t_ps=t_ps,  \
                use_time_in_trop=use_time_in_trop, \
                multiply_method=multiply_method, )
        else:
            arr = mask4troposphere( [arr], t_ps=t_ps, \
                use_time_in_trop=use_time_in_trop, \
                multiply_method=multiply_method  )[0]
    if debug and (not rtn_list):
        print [ ( i.shape, i.min(), i.max(), i.mean() ) for i in [arr ] ]    

    # Take annual mean?
    if annual_mean:
        if rtn_list:
            arr = [ i.mean( axis=-1 ) for i in arr ]
        else:
            arr = arr.mean( axis=-1 )

    # Return data and (optionally) specs
    if rtn_specs:
        return arr, specs
    else:
        return arr

# --------------
# 2.42 - Convert v/v to molec/cm3
# -------------   
def convert_v_v_2_molec_cm3( arr=None, wd=None, vol=None, a_m=None, \
            mols=None, res='4x5', trop_limit=True, debug=False):
    """ 
    Covnerts mixing ratio (v/v) into number density (molec/cm3).
    
    ARGUMENTS:
     - arr (array): arrray input
     - a_m (array): array of air mass
     - mols (array): array of molecule number density
     - trop_limit: limit output to "chemical troposphere" (level 38 )
     - res: resolution of the model input (e.g. 4x5, 2x2.5 )    

    NOTE:
     - required variables of volume (vol) and airmass (a_m) can be provided as 
    arguements or are extracted online (from provided wd )
    """

    if debug:
        print 'convert_v_v_2_molec_cm3 called for: ',
        print [ ( isinstance(i, np.ndarray), i.shape ) for i in arr, vol, a_m ]

    # Get volume ( cm^3  )
    if not isinstance(vol, np.ndarray):
        vol = get_volume_np( wd=wd, res=res, trop_limit=trop_limit, debug=debug)
        if debug:
            print 'WARNING: extracting volume online'
    # Get air mass ( kg )
    if not isinstance(a_m, np.ndarray):
        a_m = get_GC_output( wd, vars=['BXHGHT_S__AD'], trop_limit=trop_limit, 
                    dtype=np.float64)    
        if debug:
            print 'WARNING: extracting air mass online'

    # Get moles
    if not isinstance(mols, np.ndarray):
        mols = a_m*1E3/constants( 'RMM_air')

    #  Convert to molecules
    arr = ( arr * mols )*constants('AVG')

    #  convert to per unit area ( molecs / cm^3  )
    arr = arr / vol
    
    return arr

# --------------
# 2.43 - Mask non tropospheric boxes of 4D array
# -------------   
def mask4troposphere( ars=[], wd=None, t_ps=None, trop_limit=False, \
        t_lvl=None, masks4stratosphere=False, use_time_in_trop=True, \
        multiply_method=True, res='4x5', debug=False ):
    """ Mask for the troposphere using either the time in troposphere 
    diagnostic ( use_time_in_trop=True ) or troposphere level 
    ( use_time_in_trop=False ) 

    ARGUMENTS:
     - ars (list): arrays to apply masks to
     - t_ps: time in the troposphere diganostic ( float values 0 to 1 )
     - t_lvl: the model grid box level of the troposphere
     - trop_limit: limit output to "chemical troposphere" (level 38 )
     ( (boolean) options for this include using 
     - multiply_method 
     - definition of troposphere to use? 
        - use_time_in_trop (boolean): time a given box is in the troposphere
        ( if use_time_in_trop=False, the level of the troposphere is used )
     - conbine_ars (boolean): return arrays as a single array? 
     - month_eq (boolean): convert units to monthly equiivlents.
     - res: resolution of the model input (e.g. 4x5, 2x2.5 )
    
    NOTES:
     - The default option is not a mask. The time in troposphere (a fraction 
        value 0-1) is simply multipled by the provided array.
     - The tropospause is diagnosed mulitple ways. This fucntion can mask 
        for the troposphere using either time in trop ( use_time_in_trop=True ) 
        or tropopause level. Each array is done on a timestep by timestep basis. 
     - The 'chemical troposphere' in GEOS-Chem is the boxes for which 
        differential chemical equations are solved. Only these boxes (1st 38)
        are consdidered if trop_limit=True ( e.g. array shape is (72,46,38,12)
        instead of (72,46,47,12)   
    """
    logging.info( 'mask4troposphere called for arr of shape: {},'.format( \
            ars[0].shape)  )
    logging.debug('mask4troposphere - with multiply method?=', multiply_method, \
        ',use_time_in_trop=', use_time_in_trop, 'type(t_lvl): ', \
        type(t_lvl), 'type(t_ps): ', type(t_ps) )

    # --- Get time tropopause diagnostic (if not given as argument)
    if not isinstance(t_ps, np.ndarray) and use_time_in_trop: 
        t_ps = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], \
            trop_limit=trop_limit ) 
        if masks4stratosphere:
            # Extend to all full atmosphere ( 47 levels )
            a = list( get_dims4res(res) )
            a[-1] = 47-38
            a = np.zeros(  tuple( a+[t_ps.shape[-1]] ) )
            t_ps = np.ma.concatenate( (t_ps,a),  axis=-2 )

    # Get tropopause level (if not given as argument)
    if not isinstance(t_lvl, np.ndarray) and ( not use_time_in_trop): 
        t_lvl = get_GC_output( wd, vars=['TR_PAUSE__TP_LEVEL'], \
            trop_limit=False )

    # ---  Multiply by fractional time in trop. array    
    if multiply_method and use_time_in_trop:

        # Invert values if masking troposphere
        if masks4stratosphere:
            t_ps = 1 - t_ps
        # If 3D array with is given without a time dimension, average t_ps
        if len(ars[0].shape) == 3:
            t_ps = t_ps.mean(axis=-1)
        # Multiply fractional array by provided array            
        ars = [i*t_ps for i in ars ]

    # ---  Mask by fractional time in trop. array    
    elif (not multiply_method) and use_time_in_trop:
        # Mask area that are not exclusively stratospheric
        if masks4stratosphere:
            # mask troposphere
            t_ps = np.ma.masked_where( t_ps != 0, t_ps )

        # Mask area that are not exclusively tropospheric
        else:    
            t_ps = np.ma.masked_where( t_ps != 1, t_ps )

    # ---  Mask using tropopause level diagnostic values        
    else:
        # Setup dummy array with model numbers as values
        t_ps = np.zeros( ars[0].shape  )
        logging.debug( [ i.shape for i in t_ps, t_lvl, ars[0] ] )
        for i, n in enumerate( range(1,ars[0].shape[-2]+1) ):
            t_ps[:,:,i,:] =n
       
        if masks4stratosphere:
            # mask where the levels are greater that the tropopause level
            t_ps = np.ma.masked_where( t_ps<t_lvl[:,:,None,:], t_ps)
        else:
            # mask where the levels are greater that the tropopause level
            t_ps = np.ma.masked_where( t_ps>t_lvl[:,:,None,:], t_ps)
        
    # --- Set array mask to have strat mask (trop if masks4stratosphere=True)
    if (not multiply_method):
        for n, arr in enumerate( ars ):
            logging.debug( 'Using multiply_method={}, use_time_in_trop={}'.format(
                multiply_method, use_time_in_trop ) )
            try:
                ars[n] = np.ma.array( arr, mask=t_ps.mask )
            except:
                if len(arr.shape) == 3:
                    logging.debug( 'using tps averaged over time' )
                    # Apply mask
                    ars[n] = np.ma.array( arr, mask=t_ps.mean(axis=-1).mask )
                else:
                    # Log error
                    log_str = 'Using multiply_method={}, use_time_in_trop={}'.format(
                        multiply_method, use_time_in_trop )
                    logging.debug( log_str )
                    log_str = 'mask not applied for shapes', [i.shape  for i in arr, t_ps  ] 
                    logging.debug( log_str )
                    sys.exit()

    return ars 

# --------------
# 2.44 - Convert [molec/cm3/s ] to [molec/yr] 
# -------------   
def convert_molec_cm3_s2_molec_per_yr( ars=None, vol=None ):
    """ 
    Takes a list of arrays (4D) and converts their units from molec/cm3/s to     
    molec./yr
    
    ARGUMENTS:
     - vol: volumne of grid boxes
    """
    # Loop provided arrays
    for n, arr in enumerate( ars ):
        # Times by volume
        ars[n]  = arr *vol
        
        # Convert /s to /yr
        ars[n] = arr *60*60*24*365
    
    return ars


# --------------
# 2.46 - Convert molec/cm3/s to g/s
# -------------
def convert_molec_cm3_s_2_g_X_s( ars=None, specs=None, ref_spec=None, \
         months=None, years=None, vol=None, t_ps=None, trop_limit=True, \
        s_area=None, rm_strat=True, ctm_f=None, wd=None, res='4x5', \
        multiply_method=True, use_time_in_trop=True, conbine_ars=True, \
        month_eq=False, verbose=False,  debug=False ):
    """ 
    Convert molec/cm3/s to g/grid box. This is used for converting prod/loss 
    output units.

    ARGUMENTS:
     - s_area: Surface array (array of values in metres)
     - vol: volumne of grid boxes
     - specs: list of species (Prod loss variaables from input.geos)
     - integer values for months ("months") and years ("years")
     - t_ps: time in the troposphere diganostic ( float values 0 to 1 )
     - trop_limit: limit output to "chemical troposphere" (level 38 )
     - rm_strat (boolean): mask out stratosphere
     ( (boolean) options for this include using multiply_method and use_time_in_trop )
     - conbine_ars (boolean): return arrays as a single array? 
     - month_eq (boolean): convert units to monthly equiivlents.
     
    NOTES:
     - re-write of molec_cm3_s_2_Gg_Ox_np for clarity/split functionaltity
     - units of g/month can also be returned if month_eq=True
     - All functions that use "get_pl_in_Gg" should be updated to use this
     - It is most efficency to provide shared variables as arguements if this 
        function is call more that once by a single driver
    """
    logging.info( 'convert_molec_cm3_s_2_g_X_s called' )
    # --- Extract core model variables not provide
    if not isinstance(months, list):
        months = get_gc_months( ctm_f, wd=wd )
    if not isinstance(years, list):
        years  = get_gc_years( ctm_f=ctm_f, wd=wd, set_=False )
    if isinstance( s_area, np.ndarray):
        s_area = get_surface_area( res=res, debug=debug )
    if not isinstance(vol, np.ndarray):
        vol = get_volume_np( ctm_f=ctm_f, s_area=s_area, wd=wd, res=res,\
             debug=debug )
        logging.info( 'WARNING: extracting volume online - inefficent' )
    if debug:
        logging.debug( [ (i.sum(), i.shape) for i in ars ] )

    # --- loop spec ars
    for n, arr in enumerate( ars ):
        # convert from molec/cm3/s to  molec/s
        arr  = arr *vol[...,:38,:]
        # conver to to molec/s = > Gg/s
        arr =  arr / constants( 'AVG') * species_mass(ref_spec)
        # to / yr
        if month_eq:
            day_adjust = d_adjust( months, years)
            ars[n] = arr * day_adjust

    logging.debug( [ (i.sum(), i.shape) for i in ars ] )

    # only consider troposphere ( update this to use mask4troposphere )
    if rm_strat:
        ars = mask4troposphere( ars,  t_ps=t_ps, wd=wd,\
            use_time_in_trop=use_time_in_trop, multiply_method=multiply_method )

    if debug:
        logging.debug( [ (i.sum(), i.shape) for i in ars ] )
        logging.debug( [ (i.sum(), i.shape) for i in [ i[...,None] for i in ars ] ] )
        logging.debug( np.concatenate( [ i[...,None] for i in ars ], axis=-1 ).shape )
    
    if conbine_ars:
        return np.concatenate( [ i[...,None] for i in ars ], axis=-1 ) 
    else:
        return ars   
         
# --------------
# 2.47 - Print out stats on a list of 2D arrays in terms of masked areas
# -------------
def prt_2D_vals_by_region( specs=None, res='4x5', arrs=None, prt_pcent=False, \
        add_total=False, months=range(12), summate=True, \
        csv_title='regional_vals.csv', save2csv=False, \
        debug=False ):
    """ 
    Print values of a 2D (lon, lat) arry masked for regions 
    """

    # Which regions?
    m_titles = [ 'Tropics', 'Mid lats', 'Extratropics', 'Oceanic', 'NH', 'SH' ]

    # Get maskes
    masks = [ mask_all_but( i, mask2D=True, trop_limit=True, res=res)  \
        for i in m_titles ]
#    o_mask = masks[ m_titles.index( 'Oceanic' ) ]

    if debug:
        print [( m_titles[n] , i.shape ) for n, i in enumerate( masks ) ]

    # --- Average or total ?
    if add_total:
        arrs += [ np.ma.concatenate( [ i[...,None] for i in arrs ], \
            axis=-1 ).sum(axis=-1) ]
        specs += ['Total']

    # --- Print out actual values 
    pstr  = '{:<25}'+'{:<15}'*( len(m_titles)-1 )
    pstrn = '{:<25}' +'{:<15,.3f}'*( len(m_titles)-1 )
    arrsn = [ 'Species', 'Total'] + m_titles
    print m_titles, arrsn
    print pstr.format( *arrsn )    
    for n, s in enumerate( specs ):

        if summate:
            vars = [ s, np.ma.sum(arrs[n]) ]
            vars += [ np.ma.sum(arrs[n]*m) for m in masks ]        
        else:
            vars = [ s, np.ma.mean(arrs[n]) ]
            vars += [ np.ma.mean(arrs[n]*m) for m in masks ]        
        print pstrn.format( *vars )

    # --- Print out percent values 
    if prt_pcent :
        print [i.shape for i in arrs ]
        if len(arrs[0].shape) == 4:
            s_arrs = [ (i/len(months)*12).sum(axis=2) for i in arrs ]
        else:
            s_arrs = arrs
        # update titles 
        arrsn = [ 'Species', 'Run Total / Tg ','Yr. Equiv. / Tg'] + \
            [ '% '+i for i in m_titles ]
        # update numerical string to print
        pstr  = '{:<25}'+'{:<15}'*( len(arrsn)-1 )
        pstrn = '{:<25}' +'{:<15,.3f}'*( len(arrsn)-1 )
        print pstr.format( *arrsn )    
        # loop maskes and print
        vars_l = []
        for n, s in enumerate( specs ):
            vars = [ s, np.ma.sum(arrs[n]), np.ma.sum(s_arrs[n]) ]
            vars += [ np.ma.sum(s_arrs[n]*m)/np.ma.sum(s_arrs[n])*100 \
                for m in masks ]
            vars_l += [ vars ]
            print pstrn.format( *vars )

        # --- Convert to DataFrame, then save to csv
        if save2csv:
            # construct df
            d = dict( zip( specs, [i[1:] for i in vars_l ] ) )
            df = pd.DataFrame( d ).T
            # assign titles to columns
            df.columns = arrsn[1:]
            # save to csv
            df.to_csv( csv_title )




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
# (3) These will be removed when AC_tools is next re-structured.


# --------------
# 1.01 - open ctm.bpch using pygchem ( version '0.2.0' ) 
# -------------- 
def open_ctm_bpch(wd, fn='ctm.bpch', debug=False):
    """ 
    This is a vestigial programme, based on pychem version 0.2.0.
    Updates have made this incompatibile. 
        
    Update functions to use iris class. e.g. 
        
    if :
    wd = <run directory of containing ctm.bpch files/ctm.nc file to analyse>
    then:
    # setup variables 
    list_of_species =  ['NO', 'O3', 'PAN', 'CO', 'ALK4']

    # convert to Iris Cube/NetCDF naming form 
    list_of_species = ['IJ_AVG_S__'+i for i in list_of_species ]

    # this will return data as a 5D array ( species, lon, lat, alt, time)
    data = AC.get_GC_output( wd, vars=list_of_species )
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
# 1.02 - get np array (4D) of ctm.bpch (lon, lat ,alt, time) 
# --------------
def get_gc_data_np(ctm, spec='O3', category="IJ-AVG-$", debug=False):
    """ 
    This function extracts fields from diagnostics 
    
    NOTES:
     -  This was written to work with pychem version 0.2.0 and it is included 
    for back compatibility
    """

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


# ------------------ Section 6 -------------------------------------------
# -------------- Time Processing
#

# --------------
# 3.01 - Takes monthly outputs and sort to chronological order
# -------------
def ctms2chronological( ctms, debug=False ):
    """ 
    Ensure list of ctms is chronological 
    ARGUMENTS:
     - ctms is a list of ctm files 
    (from pyghcem ver <3.0 approach of passing bpch objects)   
    """
    
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
    """ 
    Ensure np array is in chronological order 

    ARGUMENTS:
     - ctms is a list of ctm files 
    (from pyghcem ver <3.0 approach of passing bpch objects)     
     - arr is a numpy array with the final dimension being time
    """

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

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# -------------- Delete/move to user specific modules 


# --------------
# 2.09 - Get tropospheric Burden - 
# -------------
# Redundent - mv'd to bottom of this module
def get_trop_burden( ctm=None, spec='O3', wd=None, a_m=None, t_p=None, \
        Iodine=False, all_data=True, total_atmos=False , res='4x5',  \
        trop_limit=True, debug=False ):
    """ REDUNDENT - Get Trosposheric burden. If requested (total_atmos=False)
        limit arrays to "chemical troposphere"  (level 38) and boxes that are in 
        the troposphere removed by multiplication of "time in troposphere"
        diagnostic
    NOTE(s):
     - use get_gc_burden instead 
    """
    logging.info('get_trop_burden called')

    # Get variables online if not provided
    if not isinstance(a_m, np.ndarray):
        a_m = get_air_mass_np( ctm_f=ctm, wd=wd, \
            trop_limit=trop_limit, debug=debug )
    if not isinstance(t_p, np.ndarray):
        # Retain PyGChem 0.2.0 approach for back compatibility
        if pygchem.__version__ == '0.2.0':
            t_p = get_gc_data_np( ctm=ctm, spec='TIMETROP', \
                category='TIME-TPS', debug=debug)
        else:
            t_p = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], \
                            trop_limit=trop_limit  ) 

    # Retain PyGChem 0.2.0 approach for back compatibility
    if pygchem.__version__ == '0.2.0':            
        ar = get_gc_data_np( ctm, spec, debug=debug )[:,:,:38,:]
    else:
        ar = get_GC_output( wd, vars=['IJ_AVG_S__'+ spec], \
                            trop_limit=trop_limit  ) 
    if debug:
        print [i.shape for i in ar, t_p, a_m ]

    # v/v * (mass total of air (kg)/ 1E3 (converted kg to g))  = moles of tracer
    ar = ar* ( a_m*1E3 / constants( 'RMM_air')) 
    if (Iodine):
        # convert moles to mass (* RMM) , then to Gg 
        ar = ar * float( species_mass('I') ) * spec_stoich(spec) /1E9 
    else:
        # convert moles to mass (* RMM) , then to Gg 
        ar = ar * float( species_mass(spec) ) /1E9 

    # Cut off at the "chemical troposphere" ( defined by GC integrator as 38th)
    if (not total_atmos):        
        ar = ar * t_p
    else:
        logging.info( 'get_trop_burden returning whole atmosphere (not troposphere)' )

    if debug:
        print 'Got buden for {} from {}'.format( spec, ctm )
    if (all_data):
        return ar
    else:
        return np.mean( ar, axis=3 )

# --------------
# 2.25 - Get O3 Burden
# -------------
# REDUNDENT - mv'd to bottom of module        
def get_O3_burden(wd=None, spec='O3', a_m=None, t_p=None, O3_arr=None, \
            ctm_f=False, trop_limit=True, all_data=False, annual_mean=True, \
            debug=False ):
    """ Get O3 burden in CTM output 

    NOTES:
     - all_data and annual_mean give the same result
     -  REDUNDENT? - this functionality is replicated in the get_burden 
     function, but it is still in use for automated output stat procsssing... 
     therefore retain for back compatibility 
    """
    # extract arrays not provided from wd
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
    if all_data or (not annual_mean):
        return ar
    else:
        return ar.mean(axis=3 )



# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# ---------------- Section X -------------------------------------------
# -------------- Non-generic Functions
# --------------------------------------------------------------------------
# 
# NOTE(s): 
# (1) These are retained, but will be migrated to a seperate non-generic module
# (2) It is not advised to use these. 


# Non generic function to move to new module

