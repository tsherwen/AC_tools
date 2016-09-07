#!/usr/bin/python

# =================================================
# --------- tms - module of Variables for re-use----------------
# -------------- 

# Section 0 - Required modules
# Section 1 - Planeflight variables
# Section 3 - GeosChem (bpch) prod loss variables
# Section 4 - GeosChem (bpch) general variables
# Section 5 - Misc
# Section 6 - Dynamic p/l processing
# Section 7 - Obervational variables
# Section 8- Drivers functions

# ---------------------------- ------------- ------------- ------------- 
# --------------  Contents
# --------------- ------------- ------------- -------------
# ---- Section 0 ----- Required modules
# --------------- ------------- ------------- -------------

# ---- Section 1 ----- Planeflight variables
# 1.01 - PF variable dictionary ***
# 1.02 - TRA_?? to Geos-Chem species name ***

# --------------- ------------- ------------- -------------
# ---- Section 2 ----- Drivers functions
# 2.01 - P/L tag to PD tag ***
# 2.02 - Get P/L dictionary for a given species ***
# 2.03 - Get prod loss (P/L) reactions for a given family.

# --------------- ------------- ------------- -------------
# ---- Section 3 ----- GeosChem (bpch) prod loss variables
# 3.01 - Spec to photolysis reaction p/l tag 
# 3.02 - Ox family for tag

# --------------- ------------- ------------- -------------
# ---- Section 4 ----- GeosChem (bpch) general variables
# 4.01 - v9-2 species in input.geos from num
# 4.02 - Get Species Mass 
# 4.03 - Get Species stoichiometry 
# 4.04 - GEOS-Chem/ctm.bpch values (current main dict ) ***
# 4.05 - latex species name
# 4.06 - converts P/L tracer mulitpler to 1
# 4.07 - Returns tracers unit and scale (if requested)
# 4.08 - Store of dirs for earth0, atmosviz1, and tms mac
# 4.09 - Ox in species (redundant now? should adapt species stoich )
# 4.10 - Get Gaw site name from GAW ID
# 4.11 - returns dict of gaw sites
# 4.12 - Return lat, lon, alt for a given resolution
# 4.13 - Get model array dimension for a given resolution
# 4.14 - Convert gamap category/species name to Iris/bpch name
# 4.15 - Get scaling for a given unit
# 4.16 - Species Class 
# 4.17 - Observation Site Class 
# 4.18 -  dictionary of category + species names from diagnostic ordering
# 4.99 - Reference data, (inc. grid data) from gchem 

# --------------- ------------- ------------- -------------
# ---- Section 5 ----- Misc
# 5.01 - dir store (standard directories on different servers )
# 5.02 - Store of  constants for use by funcs/progs

# --------------- ------------- ------------- -------------
# ---- Section 6 ----- Dynamic prod/loss dictionary processing ( For GEOS-Chem)
# 6.01 - Make rxn dict of all active reactions ***
# 6.02 - Return all reactions for a given p/l family ***
# 6.03 - Return reaction infomaiton for given tags (e.g. PD... )
# 6.04 - Create an indices list to split reaction by family (e.g. for Ox loss)
# 6.05 - Return tags for a given reaction
# 6.06 - Extract all p/l speacies in a given input.geos
# 6.07 - Extract all active tags from a given smv.log
# 6.08 - Extract all active PDs from a given smv.log
# 6.09 - get all active reaction for a given tag
# 6.10 - get all details for a given tag 
# 6.11 - get reaction coeffifecent
# 6.12 - Remove ClBrI het loss tracers during testing
# 6.13 - PD to rxn str - remake of redundant function
# 6.14 - Get OH reactants from reaction number dictionary


# --------------- ------------- ------------- -------------
# ---- Section 7 ----- Observational variables
# 7.01 - IO observation dictionary
# 7.02 - BAE flight ID dictionary
# 7.03 - CAST flight dictionary for CIMS/CIMSII
# 7.04 - Iodocarbon obs. meta data
# 7.05 - Stores locations for use by funcs/progs - LON, LAT, ALT - double up? ( with 5.02 ?? )
# 7.06 - Get Locations of observations (lats, lons, alts ) for given sites
# 7.07 - sonde station variables (list of 432 sondes)
# 7.08 - returns  (lat, lon, alt (press), timezone (UTC) ) for a given site


# --------------- ------------- ------------- -------------
# ---- Section XX ----- Redundant functions
# 4.05 - Convert species to formatted "LaTeX" form



# ------------------ Section 0 -----------------------------------
# -------------- Required modules:
#
#
# -- I/O / Low level                                                                                
import re
import os 
#import platform
import pandas as pd
from netCDF4 import Dataset
#import Scientific.IO.NetCDF as S
import sys
import glob

# - Math/Analysis                                                                                   
import numpy as np

# - tms
from AC_tools.funcs4core import *

# ------------------------------------------- Section 1 -------------------------------------------
# -------------- Planeflight variables
#

# --------------
# 1.01 - dictionary of variables used for planeflight_mod.F output
# -------------
def pf_var( input, ver='1.7', ntracers=85, JREAs=[] ):
    """ Dictionary store for planeflight names/tracers
    NOTES:
      - UPDATED NEEDED: MORE DETAILED DESCRIPT.   
      - Why is this function not in funcs4pf?
        
    """
    # planeflight variable lists
    metvars = [
    'GMAO_TEMP', 'GMAO_ABSH', 'GMAO_SURF', 'GMAO_PSFC', 'GMAO_UWND', 'GMAO_VWND'
    ]
#    species  = [
#    'O3', 'NO2', 'NO', 'NO3', 'N2O5', 'HNO4', 'HNO3', 'HNO2', 'PAN', 'PPN', 
#    'PMN', 'R4N2', 'H2O2', 'MP', 'CH2O', 'HO2', 'OH', 'RO2', 'MO2', 'ETO2', 
#    'CO', 'C2H6', 'C3H8', 'PRPE', 'ALK4', 'ACET', 'ALD2', 'MEK', 'RCHO', 'MVK', 
#    'SO2', 'DMS', 'MSA', 'SO4', 'ISOP'
#    ]
    species = [
    'A3O2', 'ATO2', 'B3O2', 'EOH', 'ETO2', 'ETP', 'GLYX', 'HO2', 'IAP', 'INO2', 'INPN', 'ISN1', 'ISNOOA', 'ISNOOB', 'ISNOHOO', 'ISNP', 'KO2', 'MAN2', 'MAO3', 'MAOP', 'MAOPO2', 'MCO3', 'MGLY', 'MO2', 'MRO2', 'MRP', 'OH', 'PO2', 'PP', 'PRN1', 'PRPN', 'R4N1', 'R4O2', 'R4P', 'RA3P', 'RB3P', 'RCO3', 'RIO2', 'ROH', 'RP', 'VRO2', 'VRP', 'LISOPOH', 'ISOPND', 'ISOPNB', 'HC5', 'DIBOO', 'HC5OO', 'DHMOB', 'MOBAOO', 'ISOPNBO2', 'ISOPNDO2', 'ETHLN', 'MACRN', 'MVKN', 'PYAC', 'IEPOXOO', 'ATOOH', 'PMNN', 'MACRNO2', 'PMNO2'
    ] 
    OH_reactivity=[
    'NO', 'ISOP', 'GMAO_TEMP', 'GMAO_PSFC', 'CO', 'ACET', 'ALD2', 'MEK', 'MVK', 'MACR', 'C3H8', 'CH2O', 'C2H6', 'SO2', 'NO2', 'ISOPNB', 'ISOPND', 'NO3', 'HNO2', 'HNO3', 'OH', 'HO2', 'H2O2', 'MP', 'ATOOH', 'HNO4', 'ALK4', 'ISN1', 'R4N2', 'RCHO', 'ROH', 'PRPE', 'PMN', 'GLYC', 'GLYX', 'MGLY', 'HAC', 'INPN', 'PRPN', 'ETP', 'RA3P', 'RB3P', 'R4P', 'RP', 'PP', 'RIP', 'IEPOX', 'IAP', 'VRP', 'MRP', 'MAOP', 'MAP', 'DMS', 'HBr', 'Br2', 'BrO', 'CHBr3', 'CH2Br2', 'CH3Br', 'HC5', 'ISOPND', 'ISOPNB', 'ISNP', 'MVKN', 'MACRN', 'DHMOB', 'MOBA', 'ETHLN', 'PROPNN'
    ] 
    OH_Extras4nic = [
    'OH', 'MCO3', 'A3O2', 'PO2', 'R4O2', 'R4O2', 'R4N1', 'ATO2', 'KO2', 'RIO2', 'VRO2', 'MRO2', 'MAN2', 'B3O2', 'INO2', 'ISNOOA', 'ISNOOB', 'ISNOHOO', 'PRN1', 'RCO3', 'MAO3', 'IEPOXOO', 'MAOPO2', 'MAOPO2', 'HC5OO', 'HC5OO', 'ISOPNDO2', 'ISOPNDO2', 'ISOPNBO2', 'ISOPNBO2', 'DIBOO', 'DIBOO', 'MOBAOO', 'MOBAOO', 'H2', 'CH4', 'HCOOH', 'MOH', 'ACTA', 'EOH', 'VRP'
    ]
    OH_rxns_17_EOH = [ 21, 
    328, 2, 387, 6, 7, 136, 9, 10, 407, 14, 15, 16, 408, 402, 403, 406, 23, 152, 25, 27, 28, 30, 415, 33, 35, 164, 167, 40, 169, 43, 172, 173, 48, 49, 178, 179, 52, 201, 58, 60, 394, 62, 319, 320, 66, 323, 196, 197, 198, 199, 200, 177, 202, 203, 204, 333, 206, 210, 211, 212, 334, 214, 215, 176, 345, 346, 347, 94, 357, 226, 101, 102, 103, 365, 113, 114, 153, 378, 379, 380, 364
    ]
    OH_rxns_17_EOH =  [ 'REA_{}'.format(i) for i in  OH_rxns_17_EOH ]

    # OH reaction species and tracers
    OH_rxn_tras = ['HBr', 'HOI', 'I2', 'CH2O', 'NO2', 'HNO3', 'CH3IT', 'NO', 'HNO3', 'HNO2', 'PRPE', 'PMN', 'HNO4', 'GLYC', 'NO3', 'C3H8', 'DMS', 'ALK4', 'SO2', 'MVK', 'ISOP', 'CHBr3', 'BrO', 'ISOP', 'CH2Br2', 'CH3Br', 'ACET', 'MEK', 'H2O2', 'CO', 'PROPNN', 'MP', 'MACR', 'HAC', 'ALD2', 'MOBA', 'RIP', 'Br2', 'IEPOX', 'MAP', 'R4N2', 'RCHO', 'HI']
    OH_rxn_tras = [ num2spec( i, ver=ver, invert=True) for i in OH_rxn_tras ]
    OH_rxn_tras = [ 'TRA_{:0>2}'.format( i)  for i in OH_rxn_tras ]
    OH_rxn_specs = ['GLYX', 'MGLY', 'HCOOH', 'MOH', 'O2', 'DHMOB', 'HO2', 'OH', 'H2', 'CH4', 'ETHLN', 'EOH', 'ATOOH', 'R4P', 'INPN', 'RA3P', 'RB3P', 'RP', 'PP', 'IAP', 'VRP', 'MRP', 'MAOP', 'ISN1', 'HC5', 'ACTA', 'ISOPNB', 'ROH', 'ISNP', 'MVKN', 'MACRN', 'ISOPND', 'ETP', 'PRPN']

    # remove inactive species 
    inactive_spec = ['ACTA', 'CH4', 'H2', 'HCOOH', 'MOH', 'O2'] # 'EOH',

    # Setup list of tracers
    if any( [(ver == i) for i in '1.6', '1.6.1', '1.6.2'] ):
        ntracers=96
    if any( [(ver == i) for i in '1.6.3', '1.6.4' ] ):
        ntracers=98
    if ver == '1.7':
        ntracers=85
    if any( [ (ver == i) for i in '2.0', '3.0' ] ):
        ntracers=103
    if ver == 'johan_br.v92':
        ntracers=87
    TRAs = ['TRA_'+ str(i) for i in range(1, ntracers+1) ] 
#    TRAs = ['TRA_{:0>2}'.format(i) for i in range(1, ntracers+1) ]

    # Setup list of reactions ( photolysis and general )
    if ver == '1.5':
        PHOT_1st, PHOT_last = 455, 533
    if ver == '1.6':
        PHOT_1st, PHOT_last = 453, 531
    if ver == '1.6.1':
        PHOT_1st, PHOT_last = 453, 530
    if ver == '1.6.2':
        PHOT_1st, PHOT_last = 452, 529
    if ver == '1.6.3':
        PHOT_1st, PHOT_last = 461, 538
    if ver == '1.7':    
        PHOT_1st, PHOT_last = 453, 529
    if ver == '2.0':    
#        PHOT_1st, PHOT_last = 413, 614 # dev
        PHOT_1st, PHOT_last = 519, 607 
    if ver == '3.0':    
#        PHOT_1st, PHOT_last = 537, 627 # coupled sim.
#        PHOT_1st, PHOT_last = 554, 644 # With IX split rxns
        PHOT_1st, PHOT_last = 548, 638 # With 0.25 IBr split rxns

    JREAs = ['REA_'+ str(i) for i in range(PHOT_1st, PHOT_last) ] 
    REAs_all = ['REA_'+ str(i) for i in range(0, 533) ] 
    
    # reduced list for high time and spatial resolution
    if any( [ input ==i for i in 'slist_v9_2_NREA_red', 
        'slist_v9_2_NREA_red_NOy'] ):
        TRAs = GC_var('active_I') + ['AERI']
        TRAs= [ num2spec( i, ver=ver, invert=True) for i in TRAs ]
        TRAs = [ 'TRA_{:0>2}'.format( i)  for i in TRAs ]

        metvars = [ i for i in metvars if not any( [ (i==ii) \
            for ii in   'GMAO_ABSH', 'GMAO_SURF', 'GMAO_PSFC' ] ) ]
        species = [ i for i in species if not any( [ (i==ii) for ii in  
        'R4N2', 'MP', 'CH2O', 'MO2', 'ETO2', 'CO', 'C2H6', 'C3H8', 'PRPE', 'ALK4', 'ACET', 'ALD2', 'MEK', 'RCHO', 'MVK', 'DMS', 'MSA', 'ISOP'  
        ]) ]

    if input =='slist_ClearFlo':
        TRAs = 'CO', 'ACET', 'ALD2', 'ISOP', 'C2H6', 'C3H8', 'CH2O', \
            'MACR', 'HNO2', 'HNO3', 'MVK', 'NO', 'NO2', 'PAN', 'O3',
        TRAs= [ num2spec( i, ver=ver, invert=True) for i in TRAs ]
        TRAs = [ 'TRA_{:0>2}'.format( i)  for i in TRAs ]
        # mannually add ethanol
        TRAs += [ 'TRA_86']  
        species = [ 'OH', 'MO2','HO2' ]

    if input =='slist_PEN_WEY_LEI':
        TRAs ='NO', 'NO2', 'PAN', 'HNO2', 'HNO3', 'O3',\
            'CO', 'ISOP', 'R4N2', 'PRPE', 'CH2O', \
            'NO3', 'N2O5', \
            'MP', 'DMS', 'SO2', \
            'BrCl', 'Cl2', 'Cl', 'ClO', 'HOCl', 'HCl', 'ClNO2', 'ICl', \
            'SO4', 'SO4s', 'MSA', 'NH3', 'NH4', 'NIT', 'NITs', \
             'SALA', 'SALC',\
#             'ACET', 'ALD2', 'C2H6', 'C3H8', 'MACR', 'MVK', 
        TRAs= [ num2spec( i, ver=ver, invert=True) for i in TRAs ]
        TRAs = [ 'TRA_{:0>2}'.format( i)  for i in TRAs ]
        # Species? ( aka ones that are not transported )
        species = [ 'OH', 'HO2' ]

        
    if input =='slist_v9_2_NREA_red_NOy':
# THIS IS NOT A GOOD APPROACH, use actual names an tranlate based on verison. 
#        missing = [  'TRA_17', 'TRA_60', 'TRA_30', 'TRA_31', 'TRA_50', \
#            'TRA_54', 'TRA_55', 'TRA_57' ]
        
#            use tracer  TRA_60, TRA_30, TRA_31, for: 'MMN' , 'NH3' , 'NH4',  'R4N2', 'BrNO2', 'BrNO3','MPN', 'PROPNN', 
        missing = 'MMN' , 'NH3' , 'NH4',  'R4N2', 'BrNO2', 'BrNO3','MPN', 'PROPNN'
        missing = [ num2spec( i, ver=ver, invert=True) for i in missing ]
        missing = [ 'TRA_{:0>2}'.format( i)  for i in missing ]
        species = species  + missing

    # use TRA_?? instead of species if species is a tracer
    species_=[]
    for s in species:
        try:
            species_ += ['TRA_{}'.format(  num2spec(s, ver=ver,\
                     invert=True) )]
        except:
            species_ = [ s ]
    species = species_

    # Construct dictionary
    d= {    
    'species' : species,
    'metvars' : metvars,
    'REAs_all' : REAs_all,
    'JREAs': JREAs,
    'TRAs' : TRAs,
    'slist' :  species +TRAs +JREAs+ metvars , 
    'slist_v9_2_NH' :   species + TRAs[:66] + metvars ,
    'slist_v9_2_NREA' :   species + TRAs + metvars ,
    'slist_v9_2_NREA_red': species + TRAs + metvars,
    'slist_REAs_all' :   species + TRAs + REAs_all + metvars,
    'slist_REAs_all_OH' :   species + TRAs  + metvars+OH_reactivity,
    'slist_REAs_all_OH_extras' :   species + TRAs  + metvars, 
    'slist_v9_2_NREA_red_NOy' : species + TRAs + metvars,
    'slist_v10_1.7_allspecs': species +TRAs+ JREAs +metvars,
    'slist_ClearFlo': species + TRAs + metvars, 
    'slist_ClearFlo_OH_rxn': species + TRAs + metvars + OH_rxns_17_EOH + OH_rxn_tras + OH_rxn_specs, 
    'slist_PEN_WEY_LEI': species + TRAs + metvars 
      } 

    # retrieve variable list from dictionary
    vars = d[input]
    
    # return unique list
    vars = sorted( list( set( vars) ) )

    # remove inactive tracers/species
    inactive_spec = [ n for n, i in enumerate( vars ) if (i in inactive_spec ) ]
    print inactive_spec
    [ vars.pop(i) for i in sorted( inactive_spec )[::-1] ]
    print vars

    return vars

# --------------
# 1.02 - Translator for planeflight species to GEOS-Chem species
# -------------
def what_species_am_i(input=None, V_9_2=True, V_9_2_C=False, \
            ver='1.7', special_case=None, invert=False, rtn_dict=False, \
            debug=False ) :
    """ Converts a GEOS-Chem (GC) species/tracer into a PF tracer (TRA_##).
            takes TRA_## & returns GC ID or other wayround
    NOTES:
     - Species have the same names in PF, but units are in molec/cm3, not 
       mixing ratio (v/v)
     -  Generic v10 and v11 need adding to this list

    """
    # select correct naming dictionary
    var ={  \
    '1.6': 'GCFP_d2TRA_all_1.6',
    '1.6.2': 'GCFP_d2TRA_all_1.6', # same as 1.6
    '1.6.3': 'GCFP_d2TRA_all_1.6.3', # 1.6 + 2
    '1.7': 'GCFP_d2TRA_all_1.7',
    '2.0': 'GCFP_d2TRA_all_2.0',    
    '3.0':  'GCFP_d2TRA_all_2.0' # Same as 2.0
    }[ver]

#    special_case = 'EOH'
#    special_case = 'EOH + actual names'    
#    if all_TRA:
#            var ={  \
#            '1.7': 'all_TRA_spec_met_1.7_EOH'
#            '1.6': 'GCFP_d2TRA_1.6'
#            }[ver]

    if not isinstance( special_case, type(None) ):
        var = {
#         'EOH':'GCFP_d2TRA_all_1.7_EOH', 
#        'EOH + actual names':'GCFP_d2TRA_all_1.7_EOH_actual_names'
        'all_TRA_spec_met_1.7_EOH.zeros' : 'TRA_spec_met_all_1.7_EOH',  #' 
       'all_TRA_spec_met_1.7_EOH':'TRA_spec_met_all_1.7_EOH_no_trailing_zeroes'
#        TRA_spec_met_all_1'
        }[special_case]

    # Get dictionary from variable store
    d =  GC_var( var ) 

    if debug:
        print d, special_case

    if invert:
        d = {v: k for k, v in d.items()}

    # return dictionary
    if rtn_dict:
        return d    
    else:
        return d[input]



# ------------------------------------------- Section 3 ----------------
# --------------  GeosChem (bpch) prod loss variables
#

# ----------------- Section 4 -------------------------------------------
# -------------- GeosChem (bpch) general variables
#

# --------------   
# 4.01 - v9-2 species in input.geos from num
# -------------    
def num2spec( num=69, rtn_dict=False, invert=False, ver = '1.7' ):
    """ Returns the tracer with a given tracer number in input.geos. Also works 
    in reverse or returns whole dictionary to remove need to remake dictionary 
    for each call.
    NOTES:
        (1) version number needed ( e.g. "Cl+Br+I" = 3.0, "1.7" = "Br+I",
        (2) UPDATE NEEDED: basecase v9-2/v10 not currently included
    """
    
    # --- Get dictionary of tracer numbers
    d = what_species_am_i( ver=ver, rtn_dict=True, special_case=None )

    # --- Slice off just numbers
    # special case for dev version?
    if any( [ (ver == i) for i in '1.6', '1.6.2', ] ): 
        d = GC_var('GCFP_d2TRA_justTRA_1.6' )
    if any( [ (ver == i) for i in'1.6.3', '1.6.4' ] ):
        d = GC_var('GCFP_d2TRA_justTRA_1.6.3' )
    # Then slice
    nums =[ int(i[4:]) for i in d.keys()]

    # --- Re-make dictionary
    d = dict( zip(nums, d.values() ) )
    
    # --- Invert to give spec for num
    if invert:
        d = { v: k for k, v in d.items() }

    if rtn_dict:
        return d
    else:
        return d[num]

# --------------
# 4.02 - RMM (Mass) (g /mol) for species 
# ------------- 
def species_mass( spec ):
    """ Function to return species mass ( in relative molecular mass ( RMM ) 
        for given species

    Note(s): 
        (1) C3H5I == C2H5I (this is a vestigle typo, left in to allow for 
        use of older model run data  )
    """
    d = {
    'HIO3': 176.0, 'OCPO': 12.0, 'Br2': 160.0, 'OCPI': 12.0, 'O3': 48.0, \
    'PAN': 121.0, 'ACET': 12.0, 'RIP': 118.0, 'BrNO3': 142.0, 'Br': 80.0, \
    'HBr': 81.0, 'HAC': 74.0, 'ALD2': 12.0, 'HNO3': 63.0, 'HNO2': 47.0, \
    'C2H5I': 168.0, 'HNO4': 79.0, 'OIO': 159.0, 'MAP': 76.0, 'PRPE': 12.0, \
    'CH2I2': 268.0, 'IONO2': 189.0, 'NIT': 62.0, 'CH3Br': 95.0, \
    'C3H7I': 170.0, 'C3H8': 12.0, 'DMS': 62.0, 'CH2O': 30.0, 'CH3IT': 142.0, \
    'NO2': 46.0, 'NO3': 62.0, 'N2O5': 105.0, 'H2O2': 34.0, 'DST4': 29.0, \
    'DST3': 29.0, 'DST2': 29.0, 'DST1': 29.0, 'MMN': 149.0, 'HOCl': 52.0, \
    'NITs': 62.0, 'RCHO': 58.0, 'C2H6': 12.0, 'MPN': 93.0, 'INO': 157.0, \
    'MP': 48.0, 'CH2Br2': 174.0, 'SALC': 31.4, 'NH3': 17.0, 'CH2ICl': 167.0, \
    'IEPOX': 118.0, 'ClO': 51.0, 'NO': 30.0, 'SALA': 31.4, 'MOBA': 114.0, \
    'R4N2': 119.0, 'BrCl': 115.0, 'OClO': 67.0, 'PMN': 147.0, 'CO': 28.0, \
    'BCPI': 12.0, 'ISOP': 12.0, 'BCPO': 12.0, 'MVK': 70.0, 'BrNO2': 126.0, \
    'IONO': 173.0, 'Cl2': 71.0, 'HOBr': 97.0, 'PROPNN': 109.0, 'Cl': 35.0, \
    'I2O2': 286.0, 'I2O3': 302.0, 'I2O4': 318.0, 'I2O5': 334.0, 'MEK': 12.0, \
    'HI': 128.0, 'ISOPN': 147.0, 'SO4s': 96.0, 'I2O': 270.0, 'ALK4': 12.0, \
    'MSA': 96.0, 'I2': 254.0, 'PPN': 135.0, 'IBr': 207.0, 'MACR': 70.0, \
    'I': 127.0, 'AERI': 127.0, 'HOI': 144.0, 'BrO': 96.0, 'NH4': 18.0, \
    'SO2': 64.0, 'SO4': 96.0, 'IO': 143.0, 'CHBr3': 253.0, 'CH2IBr': 221.0, \
    'ICl': 162.0, 'GLYC': 60.0, \
    # species, not in GEOS-Chem tracer list
    'HO2': 33.0, 'OH': 17.0,'CH4':16.0 , 'N':14.0, 'CH3I':142.0, \
    'CH2OO':46.0, 'S': 32.0, \
    # Additional 2.0 species 
    'HCl': 36.5, 'HOCl': 52.5, 'ClNO2': 81.5, 'ClNO3': 97.5 , 'ClOO': 67.5, \
    'Cl2O2': 103.0,  'CH3Cl':  50.5, 'CH2Cl2': 85.0, 'CHCl3': 119.5, \
    'BrSALA': 80., 'BrSALC': 80., 'ISALA': 127. ,  'ISALC': 127. , \
    # Additional "species" to allow for ease of  processing
    'AERI_AVG': ( (286.0+302.0+318.0)/3 )/2, 'SO4S': 96.0, 
    'IO3': 127.+(3.*16.) , 'SSBr2': 160.0, 'C': 12.0, 
    # Add families for ease of processing
    'Iodine': 127.0, 'Iy': 127., 'Bromine': 80.0, 'Bry': 80.0,'Chlorine':35.0,
    'Cly':35.0, 'NOy': 14.0, 'NOx': 14.0, 'SOx': 32.0,\
     'Sulfate': 32.0, 'sulfur': 32.0, 'VOCs':12.0, 
    }
    
    return d[spec]

# --------------
# 4.03 -  return the stoichiometry of Iodine in species
# --------------
def spec_stoich( spec, IO=False, I=False, NO=False, OH=False, N=False,
            C=False, Br=False, Cl=False, S=False, ref_spec=None, debug=False ): 
    """ Returns unit equivelent of X ( e.g. I ) for a give species. 
        
        This can be automatically set by providing a reference species
        
    Notes:
     - Update Needed: re-write to take stioch species (e.g. OH, I instead of booleans )
     - asssume I == True as default
     - C3H5I == C2H5I 
        (this is a vestigle typo, left in to allow for use of older model runs )
     - aerosol cycling specs
    # 'LO3_36' : (2.0/3.0) , 'LO3_37' : (2.0/4.0),           # aersol loss rxns... 'LO3_37' isn't true loss, as I2O4 is regen. temp
     - Aerosol loss rxns ( corrected stochio for Ox, adjsutment need for I )
    """
    # If reference species provided automatically select family
    if not isinstance( ref_spec, type(None) ):
        if any( [(ref_spec==i)  for i in 'I', 'Iy', 'Iodine' ] ):
            I=True
        if any( [(ref_spec==i) for i in 'Br', 'Bry', 'Bromine' ] ):
            Br=True
        if any( [(ref_spec==i)  for i in 'Cl', 'Cly', 'Chlorine' ] ):
            Cl=True
        if any( [(ref_spec==i) for i in  'C', 'VOC' ] ):
            C=True
        if any( [(ref_spec==i) for i in  'N', 'NOy', 'NOx' ] ):
            N=True
        if any( [(ref_spec==i) for i in  'OH', 'HO2' ] ):
            OH=True
        if ref_spec == 'IO':
            IO=True
        if ref_spec == 'NO':
            NO=True
        if any( [(ref_spec==i) for i in  'S', 'SOx', 'Sulfate' ] ):
            S=True

    if debug:
        vars = ref_spec, IO, I, NO, OH, N, C, Br, Cl
        varsn = 'ref_spec', 'IO', 'I', 'N', 'OH', 'N', 'C', 'Br', 'Cl'        
        print "'spec_stoich'  called for: ", zip( varsn, vars ) 

    # Select dictionary ( I=True is the default... )
    if IO:
        d = {
        'RD11': 2.0, 'RD10': 1.0, 'RD12': 2.0, 'LO3_36': 1./3., 'RD09': 1.0, \
        'RD66': 1.0, 'RD23': 1.0, 'RD37': 1.0, 'LO3_24': 1.0/2.0, 'RD56': 1.0, \
        'RD01': 1.0, 'RD08': 1.0, 'RD46': 2.0, 'RD30': 1.0, 'RD25': 1.0, \
        'RD27': 1.0, 'RD97':1.0
        }
    elif NO:
        d ={
        'NO2': 1.0, 'NO3': 1.0, 'N2O5': 2.0, 'NO': 1.0, 'PPN': 1.0, 'R4N2': 1.0, 
        'BrNO3': 1.0, 'INO': 1.0, 'PAN': 1.0, 'PMN': 1.0, 'HNO3': 1.0, \
        'HNO2': 1.0, 'NH3': 1.0, 'HNO4': 1.0, 'BrNO2': 1.0, \
        'IONO': 1.0, 'PROPNN': 1.0, 'NH4': 1.0, 'MPN': 1.0, 'MMN': 1.0, \
        'ISOPN': 1.0, 'IONO2': 1.0 \
        }
    elif OH:
        d = {
        'LO3_18': 2.0, 'LO3_03': 1.0,  'PO3_14': 1.0, 'RD65': 1.0,'LR25': 1.0, \
        'LOH':1.0, 'POH':1.0, 'LO3_86': 1.0, 'RD98':1.0, \
        # Redundent: 'RD95': 1.0,
        # also include HO2 and OH for HOx calculations
        'OH' :1.0, 'HO2': 1.0
        }
    elif S:
        d = {
        'S' :1.0, 'SO4': 1.0, 'SO4s': 1.0, 'SO2': 1.0
        }
    elif N:
        d= {
        'RD10': 1.0, 'LR26': 1.0, 'LR27': 1.0, 'LR20': 1.0, 'RD17': 1.0, \
        'RD16': 1.0, 'RD19': 1.0, 'RD18': 2.0, 'LR28': 1.0, 'LO3_30': 1.0, \
        'RD75': 1.0, 'LR7': 1.0, 'LR8': 1.0, 'RD56': 1.0, 'RD24': 1.0, \
        'LO3_39': 1.0, 'RD25': 1.0, 'RD81': 1.0, 'LR35': 1.0, 'LR18': 1.0, \
        'LR17': 1.0, 'LR11': 1.0, 'LR39': 1.0, 'RD20': 1.0, 'RD21': 2.0, \
        'RD22': 1.0, 'RD23': 1.0, 'RD68': 1.0, 'RD69': 1.0, \
    # NOy ( N in 'NOy')
        'NO2': 1.0, 'NO3': 1.0, 'N2O5': 2.0, 'NO': 1.0, 'PPN': 1.0, \
        'R4N2': 2.0, 'BrNO3': 1.0, 'INO': 1.0, 'PAN': 1.0, 'PMN': 1.0, \
        'HNO3': 1.0, 'HNO2': 1.0, 'NH3': 1.0, 'HNO4': 1.0, 'BrNO2': 1.0, \
        'IONO': 1.0, 'PROPNN': 1.0, 'NH4': 1.0, 'MPN': 1.0, 'MMN': 1.0, \
        'ISOPN': 1.0, 'IONO2': 1.0, 'ClNO2': 1.0, 'ClNO3':1.0, \
        }
    elif C:
        d = {
    'ACET': 3.0, 'ALD2': 2.0, 'C2H6': 2.0, 'C3H8': 3.0, 'ISOP': 5.0, \
    'PRPE': 3.0, 
        }
    elif Br:
         d= {
        'CH3Br': 1.0, 'HOBr': 1.0, 'BrO': 1.0, 'CHBr3': 3.0, 'Br2': 2.0, \
        'BrSALC': 1.0, 'CH2IBr': 1.0, 'BrCl': 1.0, 'Br': 1.0, 'CH2Br2': 2.0, \
        'IBr': 1.0, 'BrSALA': 1.0, 'BrNO2': 1.0, 'BrNO3': 1.0, 'HBr': 1.0, \
        # for ease of processing also include Seasalt Br2
         'SSBr2': 2.0, 
         # Also have reaction tracers
        'LR73' : 1.0, 
        # Note: stoichometry is for **GAS** phase Br (aka not SSA )
        # ( Aka JT03s == Br2 ( ==2 ), but one is BrSALA/BrSALC therefore =1)
        'JT03s' : 1.0, 'JT04s' :1.0, 'JT05s': 1.0
        }
    elif Cl:
        d= {
        'ClO': 1.0, 'Cl': 1.0, 'ClOO': 1.0, 'ClNO3': 1.0, 'ClNO2': 1.0, \
        'Cl2': 2.0, 'OClO': 1.0, 'HOCl': 1.0, 'HCl': 1.0, 'Cl2O2': 2.0,\
        'BrCl': 1.0, 'ICl':1.0, 'CH2Cl2':2.0, 'CHCl3': 3.0, 'CH2ICl': 1.0, \
        'CH3Cl': 1.0, 
         # Also have reaction tracers
         'LR62': 3.0, 'LR107': 3.0, 
         'LR74' : 1.0, 'LR106':1.0, 'LR103': 1.0, 
         'LR75' : 2.0, 'LR105': 2.0, 'LR104' : 2.0, 
        }
    elif Cl:
        d= {
        'ClO': 1.0, 'Cl': 1.0, 'ClOO': 1.0, 'ClNO3': 1.0, 'ClNO2': 1.0, \
        'Cl2': 2.0, 'OClO': 1.0, 'HOCl': 1.0, 'HCl': 1.0, 'Cl2O2': 2.0,
         'BrCl': 1.0, 'ICl':1.0, 
         # Also have reaction tracers
         'LR62': 3.0, 'LR107': 3.0, 
         'LR74' : 1.0, 'LR106':1.0, 'LR103': 1.0, 
         'LR75' : 2.0, 'LR105': 2.0, 'LR104' : 2.0, 
        }
    else:  # ( I=True is the default... )
        d = {
    'RD11': 1.0, 'RD10': 1.0, 'HIO3': 1.0, 'RD15': 1.0, 'RD62': 2.0, \
    'RD17': 1.0, 'RD16': 1.0, 'RD19': 1.0, 'LO3_37': 0.5, 'CH2I2': 2.0, \
    'AERII': 1.0, 'CH2ICl': 1.0, 'PIOx': 1.0, 'C3H7I': 1.0, 'RD73': 1.0, \
    'RD72': 2.0, 'RD71': 1.0, 'RD70': 1.0, 'C3H5I': 1.0, 'RD57': 1.0, \
    'CH3IT': 1.0, 'IO': 1.0, 'LO3_38': 1.0, 'RD61': 1.0, 'RD68': 1.0, \
    'I2': 2.0, 'IONO': 1.0, 'LO3_36': 0.6666666666666666, 'INO': 1.0, \
    'RD88': 1.0, 'RD89': 1.0, 'LOx': 1.0, 'RD06': 1.0, 'RD07': 1.0, \
    'RD02': 1.0, 'RD01': 1.0, 'I': 1.0,  'LO3_24': 0.5, 'AERI': 1.0, \
    'HOI': 1.0, 'RD64': 2.0, 'RD65': 1.0, 'RD66': 1.0, 'RD67': 1.0, \
    'RD60': 1.0, 'RD47': 1.0, 'C2H5I': 1.0, 'RD63': 1.0, 'RD20': 1.0, \
    'RD22': 1.0, 'RD24': 1.0, 'RD69': 1.0, 'RD27': 1.0, 'OIO': 1.0, \
    'CH2IBr': 1.0, 'LIOx': 1.0, 'L_Iy': 1.0, 'ICl': 1.0, 'IBr': 1.0, \
    'RD95': 2.0, 'I2O2': 2.0, 'I2O3': 2.0, 'I2O4': 2.0, 'I2O5': 2.0, \
    'HI': 1.0, 'I2O': 2.0, 'RD59': 1.0, 'RD93': 2.0, 'RD92': 1.0, \
    'IONO2': 1.0, 'RD58': 1.0, 'ISALA':1.0, 'ISALC':1.0, 
    # p/l for: IO, I
    'RD15': 1.0, 'RD17': 1.0, 'RD75': 1.0, 'RD72': 2.0, 'RD71': 1.0, \
    'RD70': 1.0, 'RD56': 1.0, 'RD69': 1.0, 'RD88': 1.0, 'RD89': 1.0, \
    'RD06': 1.0, 'RD07': 1.0, 'RD08': 1.0, 'RD64': 2.0, 'RD65': 1.0, \
    'RD67': 1.0, 'RD46': 2.0, 'RD47': 1.0, 'RD20': 1.0, 'RD22': 1.0, \
    'RD68': 1.0, 'RD25': 1.0, 'RD96': 1.0 ,'RD11': 1.0, 'RD12': 2.0, \
    'RD02': 1.0, 'RD16': 1.0, 'RD19': 1.0, 'RD24': 1.0, 'RD09': 1.0, \
    'RD23': 1.0, 'RD37': 1.0, 'RD97': 1.0, \
    # kludge for test analysis (HEMCO emissions )
    'ACET' : 1.0, 'ISOP': 1.0, 'CH2Br2': 1.0, 'CHBr3':1.0, 'CH3Br':1.0, \
    # Iodine in het loss/cycling reactions
    # loss to SSA/other aerosols
    # HOI
    'LR44': 1.0, 'LR45': 1.0, 'LR32': 1.0,   \
    # HI other
    'LR34': 1.0, \
    # IONO2
    'LR42': 1.0, 'LR43':1.0, 'LR35': 1.0, \
    # IONO    
    'LR46': 1.0, 'LR47': 1.0, 'LR39': 1.0
        }


    # Kludge for testing. Allow values to equal 1.0 if not defined. 
    try:
        if debug:
            print '{} (ref_spec: {}) stoichiometry : {}'.format( spec, \
                ref_spec,  d[spec] )
        return d[spec]

    except:
        print '!'*20, 'WARNING - Kludge assumming stoichiometry = 1.0, for'+ \
            ' {} (ref_spec given as: {})'.format( spec, ref_spec )
        return 1.0


# --------------
#  4.07 - Returns tracers unit and scale (if requested)
# --------------
def tra_unit(x, scale=False, adjustment=False, adjust=True, \
            global_unit=False, ClearFlo_unit=False, IUPAC_unit=False, \
            debug=False ):
    """ Get appropirate unit for Tracer. 

    NOTES:
	 - Is this redundent now with the species class?
     -  "Appropirate" unit is taken from GEOS-Chem input.geos
     -  Option to use IUPAC unit. ( set IUPAC_unit==True )
    """
    tra_unit = {
    'OCPI': 'ppbv', 'OCPO': 'ppbv', 'PPN': 'ppbv', 'HIO3': 'pptv', \
    'O3': 'ppbv', 'PAN': 'ppbv', 'ACET': 'ppbC', 'RIP': 'ppbv', \
    'BrNO3': 'pptv', 'Br': 'pptv', 'HBr': 'pptv', 'HAC': 'ppbv', \
    'ALD2': 'ppbC', 'HNO3': 'ppbv', 'HNO2': 'ppbv', 'C2H5I': 'pptv', \
    'HNO4': 'ppbv', 'OIO': 'pptv', 'MAP': 'ppbv', 'PRPE': 'ppbC', \
    'HI': 'pptv', 'CH2I2': 'pptv', 'IONO2': 'pptv', 'NIT': 'ppbv', \
    'CH3Br': 'pptv', 'C3H7I': 'pptv', 'C3H8': 'ppbC', 'DMS': 'ppbv', \
    'CH2O': 'ppbv', 'CH3IT': 'pptv', 'NO2': 'ppbv', 'NO3': 'ppbv', \
    'N2O5': 'ppbv', 'CHBr3': 'pptv', 'DST4': 'ppbv', 'DST3': 'ppbv', \
    'DST2': 'ppbv', 'DST1': 'ppbv', 'HOCl': 'ppbv', 'NITs': 'ppbv', \
    'RCHO': 'ppbv', 'C2H6': 'ppbC', 'MPN': 'ppbv', 'INO': 'pptv', \
    'MP': 'ppbv', 'CH2Br2': 'pptv', 'SALC': 'ppbv', 'NH3': 'ppbv', \
    'CH2ICl': 'pptv', 'IEPOX': 'ppbv', 'ClO': 'ppbv', 'NO': 'pptv', \
    'SALA': 'ppbv', 'MOBA': 'ppbv', 'R4N2': 'ppbv', 'BrCl': 'pptv', \
    'OClO': 'ppbv', 'PMN': 'ppbv', 'CO': 'ppbv', 'CH2IBr': 'pptv', \
    'ISOP': 'ppbC', 'BCPO': 'ppbv', 'MVK': 'ppbv', 'BrNO2': 'pptv', \
    'IONO': 'pptv', 'Cl2': 'ppbv', 'HOBr': 'pptv', 'PROPNN': 'ppbv', \
    'Cl': 'ppbv', 'I2O2': 'pptv', 'I2O3': 'pptv', 'I2O4': 'pptv', \
    'I2O5': 'pptv', 'MEK': 'ppbC', 'MMN': 'ppbv', 'ISOPN': 'ppbv', \
    'SO4s': 'ppbv', 'I2O': 'pptv', 'ALK4': 'ppbC', 'MSA': 'ppbv', \
    'I2': 'pptv', 'Br2': 'pptv', 'IBr': 'pptv', 'MACR': 'ppbv', 'I': 'pptv', \
    'AERI': 'pptv', 'HOI': 'pptv', 'BrO': 'pptv', 'NH4': 'ppbv', \
    'SO2': 'ppbv', 'SO4': 'ppbv', 'IO': 'pptv', 'H2O2': 'ppbv', \
    'BCPI': 'ppbv', 'ICl': 'pptv', 'GLYC': 'ppbv','ISALA': 'pptv', \
    'ISALC': 'pptv', 
    # Extra diagnostics to allow for simplified processing 
    'CH3I':'pptv', 'Iy':'pptv', 'PSURF': 'hPa', 'OH':'pptv', 'HO2':'pptv', \
    'MO2': 'pptv', 'NOy':'ppbv','EOH': 'ppbv' , 'CO':'ppbv', 'CH4':'ppbv', \
    'TSKIN':'K', 'GMAO_TEMP': 'K', 'GMAO_VWND' :'m/s',\
    'GMAO_UWND': 'm/s', 'RO2': 'pptv', 'U10M':'m/s','V10M': 'm/s' ,\
     'PRESS': 'hPa', 'CH2OO':'pptv', 'Bry':'ppbv', 'NOx': 'ppbv', 
    # Extra ClearFlo compounds
    u'acetylene': 'pptv', u'propene': 'pptv', u'Napthalene': 'pptv', \
    u'Styrene': 'pptv', u'1,3-butadiene': 'pptv', u'1,2-butadiene': 'pptv', \
    u'iso-butene': 'pptv', u'm+p-xylene': 'pptv', u'1-butene': 'pptv', \
    u't-2 pentene': 'pptv', u'cis-2-butene': 'pptv', u'1  pentene': 'pptv', \
    u'Trans-2-butene': 'pptv', u'o-xylene': 'pptv',\
    u'iso-pentane': 'pptv', u'n-hexane': 'pptv',  \
    u'iso-butane': 'pptv', u'Nonane, 2-methyl-': 'pptv', \
    u'Butane, 2,2,3-trimethyl-': 'pptv', u'Dodecane': 'pptv', \
    u'Pentane, 2,2,4-trimethyl-': 'pptv', u'2,3methylpentane': 'pptv', \
    u'Nonane': 'pptv', u'cyclopentane': 'pptv', u'n- heptane': 'pptv', \
    u'n-butane': 'pptv', u'n-pentane': 'pptv', u'Undecane': 'pptv', \
    u'Decane': 'pptv', u'Octane': 'pptv', u'n-octane': 'pptv',\
    # Extra Cly species 
    'ClNO2': 'pptv', 'ClNO3': 'pptv', 'HCl': 'pptv', 'ClOO': 'pptv', \
	'Cl2O2': 'pptv', 'CH2Cl2': 'pptv', 'CHCl3': 'pptv', 'CH3Cl': 'pptv', \
	'BrSALA': 'pptv', 'BrSALC':'pptv',\
    # extra tag "species" for  easy of processing 
    'PD421' : 'molec cm$^{-3}$ s$^{-1}$'
    } 
    units = tra_unit[x]

    # Adjust to appropriate scale for pf analysis 
    if adjust:
        spec_2_pptv = GC_var('spec_2_pptv')
        spec_2_pptC = GC_var('spec_2_pptC') 
        if ( x in spec_2_pptv ):
            if debug:
                print 'adjusting {} ({}) to {}'.format(x, units, 'pptv'   )
            units = 'pptv'
        if ( x in spec_2_pptC ):
            if debug:
                print 'adjusting {} ({}) to {}'.format(x, units, 'pptC' )
            units = 'pptC'

    # Over ride adjustments for globally appro. units
    if global_unit:
        spec_2_ppbv = GC_var('spec_2_ppbv')                    
        spec_2_ppbC = GC_var('spec_2_ppbC')
        if ( x in spec_2_ppbv ):
            if debug:
                print 'adjusting {} ({}) to {}'.format(x, units, 'ppbv'   )
            units = 'ppbv'
        if ( x in spec_2_ppbC ):
            if debug:
                print 'adjusting {} ({}) to {}'.format(x, units, 'ppbC'   )
            units = 'ppbC'

    if ClearFlo_unit:
        units2ppbv = [ 'NO', 'MACR', 'MVK', 'PRPE','ALK4', 'ALD2' ] 
        if any( [ x == i for i in units2ppbv  ] ):
            units = 'ppbv'
        if any( [ x == i for i in  [ 'PAN', ] ] ):
            units = 'pptv'
        if any( [ x == i for i in  [ 'ISOP' ] ] ):
            units = 'ppbC'

    if scale: 
        scaleby = get_unit_scaling( units )

    if adjustment:
        if units == 'K':
            units = 'Deg. Celcuis'
            adjustby = -273.15 
        else:
            adjustby = 0

    if IUPAC_unit:
        if units == 'ppbv':
            units = 'nmol mol$^{-1}$'
        if units == 'ppbC':
            units = 'nmol (C) mol$^{-1}$'
        if units == 'pptv':
            units = 'pmol mol$^{-1}$'
        if units == 'pptC':
            units = 'pmol (C) mol$^{-1}$'

    if scale and (not adjustment):
        return units, scaleby           
    elif (scale and adjustment):
        return units, scaleby, adjustby
    else:
        return units

# --------------
# 4.08 - Store of directories for servers ( earth0, atmosviz1, and tms MBP)
# ------------- 
# moved to AC_tools.funcs4core.py



# ----
#  4.10 - Return dictionary of gaw sites
# ----
def gaw_2_name():
    """ Returns dictionary GAW of sites
    """
    wdf = get_dir('dwd') +'ozonesurface/' + 'gaw_site_list.h5'
    df= pd.read_hdf( wdf,  'wp', mode='r' )
    names = df.values[:,1]

    # alter long name for CVO
    ind = [ n for n, i in enumerate(names) if  \
        ( i =='Cape Verde Atmospheric Observatory' ) ]
    names[ind[0]] = 'Cape Verde'

    return dict( zip( df.index, names ))

# ----
#  4.11 - Returns list of gaw sites in HDF file of O3 surface data
# ----
def get_global_GAW_sites(f='gaw_site_list_global.h5'):
    """ Get list of just GAW global sites. 
    """
    wd= get_dir('dwd') +'ozonesurface/' 
    df= pd.read_hdf( wd+f,  'wp', mode='r' )
    vars = sorted( list(df.index) )
    # Kludge: remove those not in "grouped" analysis  
    # ( Q: why are these sites not present?  - A: data control for lomb-scragle)
    [ vars.pop( vars.index(i) ) for i in ['ZUG', 'ZSF', 'ZEP', 'WLG', 'USH', 'SDK', 'PYR', 'PUY', 'PAL', 'MKN', 'IZO', 'HPB', 'DMV', 'BKT', 'AMS', 'ALT', 'ABP'] ]
#[ 'AMS', 'MKN', 'IZO' , 'WLG', 'PYR', 'USH', 'ABP', 'ALT'] ]
    return vars

# ----
#  4.12 - returns lon/lat/alt for res
# ----
# moved to AC_tools.funcs4core

# --------   
# 4.13 - Get CTM (GEOS-Chem) array dimension for a given resolution
# --------
# moved to AC_tools.funcs4core

# --------   
# 4.14 - Convert gamap category/species name to Iris/bpch name
# --------
def diagnosticname_gamap2iris( x  ):
    """ Convert ctm.bpch name into NetCDF normenclature
    """
    d={
    "IJ-AVG-$": 'IJ_AVG_S', 
    "BXHGHT-$": 'BXHEIGHT', 
    "PORL-L=$":'PORL_L_S__',
    'DAO-3D-$':'DAO_3D_S__',
    'DAO-FLDS' :'DAO_FLDS__',
    'DRYD-FLX': 'DRYD_FLX__',
    'DRYD-VEL':'DRYD_VEL__', 
    'CHEM-L=$':'CHEM_L_S__',
    'WETDCV-$':'WETDCV_S__', 
    'WETDLS-$':'WETDLS_S__',
    'WD-LSW-$':'WD_LSW_S__',
    'WD-LSR-$':'WD_LSR_S__',
    'UP-FLX-$':'UP_FLX_S__', 
    'NS-FLX-$': 'NS_FLX_S__', 
    'EW-FLX-$':'EW_FLX_S__', 
    'TURBMC-$': 'TURBMC_S__',
    'OD-MAP-$':'OD_MAP_S__',
#    'WD-FRC-$'
    'MC-FRC-$': 'MC_FRC_S__',
    }
    return d[x]

# --------   
# 4.15 - Get scaling for a given unit
# --------
def get_unit_scaling( units, scaleby=1 ):
    """ Get scaling for a given unit string 
    """
    misc = 'K', 'm/s', 'unitless', 'kg' ,'m', 'm2','kg/m2/s', \
            'molec/cm2/s', 'mol/cm3/s',  'kg/s', 'hPa', 'atoms C/cm2/s' \
            'kg S', 'mb', 'atoms C/cm2/s', 'molec/cm3', 'v/v', 'cm/s', 's-1', \
            'molec/m3'

    if any( [ (units ==  i) for i in 'pptv', 'pptC' ]):
        scaleby = 1E12
    elif any( [ (units ==  i) for i in 'ppbv', 'ppbC' ]):
        scaleby = 1E9
    elif any( [units ==i for i in misc ] ):
        scaleby = 1
    else:
        print 'WARNING: This unit is not in unit lists: ', units
    return scaleby

# --------   
# 4.16 - Species class for GEOS-Chem - Credit: Ben Newsome 
# --------
class species:
    """ Class for holding infomation about chemical speices. 

    NOTES:
        (1) code copied from Ben Newsome's function in MChem_tools
    ( https://github.com/tsherwen/MChem_tools/blob/master/MChem_tools.py )
        (2)  UPDATE Q? - should this be removed from MChem_tools and just keep 
        here? 

    """
    def __init__(self, name):
        self.name = name
        self.help = ("""This is a class to get information on species from a local CSV folder
   It might contain the following information:
   self.RMM    = The Mean Mass of the species.
   self.latex      = The latex name of the species.
   self.smiles   = The smiles string of the species.
   self.InChI    = The InChI string of the species.
   """)
        species_filename = os.path.dirname(__file__) + "/Species.csv"

        try:
            species_file = open(species_filename, 'rb')
        except IOError:
            print "Error: Species.csv does not appear to exist."
        species_csv = csv.reader(species_file)
        
        if (name == 'OH'):
            self.group = 'CHEM-L=$'
        else:
            self.group = 'IJ-AVG-$'

        for row in species_csv:
            try:
                if (str(self.name) == row[0].strip()):
                    self.formula   = row[1]
                    self.InChI     = row[2]
                    self.smiles    = row[3]
                    self.RMM       = float(row[4])
                    self.Latex     = row[5]
            except NameError:
                print "Species not found in CSV file"   

# --------   
# 4.17 -  Observational site class
# --------
class GEO_Site:
    """ Class for holding infomation about observational sites 
    
    """
    def __init__(self, name ):
        self.name = name
        self.help = (""" This class holds infomatio on obs. sites  """)
        wd = get_dir( 'tpwd' )+'/d_REF_dat_files/'
        filename = "ClBrI_ClNO2_FULL.dat" 

        # Check file exists
        if not os.path.exists( wd + filename ):
            print "ERROR. Is this file correct?: ",  wd + filename 
        # Open File and extract info on site
        df = pd.read_csv( wd+'/'+filename, skipinitialspace=True )

        #  Select Site
        df= df[ df['ID'] == name ]
        if len(df.index) > 0:
            # Select site
            df= df[ df['ID'] == name ]
            self.LAT = float( df['LAT'].values )
            self.LON = float( df['LON'].values )
            self.ALT = float( df['PRESS'].values ) # hPa
            self.UTC = float( df['UTC'].values ) # Time zone (UTC diff )
        else:
            print 'ERROR whilst reading site details', name  

# --------   
# 4.18 -  dictionary of category + species names from diagnostic ordering
# --------
def get_ctm_nc_var( variable ):
    """ Get number variable for diagnostic family where NetCDF import 
        from *.dat file has failed. This function returns a category name + a 
        number value to refer to diagnostic ording in NetCDF.    
    NOTES:
        (1) This is only called as a back up, and should only be used for test 
        ouput and not for production runs

    """
    d = {
    u'DAO_3D_S__CMFMC': u'DAO_3D_S___4',
     u'DAO_3D_S__DTRAIN': u'DAO_3D_S___3',
     u'DAO_3D_S__SPHU': u'DAO_3D_S___2',
     u'DAO_3D_S__TMPU': u'DAO_3D_S___1',
     u'DAO_3D_S__UWND': u'DAO_3D_S__',
     u'DAO_3D_S__VWND': u'DAO_3D_S___0'
     }
    return d[ variable ]




# --------------
# 5.02 - Store of  constants for use by funcs/progs
# --------------
def constants(input_x, rtn_dict=False, debug=False):
    """ Dictionary storing commonly used constants """
    con_dict ={
    'RMM_air' : ( .78*(2.*14.)+.22*(2.*16.) )  ,
    'AVG' : 6.0221413E23, 
    'mol2DU': 2.69E20
    }
    if rtn_dict:
        return con_dict
    else:
        return con_dict[input_x]


# ---------------- Section 6 -------------------------------------------
# -------------- Dynamic processing of p/l
#
    
# -------------
# 6.01 - Extract reactions to form a dictionary of active reactions 
# ------------- 
def rxn_dict_from_smvlog( wd, PHOTOPROCESS=None, ver='1.7', \
            LaTeX=False, debug=False ):
    """ build a dictionary reaction of reaction details from smv.log
          This can be used as an external call to analyse other prod/loss 
          reactions through smvgear  
    NOTES:
        (1) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained. 
    """
    
    if isinstance( PHOTOPROCESS, type(None) ):
        PHOTOPROCESS = {
        '1.6' : 457, '1.6.2': 452 , '1.6.3': 461 , '1.7' : 467, '2.0': 555,  \
        '3.0': 547
        }[ver]
    
    fn =  'smv2.log'
    if debug:
        print wd+'/'+fn
    file_ =  open( wd+'/'+fn, 'rb' )
    readrxn  = False
    for row in file_:
        row = row.split()
        if 'NMBR' in row:
            readrxn=True
        if len(row) < 1 :
            readrxn=False
        if  readrxn:
            try:
                rxns.append( row )
            except:
                rxns = [ row ]          

    # -- remove 'NMBR'
    rxns = [ i for i in rxns if (  'NMBR' not in i ) ]
    n = [int(rxn[0]) for rxn in rxns ]
    rxns = [rxn[1:] for rxn in rxns ]
    rdict = dict( zip(n, rxns) )
    
    # --- Process to Latex
    if LaTeX:
        for rxn in sorted( rdict.keys() ):
        
            # -- Use Latex formatting?
            if rxn > PHOTOPROCESS-1:
#                xarstr = r' $\xrightarrow{hv}$ '
                xarstr = r' + hv $\rightarrow$ '
            else:
                xarstr = r' + M $\rightarrow$ '     
            # --- get all print on a per tag basis the coe, rxn str 
            try:
                rxn_str = ''.join( rdict[rxn][4:] )
                if LaTeX:
                    rxn_str = rxn_str.replace('++=', xarstr)
                    rxn_str = rxn_str.replace('+', ' + ')
                    rxn_str = rxn_str.replace('+ =', r' $\rightarrow$ ' )
                else:
                    pass
                if debug:
                    print rxn_str
                try:
                    rxn_strs +=  [ rxn_str  ]
                    rxns += [ rxn   ]
                except:
                    rxn_strs = [ rxn_str ]
                    rxns = [ rxn ]
            except:
                print '!'*100, 'ERROR HERE: >{}<  >{}<'.format(  rxn, rxn_str )
        rdict = dict( zip(rxns, rxn_strs ) )
    return rdict
    
# -------------
# 6.02 - Extract reactions tracked by prod loss diag for a given p/l family
# ------------- 
def rxns_in_pl( wd, spec='LOX', debug=False ):
    """ Extract reactions tracked by p/l family in smvgear
    NOTES:
        (1) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained. 
    """
    
    fn =  'smv2.log'
    file_ =  open( wd+'/'+fn, 'rb' )
    if debug:
        print file_
    readrxn  = False
    # required strings in file line
    conditions = 'Family','coefficient' , 'rxns', spec
    for row in file_:
        row = row.split()
        if debug:
            print row, spec, all( [ i in row for i in conditions ] )
        if all( [ i in row for i in conditions ] ):
            readrxn=True
        if (len(row) < 1) or ( 'REACTANTS:' in row ):
            readrxn=False
        if  readrxn:
            try:
                rxns.append( row )
            except:
                rxns = [ row ]          

    # -- Check that rxns ahave been found? 
    if len( rxns ) < 1:
        print 'ERROR: No rxns. found for >{}<, correct family?'.format( spec )
        sys.exit(0)
        
    # -- remove 'Family' 
    rxns = [ i for i in rxns if (  'Family' not in i ) ]
    if debug:
        print 'number (len of list) of reacitons: ', len( rxns )
    n = [int(rxn[1]) for rxn in rxns ]
    rxns = [rxn[2:] for rxn in rxns ]

    rdict = dict( zip(n, rxns) )
    return rdict

# ------------- 
# 6.03 - Extract reaction infomation for given tracer tags
# ------------- 
def rxn4pl( pls, wd='example/example', rdict=None, reduce_size=True, \
            ver='1.7', debug=False ):
    """ Get information on reaction in smvgear from a provide reaction tag 
    NOTES:
        (1) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained. 
    """

    # ---  Get Dict of reaction detail
    if debug:
        print 'rxn4pl called'
    if isinstance(rdict, type(None) ):
        #   Get Dict of all reactions, Keys = #s
        rdict = rxn_dict_from_smvlog( wd, ver=ver ) 
        
    if debug:
        for i in rdict.keys():
            if any( [ (s_ in ''.join( rdict[i] ) ) for s_ in pls ] ):
                print i, 'yes'
#        print rdict#.keys()

    # --- Indices for 
    # reduce dict size
    if reduce_size:
        keys = [ i  for  i  in rdict.keys()  if \
            any( [ (s_ in ''.join( rdict[i] ) ) for s_ in pls ]) ] 

        # re-make dictionary         
        rdict = dict( zip( keys, [rdict[i] for i in keys] ))

    # loop via pl    
    keys = np.array([ [ pl, [ i  for  i  in rdict.keys()   \
        if any( [ (pl in ''.join( rdict[i] ) ) ]) ][0] ] for pl in pls ])

    # --- Return as reactions referenced by tag
    return  dict( zip( keys[:,0],  [ rdict[int(i)] for i in keys[:,1] ]) )
    
# -------------
# 6.04 - Construct a list of indicies for each fam from given tags
# ------------- 
def get_indicies_4_fam( tags, fam=False, IO_BrOx2=False, rtnspecs=False,
         NOy_as_HOx=True, Include_Chlorine=False, debug=False ):
    """ Return indicies (in list form) for in a given family """
    # assign family
#    famsn = [ 'Photolysis','HOx','NOy' ,'Bromine', 'Iodine' ]
    if Include_Chlorine:
        famsn = [ 'Photolysis','HOx' ,'Chlorine','Bromine', 'Iodine' ]
    else:
        famsn = [ 'Photolysis','HOx' ,'Bromine', 'Iodine' ]
    fams = []
    for tag in tags:
        fams.append( get_tag_fam( tag) )
    # if rm NOy family (treat as NOx)
    if NOy_as_HOx:
        fams = [x if (x!='NOy') else 'HOx' for x in fams]

    # Create dictionary from tags and fam assignment
    fd = dict(zip(tags, fams) )

#    print fams, tags
#    sys.exit()

    # Select tags with assigned family in list ("famsn")
    ll =  []
    [ ll.append([]) for i in famsn ]
    for n, tag in enumerate( tags ) :
        for fn in range(len(famsn)):
            if fd[tag] == famsn[fn] :
                ll[fn].append( n)

    # Consider Ox loss for 
    if fam:
#    if False: # Kludge for test. 
        # Kludge - to allow for counting Ox loss via XO +XO 50/50 between fams, 
        # ( add extra loss tag. )
        if IO_BrOx2:
#        if False:
            # Add extra tag for reaction ( IO + BrO )
            ll[famsn.index('Bromine')].append(  max([max(i) for i in ll])+1 )
            fams = fams +[ 'Bromine' ] 
            tags = tags +  [ 'LO3_24'] 
        if Include_Chlorine:
#        if False:
            # Add extra tag for reaction ( ClO + BrO )
            ll[famsn.index('Chlorine')].append(  max([max(i) for i in ll])+1 )
            fams = fams +[ 'Chlorine' ] 
            tags = tags +  [ 'LO3_82'] 
            # Add extra tag for reaction  ( ClO + IO )
            ll[famsn.index('Chlorine')].append(  max([max(i) for i in ll])+1 )
            fams = fams +[ 'Chlorine' ] 
            tags = tags +  [ 'LO3_87'] 

        if rtnspecs:
            return ll, fams, tags 
        else:
            return ll, fams            
    else:
        return ll

# -------------
# 6.05 - Get tags for reactions
# ------------- 
def get_p_l_tags( rxns, debug=False):
    """ get p/l tags for a given smvgear reaction """

    # (PD??, RD??, LO3_??, PO3_??, LR??)
    prefixs = 'PD', 'RD', 'PO3','LO3' , 'LR'

    for rxn in rxns:
        if debug:
            print [i for i in rxn if any( [ (x in i) for x in prefixs ]) ]
        tags = [i for i in rxn if any( [ (x in i) for x in prefixs ]) ]

        try:
            tagsl.append( tags)
        except:
            tagsl = [tags]

    return tagsl

# -------------
# 6.06 - extract reactions tracked by prod loss diag in input.geos
# ------------- 
def p_l_species_input_geos( wd, ver='1.7', 
            rm_multiple_tagged_rxs=False, debug=False ):
    """ Extract prod/loss species (input.geos) and reaction tags (globchem.dat) 
    NOTES:
        (1) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained. 
    """
    # find and open input.geos file
    fn = glob.glob(wd+'/*input.geos*')[0] 
    if  any( [ (i in fn) for i in '~', '#' ] ):
        print 'Trying next "input.geos" file - as FAIL for :', fn, 
        fn = glob.glob(wd+'/*input.geos*')[1] 

    if debug:
        print 'p_l_species_input_geos called using : ', wd, fn
    file_ =  open( fn, 'rb' )

    # Read in just the prod loss section 
    strs_in_1st_line =  'Number', 'of', 'P/L', 'families'
    section_line_divider = '------------------------+----------'+ \
        '--------------------------------------------'
    readrxn  = False
    for row in file_:
        row = row.split()

        # once at prod/loss section, start added to list
        if all( [ i in row for i in strs_in_1st_line ]):
            readrxn=True

        # if not at end of prod/loss section, add to list
        if section_line_divider in row:
            readrxn=False
        if  readrxn:
            try:
                rxns.append( row )
            except:
                rxns = [ row ]          

    # -- Only consider 'Family' ( no headers e.g. 'families' )
    rxns = [ i for i in rxns if (  'families' not in i ) ]
    rxns = [ [ i.replace(':','') for i in r ] for r in rxns ]

    # Kludge, adjust for extra space 12-99 
    # ( This is no longer required for 1.7 + )
    if ver == '1.6':
        [i.pop(0) for i in rxns if ('th' not in i[0]) ]  

    # Extract just PD (input.geos) and vars (globchem.dat vars ) 
    PD = [rxn[4] for rxn in rxns ]
    vars =  [rxn[5:] for rxn in rxns ]
    if debug:
        print rxns, PD, vars, ver

    # remove p/l with muliple values ( start from 12th input) - Kludge?
    if rm_multiple_tagged_rxs:
        PD, vars = [ i[11:] for i in PD, vars ]
        vars =  [ i[0] for i in  vars ]
        
    return PD, vars

# -------------
# 6.07 - extract all active tags from smv.log
# ------------- 
def tags_from_smvlog( wd ): #, spec='LOX' ):
    """  Get all active p/l tags in smvgear ( from smv2.log )
    NOTES:
        (1) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained.          
    """
    fn =  'smv2.log'
    file_ =  open( wd+'/'+fn, 'rb' )
    readrxn  = False
    for row in file_:
        row = row.split()
        if  all( [(i in row) for i in ['NBR', 'NAME', 'MW', 'BKGAS(VMRAT)'] ]):
            readrxn=True
        if len(row) < 1 :
            readrxn=False
        if  readrxn:
            try:
                rxns.append( row )
            except:
                rxns = [ row ]          

    # -- remove 'NMBR'
    rxns = [ i for i in rxns if (  'NBR' not in i ) ]
    rxns = [rxn[1] for rxn in rxns ]
    
    # --- only consider tags
    return [i for i in rxns if any( [x in i  \
        for x in 'PD', 'RD', 'PO3','LO3' , 'LR' ]) ]

# -------------
# 6.08 - extract all active PDs from smv.log
# ------------- 
def PDs_from_smvlog( wd, spec='LOX' ):
    """  Get all active PDs tags in smvgear ( from smv2.log ) 
    NOTES:
        (1) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained. 
    """
    fn =  'smv2.log'
    file_ =  open( wd+'/'+fn, 'rb' )
    readrxn  = False
    leniency =  0
    entries_in_title =  ['Families', 'for', 'prod', 'or', 'loss', 'output:']
    for row in file_:
        row = row.split()
        if  all( [(i in row) for i in entries_in_title ]):
            readrxn=True
            leniency = 1
        if len(row) < 1 :
            if leniency < 0:
                readrxn=False
            leniency -= 1
        if  readrxn:
            try:
                rxns.append( row )
            except:
                rxns = [ row ]          

    # -- remove 'NMBR'
    exceptions = [ 'SPECIES', '===============================================================================','Families']
    rxns = [ i for i in rxns if all(  [ (ii not in i) for ii in exceptions  ]) ]
    rxns =  [j for k in rxns for j in k ]
    return rxns
    
# -------------
# 6.09 - Give all reactions tag is active within
# ------------- 
def rxns4tag( tag, rdict=None, ver='1.7', wd=None ):
    """ get a list of all reactions with a given p/l tag 
    NOTES:
        (1) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained. 
    """
    # --- get reaction dictionary 
    if isinstance( rdict, type(None) ):
        rdict = rxn_dict_from_smvlog( wd, ver=ver )
    
    # --- Caveats - 
    # to adapt for long line errors in fortran written output
    errs = ['LO3_36'] #+ ['LO3_87']
    cerrs = ['RD95']  #+ ['LR48'] 
    # To account for reaction where not all channels result in Ox loss
    errs += ['RD48']
    cerrs += ['LO3_87'] 
    if any([ (tag == i) for i in  errs ] ):
        tag  = cerrs[  errs.index( tag) ]
    
    # -- loop reactions, if tag in reaction return reaction
    rxns = []
    for n, rxn in enumerate( rdict.values() ):

        expanded_rxn_str =  [i.split('+') for i in rxn ]
        expanded_rxn_str = [ \
            item for sublist in expanded_rxn_str for item in sublist]
                
        # ( Issue) Why endswith? Restore to use if contains any tag 
#        if any( [ (i.endswith(tag) ) for i in rxn]): 
        # This is because otherwise 'LR10' would be read as 'LR100'
#        if any( [tag in i for i in rxn]): # <= This will lead to false +ve
        # However, fortran print statment err for (  LO3_87 )
        if  any( [ i.endswith(tag) for i in expanded_rxn_str ] ):
            rxns.append(  [rdict.keys()[n] ]+ rxn )

    return rxns

# -------------
# 6.10 - get details for a given tag
# ------------- 
def get_tag_details( wd, tag=None, PDs=None,  rdict=None, \
            PHOTOPROCESS=None, ver='1.7', LaTeX=False, print_details=False, \
            debug=False ):
    """ Retriveve prod/loss tag details from smv.log
        ( rxn number + reaction description)
    NOTES:
        (1) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained. 
    """

    # what is the number of the first photolysis reaction?
    if isinstance( PHOTOPROCESS, type(None) ):
        PHOTOPROCESS = {
        '1.6' : 457,  '1.6.2': 452, '1.7' : 467, '2.0': 555, '3.0': 547
        }[ver]

    # ---  get all reactions tags are active in smv.log    
    if isinstance( rdict, type(None) ):
        rdict =  rxn_dict_from_smvlog( wd, ver=ver )
    trxns = rxns4tag( tag, wd=wd, rdict=rdict ) 
        
    # --- get all print on a per tag basis the coe, rxn str 
    try:
        rxn_str = ''.join( trxns[0][5:9]) 

        # -- Use Latex formatting?
#        if LaTeX:
#
            # setup latex arrows and replace existing arrows
#            if trxns[0][0] > PHOTOPROCESS-1:
#                xarstr = r' $\xrightarrow{hv}$ '
#            else:
#                xarstr = r' $\xrightarrow{M}$ '   
#            rxn_str = rxn_str.replace('++=', xarstr).replace('+', ' + ')
#            rxn_str = rxn_str.replace('+ =', r' $\rightarrow$ ' )

#        else:
#            pass
        if debug:
            print rxn_str
        dets = [ tag, trxns[0][0], rxn_str   ]   
    except:
        print '!'*100, 'ERROR HERE: >{}< >{}<'.format(  tag, trxns  )
    if print_details:
        print dets
    # --- return a dictionary of all active tagged reactions details by tag : PD, number, rxn str, coeeffiecn
    else:
        return dets
        
# -------------
# 6.11 - Takes a reaction number add gives the Ox Coe
# ------------- 
def get_rxn_Coe(wd, num, tag, nums=None, rxns=None, tags=None, \
            Coe=None, spec='LOX', ver='1.6', debug=False):
    """ Retrieve given reaction coefficient for smvgear (from smv2.log)

    NOTES:
        (1) This is no longer the case. However, previously if using dev. 
        Iy scheme, then the listed values from fuction
        Ox_in_species() will be used. 
        (2) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained. 
    """

    # --- get dictionaries for reactions within
    if all( [ (i == None) for i in nums, rxns, tags, Coe ] ):
        nums, rxns, tags, Coe = prod_loss_4_spec( wd,  spec, all_clean=True, \
                ver=ver )
    if debug:
        print nums, Coe
    
    # Pull reaction coefficient  from dictionary
    Coe_dict = dict( zip(nums, Coe) )             
    Coe = float(Coe_dict[ num ])
    # Consider all change positive - Kludge 
    # ( This is due to the assignment approach, where P=prod, L=loss )
    if ('P' not in tag ): 
        Coe =  Coe*-1.0

    # --- over write Ox in species with prod_mod_tms values
    # ( This is only used for the Ox tracers from dev. Iy simulation )
    # This has been removed to avoid possibility of double up. 
#    try:
#        Coe = Ox_in_species(tag, rxns=True)
#        if debug:
#            print 'using values from Ox_in_species'
#    except:
#        if debug:
#            print 'using values from smv.log @: {}'.format(wd), 'for ', tag, Coe
    return Coe


# --------------
# 6.12 - Remove ClBrI het loss tracers during testing
# -------------
def rm_ClBrI_het_loss( spec_l=None, r_=None, fam=None, debug=False):
    """ Allow for remove of het loss routes during testing 
            Can return species list (spec_l) + optionally 
    """

    # Print argument variables
    if debug:
        print 'before ind removal', spec_l, fam
        print [ len(i) for i in spec_l, fam ], \
    

    # --- Local variables
    rm_tracers = [ \
        'LR44', 'LR45', 'LR42', 'LR43', 'LR33', 'LR35', 'LR39', 'LR32', \
        'LR47', 'LR46']
    # -- get indices of tracers to rm, then pop from lists
    ind = [ n for n,i in enumerate( spec_l ) if ( i in rm_tracers) ]
    # remove species from list
    [ spec_l.pop(i) for i in sorted( ind )[::-1] ]         
    rtn_list = [ spec_l]
    # remove ind from fam list
    if not isinstance( fam, type(None) ):
        [ fam.pop(i) for i in sorted( ind )[::-1] ] 
        rtn_list += [ fam ]

    # remove ind from "r_" list
    if not isinstance( r_, type(None) ):
        if debug:
            print len( [item for sublist in r_ for item in sublist] ), len(r_)
        count = len( spec_l )        
        for list_ in r_[::-1]:
            for element in list_[::-1]:
                if count in ind:
                    list_.pop( list_[::-1].index( element ) )
                # reduce count
                count = count - 1
        if debug:
            print len( [item for sublist in r_ for item in sublist] )
        rtn_list += [ r_ ]

    if debug:
        print 'after ind removal', spec_l, fam, ind, sorted( ind )[::-1]
        print [ len(i) for i in spec_l, fam ], \

    return rtn_list

# --------------
# 6.13 - PD to rxn str - remake of redundant function
# -------------
#def PD_to_rxn_string( tag, ver='1.6', rdict=None, wd=None, \
#        verbose=False, debug=False ):
#    """ 
#        NOTE
#            - This is a remake of the function, the last version was lost in a 
#            reconstruction of AC_tools/PhD_progs
#    """
#
#    # Get rxn/wd dictionary if not provided
#    if isinstance( wd, type(None) ):
#        wd =  MUTD_runs( ver=ver )[0]
#    if isinstance( rdict, type(None) ):
#        rdict =  rxn_dict_from_smvlog( wd, ver=ver )
#
#    # Extract rxn details, then only use the reaction string
#    rxn_ = rxns4tag( tag=tag, rdict=rdict, ver=ver, wd=wd )
#    rxn_str =  ''.join( rxn_[0][5:9])
#    
#    return rxn_str

# --------------
# 6.14 - Get OH reactants from reaction number dictionary
# -------------
def get_pldict_reactants( pl_dict=None, only_rtn_tracers=True, \
            rm_OH=True, rm_Cl=True, tags=None, \
            debug=False ):
    """ Get reactant from smv2.log dictionary 
    NOTE: 
        (1) some reactants are not tracers ( see list: non_TRAs )  
        ( to remove these set only_rtn_tracers=True )
        (2) to remove OH from reactant list set rm_OH=True
        (3) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained. 

    """
#    non_TRAs = ['CH4', '', 'ETHLN', 'ISOPND', 'E', 'M', 'HCO', 'MVKN', 'ACTA']

    # Get reaction strings
    if isinstance( tags, type(None) ):
        tags = pl_dict.keys() 
    strs = [pl_dict[i][1] for i in tags ]    
    if debug:
        print 'rxn strs: ',  strs, len(strs)
    # remove arrows from reactions 
    strs = [ i.replace('+M=','+=').replace('+O2=','+=').replace('+N2=','+=') \
        for i in strs ]
    if debug:
        print strs, len(strs)
    # select reactants
    strs = [ i.split('+=')[0] for i in strs ]
    if debug:
        print strs, len(strs)
    if rm_OH:     # remove OH from reaction strings
        strs = [ i.replace('+OH','+') for i in strs ] 
        for n, str in enumerate( strs ):
            if str.startswith('OH+'):
                strs[n] = str.replace('OH','+')
    if rm_Cl:     # remove Cl from reaction strings
        strs = [ i.replace('+Cl','+') for i in strs ] 
        for n, str in enumerate( strs ):
            if str.startswith('Cl+'):
                strs[n] = str.replace('Cl+','+')

    # remove "+" punctuation 
    strs = [ i.replace('+','').strip() for i in strs ] 
    if debug:
        print strs, len(strs)

    return strs

# --------------
# 6.14 - 
# -------------
def get_adjustment4tags( tags, PDs=None, pl_dict=None, ver='1.6', \
            verbose=False, wd=None, Include_Chlorine=False, IO_BrOx2=False, \
            debug=False):
    """  Get coefficent for rxn tag from smv2.log using the provided family
    and adjusts the reaction tag to unity.

    The adjust to uniuty only occurs for older LO3_??/PO3_?? tags, which 
    included coefficents for the reactions in globchem.dat

    This function is a cousin to "get_rxn_Coe", but takes rxn tags as arguements 
    and adjusts a tag to unity. 
    NOTE(s):
        (1) This function is useful, but update to GEOS-Chem flexchem ( in >v11) 
        will make it redundent and therefore this is not being maintained.     
    """
        
    # --- get dictionaries for reactions + PDs, if not provided
    if isinstance( pl_dict, type(None) ):
        pl_dict = get_pl_dict( wd, spec='LOX', rmx2=True, ver=ver, debug=debug )
    if isinstance( PDs, type(None) ):
        PDs = [ PLO3_to_PD(i, wd=wd, ver=ver) for i in spec_l ]

    # Extract smv2.log change in rxn for family 
    # ( note: this won't include tag coefficients as these tag rxn not family )    
    Coes = [ pl_dict[i][-1] for i in tags ]

    # --- If species is manually known or non zero, then times by tag Coe
    # ( This is only used for the Ox tracers from dev. Iy...  simulation )
    for n, tag in enumerate( tags ):

        # Consider all change positive - Kludge 
        # ( This is due to the assignment approach, where P=prod, L=loss )
#        if ('P' not in PDs[n] ): 
        if Coes[n] < 0:
            Coes[n] =  Coes[n]*-1.0

        # Ensure all tag coefficients start from unity
        # times by tag Coefficient (if not ==1) all -Coes start from unity
        try:            
            if debug:
                print 'Accounting for (non unity) Coe in globchem.dat for:' + \
                    '{}, PD:{}, to Coe:{} (from {})'.format(  tag, PDs[n],  \
                    Coes[n]*p_l_unity(tag), Coes[n]  )
            Coes[n] = Coes[n]*p_l_unity(tag)

        except:
            if debug:
                print 'Just using Coe from smv.log @:'+ \
                    '{} for {}, PD:{}, Coe:{}'.format( wd, tag, PDs[n], Coes[n])

    # Reduce route by half if considered twice (e.g. for two families )
    for n, tag in enumerate( tags ):    
        
        # If Br + I
        if ( tag == 'LO3_24' ) and IO_BrOx2:
            if debug:
                print 'before: ', tag, Coe[n] 
            Coes[n] = Coes[n] * adjust2half4crossover( tag )
            if debug:
                print 'after: ', tag, Coe[n] 

        # If  Cl+Br+I ("Include_Chlorine")
        if any( [ (tag == i ) for i in 'LO3_87' , 'LO3_82' ] ) and \
             Include_Chlorine:
            if debug:
                print 'before: ', tag, Coes[n] 
            Coes[n] = Coes[n] * adjust2half4crossover( tag )
            if debug:
                print 'after: ', tag, Coes[n] 
    
    return Coes

def adjust2half4crossover( tag='LO3_24', ):
    """ Consider half the value for halogen cross-over reaction tags. 
        This allows for the tags to be included once per family 
    NOTE(s):
        (1) This function is redundent (and no longer in use?)

    """
    
    d = {
    'LO3_24' :  0.5,  # IO + BrO
    'LO3_87' : 0.5, # IO + ClO  
    'LO3_82' : 0.5 # BrO + ClO
    }
    
    return d[tag]

# -------------- Section 7 -------------------------------------------
# -------------- Observational Variables
#

# 5.02 - obs sites (e.g. deposition locs, lons, alts )
# 5.04 - Store site locations ( LAT,  LON,  ALT)




# ---------------- Section 2 -------------------------------------------
# -------------- Drivers
#

# --------------
# 2.01 - Convert Production/Loss RD IDs for O3 to PD## for input.geos/tracer.dat linked files
# -------------
def PLO3_to_PD(PL, fp=True, wd=None, ver='1.6', res='4x5',  \
            verbose=False, debug=False): 
    """ Converts globchem.dat tracer to PD/LD from prod/loss diag in 
    input.geos
    NOTES
        (1) 'fp' option is now obselete. 
        (2) UPDATED NEEDED: MORE DETAILED DESCRIPT.    
    """

    if verbose:
        print 'PLO3_to_PD called for wd = ', wd

    versions = [ \
    '1.3' ,'1.4' ,'1.5' , '1.6', '1.6.1','1.6.2', '1.6.3', '1.7', '2.0', '3.0' ]
    if any( [(ver ==i) for i in versions ]):

        if isinstance( wd, type(None) ):
            print 'WARNING: Using MUTD wd'
            wd = MUTD_runs(ver=ver, res=res, debug=debug)[0]
            
        # Get list of assigned PDs for vars
        PDs, vars = p_l_species_input_geos( wd, ver=ver,
             rm_multiple_tagged_rxs=True, debug=debug )

        # Add other (non 'PD') vars for ease of processing
        non_PDs = [\
        'PIOx', 'iLOX', 'LIOx', 'iPOX', 'POX', 'LOX', 'LOx', 'L_Iy', 'LOH', \
        'LCl', 'POH', 'PCl', 'P_Iy', 'L_Bry', 'P_Bry','L_Cly','P_Cly',\
        'LSBrA', 'LSBrC', 'PSBrA', 'PSBrC'
        ]
        vars += non_PDs
        PDs += non_PDs
        
        if debug:
            print vars, PDs
    
        return dict( zip(vars, PDs))[PL ]
    else:
        print 'update programme - manual PLdict now obsolete. '

# -------------
# 2.02 - Uses functions to build a dictionary for a given family of loss
# ------------- 
def get_pl_dict( wd, spec='LOX' , rmx2=False, ver='1.7', \
        rm_redundent_ClBrI_tags=False, debug=False):
    """ Get reaction IDs for each rxn. in spec (p/l, e.g. LOX) 
        This is the driver for the prod/loss programmes 
    NOTES:
        -  UPDATED NEEDED: MORE DETAILED DESCRIPT.
                
    """
    if debug:
        print 'get_pl_dict called for ', ver, spec, wd
        
    # Extract details on reactio in p/l family
    nums, rxns, tags, Coe = prod_loss_4_spec( wd,  spec, all_clean=True, \
        ver=ver, debug=debug )

    # Make a dictionary of coeffiecnts of reaction
    Coe_dict = dict( zip(nums, Coe) )

    # unpack for mulutple tags of same reactions, then get details
    unpacked_tags = [ j for k in tags for j in k ]
    
    # Kludge - remove 'PO3_10' temporarily from dictionary as these contain
    # partial tag names (  Fortran print statment cut off )
    unpacked_tags = [ i for i in unpacked_tags if ( 'PO3_10' not in i) ]
    
    if debug:
        print unpacked_tags
    
    details  = [ get_tag_details( wd, tag, ver=ver ) for tag in unpacked_tags ]
    
    # Kludge - 1 rxn missing from POx tracking? - 999  + "'ISOPND+OH', '+', '=1.0ISOPND'"
    [ details.pop(n) for n, i in enumerate( details ) if i[1]==364 ]
    ind =  [n for n, i  in enumerate( nums ) if i ==354 ]
    
    # Get coefficients for change in reaction family
    # NOTE: This does not exclude adjustment for non unity globchem.dat tags
    Coes = [ get_rxn_Coe( wd, d[1], unpacked_tags[n], nums=nums, \
                    rxns=rxns, tags=tags, Coe=Coe, spec=spec, debug=debug ) \
                    for  n, d in enumerate( details ) ]

    # Remove double ups, which are present due to Loss (LO3_??) and 
    # rate tagging (RD??) originally performed separately  
    if rmx2:
        d = [ 
        # I2O2 => AERI (LOSS) ( either RD62 or LO3_38 fine. )
        ['RD62', 'LO3_38'], 
        # IONO2 ( use LO3_30 or RD59, LR42, and LR43 )
        ['RD59', 'LO3_30'],  \
        # HOI +hv =>  ( either RD65 or LO3_34 fine. )
        ['RD65', 'LO3_34'],  \
        # I2O4 => AERI (LOSS) ( either RD93 or LO3_55 fine. )
        ['RD93', 'LO3_55'], 
        # IONO => IX/AERI (LOSS) ( use LO3_39 ( RD92 is not equiv) )
        ['RD92', 'LO3_39'], \
        # I2O3 => AERI (LOSS) ( either RD95 or LO3_36 fine. )
        [ 'RD95', 'LO3_36'],   \
        # OIO +hv =>  => AERI (LOSS) ( either RD67 or LO3_35 fine. )
        ['RD67', 'LO3_35'],  
        # HOBr + hv=>  ( either LR25 or LO3_84 fine. )
        ['LR25', 'LO3_84' ] 
        ]
        # Use Justin's 'LR??'/Johan's JTO1 tags in preference to 'LO3??' tags
        # This is an Kludge that only is necessary for "NOHAL" runs
        if rm_redundent_ClBrI_tags:
            d += [ \
        ['LR6', 'LO3_73'], ['LR5',  'LO3_73'], ['LR10', 'LO3_74'], \
        ['LR25', 'LO3_84'], ['LO3_82', 'LO3_82'], ['LO3_76','LO3_76'], \
        ['LO3_75','LO3_75'] \
            ]
            # make sure 'LO3_73' is dropped regardless of index selection
            d += [  ['LO3_73', 'LO3_73'],  ['LO3_73', 'LO3_73'], 
             ['LO3_74', 'LO3_74'], ['LO3_84', 'LO3_84']
             ]

        # Either one can be removed as currently two will be present and both # 
        # are  equally weighted by use of the Ox_in_species diagnostic ... 
        # ( However, this  only if all are equiv. )
        if rm_redundent_ClBrI_tags:
            # Drop the 2nd element in "d" list # 
            d = [i[1] for i in d ]  
        else:
            # Drop the 1st element in "d" list # <= MUST USE for Halogen sims. 
            d = [i[0] for i in d ]  
        
        ind = [ n for n, i in enumerate(details) if any( [ i[0] == ii \
            for ii in d ] )  ]
        # make sure indices ('ind') are only removed once
        ind = list( sorted( set( ind ) ) )

        if debug:
            print d, ind, [len(i) for i in details, Coes ],  \
                [ [i[0] for i in details][ii] for ii in ind ][::-1]
        # If cases have been found, remove these
        if len( ind ) > 0:
            [ l.pop(i) for i in ind[::-1] for l in details, Coes ]
        if debug:
            print  [len(i) for i in details, Coes ]
            
    # return a dictionary indexed by p/l tracer, with rxn #, 
    # reaction str and Coe of rxn.
    return dict( zip( [i[0] for i in details], [ i[1:] + [ Coes[n] ]  \
            for n, i in enumerate( details) ] ) )

# -------------
# 2.03 - Get prod loss reactions for a given family.
# ------------- 
def prod_loss_4_spec( wd, fam, all_clean=True, \
        ver='1.7', debug=False ):
    """ Retrieve reaction numbers for family of tags
    
    NOTES
        - coefficecents ("Coe") returned are for the family (e.g LOX)
        within a reaciton. ( aka not for the tag )
        -  UPDATED NEEDED: MORE DETAILED DESCRIPT.
     """

    # ---  Get Dict of all reactions, Keys = #s
    rdict = rxn_dict_from_smvlog( wd, ver=ver )

    # ---  Get reaction # tracked by p/l diag for spec and coefficient.
    rxns = rxns_in_pl( wd, fam )
    nums = rxns.keys() 
    Coe = [ rxn[-1] for rxn in rxns.values() ]

    # --- get all details from full reaction dictionary
    rxns =  [ rdict[i] for i in nums ]

    # --- get tags for tracked reactions, state where reactions are un tracked
    tags = get_p_l_tags( rxns )

    # --- cleaned tags
    if all_clean:
        tags = [ [re.sub('\+\d.\d', '',  i) for i in u ] for u in tags ]  
        tags = [ [re.sub('\=\d.\d', '',  i) for i in u ] for u in tags ]  
    
        # -- remove erroneous read/  Kludge on val
        # ---  Fortran write error leads to combination of species at the
        #  end of long line of a chemical reaction in globchem.dat
        if debug:
            print [  i[:3] for i in nums, rxns, tags, Coe]
            print [ len(i) for i in nums, rxns, tags, Coe]

        # LO3_36RD95 is present twice as this reaction is present 2 times in the code
        # update 16 01 11: LO3_36RD95 is now present x3 in version 3.0 
        #( due split uptake in iodine to aerosol )
        errs = [ \
        'LO3_36RD95', 'LO3_36RD95','LO3_36RD95', 
        'ISOPNDPO3_50', 
        'ISOPNDLR40', 
        'LO3_30LR42', 
        'LO3_30LR43', 
        'LO3_39LR46', 
        'LO3_39LR47', 
        'LO3_87LR48', 'LO3_87LR48', 
        'LISOPOLR86',  
        'LO3_50LR113', 
        'LO3_50LR114'
        ]
        cerrs = [ \
        ['LO3_36', 'RD95'], ['LO3_36', 'RD95'], ['LO3_36', 'RD95'], \
        ['PO3_50'], \
        ['LR40'], \
        ['LO3_30', 'LR42'], \
        [ 'LO3_30','LR43'], \
        ['LO3_39', 'LR46'], \
        ['LO3_39', 'LR47'], \
        ['LO3_87'], ['LO3_87'], \
        ['LR86' ], \
#        ['LO3_50', 'LR113' ], \
#        ['LO3_50', 'LR114' ]
        # revserse LR tags also ( LO3_50 now redundant? ) - NO
        ['LR113', 'LO3_50' ], \
        ['LR114', 'LO3_50' ]
#        ['LR113'], 
#        ['LR114'], 
        ]
#        errs = ['LO3_36RD95' , 'ISOPNDPO3_50', 'ISOPNDLR40']
#        cerrs = [ ['RD95'], ['PO3_50'], ['LR40'] ]
        for n, e in enumerate( errs ):
            try:
                # Get index of erroneous rxn tag
                ind = [ nn  for nn, i in enumerate( tags) if \
                                any([ ( e in ii) for ii in i ]) ] [0]
                # Extract vars for a given index
                vars =  [ i[ind] for i in nums, rxns, tags, Coe]
                if debug:
                    print 3, [ i[-1] for i in nums, rxns, tags, Coe], vars,  \
                            [ len(i) for i in nums, rxns, tags, Coe]
                # remove index ( "ind" ) value from nums, rxns, tags, and Coe 
                [ i.pop(ind) for i in nums, rxns, tags, Coe ]

                # Add the cerrs values on the end
                if debug:
                    print 4, [ i[-1] for i in nums, rxns, tags, Coe],  \
                            [ len(i) for i in nums, rxns, tags, Coe]
                nums +=  [ vars[0] ]
                rxns +=  [ vars[1] ]
                tags += [ cerrs[n] ]
                Coe  +=  [ vars[-1] ]

                if debug:
                    print 6, [ i[-1] for i in nums, rxns, tags, Coe], \
                         [ len(i) for i in nums, rxns, tags, Coe]
                    print '->'*30,  'SUCCESS >{}<  >{}<'.format( n, e )
            except:
                print '>'*50, 'FAIL (NOT REPLACED) >{}< >{}<'.format( n, e )

    # KLUDGE! - rm empty list values of ones that contain errs
#    ind = [ n for n,i in enumerate(tags) if ( (len(i)==0) or (i[0] in errs) ) ] 
#    [ [ l.pop(i) for i in sorted(ind)[::-1] ] for  l in nums, rxns, tags, Coe ]
    if debug:
        print tags 

#    print '1'*300, tags 

    return nums, rxns, tags, Coe

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

# --------------
# 4.05 -  GEOS-Chem/ctm.bpch values
# --------------
def latex_spec_name(input_x, debug=False):
    """ Formatted ( Latex ) strings for species and analysis  
        REDUNDENT: now using class structure ( see MChem_tools ) """
    spec_dict = {
    'OIO': 'OIO', 'C3H7I': 'C$_{3}$H$_{7}$I', 'IO': 'IO', 'I': 'I', \
    'I2': 'I$_{2}$', 'CH2ICl': 'CH$_{2}$ICl', 'HOI': 'HOI', \
    'CH2IBr': 'CH$_{2}$IBr', 'C2H5I': 'C$_{2}$H$_{5}$I', \
    'CH2I2': 'CH$_{2}$I$_{2}$', 'CH3IT': 'CH$_{3}$I', \
    'IONO': 'INO$_{2}$','HIO3': 'HIO$_{3}$', 'ICl': 'ICl', \
    'I2O3': 'I$_{2}$O$_{3}$', 'I2O4': 'I$_{2}$O$_{4}$', \
    'I2O5': 'I$_{2}$O$_{5}$', 'INO': 'INO', 'I2O': 'I$_{2}$O', \
    'IBr': 'IBr','I2O2': 'I$_{2}$O$_{2}$', 'IONO2': 'INO$_{3}$', 'HI':'HI', \
    'BrO':'BrO', 'Br':'Br', 'HOBr':'HOBr', 'Br2':'Br$_{2}$', \
    'CH3Br':'CH$_{3}$Br', 'CH2Br2':'CH$_{2}$Br$_{2}$', \
    'CHBr3':'CHBr$_{3}$','O3':'O$_{3}$', 'CO':'CO' , 'DMS':'DMS', \
    'NO':'NO', 'NO2':'NO$_{2}$',\
    'NO3':'NO$_{3}$','HNO3':'HNO$_{3}$', 'HNO4':'HNO$_{4}$',\
    'PAN':'PAN', 'HNO2':'HNO$_{2}$', 'N2O5':'N$_{2}$O$_{5}$',\
    'ALK4':'$\geq$C4 alkanes','ISOP':'Isoprene', 'H2O2':'H$_{2}$O$_{2}$', \
    'ACET':'CH$_{3}$C(O)CH$_{3}$', 'MEK':'>C3 ketones', \
    'RCHO': 'CH$_{3}$CH$_{2}$CHO', \
    'MVK':'CH$_{2}$=CHC(O)CH$_{3}$', 'MACR':'Methacrolein', \
    'PMN':'PMN', 'PPN':'PPN', \
    'R4N2':'$\geq$C4 alkylnitrates','PRPE':'$\geq$C3 alkenes', \
    'C3H8':'C$_{3}$H$_{8}$','CH2O':'CH$_{2}$O', \
    'C2H6':'C$_{2}$H$_{6}$', 'MP':'CH$_{3}$OOH', 'SO2':'SO$_{2}$',\
    'SO4':'SO$_{4}$','SO4s':'SO$_{4}$ on SSA', \
    'MSA':'CH$_{4}$SO$_{3}$','NH3':'NH$_{3}$', 'NH4': 'NH$_{4}$', \
    'NIT': 'InOrg N', 'NITs': 'InOrg N on SSA', 'BCPI':'BCPI', \
    'OCPI':'OCPI', 'BCPO':'BCPO','OCPO':'OCPO', 'DST1':'DST1', \
    'DST2':'DST2','DST3':'DST3','DST4':'DST4','SALA':'SALA', \
    'SALC':'SALC',  'HBr':'HBr', 'BrNO2': 'BrNO$_{2}$', \
    'BrNO3': 'BrNO$_{3}$', 'MPN':'CH$_{3}$O$_{2}$NO$_{2}$', \
    'ISOPN':'ISOPN', 'MOBA':'MOBA', 'PROPNN':'PROPNN', \
    'HAC':'HAC', 'GLYC':'GLYC', 'MMN':'MMN', 'RIP':'RIP', \
    'IEPOX':'IEPOX','MAP':'MAP', 'AERI':'Aerosol Iodine', 'Cl2':'Cl$_{2}$', \
    'Cl':'Cl','HOCl':'HOCl','ClO':'ClO','OClO':'OClO','BrCl':'BrCl', \
    'HI+OIO+IONO+INO':'HI+OIO+INO$_{2}$+INO', \
    'CH2IX':'CH$_{2}$IX (X=Cl, Br, I)', \
    'IxOy':u'I$_{2}$O$_{X}$ ($_{X}$=2,3,4)',\
    'CH3I':'CH$_{3}$I', 'OH':'OH', 'HO2':'HO$_{2}$', 'MO2':'MO$_{2}$', \
    'RO2':'RO$_{2}$' , 'ISALA': 'Iodine on SALA',  \
    'ISALC': 'Iodine on SALC', 'CH4': 'CH$_{4}$', 'MOH': 'Methanol', \
    'RD01':r'I + O$_{3}$ $\rightarrow$ IO + O$_{2}$', 
    # Adjusted names
    'ALD2':'Acetaldehyde', 
    # Analysis names 
    'iodine_all':'All Iodine', 'Iy': u'I$_{\\rm y}$',\
    'IOy': u'IO$_{\\rm y}$', \
    'IyOx': u'I$_{y}$O$_{x}$', 
    'IOx': u'IO$_{\\rm x}$', \
    'iodine_all_A':'All Iodine (Inc. AERI)',  \
    'I2Ox': u'I$_{2}$O$_{\\rm X}$' , 'AERI/SO4': 'AERI/SO4', \
    'EOH':'Ethanol','OH reactivity / s-1': u'OH reactivity / s$^{-1}$', \
    'PSURF': 'Pressure at the bottom of level', \
    'GMAO_TEMP' : 'Temperature', 'TSKIN' : 'Temperature at 2m', \
    'GMAO_UWND':'Zonal Wind', 'GMAO_VWND':'Meridional Wind', \
    'U10M':'10m Meridional Wind', 'V10M': '10m Zonal Wind', \
    'CH2OO':'CH$_{2}$OO', 'Sulfate': 'Sulfate', 'VOCs': 'VOCs', \
    'GMAO_ABSH' : 'Absolute humidity', 'GMAO_SURF': 'Aerosol surface area', \
    'GMAO_PSFC': 'Surface pressure', 
    # Family Names
    'N_specs':u'NO$_{\\rm y}$', 'NOy':u'NO$_{\\rm y}$', 
     'Bry':u'Br$_{\\rm y}$', 'Cly':u'Cl$_{\\rm y}$',  \
    'N_specs_no_I': u'NO$_{\\rm y}$ exc. iodine', 
    'NOx':u'NO$_{\\rm x}$', 'HOx':u'HO$_{\\rm x}$',\
    'SOx':u'SO$_{\\rm x}$', \
    # typos
    'CH2BR2':'CH$_{2}$Br$_{2}$',\
    # Cly names
    'ClOO':'ClOO', 'Cl2':'Cl$_{2}$', \
    'BrCl': 'BrCl','ICl': 'ICl', 'HOCl': 'HOCl', 'ClO':'ClO', 'ClOO':'ClOO', \
    'OClO':'OClO', 'Cl2O2':'Cl$_{2}$O$_{2}$', 'HCl':'HCl', \
    'ClNO2': 'ClNO$_{2}$','ClNO3':'ClNO$_{3}$', 'Cl':'Cl',\
    'CH3Cl': 'CH$_{3}$Cl',  'CH2Cl2': 'CH$_{2}$Cl$_{2}$', \
    'CHCl3': 'CHCl$_{3}$', 
    # Bry names 
    'BrSALC': 'Br- on SALC', 'BrSALA': 'Br- on SALA',
            }
    return spec_dict[input_x]
    

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


# --------------
#  7.01 - open ocean IO data  
# --------------
def IO_obs_data(just_IO=False):
    """ Dictionary of open ocean IO observations for automated comparisons

    key = ref (e.g. name_year )     
    values =  alt, lon, lat, times, full_name  , IO (avg), BrO (avg), CAM-Chem IO, CAM-Chem BrO, group
    NOTE(s):
        (1) This function is redundent (and no longer in use?)s
        (2) UPDATE NEEDED: move this to func_vars4obs
    """

    IO_obs_dict={
    'Read_2008': [0.0, -24.87, 16.85, 6.0, 'Read, 2008', '1', '2', '1', '2', 'Leeds'], 
    'Weddel_Sea': [0.03, -50.0, 75.0, 10.0, 'Weddel Sea', '-', '-', '-', '-', '-'], 
    'Oetjen_2009': [0.0, 73.5, 5.0, 2.0, 'Oetjen, 2009', '2.4', '-', '1', '-', 'Leeds'], 
    'Jones_2010_II': [0.0, -19.0, 30.0, 7.0, 'Jones, 2010, RHaMBLe II (26-36N)', '-', '-', '-', '-', '-'], 
    'Leser_2003': [0.0, -24.87, 16.85, 10.0, 'Leser, 2003', '-', '3.6', '-', '0.8', 'Heidelberg'], 
    'Jones_2010_I': [0.0, -19.0, 20.0, 7.0, 'Jones, 2010, RHaMBLe I (15-25N)', '-', '-', '-', '-', '-'], 
    'Halley_BAS': [0.03, -26.34, 75.35, 10.0, 'Halley, BAS', '-', '-', '-', '-', '-'], 
    'Theys_2007': [0.0, -20.9, -55.5, 8.0, 'Theys, 2007', '-', '<0.5', '-', '0.8', 'Belgium'], 
    'Dix_2013': [9.0, -160.0, 10.0, 1.0, 'Dix, 2013', '0.1', '-', '', '-', 'NCAR'], 
    'Allan_2000': [0.162, -16.6, 28.4, 6.0, 'Allan, 2000', '1.2', '-', '0.4', '-', 'Leeds'], 
    'Jones_2010_MAP': [0.0, -10.0, 55.0, 6.0, 'Jones, 2010, MAP', '-', '-', '-', '-', '-'], 
    'Grobmann_2013': [0.0, 150.0, 15.0, 1.0, 'Grobmann, 2013', '0.72-1.8', '-', '-', '-', 'Heidelberg'], 
    'Dix_2013_j_comp': [9.0, -160.0, 10.0, 5.0, 'Dix, 2013 (wrong month to allow with jones runs comparisons)', '-', '-', '-', '-', '-'], 
    'schonart_2009_II': [0.0, -80.0, -20.0, 10.0, 'Schonart, 2009 (satilitte)', '3.3', '-', '1', '-', 'Bremen'], 
    'Martin_2009': [0.0, -24.87, 16.85, 2.0, 'Martin, 2009', '-', '<3.0', '-', '1.2', 'Heidelberg'], 
    'schonart_2009_I': [0.0, -90.0, -5.0, 10.0, 'Schonart, 2009 (satilitte)', '3.3', '-', '1', '-', 'Bremen'], 
    'Yokouchi_2013': [0.0, 120.0, 15.0, 2.0, 'Yokouchi, 2013', '-', '-', '-', '-', '-'], 'Butz_2009': [0.045, -44.4, -2.4, 12.0, 'Butz, 2009', '0.1', '~1.0', '0.02', '0.5', 'Leeds']}
    IO_obs_dict_just = {
    'Read_2008': [0.0, -24.87, 16.85, 6.0, 'Read, 2008', '1', '2', '1', '2', 'Leeds'], 
    'schonart_2009_II': [0.0, -80.0, -20.0, 10.0, 'Schonart, 2009 (satilitte)', '3.3', '-', '1', '-', 'Bremen'], 
    'Oetjen_2009': [0.0, 73.5, 5.0, 2.0, 'Oetjen, 2009', '2.4', '-', '1', '-', 'Leeds'], 
    'Allan_2000': [0.162, -16.6, 28.4, 6.0, 'Allan, 2000', '1.2', '-', '0.4', '-', 'Leeds'], 
    'Leser_2003': [0.0, -24.87, 16.85, 10.0, 'Leser, 2003', '-', '3.6', '-', '0.8', 'Heidelberg'], 
    'Theys_2007': [0.0, -20.9, -55.5, 8.0, 'Theys, 2007', '-', '<0.5', '-', '0.8', 'Belgium'], 
    'Dix_2013': [9.0, -160.0, 10.0, 1.0, 'Dix, 2013', '0.1', '-', '', '-', 'NCAR'], 
    'Grobmann_2013': [0.0, 135.0, 15.0, 1.0, 'Grobmann, 2013', '1.095', '-', '-', '-', 'Heidelberg'], 
    'schonart_2009_I': [0.0, -90.0, -5.0, 10.0, 'Schonart, 2009 (satilitte)', '3.3', '-', '1', '-', 'Bremen'], 
    'Martin_2009': [0.0, -24.87, 16.85, 2.0, 'Martin, 2009', '-', '<3.0', '-', '1.2', 'Heidelberg'], 
    'Butz_2009': [0.045, -44.4, -2.4, 12.0, 'Butz, 2009', '0.1', '~1.0', '0.02', '0.5', 'Leeds'], 
    'Mahajan_2010': [0., -100, -10, 4.0, 'Mahajan, 2010', '0.58', '-', '-', '-', 'Leeds'] }

    if (just_IO):
        return IO_obs_dict_just
    else:
        return IO_obs_dict

# --------------
# 7.02 - BAE  flight ID to date
# -------------
def bae_flight_ID_2_date( f_ID, debug=False):
    """ BAE flight flight ID to date for CAST campaign
    NOTES:
        - UPDATE NEEDED: move this to func_vars4obs
    """

    if debug:
        print 'bae_flight_ID_2_date called'
    f_ID_dict = {'847':'2014-02-18','846':'2014-02-17','845':'2014-02-17', '844':'2014-02-16','843':'2014-02-15','842':'2014-02-17','840':'2014-02-13','839':'2014-02-12','838':'2014-02-05','837':'2014-02-04','836':'2014-02-04','835':'2014-02-03','834':'2014-02-01','833':'2014-02-01','832':'2014-01-30','831':'2014-01-30','830':'2014-01-30','829':'2014-01-28','828':'2014-01-26','827':'2014-01-26','826':'2014-01-25','825':'2014-01-24','824':'2014-01-21','823':'2014-01-18'}
    return f_ID_dict[f_ID]

# --------------
# 7.03 - details on flights from FAAM's NETCDF core data files, 
# --------------
def CAST_flight(all=True, CIMS=False, CIMSII=False):
    """ Callable dictionary of CAST flight details 
        removed 'CAST_flight, 'b822': '20140108'' 
    NOTES:
        - UPDATE NEEDED: move this to func_vars4obs
    """
    flight_dict={
    'b828': ['20140126', '0559', '20140126', '0930', 1390715993, 1390728623], 
    'b829': ['20140128', '1932', '20140129', '0224', 1390937523, 1390962280], 
    'b824': ['20140121', '1956', '20140122', '0413', 1390334165, 1390363990], 
    'b825': ['20140124', '1940', '20140125', '0140', 1390592404, 1390614058],
     'b826': ['20140125', '0128', '20140125', '0637', 1390613301, 1390631839], 
    'b827': ['20140126', '0041', '20140126', '0445', 1390696903, 1390711511], 
    'b823': ['20140118', '1713', '20140119', '0147', 1390065184, 1390096069], 
    'b839': ['20140212', '0320', '20140212', '0824', 1392175209, 1392193471], 
    'b838': ['20140205', '2246', '20140206', '0545', 1391640363, 1391665549], 
    'b837': ['20140204', '1901', '20140205', '0318', 1391540463, 1391570305], 
    'b836': ['20140204', '0134', '20140204', '0627', 1391477699, 1391495264], 
    'b835': ['20140203', '2101', '20140204', '0104', 1391461263, 1391475894], 
    'b834': ['20140201', '0610', '20140201', '1217', 1391235033, 1391257044], 
    'b833': ['20140201', '0045', '20140201', '0606', 1391215547, 1391234764], 
    'b832': ['20140130', '0622', '20140130', '1131', 1391062943, 1391081485], 
    'b831': ['20140130', '0101', '20140130', '0614', 1391043666, 1391062492], 
    'b830': ['20140129', '0239', '20140129', '0830', 1390963178, 1390984251], 
    'b846': ['20140217', '2058', '20140218', '0510', 1392670683, 1392700201], 
    'b847': ['20140218', '0503', '20140218', '0857', 1392699809, 1392713860], 
    'b844': ['20140216', '1857', '20140217', '0324', 1392577026, 1392607467], 
    'b845': ['20140217', '0323', '20140217', '0753', 1392607430, 1392623590], 
    'b842': ['20140214', '0508', '20140214', '1035', 1392354530, 1392374130], 
    'b843': ['20140215', '1903', '20140216', '0400', 1392490986, 1392523256], 
    'b840': ['20140213', '0005', '20140213', '0658', 1392249904, 1392274686],
     'b841': ['20140214', '0004', '20140214', '0458', 1392336265, 1392353901]}

    if (CIMS):
        flight_dict={'b828': ['20140126', '0559', '20140126', '0930', 1390715993, 1390728623], 'b829': ['20140128', '2322', '20140129', '0224', 1390937523, 1390962280], 'b824': ['20140122', '0155', '20140122', '0411', 1390334165, 1390363990], 'b825': ['20140124', '2230', '20140125', '0034', 1390592404, 1390614058], 'b826': ['20140125', '0128', '20140125', '0637', 1390613301, 1390631839], 'b827': ['20140126', '0041', '20140126', '0445', 1390696903, 1390711511], 'b823': ['20140118', '2006', '20140119', '0055', 1390065184, 1390096069], 'b839': ['20140212', '0320', '20140212', '0824', 1392175209, 1392193471], 'b838': ['20140205', '2246', '20140206', '0545', 1391640363, 1391665549], 'b837': ['20140204', '2333', '20140205', '0315', 1391540463, 1391570305], 'b836': ['20140204', '0134', '20140204', '0608', 1391477699, 1391495264], 'b835': ['20140203', '2316', '20140204', '0104', 1391461263, 1391475894], 'b834': ['20140201', '0610', '20140201', '1220', 1391235033, 1391257044], 'b833': ['20140201', '0145', '20140201', '0606', 1391215547, 1391234764], 'b832': ['20140130', '0622', '20140130', '1123', 1391062943, 1391081485], 'b831': ['20140130', '0201', '20140130', '0614', 1391043666, 1391062492], 'b830': ['20140129', '0249', '20140129', '0825', 1390963178, 1390984251], 'b846': ['20140217', '2058', '20140218', '0510', 1392670683, 1392700201], 'b847': ['20140218', '0503', '20140218', '0857', 1392699809, 1392713860], 'b844': ['20140216', '1857', '20140217', '0324', 1392577026, 1392607467], 'b845': ['20140217', '0323', '20140217', '0750', 1392607430, 1392623590], 'b842': ['20140214', '0611', '20140214', '1011', 1392354530, 1392374130], 'b843': ['20140215', '2227', '20140216', '0324', 1392490986, 1392523256], 'b840': ['20140213', '0304', '20140213', '0648', 1392249904, 1392274686], 'b841': ['20140214', '0004', '20140214', '0456', 1392336265, 1392353901]}
    if (CIMSII):
        flight_dict={'b828': ['20140126', '0559', '20140126', '0930', 1390715993, 1390728623], 'b829': ['20140128', '2322', '20140129', '0224', 1390937523, 1390962280], 'b824': ['20140122', '0155', '20140122', '0411', 1390334165, 1390363990], 'b825': ['20140124', '2230', '20140125', '0034', 1390592404, 1390614058], 'b826': ['20140125', '0141', '20140125', '0619', 1390613301, 1390631839], 'b827': ['20140126', '0100', '20140126', '0442', 1390696903, 1390711511], 'b823': ['20140118', '2006', '20140119', '0055', 1390065184, 1390096069], 'b839': ['20140212', '0320', '20140212', '0824', 1392175209, 1392193471], 'b838': ['20140205', '2246', '20140206', '0545', 1391640363, 1391665549], 'b837': ['20140204', '2333', '20140205', '0315', 1391540463, 1391570305], 'b836': ['20140204', '0134', '20140204', '0608', 1391477699, 1391495264], 'b835': ['20140203', '2338', '20140204', '0104', 1391461263, 1391475894], 'b834': ['20140201', '0619', '20140201', '1210', 1391235033, 1391257044], 'b833': ['20140201', '0223', '20140201', '0558', 1391215547, 1391234764], 'b832': ['20140130', '0622', '20140130', '1123', 1391062943, 1391081485], 'b831': ['20140130', '0201', '20140130', '0614', 1391043666, 1391062492], 'b830': ['20140129', '0249', '20140129', '0825', 1390963178, 1390984251], 'b846': ['20140218', '0035', '20140218', '0507', 1392670683, 1392700201], 'b847': ['20140218', '0503', '20140218', '0857', 1392699809, 1392713860], 'b844': ['20140216', '2224', '20140217', '0324', 1392577026, 1392607467], 'b845': ['20140217', '0323', '20140217', '0750', 1392607430, 1392623590], 'b842': ['20140214', '0611', '20140214', '1011', 1392354530, 1392374130], 'b843': ['20140215', '2227', '20140216', '0324', 1392490986, 1392523256], 'b840': ['20140213', '0342', '20140213', '0648', 1392249904, 1392274686], 'b841': ['20140214', '0050', '20140214', '0456', 1392336265, 1392353901]}
    if (all) :
        return flight_dict

# --------------
# 7.04 - Iodocarbon obs. meta data
# --------------
def iodocarbon_obs():
    """ dictionary of Iodocarbon observations for automated 
        analysis/comparions with observations 
    NOTES:
        - UPDATE NEEDED: move this to func_vars4obs    
    """
    Org_obs= { 
    'CH3IT'  : [  \
    ['Chuck et al (2005)' ,  list(np.linspace( -36, -49, 3) ) + \
    list( np.linspace( -20, 28, 3) ), [-20] *6,  [0.71]*3 + [ 1.94]*3  ]  ] , 
    'CH2ICl' : [  \
    ['Chuck et al (2005)' , list(  np.linspace(  -36, -49, 3 )  )+  \
    list(np.linspace( -20, 28, 3) ), [-20] *6, [0.23]*3+ [0.32]*3 ],\
    ['Jones et al (2010)' , np.linspace( 60, 15, 5), [-15] *5, [0.12]*5 ] ],
    'CH2I2' : [  \
    ['Jones et al (2010)' , np.linspace( 60,15, 5), [-15] *5, [0.01]*5 ]  ],
    'CH2IBr' : [  \
    ['Jones et al (2010)' , np.linspace( 60, 15, 5), [-15] *5, [0.01]*5 ]  ],
    'C2H5I' : [  \
    ['Jones et al (2010)' , np.linspace( 26, 36, 3), [120] *3, [0.09]*3 ]  ],
    'I2' : [  ['Lawler et al (2014)' , [16.51]   , [-24], [0.2] ]  ],
               }
    return Org_obs

# --------------
# 7.05 - Stores locations for use by funcs/progs - 
# --------------
def get_loc( loc=None, rtn_dict=False, debug=False ):
    """
        Dictionary to store locations for automated analysis

        Data arranged: LON, LAT, ALT 
        ( LON in deg E, LAT in deg N, ALT in metres a.s.l. )

        Notes:
            - double up? ( with 5.02 ?? )
            - Now use Class of GEO_site in preference to this func. 
    NOTES:
        - UPDATE NEEDED: move this to func_vars4obs    
    """

    loc_dict ={    
    # --- CAST/CONTRAST
    'GUAM' :  ( 144.800, 13.500, 0  ),
    'CHUUK' : ( 151.7833, 7.4167, 0 ),
    'PILAU' : ( 134.4667,7.3500, 0 ),      
    'London': ( -0.1275, 51.5072, 0 ),
    'Weyborne' : (1.1380, 52.9420,  0 ), 
    'WEY' : (1.1380, 52.9420,  0 ),  # Weyboure ID 
    'Cape Verde': ( -24.871, 16.848, 0 ), 
    'CVO': (-24.871,16.848,  0 ),  # Cape Verde ID
    # --- ClearFlo
    'North Ken' :  (-0.214174, 51.520718, 0),
    'KEN' :  (-0.214174, 51.520718, 0),
    'BT tower': (-0.139055, 51.521556, 190),
    'BTT': (-0.139055, 51.521556, 190),
    # --- ClNO2 sites
    'HOU' : (-95.22, 29.45, 0  ) ,
    'BOL' : ( -105.27, 40.0, 1655+ 150    ) ,
    'LAC' : ( -118.23, 34.05, 	0 ),
    'HES' : ( 8.45, 50.22, 	825 ),
    'SCH': ( 114.25, 22.22, 60 ),
    'TEX': (-95.425000, 30.350278,  60 ), 
    'CAL' : ( -114.12950, 51.07933,  1100), 
    'PAS':  ( -118.20, 34.23, 246  ), 
    # --- ClNO2 (UK) sites
    'PEN':  ( -4.1858, 50.3214 , 0  ), 
    'LEI_AUG' :  ( -1.127311, 52.619823, 0 ),
    'LEI_MAR' :  ( -1.127311, 52.619823, 0 ),
    'LEI' :  ( -1.127311, 52.619823, 0 ),
    # --- O3 preindustrial
    'MON' :  ( 2.338333, 48.822222,  75+5 ), 
#    'MON' : (2.3, 48.8, 80), # Monsoursis 
    # Pavelin  et al. (1999) / Mickely  et al. (2001)
    # 0 m. a. l. assumed following ( "<500m") in Mickely  et al. (2001)
    'ADE' : (138.0, -35.0, 0), # Adelaide 
    'COI' : (-8.0, 40.0, 0), # Coimbra 
    'HIR' : (132.0, 34.0, 0), # Hiroshima 
    'HOB' : (147.0, -43.0, 0), # Hobart 
    'HOK' : (114.0, 22.0, 0), # Hong Kong 
    'LUA' : (14.0, -9.0, 0), # Luanda 
    'MAU' : (57.0, -20.0, 0), # Mauritius 
    'MOV' : (-56.0, -35.0, 0), # Montevideo 
    'MVT' : (4.0, 44.0, 1900), # Mont-Ventoux 
    'NEM' : (145.0, 43.0, 0), # Nemuro 
    'TOK' : (139.0, 35.0, 0), # Tokyo 
    'VIE' : (16.0, 48.0, 0), # Vienna 
    'PDM' : (0.0, 43.0, 1000), # Pic du midi 
    #
    }
    if rtn_dict:
        return loc_dict
    else:
        return loc_dict[loc]

# --------------
# 7.06 - Get Locations of observations (lats, lons, alts ) for given sites
# --------------
def get_obs_loc(loc, debug=False):
    """ Dictionary to store groups of locations for automated analysis 
    NOTES:
        - UPDATE NEEDED: move this to func_vars4obs    
    """
    d = {  
    'Denmark' :[ 
    [ 68.35, 59.85, 56.13, 55.69 ], [ 18.81, 17.63, 13.05, 12.10  ] ], 
    # split by obs site...
    'Denmark1': [[68.35], [18.81]], 
    'Denmark2': [[59.85], [17.63]], 
    'Denmark3': [[56.13], [13.05]], 
    'Denmark4': [[55.69], [12.1]] ,
    # Berlin, bonn, hamburg, Westerland, Munich, Brotjacklriegel,  Deuselbach, Schauinsland
    'Germany' :[ 
    [52.5167, 50.7340, 53.5653,54.9100, 52.4231, 48.1333, \
    48.491, 49.4508, 47.91111  ] , 
    [ 13.3833 , 7.0998, 10.0014, 8.3075, 10.7872, 11.5667, 13.133, 7.302, \
     7.8894] ],  \
    'Weyborne' : [  [ 52.9420],  [1.1380 ]   ] , \
    'Penlee' : [  [ 50.3214],  [-4.1858 ]   ] ,\
    'Penlee_M2' :[  [ 49.7795272],  [-2.0229414 ]   ] , \
    'Penlee_M3' :[  [ 49.8370764,],  [-5.3652425 ]   ] , \
    'Penlee_M4' :[  [ 50.25],  [-4.15 ]   ] , \
    'Penlee_M5' :[  [ 50.25],  [-0.85  ]   ] , \
    'Penlee_M6' :[  [ 50.25],  [-7.05 ]   ] , \
    'Mace_head_M3' :[  [ 53.209003], [  -10.846408 ]   ] , \
    'Leicester' :[  [ 52.619823],  [-1.127311 ]   ]    }

    return d[loc] 
    
# --------------
# 7.07 - sonde station variables (list of 432 sondes)
# -------------
def sonde_STNs():
    """ Dictionary of WOUDC sonde location variables 
    NOTES:
        - redundent. Now using NetCDF meta data online
        - UPDATE NEEDED: move this to func_vars4obs    
    """
    
    sonde_dict = {
    101: ['SYOWA', 101.0, -69.0, 39.58, 22.0, 'JPN', 'ANTARCTICA'], \
    104: ['BEDFORD', 104.0, 42.45, -71.267, 80.0, 'USA', 'IV'], \
    105: ['FAIRBANKS (COLLEGE)', 105.0, 64.817, -147.867, 138.0, 'USA', 'IV'], \
    107: ['WALLOPS ISLAND', 107.0, 37.898, -75.483, 13.0, 'USA', 'IV'], \
    108: ['CANTON ISLAND', 108.0, -2.76, -171.7, 3.0, 'USA', 'V'], \
    109: ['HILO', 109.0, 19.5735, -155.0485, 11.0, 'USA', 'V'], \
    111: ['AMUNDSEN-SCOTT (SOUTH POLE)', 111.0, -89.983, 0.0, 2820.0, 'ATA', \
    'ANTARCTICA'], 
    131: ['PUERTO MONTT', 131.0, -41.45, -72.833, 5.0, 'CHL', 'III'], \
    132: ['SOFIA', 132.0, 42.817, 23.383, 588.0, 'BGR', 'VI'], \
    137: ['TOPEKA', 137.0, 39.067, -95.633, 270.0, 'USA', 'IV'], \
    138: ['CHRISTCHURCH', 138.0, -43.483, 172.55, 34.0, 'NZL', 'V'], \
    149: ['OVEJUYO (LA PAZ)', 149.0, -16.517, -68.033, 3420.0, 'BOL', 'III'], \
    156: ['PAYERNE', 156.0, 46.49, 6.57, 491.0, 'CHE', 'VI'], \
    157: ['THALWIL', 157.0, 46.817, 8.455, 515.0, 'CHE', 'VI'], \
    163: ['WILKES', 163.0, -66.25, 110.517, 12.0, 'USA', 'ANTARCTICA'], \
    174: ['LINDENBERG', 174.0, 52.21, 14.12, 112.0, 'DEU', 'VI'], \
    175: ['NAIROBI', 175.0, -1.267, 36.8, 1745.0, 'KEN', 'I'], \
    181: ['BERLIN/TEMPLEHOF', 181.0, 52.467, 13.433, 50.0, 'DEU', 'VI'], \
    187: ['PUNE', 187.0, 18.553, 73.86, 559.0, 'IND', 'II'], \
    190: ['NAHA', 190.0, 26.2, 127.683, 27.0, 'JPN', 'II'], \
    191: ['SAMOA', 191.0, -14.25, -170.56, 82.0, 'ASM', 'V'], \
    194: ['YORKTON', 194.0, 51.263, -102.467, 504.0, 'CAN', 'IV'], \
    197: ['BISCARROSSE/SMS', 197.0, 44.367, -1.233, 18.0, 'FRA', 'VI'], \
    198: ['COLD LAKE', 198.0, 54.783, -110.05, 702.0, 'CAN', 'IV'], \
    199: ['BARROW', 199.0, 71.317, -156.635, 11.0, 'USA', 'IV'], \
    203: ['FT. SHERMAN', 203.0, 9.33, -79.983, 57.0, 'PAN', 'IV'], \
    205: ['THIRUVANANTHAPURAM', 205.0, 8.483, 76.97, 60.0, 'IND', 'II'], \
    206: ['BOMBAY', 206.0, 19.117, 72.85, 145.0, 'IND', 'II'], \
    210: ['PALESTINE', 210.0, 31.8, -95.717, 121.0, 'USA', 'IV'], \
    213: ['EL ARENOSILLO', 213.0, 37.1, -6.733, 41.0, 'ESP', 'VI'], \
    217: ['POKER FLAT', 217.0, 65.133, -147.45, 357.5, 'USA', 'IV'], \
    219: ['NATAL', 219.0, -5.71, -35.21, 30.5, 'BRA', 'III'], \
    221: ['LEGIONOWO', 221.0, 52.4, 20.967, 96.0, 'POL', 'VI'], \
    224: ['CHILCA', 224.0, -12.5, -76.8, -1.0, 'PER', 'III'], \
    225: ['KOUROU', 225.0, 5.333, -52.65, 4.0, 'GUF', 'III'], \
    227: ['MCDONALD OBSERVATORY', 227.0, 30.666, -90.933, 2081.0, 'USA', 'IV'],\
    228: ['GIMLI', 228.0, 50.633, -97.05, 228.0, 'CAN', 'IV'], \
    229: ['ALBROOK', 229.0, 8.983, -79.55, 66.0, 'PAN', 'IV'], \
    231: ['SPOKANE', 231.0, 47.667, -117.417, 576.0, 'USA', 'IV'], \
    233: ['MARAMBIO', 233.0, -64.233, -56.623, 196.0, 'ATA', 'ANTARCTICA'], \
    234: ['SAN JUAN', 234.0, 18.483, -66.133, 17.0, 'PRI', 'IV'], \
    235: ['LONG VIEW', 235.0, 32.5, -94.75, 103.0, 'USA', 'IV'], \
    236: ['COOLIDGE FIELD', 236.0, 17.283, -61.783, 10.0, 'ATG', 'IV'], \
    237: ['GREAT FALLS', 237.0, 47.483, -111.35, 1118.0, 'USA', 'IV'], \
    238: ['DENVER', 238.0, 39.767, -104.883, 1611.0, 'USA', 'IV'], \
    239: ['SAN DIEGO', 239.0, 32.76, -117.19, 72.5, 'USA', 'IV'], \
    242: ['PRAHA', 242.0, 50.02, 14.45, 304.0, 'CZE', 'VI'], \
    254: ['LAVERTON', 254.0, -37.867, 144.75, 21.0, 'AUS', 'V'], \
    255: ['AINSWORTH (AIRPORT)', 255.0, 42.583, -100.0, 789.0, 'USA', 'IV'], \
    256: ['LAUDER', 256.0, -45.03, 169.683, 370.0, 'NZL', 'V'], \
    257: ['VANSCOY', 257.0, 52.115, -107.165, 510.0, 'CAN', 'IV'], \
    260: ['TABLE MOUNTAIN (CA)', 260.0, 34.4, -117.7, 2286.0, 'USA', 'IV'], 
    262: ['SODANKYLA', 262.0, 67.335, 26.505, 179.0, 'FIN', 'VI'], 
    265: ['IRENE', 265.0, -25.91, 28.211, 1524.0, 'ZAF', 'I'], 
    280: ['NOVOLASAREVSKAYA / FORSTER', 280.0, -70.767, 11.867, 110.0, 'ATA',\
     'ANTARCTICA'], 
    297: ['S.PIETRO CAPOFIUME', 297.0, 44.65, 11.617, 11.0, 'ITA', 'VI'], \
    303: ['IQALUIT', 303.0, 63.75, -68.55, 20.0, 'CAN', 'IV'], \
    308: ['MADRID / BARAJAS', 308.0, 40.46, -3.65, 650.0, 'ESP', 'VI'], \
    315: ['EUREKA / EUREKA LAB', 315.0, 80.04, -86.175, 310.0, 'CAN', 'IV'], \
    316: ['DE BILT', 316.0, 52.1, 5.18, 4.0, 'NLD', 'VI'], \
    318: ['VALENTIA OBSERVATORY', 318.0, 51.93, -10.25, 14.0, 'IRL', 'VI'], \
    323: ['NEUMAYER', 323.0, -70.65, -8.25, 42.0, 'ATA', 'ANTARCTICA'], \
    328: ['ASCENSION ISLAND', 328.0, -7.98, -14.42, 91.0, 'SHN', 'I'], \
    329: ['BRAZZAVILLE', 329.0, -4.28, 15.25, 314.0, 'COG', 'I'], \
    330: ['HANOI', 330.0, 21.033, 105.84, 5.0, 'VNM', 'II'], \
    333: ['PORTO NACIONAL', 333.0, -10.8, -48.4, 240.0, 'BRA', 'III'], \
    334: ['CUIABA', 334.0, -15.6, -56.1, 990.0, 'BRA', 'III'], \
    335: ['ETOSHA PAN', 335.0, -19.2, 15.9, 1100.0, 'NAM', 'I'], \
    336: ['ISFAHAN', 336.0, 32.477, 51.425, 1550.0, 'IRN', 'II'], \
    338: ['BRATTS LAKE (REGINA)', 338.0, 50.205, -104.705, 592.0, 'CAN', 'IV'],\
    339: ['USHUAIA', 339.0, -54.85, -68.308, 15.0, 'ARG', 'III'], \
    344: ['HONG KONG OBSERVATORY', 344.0, 22.31, 114.17, 66.0, 'HKG', 'II'], \
    348: ['ANKARA', 348.0, 39.95, 32.883, 896.0, 'TUR', 'VI'], \
    360: ['PELLSTON (MI)', 360.0, 45.56, -84.67, 238.0, 'USA', 'IV'], \
    361: ['HOLTVILLE (CA)', 361.0, 32.81, -115.42, -18.0, 'USA', 'IV'], \
    394: ['BROADMEADOWS', 394.0, -37.6914, 144.9467, 108.0, 'AUS', 'V'], \
    400: ['MAITRI', 400.0, -70.46, 11.45, 223.5, 'ATA', 'ANTARCTICA'], \
    401: ['SANTA CRUZ', 401.0, 28.42, -16.26, 36.0, 'ESP', 'I'], \
    404: ['JOKIOINEN', 404.0, 60.81, 23.5, 103.0, 'FIN', 'VI'], \
    406: ['SCORESBYSUND', 406.0, 70.49, -21.98, 50.0, 'GRL', 'VI'], \
    418: ['HUNTSVILLE', 418.0, 34.72, -86.64, 196.0, 'USA', 'IV'], \
    420: ['BELTSVILLE (MD)', 420.0, 39.02, -76.74, 64.0, 'USA', 'IV'], \
    432: ['PAPEETE (TAHITI)', 432.0, -18.0, -149.0, 2.0, 'PYF', 'V'], \
    434: ['SAN CRISTOBAL', 434.0, -0.92, -89.6, 8.0, 'ECU', 'III'], \
    435: ['PARAMARIBO', 435.0, 5.81, -55.21, 22.5, 'SUR', 'III'], \
    436: ['LA REUNION ISLAND', 436.0, -20.99, 55.48, 61.5, 'REU', 'I'], \
    437: ['WATUKOSEK (JAVA)', 437.0, -7.57, 112.65, 50.0, 'IDN', 'V'], \
    438: ['SUVA (FIJI)', 438.0, -18.13, 178.315, 6.0, 'FJI', 'V'], \
    439: ['KAASHIDHOO', 439.0, 5.0, 73.5, 1.0, 'MDV', 'V'], \
    441: ['EASTER ISLAND', 441.0, -27.17, -109.42, 62.0, 'CHL', 'III'], \
    443: ['SEPANG AIRPORT', 443.0, 2.73, 101.7, 17.0, 'MYS', 'V'], \
    444: ['CHEJU', 444.0, 33.5, 126.5, 300.0, 'KOR', 'II'], \
    445: ['TRINIDAD HEAD', 445.0, 40.8, -124.16, 55.0, 'USA', 'IV'], \
    448: ['MALINDI', 448.0, -2.99, 40.19, -6.0, 'KEN', 'I'], \
    450: ['DAVIS', 450.0, -68.577, 77.973, 16.0, 'ATA', 'ANTARCTICA'], \
    456: ['EGBERT', 456.0, 44.23, -79.78, 253.0, 'CAN', 'IV'], \
    457: ['KELOWNA', 457.0, 49.93, -119.4, 456.0, 'CAN', 'IV'], \
    458: ['YARMOUTH', 458.0, 43.87, -66.1, 9.0, 'CAN', 'IV'], \
    459: ['TBD', 459.0, 0.0, 0.0, 0.0, '', 'VI'], \
    460: ['THULE', 460.0, 76.53, -68.74, 57.0, 'GRL', 'VI'], \
    466: ['MAXARANGUAPE (SHADOZ-NATAL)', 466.0, -5.445, -35.33, 32.0, \
    'BRA', 'III'], \
    472: ['COTONOU', 472.0, 6.21, 2.23, 10.0, 'BEN', 'I'], \
    477: ['HEREDIA', 477.0, 10.0, -84.11, 1176.0, 'CRI', 'IV'], \
    480: ['SABLE ISLAND', 480.0, 43.93, -60.02, 4.0, 'CAN', 'IV'], \
    482: ['WALSINGHAM', 482.0, 42.6, -80.6, 200.0, 'CAN', 'IV'], \
    483: ['BARBADOS', 483.0, 13.16, -59.43, 32.0, 'BRB', 'III'], \
    484: ['HOUSTON (TX)', 484.0, 29.72, -95.4, 19.0, 'USA', 'IV'], \
    485: ['TECAMEC (UNAM)', 485.0, 19.33, -99.18, 2272.0, 'MEX', 'IV'], \
    487: ['NARRAGANSETT', 487.0, 41.49, -71.42, 21.0, 'USA', 'IV'], \
    488: ['PARADOX', 488.0, 43.92, -73.64, 284.0, 'USA', 'IV'], \
    489: ['RICHLAND', 489.0, 46.2, -119.16, 123.0, 'USA', 'IV'], \
    490: ['VALPARAISO (IN)', 490.0, 41.5, -87.0, 240.0, 'USA', 'IV'], \
    494: ['ALAJUELA', 494.0, 9.98, -84.21, 899.0, 'CRI', 'IV']
    }
    return sonde_dict 

# ----
#  7.08 - returns  (lat, lon, alt (press), timezone (UTC) ) for a given site
# ----
def gaw_2_loc(site,  f =  'GLOBAL_SURFACE_O3_2006_2012.nc' ):
    """ Extract GAW site locations for a given site 
     Another file is availible with just GAW sites:
          'GAW_SURFACE_O3_2006_2012.nc'  
    NOTE:
        - Also stores non GAW sites ( obs)
        - UPDATE NEEDED: move this to func_vars4obs    
    """
    from AC_tools.funcs4generic import hPa_to_Km

    # Use simple dictionary if site listed
    try:
        gaw_sites= {
        'SMO': (-14.247,  -170.565,1002.7885270480558, -11), \
        'MNM':(24.285, 153.981, 1011.9342452324959, 9), \
        'BMW':(32.27, -64.88, 1008.6109830510485, -4 ) ,\
        'CVO': (16.848, -24.871, 1011.6679817831093, -1), \
        'RPB':(13.17000, -59.43000, 1007.0196960034474, -4 ), \
        'ogasawara': (26.38, 142.10,996.08181619552602, 9 ), \
        'OGA': (26.38, 142.10,996.08181619552602, 9 ), \
        # Add extras for ease of analysis (e.g. Roscoff ... )
        'ROS': (48.433, -3.5904, 1011.6679817831093, +1)       
        }
        return gaw_sites[ site ] 

    # If not in list then extract details from NetCDF
    except:
        wd= get_dir('dwd') +'ozonesurface/' 
        with  Dataset(wd+f, 'r', format='NETCDF4') as f:
            lon= f.groups[site].longitude
            alt =  f.groups[site].altitude /1E3
            lat =  f.groups[site].latitude
            print [ (i, type(i) ) for i in lat, lon, alt ]
        return (lat, lon, float( hPa_to_Km([alt], reverse=True)[0] ), -9999 )