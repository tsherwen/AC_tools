#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Variable store/dictionarys for use in AC_tools

Use help(<name of function>) to get details on a particular function.


Notes
 - This code will be updated to use a user configuration approach (*.rc file) shortly
 - This module is underdevelopment vestigial/inefficient code is being removed/updated
 - Where external code is used credit is given
"""
# - Required modules:
# I/O / Low level
import re
import os
#import platform
import pandas as pd
from netCDF4 import Dataset
#import Scientific.IO.NetCDF as S
import sys
import glob
# Math/Analysis
import numpy as np

# the below import needs to be updated,
# imports should be specific and in individual functions
# import tms modules with shared functions
from .. core import *
from .. variables import get_unit_scaling


def pf_var(input, ver='3.0', ntracers=85, fill_var_with_zeroes=False):
    """
    Dictionary store for planeflight ("pf") names/tracers

    Parameters
    ----------
    input (string): tracer name to convert from GEOS-Chem to pf nomenclature
    ver (string): version number of halogen code (atop of GEOS-Chem)
    ntracers (integer): number of tracers in a given version of GEOSchem
    fill_var_with_zeroes (bool): fill tracer name with zeroes (e.g.'1'=='001')

    Returns
    -------
    (string) of tracer/output variable in pf normenclature

    Notes
    -----
     - UPDATED NEEDED: MORE DETAILED DESCRIPT.
     - Why is this function not in planeflight?

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
    OH_reactivity = [
        'NO', 'ISOP', 'GMAO_TEMP', 'GMAO_PSFC', 'CO', 'ACET', 'ALD2', 'MEK', 'MVK', 'MACR', 'C3H8', 'CH2O', 'C2H6', 'SO2', 'NO2', 'ISOPNB', 'ISOPND', 'NO3', 'HNO2', 'HNO3', 'OH', 'HO2', 'H2O2', 'MP', 'ATOOH', 'HNO4', 'ALK4', 'ISN1', 'R4N2', 'RCHO', 'ROH', 'PRPE', 'PMN', 'GLYC', 'GLYX', 'MGLY', 'HAC', 'INPN', 'PRPN', 'ETP', 'RA3P', 'RB3P', 'R4P', 'RP', 'PP', 'RIP', 'IEPOX', 'IAP', 'VRP', 'MRP', 'MAOP', 'MAP', 'DMS', 'HBr', 'Br2', 'BrO', 'CHBr3', 'CH2Br2', 'CH3Br', 'HC5', 'ISOPND', 'ISOPNB', 'ISNP', 'MVKN', 'MACRN', 'DHMOB', 'MOBA', 'ETHLN', 'PROPNN'
    ]
    OH_Extras4nic = [
        'OH', 'MCO3', 'A3O2', 'PO2', 'R4O2', 'R4O2', 'R4N1', 'ATO2', 'KO2', 'RIO2', 'VRO2', 'MRO2', 'MAN2', 'B3O2', 'INO2', 'ISNOOA', 'ISNOOB', 'ISNOHOO', 'PRN1', 'RCO3', 'MAO3', 'IEPOXOO', 'MAOPO2', 'MAOPO2', 'HC5OO', 'HC5OO', 'ISOPNDO2', 'ISOPNDO2', 'ISOPNBO2', 'ISOPNBO2', 'DIBOO', 'DIBOO', 'MOBAOO', 'MOBAOO', 'H2', 'CH4', 'HCOOH', 'MOH', 'ACTA', 'EOH', 'VRP'
    ]
    OH_rxns_17_EOH = [21,
                      328, 2, 387, 6, 7, 136, 9, 10, 407, 14, 15, 16, 408, 402, 403, 406, 23, 152, 25, 27, 28, 30, 415, 33, 35, 164, 167, 40, 169, 43, 172, 173, 48, 49, 178, 179, 52, 201, 58, 60, 394, 62, 319, 320, 66, 323, 196, 197, 198, 199, 200, 177, 202, 203, 204, 333, 206, 210, 211, 212, 334, 214, 215, 176, 345, 346, 347, 94, 357, 226, 101, 102, 103, 365, 113, 114, 153, 378, 379, 380, 364
                      ]
    OH_rxns_17_EOH = ['REA_{}'.format(i) for i in OH_rxns_17_EOH]

    # OH reaction species and tracers
    OH_rxn_tras = ['HBr', 'HOI', 'I2', 'CH2O', 'NO2', 'HNO3', 'CH3IT', 'NO', 'HNO3', 'HNO2', 'PRPE', 'PMN', 'HNO4', 'GLYC', 'NO3', 'C3H8', 'DMS', 'ALK4', 'SO2', 'MVK', 'ISOP',
                   'CHBr3', 'BrO', 'ISOP', 'CH2Br2', 'CH3Br', 'ACET', 'MEK', 'H2O2', 'CO', 'PROPNN', 'MP', 'MACR', 'HAC', 'ALD2', 'MOBA', 'RIP', 'Br2', 'IEPOX', 'MAP', 'R4N2', 'RCHO', 'HI']
    OH_rxn_tras = [num2spec(i, ver=ver, invert=True) for i in OH_rxn_tras]
    OH_rxn_tras = ['TRA_{:0>2}'.format(i) for i in OH_rxn_tras]
    OH_rxn_specs = ['GLYX', 'MGLY', 'HCOOH', 'MOH', 'O2', 'DHMOB', 'HO2', 'OH', 'H2', 'CH4', 'ETHLN', 'EOH', 'ATOOH', 'R4P', 'INPN', 'RA3P',
                    'RB3P', 'RP', 'PP', 'IAP', 'VRP', 'MRP', 'MAOP', 'ISN1', 'HC5', 'ACTA', 'ISOPNB', 'ROH', 'ISNP', 'MVKN', 'MACRN', 'ISOPND', 'ETP', 'PRPN']

    # remove inactive species
    inactive_spec = ['ACTA', 'CH4', 'H2', 'HCOOH', 'MOH', 'O2']  # 'EOH',

    # Setup list of tracers
    if any([(ver == i) for i in ('1.6', '1.6.1', '1.6.2')]):
        ntracers = 96
    elif any([(ver == i) for i in ('1.6.3', '1.6.4')]):
        ntracers = 98
    elif ver == '1.7':
        ntracers = 85
    elif any([(ver == i) for i in ('2.0', '3.0')]):
        ntracers = 103
    if ver == 'johan_br.v92':
        ntracers = 87
    # fill with zeros?
    TRAs = ['TRA_' + str(i) for i in range(1, ntracers+1)]
    if fill_var_with_zeroes:
        #    TRAs = ['TRA_{:0>2}'.format(i) for i in range(1, ntracers+1) ]
        TRAs = ['TRA_{:0>3}'.format(i) for i in range(1, ntracers+1)]

    # Setup list of reactions ( photolysis and general )
    if ver == '1.5':
        PHOT_1st, PHOT_last = 455, 533
    elif ver == '1.6':
        PHOT_1st, PHOT_last = 453, 531
    elif ver == '1.6.1':
        PHOT_1st, PHOT_last = 453, 530
    elif ver == '1.6.2':
        PHOT_1st, PHOT_last = 452, 529
    elif ver == '1.6.3':
        PHOT_1st, PHOT_last = 461, 538
    elif ver == '1.7':
        PHOT_1st, PHOT_last = 453, 529
    elif ver == '2.0':
        #        PHOT_1st, PHOT_last = 413, 614 # dev
        PHOT_1st, PHOT_last = 519, 607
    elif ver == '3.0':
        #        PHOT_1st, PHOT_last = 537, 627 # coupled sim.
        #        PHOT_1st, PHOT_last = 554, 644 # With IX split rxns
        PHOT_1st, PHOT_last = 548, 638  # With 0.25 IBr split rxns
    else:
        print(('pf_var not setup for ver: ', ver))

    JREAs = ['REA_' + str(i) for i in range(PHOT_1st, PHOT_last)]
    REAs_all = ['REA_' + str(i) for i in range(0, 533)]

    # reduced list for high time and spatial resolution
    if any([input == i for i in ('slist_v9_2_NREA_red',
                                 'slist_v9_2_NREA_red_NOy')]):
        TRAs = GC_var('active_I') + ['AERI']
        TRAs = [num2spec(i, ver=ver, invert=True) for i in TRAs]
        TRAs = ['TRA_{:0>2}'.format(i) for i in TRAs]

        metvars = [i for i in metvars if not any([(i == ii)
                                                  for ii in ('GMAO_ABSH', 'GMAO_SURF', 'GMAO_PSFC')])]
        species = [i for i in species if not any([(i == ii) for ii in
                                                  ('R4N2', 'MP', 'CH2O', 'MO2', 'ETO2', 'CO', 'C2H6', 'C3H8', 'PRPE',
                                                   'ALK4', 'ACET', 'ALD2', 'MEK', 'RCHO', 'MVK', 'DMS', 'MSA', 'ISOP')
                                                  ])]

    if input == 'slist_ClearFlo':
        TRAs = 'CO', 'ACET', 'ALD2', 'ISOP', 'C2H6', 'C3H8', 'CH2O', \
            'MACR', 'HNO2', 'HNO3', 'MVK', 'NO', 'NO2', 'PAN', 'O3',
        TRAs = [num2spec(i, ver=ver, invert=True) for i in TRAs]
        TRAs = ['TRA_{:0>2}'.format(i) for i in TRAs]
        # mannually add ethanol
        TRAs += ['TRA_86']
        species = ['OH', 'MO2', 'HO2']

    if input == 'slist_PEN_WEY_LEI':
        TRAs = 'NO', 'NO2', 'PAN', 'HNO2', 'HNO3', 'O3',\
            'CO', 'ISOP', 'R4N2', 'PRPE', 'CH2O', \
            'NO3', 'N2O5', \
            'MP', 'DMS', 'SO2', \
            'BrCl', 'Cl2', 'Cl', 'ClO', 'HOCl', 'HCl', 'ClNO2', 'ICl', \
            'SO4', 'SO4s', 'MSA', 'NH3', 'NH4', 'NIT', 'NITs', \
            'SALA', 'SALC',\
            #             'ACET', 'ALD2', 'C2H6', 'C3H8', 'MACR', 'MVK',
        TRAs = [num2spec(i, ver=ver, invert=True) for i in TRAs]
        TRAs = ['TRA_{:0>2}'.format(i) for i in TRAs]
        # Species? ( aka ones that are not transported )
        species = ['OH', 'HO2']

    if input == 'slist_v9_2_NREA_red_NOy':
        # THIS IS NOT A GOOD APPROACH, use actual names an tranlate based on verison.
        #        missing = [  'TRA_17', 'TRA_60', 'TRA_30', 'TRA_31', 'TRA_50', \
        #            'TRA_54', 'TRA_55', 'TRA_57' ]

        #            use tracer  TRA_60, TRA_30, TRA_31, for: 'MMN' , 'NH3' , 'NH4',  'R4N2', 'BrNO2', 'BrNO3','MPN', 'PROPNN',
        missing = 'MMN', 'NH3', 'NH4',  'R4N2', 'BrNO2', 'BrNO3', 'MPN', 'PROPNN'
        missing = [num2spec(i, ver=ver, invert=True) for i in missing]
        missing = ['TRA_{:0>2}'.format(i) for i in missing]
        species = species + missing

    # use TRA_?? instead of species if species is a tracer
    species_ = []
    for s in species:
        try:
            species_ += ['TRA_{}'.format(num2spec(s, ver=ver,
                                                  invert=True))]
        except:
            species_ = [s]
    species = species_

    # Construct dictionary
    d = {
        'species': species,
        'metvars': metvars,
        'REAs_all': REAs_all,
        'JREAs': JREAs,
        'TRAs': TRAs,
        'slist':  species + TRAs + JREAs + metvars,
        'slist_v9_2_NH':   species + TRAs[:66] + metvars,
        'slist_v9_2_NREA':   species + TRAs + metvars,
        'slist_v9_2_NREA_red': species + TRAs + metvars,
        'slist_REAs_all':   species + TRAs + REAs_all + metvars,
        'slist_REAs_all_OH':   species + TRAs + metvars+OH_reactivity,
        'slist_REAs_all_OH_extras':   species + TRAs + metvars,
        'slist_v9_2_NREA_red_NOy': species + TRAs + metvars,
        'slist_v10_1.7_allspecs': species + TRAs + JREAs + metvars,
        'slist_ClearFlo': species + TRAs + metvars,
        'slist_ClearFlo_OH_rxn': species + TRAs + metvars + OH_rxns_17_EOH + OH_rxn_tras + OH_rxn_specs,
        'slist_PEN_WEY_LEI': species + TRAs + metvars
    }

    # retrieve variable list from dictionary
    vars = d[input]

    # return unique list
    vars = sorted(list(set(vars)))

    # remove inactive tracers/species
    inactive_spec = [n for n, i in enumerate(vars) if (i in inactive_spec)]
    print(inactive_spec)
    [vars.pop(i) for i in sorted(inactive_spec)[::-1]]
    print(vars)

    return vars


def what_species_am_i(input=None, V_9_2=True, V_9_2_C=False, ver='1.7',
                      special_case=None, invert=False, rtn_dict=False, debug=False):
    """
    Converts a GEOS-Chem (GC) species/tracer into a PF tracer (TRA_##).
            takes TRA_## & returns GC ID or other wayround

    Parameters
    ----------
    wd (string): the wd to get the results from a run
    res (string): the resolution if wd not given ( e.g. '4x5')
    invert (bool): invert dictionary keys and values
    debug = False (legacy debug, replaced by logging)
    V_9_2, V_9_2_C = redundent oiption swicthes for previous GEOS-Chem versions
    special_case = overide seclected species dictionary

    Returns
    -------
    (string) species name in GEOS-Chem tracer naming nomenclature (or entire directory
    if rtn_dict=True)

    Notes
    -----
     - Species have the same names in PF, but units are in molec/cm3, not
       mixing ratio (v/v)
     -  Generic v10 and v11 need adding to this list

    """
    # select correct naming dictionary
    var = {
        '1.6': 'GCFP_d2TRA_all_1.6',
        '1.6.2': 'GCFP_d2TRA_all_1.6',  # same as 1.6
        '1.6.3': 'GCFP_d2TRA_all_1.6.3',  # 1.6 + 2
        '1.7': 'GCFP_d2TRA_all_1.7',
        '2.0': 'GCFP_d2TRA_all_2.0',
        '3.0':  'GCFP_d2TRA_all_2.0',  # Same as 2.0
        '4.0':  'GCFP_d2TRA_all_2.0',  # Same as 2.0
        '5.0':  'GCFP_d2TRA_all_2.0',  # Same as 2.0? - Kludge for now
    }[ver]

#    special_case = 'EOH'
#    special_case = 'EOH + actual names'
#    if all_TRA:
#            var ={  \
#            '1.7': 'all_TRA_spec_met_1.7_EOH'
#            '1.6': 'GCFP_d2TRA_1.6'
#            }[ver]

    if not isinstance(special_case, type(None)):
        var = {
            #         'EOH':'GCFP_d2TRA_all_1.7_EOH',
            #        'EOH + actual names':'GCFP_d2TRA_all_1.7_EOH_actual_names'
            'all_TRA_spec_met_1.7_EOH.zeros': 'TRA_spec_met_all_1.7_EOH',  # '
            'all_TRA_spec_met_1.7_EOH': 'TRA_spec_met_all_1.7_EOH_no_trailing_zeroes'
            #        TRA_spec_met_all_1'
        }[special_case]

    # Get dictionary from variable store
    d = GC_var(var)

    if debug:
        print((d, special_case))

    if invert:
        d = {v: k for k, v in list(d.items())}

    # return dictionary
    if rtn_dict:
        return d
    else:
        return d[input]


def num2spec(num=69, rtn_dict=False, invert=False, ver='1.7'):
    """
    Parameters
    ----------
    ver (str): Specify the wd to get the results from a run.
    num (int): the resolution if wd not given (e.g. '4x5' )
    rtn_dict (bool): legacy debug option, replaced by python logging
    invert (bool): uses spec as the keys in the dictionary

    Returns
    -------
    (int) number of tracer in GEOS-Chem or dictionary of tracer names and indicies (dict)

    Notes
    -----
     - Also works in reverse or returns whole dictionary to remove need to remake
    dictionary for each call.
     - Version number needed ( e.g. "Cl+Br+I" = 3.0, "1.7" = "Br+I",
     - UPDATE NEEDED: basecase v9-2/v10 not currently included
    """

    # --- Get dictionary of tracer numbers
    d = what_species_am_i(ver=ver, rtn_dict=True, special_case=None)

    # --- Slice off just numbers
    # special case for dev version?
    if any([(ver == i) for i in ('1.6', '1.6.2',)]):
        d = GC_var('GCFP_d2TRA_justTRA_1.6')
    if any([(ver == i) for i in('1.6.3', '1.6.4')]):
        d = GC_var('GCFP_d2TRA_justTRA_1.6.3')
    # Then slice
    nums = [int(i[4:]) for i in list(d.keys())]

    # --- Re-make dictionary
    d = dict(list(zip(nums, list(d.values()))))

    # --- Invert to give spec for num
    if invert:
        d = {v: k for k, v in list(d.items())}

    if rtn_dict:
        return d
    else:
        return d[num]


def diagnosticname_gamap2iris(x):
    """
    Convert ctm.bpch name into NetCDF nomenclature
    """
    d = {
        "IJ-AVG-$": 'IJ_AVG_S',
        "BXHGHT-$": 'BXHEIGHT',
        "PORL-L=$": 'PORL_L_S__',
        'DAO-3D-$': 'DAO_3D_S__',
        'DAO-FLDS': 'DAO_FLDS__',
        'DRYD-FLX': 'DRYD_FLX__',
        'DRYD-VEL': 'DRYD_VEL__',
        'CHEM-L=$': 'CHEM_L_S__',
        'WETDCV-$': 'WETDCV_S__',
        'WETDLS-$': 'WETDLS_S__',
        'WD-LSW-$': 'WD_LSW_S__',
        'WD-LSR-$': 'WD_LSR_S__',
        'UP-FLX-$': 'UP_FLX_S__',
        'NS-FLX-$': 'NS_FLX_S__',
        'EW-FLX-$': 'EW_FLX_S__',
        'TURBMC-$': 'TURBMC_S__',
        'OD-MAP-$': 'OD_MAP_S__',
        #    'WD-FRC-$'
        'MC-FRC-$': 'MC_FRC_S__',
    }
    return d[x]


def get_ctm_nc_var(variable):
    """
    Get number variable for diagnostic family where NetCDF import
    from *.dat file has failed. This function returns a category name + a
    number value to refer to diagnostic ording in NetCDF.

    Parameters
    ----------
    variable (str): NetCDF variable within number suffix

    Returns
    -------
    variable (str): NetCDF variable in GEOS-Chem output

    Notes
    -----
     - This is only called as a back up, and should only be used for test
    ouput and not for production runs

    """
    d = {
        'DAO_3D_S__CMFMC': 'DAO_3D_S___4',
        'DAO_3D_S__DTRAIN': 'DAO_3D_S___3',
        'DAO_3D_S__SPHU': 'DAO_3D_S___2',
        'DAO_3D_S__TMPU': 'DAO_3D_S___1',
        'DAO_3D_S__UWND': 'DAO_3D_S__',
        'DAO_3D_S__VWND': 'DAO_3D_S___0'
    }
    return d[variable]




# -------------- Non-generic Functions
#
# NOTE(s):
# (1) These are retained, but will be migrated to a seperate non-generic module
# (2) It is not advised to use these.


def gaw_2_name():
    """
    Returns dictionary GAW of sites
    """
    wdf = get_dir('dwd') + 'ozonesurface/' + 'gaw_site_list.h5'
    df = pd.read_hdf(wdf,  'wp', mode='r')
    names = df.values[:, 1]
    # alter long name for CVO
    ind = [n for n, i in enumerate(names) if
           (i == 'Cape Verde Atmospheric Observatory')]
    names[ind[0]] = 'Cape Verde'

    return dict(list(zip(df.index, names)))


def get_global_GAW_sites(f='gaw_site_list_global.h5'):
    """
    Get list of (just) GAW global sites

    Notes
    -----
     - This data is from Sofen et al 2016 and freely availible.
     (BADC doi:10.5285/08fbe63d-fa6d-4a7a-b952-5932e3ab0452 )
    """
    wd = get_dir('dwd') + 'ozonesurface/'
    df = pd.read_hdf(wd+f,  'wp', mode='r')
    vars = sorted(list(df.index))
    # Kludge: remove those not in "grouped" analysis
    # These sites are not included due to a data control for lomb-scragle work
    sites2exclude = [
        'ZUG', 'ZSF', 'ZEP', 'WLG', 'USH', 'SDK', 'PYR', 'PUY', 'PAL',
        'MKN', 'IZO', 'HPB', 'DMV', 'BKT', 'AMS', 'ALT', 'ABP'
    ]
    [vars.pop(vars.index(i)) for i in sites2exclude]
    return vars




