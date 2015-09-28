# =================================================
# --------- tms - module of Variables for re-use----------------
# -------------- 

# Section 0 - Required modules
# Section 1 - Planeflight variables
# Section 2 - Drivers functions
# Section 3 - GeosChem (bpch) prod loss variables
# Section 4 - GeosChem (bpch) general variables
# Section 5 - Misc
# Section 6 - Dynamic p/l processing
# Section 7 - Obervational variables

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
# 2.02 - Get P/L dictionary for a given family ***
# 2.03 - Get P/L dictionary for a given species ***

# --------------- ------------- ------------- -------------
# ---- Section 3 ----- GeosChem (bpch) prod loss variables
# 3.01 - Spec to photolysis reaction p/l tag 
# 3.02 - Ox family for tag

# --------------- ------------- ------------- -------------
# ---- Section 4 ----- GeosChem (bpch) general variables
 # 4.01 - v9-2 species in input.geos from num
# 4.02 - Get Species Mass 
# 4.03 - Get Species stioch ***
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
# 4.99 - Reference data, (inc. grid data) from gchem - credit: GK (Gerrit Kuhlmann )

# --------------- ------------- ------------- -------------
# ---- Section 5 ----- Misc
# 5.01 - dir store (for standard dirs on different computers )
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

# --------------- ------------- ------------- -------------
# ---- Section 7 ----- Observational variables
# 7.01 - IO observation dictionary
# 7.02 - BAE flight ID dictionary
# 7.03 - CAST flight dictionary for CIMS/CIMSII
# 7.04 - Iodocarbon obs. meta data
# 7.05 - Stores locations for use by funcs/progs - LON, LAT,   ALT - double up? ( with 5.02 ?? )
# 7.06 - Get Locations of observations (lats, lons, alts ) for given sites
# 7.07 - sonde station variables (list of 432 sondes)
# 7.08 - returns  (lat, lon, alt (press), timezone (UTC) ) for a given site


# ------------------------------------------- Section 0 -------------------------------------------
# -------------- Required modules:
#
#!/usr/bin/python
#
# -- I/O / Low level                                                                                
import re
import platform
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
def pf_var( input, ver='1.7', ntracers=85 ):

    # planeflight variable lists
    metvars = [
    'GMAO_TEMP', 'GMAO_ABSH', 'GMAO_SURF', 'GMAO_PSFC', 'GMAO_UWND', 'GMAO_VWND'
    ]
    species  = [
    'O3', 'NO2', 'NO', 'NO3', 'N2O5', 'HNO4', 'HNO3', 'HNO2', 'PAN', 'PPN', 
    'PMN', 'R4N2', 'H2O2', 'MP', 'CH2O', 'HO2', 'OH', 'RO2', 'MO2', 'ETO2', 
    'CO', 'C2H6', 'C3H8', 'PRPE', 'ALK4', 'ACET', 'ALD2', 'MEK', 'RCHO', 'MVK', 
    'SO2', 'DMS', 'MSA', 'SO4', 'ISOP'
    ]
    all_species_not_TRA = [
    'A3O2', 'ATO2', 'B3O2', 'EOH', 'ETO2', 'ETP', 'GLYX', 'HO2', 'IAP', 'INO2', 'INPN', 'ISN1', 'ISNOOA', 'ISNOOB', 'ISNOHOO', 'ISNP', 'KO2', 'MAN2', 'MAO3', 'MAOP', 'MAOPO2', 'MCO3', 'MGLY', 'MO2', 'MRO2', 'MRP', 'OH', 'PO2', 'PP', 'PRN1', 'PRPN', 'R4N1', 'R4O2', 'R4P', 'RA3P', 'RB3P', 'RCO3', 'RIO2', 'ROH', 'RP', 'VRO2', 'VRP', 'LISOPOH', 'ISOPND', 'ISOPNB', 'HC5', 'DIBOO', 'HC5OO', 'DHMOB', 'MOBAOO', 'ISOPNBO2', 'ISOPNDO2', 'ETHLN', 'MACRN', 'MVKN', 'PYAC', 'IEPOXOO', 'ATOOH', 'PMNN', 'MACRNO2', 'PMNO2'
    ] 
    OH_reactivity=[
    'NO', 'ISOP', 'GMAO_TEMP', 'GMAO_PSFC', 'CO', 'ACET', 'ALD2', 'MEK', 'MVK', 'MACR', 'C3H8', 'CH2O', 'C2H6', 'SO2', 'NO2', 'ISOPNB', 'ISOPND', 'NO3', 'HNO2', 'HNO3', 'OH', 'HO2', 'H2O2', 'MP', 'ATOOH', 'HNO4', 'ALK4', 'ISN1', 'R4N2', 'RCHO', 'ROH', 'PRPE', 'PMN', 'GLYC', 'GLYX', 'MGLY', 'HAC', 'INPN', 'PRPN', 'ETP', 'RA3P', 'RB3P', 'R4P', 'RP', 'PP', 'RIP', 'IEPOX', 'IAP', 'VRP', 'MRP', 'MAOP', 'MAP', 'DMS', 'HBr', 'Br2', 'BrO', 'CHBr3', 'CH2Br2', 'CH3Br', 'HC5', 'ISOPND', 'ISOPNB', 'ISNP', 'MVKN', 'MACRN', 'DHMOB', 'MOBA', 'ETHLN', 'PROPNN'
    ] 
    OH_Extras4nic = [
    'OH', 'MCO3', 'A3O2', 'PO2', 'R4O2', 'R4O2', 'R4N1', 'ATO2', 'KO2', 'RIO2', 'VRO2', 'MRO2', 'MAN2', 'B3O2', 'INO2', 'ISNOOA', 'ISNOOB', 'ISNOHOO', 'PRN1', 'RCO3', 'MAO3', 'IEPOXOO', 'MAOPO2', 'MAOPO2', 'HC5OO', 'HC5OO', 'ISOPNDO2', 'ISOPNDO2', 'ISOPNBO2', 'ISOPNBO2', 'DIBOO', 'DIBOO', 'MOBAOO', 'MOBAOO', 'H2', 'CH4', 'HCOOH', 'MOH', 'ACTA', 'EOH', 'VRP'
    ]

    # remove inactive species 
    inactive_spec = ['ACTA', 'CH4', 'EOH', 'H2', 'HCOOH', 'MOH']
    [ OH_Extras4nic.pop(ii) for ii in sorted([ OH_Extras4nic.index(i)  \
        for i in inactive_spec ])[::-1] ]

    # Deal 
    if ver == '1.7':
        ntracers=85
    if ver == 'johan_br.v92':
        ntracers=87
        JREAs=[]
    TRAs = ['TRA_'+ str(i) for i in range(1, ntracers+1) ] 
#    TRAs = ['TRA_{:0>2}'.format(i) for i in range(1, ntracers+1) ]
    if ver == '1.7':
        JREAs = ['REA_'+ str(i) for i in range(453, 529) ] 
    if ver == '1.6':
        JREAs = ['REA_'+ str(i) for i in range(453, 531) ] 
    if ver == '1.6.1':
        JREAs = ['REA_'+ str(i) for i in range(453, 530) ] 
    if ver == '1.5':
        JREAs = ['REA_'+ str(i) for i in range(455, 533) ] 
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

    if input =='slist_v9_2_NREA_red_NOy':
# THIS IS NOT A GOOD APPROACH, use actual names an tranlate based on verison. 
#        missing = [  'TRA_17', 'TRA_60', 'TRA_30', 'TRA_31', 'TRA_50', \
#            'TRA_54', 'TRA_55', 'TRA_57' ]
        
#            use tracer  TRA_60, TRA_30, TRA_31, for: 'MMN' , 'NH3' , 'NH4',  'R4N2', 'BrNO2', 'BrNO3','MPN', 'PROPNN', 
        missing = 'MMN' , 'NH3' , 'NH4',  'R4N2', 'BrNO2', 'BrNO3','MPN', 'PROPNN'
        missing = [ num2spec( i, ver=ver, invert=True) for i in missing ]
        missing = [ 'TRA_{:0>2}'.format( i)  for i in missing ]
        species = species  + missing


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
     'slist_REAs_all_OH_extras' :   all_species_not_TRA + TRAs  + metvars, 
    'slist_v9_2_NREA_red_NOy' : species + TRAs + metvars,
    'slist_v10_1.7_allspecs': all_species_not_TRA +TRAs+ JREAs +metvars
      } 

    # 
    vars = d[input]
    vars = sorted( list( set( vars) ) )
    print vars

    return vars

# --------------
# 1.02 - What GEOS-Chem (GC) Species am i? takes TRA_## & returns GC ID or other wayround 
# -------------
def what_species_am_i(input=None, V_9_2=True, V_9_2_C=False, \
            ver='1.7', special_case=None, invert=False, rtn_dict=False, \
            debug=False ) :

    # select correct naming dictionary
    var ={  \
    '1.7': 'GCFP_d2TRA_all_1.7',
    '1.6': 'GCFP_d2TRA_all_1.6'
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
#        'all_TRA_spec_met_1.7_EOH' : 'TRA_spec_met_all_1.7_EOH' #' 
       'all_TRA_spec_met_1.7_EOH':'TRA_spec_met_all_1.7_EOH_no_trailing_zeroes'
#        TRA_spec_met_all_1'
        }[special_case]

    # Get dictionary from variable store
    d =  GC_var( var ) 

    if debug:
        print d

    if invert:
        d = {v: k for k, v in d.items()}

    # return dictionary
    if rtn_dict:
        return d    
    else:
        return d[input]


# ---------------- Section 2 -------------------------------------------
# -------------- Drivers
#

# --------------
# 2.01 - Convert Production/Loss RD IDs for O3 to PD## for input.geos/tracer.dat linked files
# -------------
def PLO3_to_PD(PL, fp=True, wd=None, ver='1.6', res='4x5',debug=False): 

    """ Converts """

    if any( [(ver ==i) for i in  '1.3' ,'1.4' ,'1.5' , '1.6', '1.7' ]):

        if wd==None:
            if debug:
                print 'WARNING: Using MUTD wd'
            wd = MUTD_runs(ver=ver, res=res, debug=debug)[0]
        PDs, vars = p_l_species_input_geos( wd, ver=ver,
             rm_multiple_tagged_rxs=True)

        # Add other vars for ease of processing
        vars += ['PIOx', 'iLOX', 'LIOx', 'iPOX', 'POX', 'LOX', 'LOx', 'L_Iy']
        PDs += ['PIOx', 'iLOX', 'LIOx', 'iPOX', 'POX', 'LOX', 'LOx', 'L_Iy']
    
        return dict( zip(vars, PDs))[PL ]
    else:
        print 'update programme - manual PLdict now obsolete. '

# -------------
# 2.02 - DRIVER - uses functions to build a dictionary for a given family of loss
# ------------- 
def get_pl_dict( wd, spec='LOX' , rmx2=False, debug=False):
    # Get reaction IDs for each rxn. in spec (p/l, e.g. LOX)
    nums, rxns, tags, Coe = prod_loss_4_spec( wd,  spec, all_clean=True, debug=debug )

    # Make a dictionary of coeffiecnts of reaction
    Coe_dict = dict(zip(nums, Coe) )

    # unpack for mulutple tags of same reactions, then get details
    unpacked_tags = [j for k in tags for j in k ]
    details  = [  get_tag_details( wd, tag  ) for  tag in unpacked_tags ]

    # Kludge - 1 rxn missing from POx tracking? - 999  + "'ISOPND+OH', '+', '=1.0ISOPND'"
    [ details.pop(n) for n, i in enumerate( details ) if i[1]==364]
    ind =  [n for n, i  in enumerate( nums ) if i ==354 ]
    
    # Get Coes and overwrite where prog_mod_tms has values
    Coes = [  get_rxn_Coe( wd, d[1], unpacked_tags[n], nums=nums, rxns=rxns, tags=tags, Coe=Coe, spec=spec, debug=debug ) 
                    for  n, d in enumerate( details ) ]

    # Remove double ups, which are present due to Loss (LO3_??) and rate tagging (RD??) originally performed separately  
    if rmx2:
        d = [ 
        ['RD62', 'LO3_38'], ['RD59', 'LO3_30'], ['RD65', 'LO3_34'],  \
        ['RD93', 'LO3_55'], ['RD92', 'LO3_39'], [ 'RD95', 'LO3_36'],   \
        ['RD67', 'LO3_35'] 
        ]
        d = [i[0] for i in d ]
        ind = [ n for n, i in enumerate(details) if any( [ i[0] == ii \
            for ii in d ] )  ]
        if debug:
            print d, ind, [len(i) for i in details, Coes ] ,  [ [i[0] \
                for i in details][ii] for ii in ind ][::-1]
        [ l.pop(i) for i in ind[::-1] for l in details, Coes ]
        if debug:
            print  [len(i) for i in details, Coes ]
            
    # return a dictionary indexed by p/l tracer, with rxn #, 
    # reaction str and Coe of rxn.
    return dict( zip( [i[0] for i in details], [ i[1:] + [ Coes[n] ]  for n, i in enumerate( details) ] ) )

# -------------
# 2.03 - Get prod loss reactions for a given family.
# ------------- 
def prod_loss_4_spec( wd, fam, all_clean=True, debug=False ):
    # ---  Get Dict of all reactions, Keys = #s
    rdict = rxn_dict_from_smvlog(wd)

    # ---  Get reaction # tracked by p/l diag for spec and coefficient.
    rxns = rxns_in_pl(wd, fam)
    nums =  rxns.keys() 
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
        #  end of long line 
        if (debug):
            print [  i[:3] for i in nums, rxns, tags, Coe]
            print [ len(i) for i in nums, rxns, tags, Coe]

        errs = ['LO3_36RD95' , 'ISOPNDPO3_50']
        cerrs = [['LO3_36', 'RD95'], ['PO3_50'] ]
        for n, e in enumerate( errs ):
            try:
                ind = [ nn  for nn, i in enumerate( tags) if any([ ( e in ii) for ii in i ]) ] [0]
                vars =  [ i[ind] for i in nums, rxns, tags, Coe]
                if (debug):
                    print 3, [ i[-1] for i in nums, rxns, tags, Coe], vars,  [ len(i) for i in nums, rxns, tags, Coe]
                [i.pop(ind) for i in nums, rxns, tags, Coe ]

                # add the cerrs values on the end
                if (debug):
                    print 4, [ i[-1] for i in nums, rxns, tags, Coe],  [ len(i) for i in nums, rxns, tags, Coe]
                nums +=  [ vars[0] ]
                rxns   +=  [vars[1] ]
                tags    += [cerrs[n]]
                Coe    +=  [vars[-1] ]

                if (debug):
                    print 6, [ i[-1] for i in nums, rxns, tags, Coe], \
                         [ len(i) for i in nums, rxns, tags, Coe]
                    print '->'*30,  'SUCCESS', n, e
            except:
                print '>'*100, 'FAIL' , n, e     

    return nums, rxns, tags, Coe

# ------------------------------------------- Section 3 -------------------------------------------
# --------------  GeosChem (bpch) prod loss variables
#

# --------------
# 3.01 - Spec to photolysis reaction p/l tag 
# -------------
def spec_phot_2_RD(spec):
    d = {
    'OIO': 'RD67', 'ICl': 'RD74', 'I2O2': 'RD70', 'I2': 'RD64', \
    'CH2ICl': 'RD88', 'HOI': 'RD65', 'CH2IBr': 'RD89', 'INO': 'RD75',  \
    'IO': 'RD66', 'CH2I2': 'RD72','CH3IT': 'RD71', 'IONO2': 'RD69',  \
    'IONO': 'RD68', 'IBr': 'RD73' \
    }
    return d[spec]

# -------------
# 3.02 - Get families for reactions 
# ------------- 
def get_tag_fam( tag ):
    """
        dictionary of manually constructed assignment list
            - Ox loss familes
            - addition for paranox Kludge 
            (in v9-2 (patched), but removed in v10? )

    """
    # Ox family dictionary
    fam_d = {
    'LO3_18': 'Photolysis', 'LR25': 'Bromine', 'LR21': 'Bromine', 'LO3_38': 'Iodine', 'LO3_63': 'NOy', 'LO3_10': 'HOx', 'LO3_34': 'Iodine', 'LO3_35': 'Iodine', 'LO3_30': 'Iodine', 'LR5': 'Bromine', 'LR6': 'Bromine', 'LO3_61': 'NOy', 'LO3_60': 'NOy', 'LO3_39': 'Iodine', 'LO3_05': 'HOx', 'LO3_07': 'NOy', 'LO3_06': 'HOx', 'LO3_49': 'NOy', 'LO3_62': 'NOy', 'LO3_03': 'HOx', 'LO3_02': 'HOx', 'LO3_67': 'NOy', 'LO3_66': 'NOy', 'LO3_69': 'NOy', 'LO3_42': 'NOy', 'LO3_41': 'NOy', 'LO3_40': 'NOy', 'LO3_47': 'HOx', 'LO3_46': 'NOy', 'LO3_09': 'HOx', 'LO3_44': 'NOy', 'LR37': 'HOx', 'LR36': 'NOy', 'LO3_65': 'NOy', 'LR30': 'Bromine', 'LO3_24': 'Iodine', 'LR10': 'Bromine', 'LR38': 'NOy', 'LO3_68': 'NOy', 'LO3_64': 'NOy', 'LO3_36': 'Iodine', 'LO3_57': 'NOy', 'LO3_72': 'NOy', 'RD98': 'Photolysis', 'LO3_71': 'NOy', 'LO3_58': 'NOy', 'LO3_54': 'Photolysis', 'LO3_55': 'Iodine', 'LO3_56': 'HOx', 'LO3_08': 'HOx', 'LO3_50': 'NOy', 'LO3_51': 'NOy', 'LO3_52': 'NOy', 'LO3_53': 'HOx'
    # added
    , 'RD63': 'Iodine', 'RD62':'Iodine',  'LO3_38': 'Iodine', 'RD59': 'Iodine', 'LO3_30' : 'Iodine', 'RD65': 'Iodine', 'LO3_34': 'Iodine',  'RD93': 'Iodine', 'LO3_55': 'Iodine', 'RD92': 'Iodine', 'LO3_39': 'Iodine' , 'LO3_36': 'Iodine','RD95': 'Iodine' , 'RD67': 'Iodine', 'LO3_35': 'Iodine' 
#    , 'RD36': 'Bromine' # Kludge to allow combination reactions 
    # Kludge - from Chris Holmes (paranox deposition, goes through p/l as Ox losss )
    ,'LO3_70' : 'Photolysis'

    # Extra tags not in list? 
    #  (these reactions are appearing due to lack of inclusion of iodine 
    # species in Ox family... )
#    ,'RD19': 'iodine', 'RD37': 'iodine', 'RD01': 'iodine' 
    }
    
    # Creigee reaction class/"family" dictionary
#    if cregiee:
#        fam_d ={
#        }
        
    return fam_d[tag]

# ----------------- Section 4 -------------------------------------------
# -------------- GeosChem (bpch) general variables
#
# 4.01 - v9-2 species in input.geos from num
# 4.02 - Get Species Mass 
# 4.03 - Get Species stioch
# 4.04 -  GEOS-Chem/ctm.bpch values (current main dict )
# 4.05 - latex species name
# 4.06 - converts P/L tracer mulitpler to 1
# 4.07 - Returns tracers unit and scale (if requested)
# 4.08 - Store of dirs for earth0, atmosviz1, and tms mac
# 4.09 - Ox in species (redundant now? should adapt species stoich )
# 4.10 - Get GAW site info  (lat, lon, alt (press), timezone (UTC) ) 

# 4.99 - Reference data, (inc. grid data) from gchem - credit: GK (Gerrit Kuhlmann )

# --------------   
# 4.01 - v9-2 species in input.geos from num
# -------------    
def num2spec( num=69, rtn_dict=False, invert=False, ver = '1.7' ):

    # get dictionary of tracer numbers
    d= what_species_am_i( ver=ver, rtn_dict=True, special_case=None )

    # slice off just numbers
    nums =[ int(i[4:]) for i in d.keys()]

    # re-make dictionary
    d = dict( zip(nums, d.values() ) )
    
    # inver to give spec for num
    if invert:
        d = { v: k for k, v in d.items() }

    if rtn_dict:
        return d
    else:
        return d[num]

# --------------
# 4.02 - RMM (Mass) (g /mol) for species 
# ------------- 
# C3H5I == C2H5I (this is a vestigle typo, left in to allow for use of older model runs
def species_mass(spec):
    d = {
    'HIO3': 176.0, 'OCPO': 12.0, 'Br2': 160.0, 'OCPI': 12.0, 'O3': 48.0, 'PAN': 121.0, 'ACET': 12.0, 'RIP': 118.0, 'BrNO3': 142.0, 'Br': 80.0, 'HBr': 81.0, 'HAC': 74.0, 'ALD2': 12.0, 'HNO3': 63.0, 'HNO2': 47.0, 'C2H5I': 168.0, 'HNO4': 79.0, 'OIO': 159.0, 'MAP': 76.0, 'PRPE': 12.0, 'CH2I2': 268.0, 'IONO2': 189.0, 'NIT': 62.0, 'CH3Br': 95.0, 'C3H7I': 170.0, 'C3H8': 12.0, 'DMS': 62.0, 'CH2O': 30.0, 'CH3IT': 142.0, 'NO2': 46.0, 'NO3': 62.0, 'N2O5': 105.0, 'H2O2': 34.0, 'DST4': 29.0, 'DST3': 29.0, 'DST2': 29.0, 'DST1': 29.0, 'MMN': 149.0, 'HOCl': 52.0, 'NITs': 62.0, 'RCHO': 58.0, 'C2H6': 12.0, 'MPN': 93.0, 'INO': 157.0, 'MP': 48.0, 'CH2Br2': 174.0, 'SALC': 31.4, 'NH3': 17.0, 'CH2ICl': 167.0, 'IEPOX': 118.0, 'ClO': 51.0, 'NO': 30.0, 'SALA': 31.4, 'MOBA': 114.0, 'R4N2': 119.0, 'BrCl': 115.0, 'OClO': 67.0, 'PMN': 147.0, 'CO': 28.0, 'BCPI': 12.0, 'ISOP': 12.0, 'BCPO': 12.0, 'MVK': 70.0, 'BrNO2': 126.0, 'IONO': 173.0, 'Cl2': 71.0, 'HOBr': 97.0, 'PROPNN': 109.0, 'Cl': 35.0, 'I2O2': 286.0, 'I2O3': 302.0, 'I2O4': 318.0, 'I2O5': 338.0, 'MEK': 12.0, 'HI': 128.0, 'ISOPN': 147.0, 'SO4s': 96.0, 'I2O': 270.0, 'ALK4': 12.0, 'MSA': 96.0, 'I2': 254.0, 'PPN': 135.0, 'IBr': 207.0, 'MACR': 70.0, 'I': 127.0, 'AERI': 127.0, 'HOI': 144.0, 'BrO': 96.0, 'NH4': 18.0, 'SO2': 64.0, 'SO4': 96.0, 'IO': 143.0, 'CHBr3': 253.0, 'CH2IBr': 221.0, 'ICl': 162.0, 'GLYC': 60.0
    # species, not in tracer list
    , 'HO2': 33.0, 'OH': 17.0,'CH4':16.0 , 'N':14.0, 'CH3I':142.0, 'CH2OO':46.0
    }
    
    return d[spec]

# --------------
# 4.03 -  return the stiochometry of Iodine in species
# --------------
def spec_stoich( spec, IO=False, I=False, NO=False, OH=False, N=False,
            C=False ): 
#    if I:  # note -  re-write to take stioch species (e.g. OH, I instead of booleans ) - asssume I == True as default
    #  C3H5I == C2H5I (this is a vestigle typo, left in to allow for use of older model runs    
    # aerosol cycling specs
    # 'LO3_36' : (2.0/3.0) , 'LO3_37' : (2.0/4.0),           # aersol loss rxns... 'LO3_37' isn't true loss, as I2O4 is regen. temp
    # aerosol loss rxns - corrected stochio for Ox, adjsutment need for I
    d = {
    'RD11': 1.0, 'RD10': 1.0, 'HIO3': 1.0, 'RD15': 1.0, 'RD62': 2.0, 'RD17': 1.0, 'RD16': 1.0, 'RD19': 1.0, 'LO3_37': 0.5, 'CH2I2': 2.0, 'AERII': 1.0, 'CH2ICl': 1.0, 'PIOx': 1.0, 'C3H7I': 1.0, 'RD73': 1.0, 'RD72': 2.0, 'RD71': 1.0, 'RD70': 1.0, 'C3H5I': 1.0, 'RD57': 1.0, 'CH3IT': 1.0, 'IO': 1.0, 'LO3_38': 1.0, 'RD61': 1.0, 'RD68': 1.0, 'I2': 2.0, 'IONO': 1.0, 'LO3_36': 0.6666666666666666, 'INO': 1.0, 'RD88': 1.0, 'RD89': 1.0, 'LOx': 1.0, 'RD06': 1.0, 'RD07': 1.0, 'RD02': 1.0, 'RD01': 1.0, 'I': 1.0, 'LO3_24': 0.5, 'AERI': 1.0, 'HOI': 1.0, 'RD64': 2.0, 'RD65': 1.0, 'RD66': 1.0, 'RD67': 1.0, 'RD60': 1.0, 'RD47': 1.0, 'C2H5I': 1.0, 'RD63': 1.0, 'RD20': 1.0, 'RD22': 1.0, 'RD24': 1.0, 'RD69': 1.0, 'RD27': 1.0, 'OIO': 1.0, 'CH2IBr': 1.0, 'LIOx': 1.0, 'L_Iy': 1.0, 'ICl': 1.0, 'IBr': 1.0, 'RD95': 2.0, 'I2O2': 2.0, 'I2O3': 2.0, 'I2O4': 2.0, 'I2O5': 2.0, 'HI': 1.0, 'I2O': 2.0, 'RD59': 1.0, 'RD93': 2.0, 'RD92': 1.0, 'IONO2': 1.0, 'RD58': 1.0,
    # p/l for: IO, I
    'RD15': 1.0, 'RD17': 1.0, 'RD75': 1.0, 'RD72': 2.0, 'RD71': 1.0, 'RD70': 1.0, 'RD56': 1.0, 'RD69': 1.0, 'RD88': 1.0, 'RD89': 1.0, 'RD06': 1.0, 'RD07': 1.0, 'RD08': 1.0, 'RD64': 2.0, 'RD65': 1.0, 'RD67': 1.0, 'RD46': 2.0, 'RD47': 1.0, 'RD20': 1.0, 'RD22': 1.0, 'RD68': 1.0, 'RD25': 1.0, 'RD96': 1.0 ,
    'RD11': 1.0, 'RD12': 2.0, 'RD02': 1.0, 'RD16': 1.0, 'RD19': 1.0, 'RD24': 1.0, 'RD09': 1.0, 'RD23': 1.0, 'RD37': 1.0, 'RD97': 1.0, 
    # kludge for test analysis (HEMCO emissions )
    'ACET' : 1.0, 'ISOP': 1.0, 'CH2Br2': 1.0, 'CHBr3':1.0, 'CH3Br':1.0

    }

    if IO:
        d = {
        'RD11': 2.0, 'RD10': 1.0, 'RD12': 2.0, 'LO3_36': 1./3., 'RD09': 1.0, 'RD66': 1.0, 'RD23': 1.0, 'RD37': 1.0, 'LO3_24': 1.0/2.0, 'RD56': 1.0, 'RD01': 1.0, 'RD08': 1.0, 'RD46': 2.0, 'RD30': 1.0, 'RD25': 1.0, 'RD27': 1.0, 'RD97':1.0
        }
    if NO:
        d ={
        'NO2': 1.0, 'NO3': 1.0, 'N2O5': 2.0, 'NO': 1.0, 'PPN': 1.0, 'R4N2': 1.0, 'BrNO3': 1.0, 'INO': 1.0, 'PAN': 1.0, 'PMN': 1.0, 'HNO3': 1.0, 'HNO2': 1.0, 'NH3': 1.0, 'HNO4': 1.0, 'BrNO2': 1.0, 'IONO': 1.0, 'PROPNN': 1.0, 'NH4': 1.0, 'MPN': 1.0, 'MMN': 1.0, 'ISOPN': 1.0, 'IONO2': 1.0
        }
    if OH:
        d = {
        'LO3_18': 2.0, 'LO3_03': 1.0, 'RD95': 1.0, 'PO3_14': 1.0
        }
    if N:
        d= {
        'RD10': 1.0, 'LR26': 1.0, 'LR27': 1.0, 'LR20': 1.0, 'RD17': 1.0, 'RD16': 1.0, 'RD19': 1.0, 'RD18': 2.0, 'LR28': 1.0, 'LO3_30': 1.0, 'RD75': 1.0, 'LR7': 1.0, 'LR8': 1.0, 'RD56': 1.0, 'RD24': 1.0, 'LO3_39': 1.0, 'RD25': 1.0, 'RD81': 1.0, 'LR35': 1.0, 'LR18': 1.0, 'LR17': 1.0, 'LR11': 1.0, 'LR39': 1.0, 'RD20': 1.0, 'RD21': 2.0, 'RD22': 1.0, 'RD23': 1.0, 'RD68': 1.0, 'RD69': 1.0
    # NOy ( N in 'NOy')
, 'NO2': 1.0, 'NO3': 1.0, 'N2O5': 2.0, 'NO': 1.0, 'PPN': 1.0, 'R4N2': 2.0, 'BrNO3': 1.0, 'INO': 1.0, 'PAN': 1.0, 'PMN': 1.0, 'HNO3': 1.0, 'HNO2': 1.0, 'NH3': 1.0, 'HNO4': 1.0, 'BrNO2': 1.0, 'IONO': 1.0, 'PROPNN': 1.0, 'NH4': 1.0, 'MPN': 1.0, 'MMN': 1.0, 'ISOPN': 1.0, 'IONO2': 1.0
        }
    if C:
        d = {
    'ACET': 3.0, 'ALD2': 2.0, 'C2H6': 2.0, 'C3H8': 3.0, 'ISOP': 5.0
        }

    return d[spec]

# --------------
# 4.04 -  GEOS-Chem/ctm.bpch values
# --------------
def GC_var(input_x=None, rtn_dict=False, debug=False):
    """
        Note: Most of this dictionary is vestigial. 
        <= ACTION NEEDED  ( remove redundant variables )
    f_var     = GC flux (EW, NS , UP) variables
    Ox        = 'Ox', 'POX', 'LOX' + list of drydep species # not inc. 'NO3df', 'HNO4df', 'BrOdf' , 'BrNO2', 'IO', 'IONO', 'OIO', 
    Ox_p      = Ox prod list
    Ox-l      = Ox loss list
    d_dep  = dry dep (category, name = species)
    w_dep  = wet dep ( 3x categories (WETDCV = rain out loss in convective updrafts (kg/s), WETDLS = rainout in large scale precip (kg/s), CV-FLX = Mass change due to cloud convection (kg/s); name = species)
    BL_m    = UPWARD MASS FLUX FROM BOUNDARY-LAYER MIXING,  (category, name = species) 
    f_strat  = strat flux (to tropsosphere) (category, name = species) 
    """
    if (debug):
        print 'GC_var called'
    GC_var_dict = { 
                    # Ox budget analysis
                    'f_var'      : ['EW-FLX-$', 'NS-FLX-$', 'UP-FLX-$' ], 
#                    'r_t'         : [ 'Photolysis','HOx','NOy'    ,'Bromine', 'Iodine' ], 
#                    'r_tn'       : ['NOy'  ,'Photolysis','HOx' ,'Bromine', 'Iodine' ],                     
                    'r_t'         : [ 'Photolysis','HOx','Bromine', 'Iodine' ], 
                    'r_tn'       : ['Photolysis','HOx' ,'Bromine', 'Iodine' ],                     
                    'r_tn_lc'  : ['photolysis','HOx' ,'bromine', 'iodine' ],                     

#                    'Ox'         : ['Ox','POX','LOX','O3df','NO2df', 'PANdf', 'PMNdf', 'PPNdf', 'N2O5df','HNO3df', 'HOBrdf','BrNO3df','HOIdf','IONO2df', 'I2O2df', 'I2O4df','I2O3df'], 
#                    'Ox1.1'         : ['Ox','POX','LOX','O3df','NO2df', 'PANdf', 'PMNdf', 'PPNdf', 'N2O5df','HNO3df', 'HOBrdf','BrNO3df','HOIdf','IONO2df', 'I2O2df'], 
#                    'Ox_spec'    : ['O3', 'NO2', 'NO3', 'PAN', 'PMN', 'PPN', 'HNO4', 'N2O5', 'HNO3', 'BrO', 'HOBr', 'BrNO2', 'BrNO3', 'MPN', 'IO', 'HOI', 'IONO', 'IONO2', 'OIO', 'I2O2', 'I2O4', 'I2O3'],
#                    'Ox_spec1.1'    : ['O3', 'NO2', 'NO3', 'PAN', 'PMN', 'PPN', 'HNO4', 'N2O5', 'HNO3', 'BrO', 'HOBr', 'BrNO2', 'BrNO3', 'MPN', 'IO', 'HOI', 'IONO', 'IONO2', 'OIO', 'I2O2', 'I2O4'],                    
#                    'Ox_p_1.3'  : ['PO3_85', 'PO3_76', 'PO3_62', 'PO3_63', 'RD06', 'PO3_01', 'PO3_77', 'PO3_86', 'PO3_79', 'PO3_14', 'PO3_87', 'PO3_66', 'PO3_15', 'PO3_67', 'PO3_51', 'PO3_88', 'PO3_56', 'PO3_89', 'PO3_90', 'PO3_72', 'PO3_91', 'PO3_02', 'PO3_03', 'PO3_92', 'PO3_64', 'PO3_68', 'LR9', 'PO3_58', 'PO3_57', 'PO3_70', 'PO3_60', 'PO3_65', 'PO3_18', 'PO3_73', 'PO3_19', 'PO3_20', 'PO3_21', 'PO3_22', 'PO3_24', 'PO3_25', 'PO3_26', 'PO3_27', 'PO3_34', 'PO3_35', 'PO3_37', 'PO3_38', 'PO3_39', 'PO3_05', 'PO3_45', 'PO3_69', 'PO3_46', 'PO3_97', 'PO3_47', 'PO3_48', 'PO3_52', 'PO3_53', 'PO3_80', 'PO3_40', 'PO3_54', 'PO3_81', 'PO3_55', 'PO3_74', 'PO3_41', 'PO3_43', 'PO3_93', 'PO3_94', 'PO3_95', 'PO3_96', 'PO3_84', 'PO3_61', 'PO3_83', 'PO3_59', 'PO3_71', 'PO3_98', 'PO3_99', 'PO3100', 'PO3101', 'PO3_75',  'PO3_50'], # 'ISOPND',
#                    'Ox_l_fp'    : ['LO3_18','RD98','LO3_54', # hv
#                                    'LO3_03', 'LO3_02', 'LO3_08', 'LO3_09', 'LO3_06', 'LO3_05', 'LO3_10', 'LO3_53', 'LO3_47', # HOx (ROx route) # HOx
#                                    'LO3_52', 'LO3_51','LO3_50', 'LO3_49', 'LO3_46', 'LO3_44', 'LO3_42', 'LO3_41', 'LO3_40',#'LO3_48','LO3_04', # NOy route
#                                    'LR25','LR21','LR5','LR6' ,'LR30',  'LR10', 'LO3_24', # Bry route
#                                    'LO3_34','LO3_24','LO3_35','LO3_30', 'LO3_39', 'LO3_38','LO3_55' , 'LO3_36', 'RD63'], # Iy route
#                    'Ox_l_fp1.1'    : ['LO3_18','LO3_54', # hv
#                                    'LO3_03', 'LO3_02', 'LO3_08', 'LO3_09', 'LO3_06', 'LO3_05', 'LO3_10', 'LO3_53', 'LO3_47', # HOx (ROx route) # HOx
#                                    'LO3_52', 'LO3_51','LO3_50', 'LO3_49', 'LO3_46', 'LO3_44', 'LO3_42', 'LO3_41', 'LO3_40',#'LO3_48','LO3_04', # NOy route
#                                    'LR25','LR21','LR5','LR6' ,'LR30', 'LR10','LO3_24', # Bry route
#                                    'LO3_34','LO3_24','LO3_35','LO3_30','LO3_39',  'LO3_38','LO3_55' , 'LO3_36', 'RD63'], # Iy route
#                    'Ox_l_fp1.3'    : ['LO3_18','RD98','LO3_54', # hv
#                                    'LO3_03', 'LO3_02', 'LO3_08', 'LO3_09', 'LO3_06', 'LO3_05', 'LO3_10', 'LO3_53', 'LO3_47', 'LO3_56','LR37', 'LR38',  # HOx (ROx route) # HOx
#                                    'LO3_52', 'LO3_51','LO3_50', 'LO3_49', 'LO3_46', 'LO3_44', 'LO3_42', 'LO3_41', 'LO3_40','LR36', 'LO3_58', #'LO3_48','LO3_04',,# NOy route
#                                    'LO3_57', 'LO3_07', 'LO3_60', 'LO3_61', 'LO3_62', 'LO3_63', 'LO3_64', 'LO3_65', 'LO3_66', 'LO3_67', 'LO3_68', 'LO3_69', 'LO3_72','LO3_71', # NOy route
#                                    'LR25','LR21','LR5','LR6' ,'LR30',  'LR10', 'LO3_24', # Bry route
#                                    'LO3_34','LO3_24','LO3_35','LO3_30', 'LO3_39', 'LO3_38','LO3_55' , 'LO3_36'],  # Iy route

#                    'Ox_l_fp_r_' :  [(0, 3), (3, 11), (11, 21), (21, 28), (28, None)], #, (-2, -1)]
#                    'Ox_l_fp_r_1.1' :  [(0, 2), (2, 11), (11, 21), (21, 28), (28, None)], #, (-2, -1)]
#                    'Ox_l_fp_r_1.3' :  [(0, 3), (3, 14), (14, 40), (40, 47), (47, None)], #, (-2, -1)]                    
                    'fams'    :  ['I2','HOI','IO', 'I', 'HI+OIO+IONO+INO', 'IONO2','IxOy', 'CH3I', 'CH2IX'],   # Iy families
                    'fams_A'    :  ['I2','HOI','IO', 'I', 'HI+OIO+IONO+INO', 'IONO2','IxOy', 'CH3I', 'CH2IX', 'AERI'],   # Iy families            
                    'fam_slice' : [(0, 1), (1, 2), (2, 3), (3,4 ),(4, 8), (8, 9), (9, 12), (12, 13), (13, None)],   # slice 
                    'fam_slice_A' : [(0, 1), (1, 2), (2, 3), (3,4 ),(4, 8), (8, 9), (9, 12), (12, 13), (13, 16),(16, None)],   # slice 
#                    'POx_l_fp'   : ['PO3_01', 'PO3_03', 'PO3_02', 'PO3_05','PO3_14', 'PO3_15', 'PO3_18', 'PO3_19', 'PO3_20', 'PO3_21', 'PO3_22', 'PO3_24', 'PO3_25', 'PO3_26', 'PO3_27', 'PO3_30', 'PO3_31', 'PO3_32', 'PO3_33', 'PO3_34', 'PO3_35', 'PO3_37', 'PO3_38', 'PO3_39', 'PO3_40', 'PO3_41', 'PO3_43'],
                    'Ox_key'     : ['POX', 'PO3_14', 'PO3_15',   'LOX'],#, 'LO3_18', 'LO3_03', 'LO3_02','LR25', 'LR21', 'LR5','LR6','LO3_34', 'LO3_33','LO3_24', 'LO3_35'  ],
                    'POxLOx'     : ['POX', 'LOX'],
                    'iPOxiLOx'   : ['POX', 'LOX', 'iPOX', 'iLOX'],

                    # Iy/ Iodine budget analysis
                    'BL_FT_UT'   : [(0, 6), (6, 26), (26, 38)] ,            
                    'n_order'  :['CH2IX','CH3I', 'I2', 'HOI','IO', 'I', 'IONO2','HI+OIO+IONO+INO','IxOy' ] ,
                    'n_order_A'  :['CH2IX','CH3I', 'I2', 'HOI','IO', 'I', 'IONO2','HI+OIO+IONO+INO','IxOy', 'AERI' ] ,
                    'I_l'        : ['RD01', 'RD02', 'RD16', 'RD19', 'RD24', 'RD27'],
                    'IO_l'       : ['RD09', 'RD10', 'RD11', 'RD12', 'RD23', 'LO3_24', 'RD37', 'RD97', 'RD66'], # LO37 swaped for RD97 as LO37 assigned to loss point of I2O3 uptake
                    'I_p'        : ['RD06','RD07','RD10','RD11','RD47','RD15','RD17','RD20','RD22', 'LO3_24', 'RD64', 'RD65', 'RD66', 'RD67','RD68','RD69', 'RD70', 'RD71','RD72', 'RD73', 'RD88', 'RD89'],
                    'IO_p'       : [ 'RD01', 'RD08', 'RD46', 'RD25', 'RD27','RD56'],
                     'sOH'            :  ['LO3_18'],
                    'd_dep'      : ['DRYD-FLX'],
                    'w_dep'      : ['WETDCV-$','WETDLS-$'], 
                    'BL_m'       : ['TURBMC-$'],
                    'f_strat'    : ['STRT-FL'],
                    'p_l'        : ['PORL-L=$'],
                    'Cld_flx'    : ['CV-FLX-$'],
                    'I_Br_O3'    : ['IO', 'OIO','HOI','I2','I','CH3IT','CH2I2','CH2ICl', 'CH2IBr', 'C3H7I','C2H5I', 'BrO', 'Br', 'HOBr','Br2','CH3Br', 'CH2Br2', 'CHBr3', 'O3', 'CO'],
                    'IOrg_RIS'   : ['CH3IT','CH2ICl','CH2I2', 'CH2IBr', 'I2','HOI','I','IO', 'OIO', 'HI','IONO','IONO2'],
                    'I_specs'    : ['I2','HOI','IO', 'OIO', 'HI','IONO', 'IONO2','I2O2', 'I2O3','I2O4''CH3IT','CH2I2','I','INO'] ,
                    'Iy'         : ['I2','HOI','IO', 'OIO', 'HI','INO','IONO', 'IONO2','I2O2', 'I2O3','I2O4','I'],
                    'Iy1.1'      : ['I2','HOI','IO', 'OIO', 'HI','IONO', 'IONO2','I2O2','I2O4','I','INO'],
                    'IOy'        : ['HOI','IO', 'OIO','IONO','IONO2','INO','I2O2','I2O4', 'I2O3'],
                    'IOy1.1'     : ['HOI','IO', 'OIO','IONO','IONO2','INO','I2O2','I2O4'],
                    'I2Ox'       : ['I2O2','I2O4','I2O3'],
                    'I2Ox'       : ['I2O2','I2O4','I2O3'],
                    'IyOx1.1'       : ['I2O2','I2O4'],
                    'Iy_no_i2o4' : ['I2','HOI','IO', 'OIO', 'HI','IONO', 'IONO2','I2O2','I','INO', 'I2O3'],
                    'Iy_no_i2o41.1' : ['I2','HOI','IO', 'OIO', 'HI','IONO', 'IONO2','I2O2','I','INO'],
                    'Phot_s_Iy'  : ['CH3IT','CH2ICl','CH2I2', 'CH2IBr'],#['RD89', 'RD88', 'RD71', 'RD72'],
                    'HOI'        : ['HOI'],
                    'IOx'        : ['IO','I',],
                    'IO'         : ['IO'],
                    'I'          : ['I',],
                    'OIO'         : ['OIO'],
                    'LIOx'       : ['LIOx'],  # LOx is p/l tracer name, for Loss of IOx
                    'PIOx'       : ['PIOx'],  # LOx is p/l tracer name, for Loss of IOx
                    'iodine_all'  : ['I2','HOI','IO', 'I', 'HI', 'OIO', 'INO',  'IONO','IONO2','I2O2', 'I2O4', 'I2O3', 'I2O5', 'CH3IT', 'CH2I2', 'CH2ICl', 'CH2IBr', 'C3H7I','C2H5I','ICl', 'I2O', 'IBr', 'HIO3', ],
                    'iodine_all_A': ['I2','HOI','IO', 'I', 'HI', 'OIO', 'INO',  'IONO','IONO2','I2O2', 'I2O4', 'I2O3', 'I2O5', 'CH3IT', 'CH2I2', 'CH2ICl', 'CH2IBr', 'C3H7I','C2H5I','ICl', 'I2O', 'IBr', 'HIO3','AERI' ],

                    # Misc analysis
                    'LHOI'          : ['RD65', 'RD63', 'RD08'],
                     'LHOBr'     : ['LR25', 'LR30','LR21'],
                    'LI2'           : ['RD64', 'RD06', 'RD22'],
                    'LCH2I2'     : ['RD72' ],
                    'LCH2Cl'     : ['RD88' ] ,                       
                    'LCH2Br'     : ['RD89' ],                        
                    'LCH3IT'     : ['RD15' , 'RD71'],                        
                     'sHOX'       : ['HOI', 'HOBr'],
                    'HO2_loss'   : ['PO3_14','RD09','RD02', 'LR2', 'LR3', 'PO3_46','PO3_02', 'PO3_03', 'PO3_05'],
                    'CAST_int'   : ['IO', 'OIO','HOI','I2','I','HOI','CH3I','CH2I2','CH2ICl', 'CH2IBr', 'C3H7I','C3H5I', 'BrO', 'Br', 'HOBr','Br2','CH3Br', 'CH2Br2', 'CHBr3', 'O3', 'CO', 'OH', 'HO2','NO','NO2'],
                    'CAST_intn'  : ['IO', 'OIO','HOI','I2','I','HOI','CH3IT','CH2I2','CH2ICl', 'CH2IBr', 'C3H7I','C2H5I', 'BrO', 'Br', 'HOBr','Br2','CH3Br', 'CH2Br2', 'CHBr3', 'O3', 'CO', 'DMS', 'NO','HNO3','HNO4', 'NO2','NO3' , 'PAN' , 'HNO2', 'N2O5'],
                    'CAST_int_n' : ['IO', 'OIO','HOI','I2','I','HOI','CH3I','CH2I2','CH2ICl', 'CH2IBr', 'C3H7I','C2H5I', 'BrO', 'Br', 'HOBr','Br2','CH3Br', 'CH2Br2', 'CHBr3', 'O3', 'CO', 'OH', 'HO2','NO','NO2'],
                    'diurnal_sp' : ['IO','I2', 'CH2I2', 'BrO' ]  ,
                    'obs_comp'   : ['CH3IT','CH2I2','CH2ICl','CH2IBr','C2H5I','C3H7I','I2','IO'] ,                    
                    'emiss_specs': ['CH3IT','CH2I2','CH2ICl','CH2IBr','I2','HOI'] ,
                    'w_dep_specs': ['I2' ,'HI'  ,'HOI'  ,'IONO', 'IONO2','I2O2', 'I2O4', 'I2O3' ,'AERI'],  #, 'IBr', 'ICl']
  #                  'd_dep_specsl' : ['I2', 'HI', 'HOI', 'IONO', 'IONO2',  'I2O2', 'I2O4', 'I2O3','AERI'], #, 'IO', 'OIO'] ,                      
                    'd_dep_specsl1.1' : ['I2', 'HI', 'HOI', 'IONO', 'IONO2',  'I2O2', 'I2O4', 'AERI'], #, 'IO', 'OIO'] ,                      
                    'd_dep_specs': ['I2df', 'HIdf', 'HOIdf', 'IONOdf', 'IONO2df',  'I2O2df', 'I2O4df', 'I2O3df', 'AERIdf',], #, 'IOdf', 'OIOdf'], #
#                    'd_dep_specs1.1': ['I2df', 'HIdf', 'HOIdf', 'IONOdf', 'IONO2df',  'I2O2df', 'I2O4df','AERIdf',], #, 'IOdf', 'OIOdf'], #
                    'I2_het_cyc'  : ['RD59','RD92','RD63'],  # IONO2, IONO, HOI
                    'I_het_loss'  : [ 'RD58', 'RD62', 'RD93' ,'RD95'], # HI, I2O2, I2O4, I2O3 uptake (prev: 2OIO excuded as I2Ox formaed, IO+OIO included as I2O3 not treated )
 #['RD60','RD61','RD62','RD52','RD53','RD54','RD55','RD13'],  # RD13 = OIO + OH => HIO3  86 => AERI loss
                    'NOx'          : ['NO', 'NO2' ],
                    'N_specs'  :  ['NO', 'NO2', 'PAN', 'HNO3', 'PMN', 'PPN', 'R4N2', 'N2O5', 'HNO4', 'NH3', 'NH4', 'BrNO2', 'BrNO3', 'MPN', 'ISOPN', 'PROPNN', 'MMN', 'NO3', 'HNO2', 'IONO', 'IONO2', 'INO'],
                    'N_specs_no_I'  :  ['NO', 'NO2', 'PAN', 'HNO3', 'PMN', 'PPN', 'R4N2', 'N2O5', 'HNO4', 'NH3', 'NH4', 'BrNO2', 'BrNO3', 'MPN', 'ISOPN', 'PROPNN', 'MMN', 'NO3', 'HNO2'],
                    'I_N_tags' : ['RD10', 'RD23', 'RD19', 'RD16', 'RD22', 'RD56', 'RD24', 'LO3_30', 'RD69', 'RD68', 'RD20', 'RD21', 'RD25', 'LO3_39', 'RD17', 'RD18', 'RD75'],
                    'Br_N_tags' : ['LR7', 'LR18', 'LR17', 'LR11', 'LR8', 'LR20', 'LR26', 'LR28', 'LR27'],
                    'inactive_I'  : ['BrCl', 'OClO', 'ClO', 'HOCl', 'Cl', 'Cl2', 'I2O5', 'I2O', 'HIO3', 'IBr', 'ICl', 'C2H5I','C3H7I'], # I2O3 now active.
                    'active_I' : ['I2', 'HOI', 'IO', 'I', 'HI', 'OIO', 'INO', 'IONO', 'IONO2', 'I2O2', 'I2O4', 'I2O3', 'CH3IT', 'CH2I2', 'CH2ICl', 'CH2IBr'], 
                    'surface_specs' : ['O3', 'NO', 'NO2', 'NO3' ,'N2O5', 'IO', 'IONO2' ],
                
                    # Model run title dictionaries
                    'run_name_dict': {'run': 'Halogens (I+,Br+)', 'Br_2ppt': 'Halogens (I+,Br+) + fixed 2 pptv BrO', 'just_I': 'Just Iodine (I+,Br-)', 'no_hal': 'No Halogens', 'just_Br': 'Just Bromine (I-,Br+)', 'Br_1ppt': 'Halogens (I+,Br+) + fixed 1 pptv BrO', 'obs': 'Observations'}   ,
                    'latex_run_names': {'I2Ox_half': 'I$_{2}$Ox loss ($\\gamma$) /2', 'run': 'Iodine simulation.', 'MacDonald_iodide': 'Ocean iodide', 'Sulfate_up': 'Sulfate Uptake', 'I2Ox_phot_exp': 'I$_{2}$Ox exp. X-sections',  'het_double': 'het. cycle ($\\gamma$) x2', 'I2Ox_phot_x2': 'I$_{2}$Ox X-sections x2', 'no_het': 'no het. cycle ', 'I2Ox_double': 'I$_{2}$Ox loss ($\\gamma$) x2', 'just_I': '(I+,Br-)', 'BrO1pptv': 'MBL BrO 1 pptv', 'het_half': 'het. cycle ($\\gamma$) /2', 'Just_I_org': 'Just org. I', 'no_I2Ox': 'No I$_{2}$Ox Photolysis', 'BrO1pptv_ALL' : 'BrO 1 pptv in Trop.', 'BrO2pptv' : 'MBL BrO 2 pptv',
                    # adjust from GBC to ACP names
#                    'no_hal': '(I-,Br-)', 'Just_Br': '(I-,Br+)', 
                    'no_hal': 'No Halogens', 'Just_Br': 'GEOS-Chem (v9-2)', 
                   # kludge for diurnal plot
                   'Iodine simulation.':'Iodine simulation.', '(I+,Br+)': 'Iodine simulation.','(I+,Br-)': 'Just Iodine', '(I-,Br+)': 'GEOS-Chem (v9-2)', '(I-,Br-)': 'No Halogens'},
                    
                    # tracer unit handling
                    'spec_2_pptv' : ['I2', 'HOI', 'IO', 'OIO', 'HI', 'IONO', 'IONO2', 'I2O2', 'CH3IT', 'CH2I2', 'IBr', 'ICl', 'I', 'HIO3', 'I2O', 'INO', 'I2O3', 'I2O4', 'I2O5', 'AERI', 'Cl2', 'Cl', 'HOCl', 'ClO', 'OClO', 'BrCl', 'CH2ICl', 'CH2IBr', 'C3H7I', 'C2H5I', 'Br2', 'Br', 'BrO', 'HOBr', 'HBr', 'BrNO2', 'BrNO3', 'CHBr3', 'CH2Br2', 'CH3Br','RCHO', 'MVK', 'MACR', 'PMN', 'PPN', 'R4N2', 'DMS', 'SO4s', 'MSA', 'NITs', 'BCPO', 'DST4', 'ISOPN', 'MOBA', 'PROPNN', 'HAC', 'GLYC', 'MMN', 'RIP', 'IEPOX', 'MAP' ,'N2O5','NO3'], # 'HNO4',  'HNO2'],
                    'spec_2_pptC' : ['PRPE', 'ISOP'],
                    # global 
                    'spec_2_ppbv': ['NO','DMS',  'RIP', 'IEPOX','BCPO', 'DST4', 'HAC', 'GLYC','MACR', 'ISOP'],
                    'spec_2_ppbC' : ['ALK4'],
                    # pf dictionaries
                    # WARNING - remove non forwards combatible dicts: 
                    # (GCFP_TRA_d ... GCFP_d2TRA  ... GCFP_d2TRA_justTRA . etc)
# GCFP_TRA_d is in use by CVO plotters -  use what_species_am_i instead!
                    'GCFP_TRA_d' : {'TRA_17': 'R4N2', 'TRA_16': 'PPN', 'TRA_15': 'PMN', 'TRA_14': 'MACR', 'TRA_13': 'MVK', 'TRA_12': 'RCHO', 'TRA_11': 'ALD2', 'TRA_19': 'C3H8', 'TRA_18': 'PRPE', 'TRA_96': 'C2H5I', 'TRA_95': 'C3H7I', 'TRA_94': 'CH2IBr', 'TRA_93': 'CH2ICl', 'TRA_92': 'BrCl', 'TRA_91': 'OClO', 'TRA_90': 'ClO', 'TRA_62': 'IEPOX', 'TRA_63': 'MAP', 'TRA_60': 'MMN', 'TRA_61': 'RIP', 'TRA_66': 'HNO2', 'TRA_67': 'I2', 'TRA_64': 'NO2', 'TRA_65': 'NO3', 'TRA_68': 'HOI', 'TRA_69': 'IO', 'TRA_71': 'HI', 'TRA_70': 'OIO', 'TRA_73': 'IONO2', 'TRA_72': 'IONO', 'TRA_75': 'CH3IT', 'TRA_74': 'I2O2', 'TRA_77': 'IBr', 'TRA_76': 'CH2I2', 'TRA_79': 'I', 'TRA_78': 'ICl', 'TRA_48': 'HBr', 'TRA_49': 'BrNO2', 'TRA_44': 'Br2', 'TRA_45': 'Br', 'TRA_46': 'BrO', 'TRA_47': 'HOBr', 'TRA_40': 'DST3', 'TRA_41': 'DST4', 'TRA_42': 'SALA', 'TRA_43': 'SALC', 'TRA_59': 'GLYC', 'TRA_58': 'HAC', 'TRA_53': 'CH3Br', 'TRA_52': 'CH2Br2', 'TRA_51': 'CHBr3', 'TRA_50': 'BrNO3', 'TRA_57': 'PROPNN', 'TRA_56': 'MOBA', 'TRA_55': 'ISOPN', 'TRA_54': 'MPN', 'TRA_28': 'SO4s', 'TRA_29': 'MSA', 'TRA_26': 'SO2', 'TRA_27': 'SO4', 'TRA_24': 'MP', 'TRA_25': 'DMS', 'TRA_22': 'N2O5', 'TRA_23': 'HNO4', 'TRA_20': 'CH2O', 'TRA_21': 'C2H6', 'TRA_39': 'DST2', 'TRA_38': 'DST1', 'TRA_35': 'OCPI', 'TRA_34': 'BCPI', 'TRA_37': 'OCPO', 'TRA_36': 'BCPO', 'TRA_31': 'NH4', 'TRA_30': 'NH3', 'TRA_33': 'NITs', 'TRA_32': 'NIT', 'TRA_08': 'H2O2', 'TRA_09': 'ACET', 'TRA_01': 'NO', 'TRA_02': 'O3', 'TRA_03': 'PAN', 'TRA_04': 'CO', 'TRA_05': 'ALK4', 'TRA_06': 'ISOP', 'TRA_07': 'HNO3', 'TRA_80': 'HIO3', 'TRA_81': 'I2O', 'TRA_82': 'INO', 'TRA_83': 'I2O3', 'TRA_84': 'I2O4', 'TRA_85': 'I2O5', 'TRA_86': 'AERI', 'TRA_87': 'Cl2', 'TRA_88': 'Cl', 'TRA_89': 'HOCl','O3':'O3', 'CO':'CO'} ,
# GCFP_d2TRA is in use by IO plotters -  use what_species_am_i instead!
                    'GCFP_d2TRA' : {'HIO3': 'TRA_80', 'OCPO': 'TRA_37', 'PPN': 'TRA_16', 'OCPI': 'TRA_35', 'O3': 'TRA_2', 'PAN': 'TRA_3', 'ACET': 'TRA_9', 'IEPOX': 'TRA_62', 'BrNO3': 'TRA_50', 'Br': 'TRA_45', 'HBr': 'TRA_48', 'HAC': 'TRA_58', 'ALD2': 'TRA_11', 'HNO3': 'TRA_7', 'HNO2': 'TRA_66', 'C2H5I': 'TRA_96', 'HNO4': 'TRA_23', 'OIO': 'TRA_70', 'MAP': 'TRA_63', 'PRPE': 'TRA_18', 'HI': 'TRA_71', 'CH2I2': 'TRA_76', 'IONO2': 'TRA_73', 'NIT': 'TRA_32', 'CH3Br': 'TRA_53', 'C3H7I': 'TRA_95', 'C3H8': 'TRA_19', 'DMS': 'TRA_25', 'CH2O': 'TRA_20', 'CH3IT': 'TRA_75','CH3I': 'TRA_75', 'NO2': 'TRA_64', 'NO3': 'TRA_65', 'N2O5': 'TRA_22', 'CHBr3': 'TRA_51', 'DST4': 'TRA_41', 'DST3': 'TRA_40', 'DST2': 'TRA_39', 'DST1': 'TRA_38', 'HOCl': 'TRA_89', 'NITs': 'TRA_33', 'RCHO': 'TRA_12', 'C2H6': 'TRA_21', 'MPN': 'TRA_54', 'INO': 'TRA_82', 'MP': 'TRA_24', 'CH2Br2': 'TRA_52', 'SALC': 'TRA_43', 'NH3': 'TRA_30', 'CH2ICl': 'TRA_93', 'RIP': 'TRA_61', 'ClO': 'TRA_90', 'NO': 'TRA_1', 'SALA': 'TRA_42', 'MOBA': 'TRA_56', 'R4N2': 'TRA_17', 'BrCl': 'TRA_92', 'OClO': 'TRA_91', 'PMN': 'TRA_15', 'CO': 'TRA_4', 'CH2IBr': 'TRA_94', 'ISOP': 'TRA_6', 'BCPO': 'TRA_36', 'MVK': 'TRA_13', 'BrNO2': 'TRA_49', 'IONO': 'TRA_72', 'Cl2': 'TRA_87', 'HOBr': 'TRA_47', 'PROPNN': 'TRA_57', 'Cl': 'TRA_88', 'I2O2': 'TRA_74', 'I2O3': 'TRA_83', 'I2O4': 'TRA_84', 'I2O5': 'TRA_85', 'MEK': 'TRA_10', 'MMN': 'TRA_60', 'ISOPN': 'TRA_55', 'SO4s': 'TRA_28', 'I2O': 'TRA_81', 'ALK4': 'TRA_5', 'MSA': 'TRA_29', 'I2': 'TRA_67', 'Br2': 'TRA_44', 'IBr': 'TRA_77', 'MACR': 'TRA_14', 'I': 'TRA_79', 'AERI': 'TRA_86', 'HOI': 'TRA_68', 'BrO': 'TRA_46', 'NH4': 'TRA_31', 'SO2': 'TRA_26', 'SO4': 'TRA_27', 'IO': 'TRA_69', 'H2O2': 'TRA_8', 'BCPI': 'TRA_34', 'ICl': 'TRA_78', 'GLYC': 'TRA_59','ALK4': 'ALK4', 'MSA': 'MSA', 'MO2': 'MO2', 'C3H8': 'C3H8', 'ISOP': 'ISOP', 'DMS': 'DMS', 'CH2O': 'CH2O', 'O3': 'O3', 'PAN': 'PAN', 'NO3': 'NO3', 'N2O5': 'N2O5', 'H2O2': 'H2O2', 'NO': 'NO', 'PPN': 'PPN', 'R4N2': 'R4N2', 'HO2': 'HO2', 'NO2': 'NO2', 'PMN': 'PMN', 'ACET': 'ACET', 'CO': 'CO', 'ALD2': 'ALD2', 'RCHO': 'RCHO', 'HNO3': 'HNO3', 'HNO2': 'HNO2', 'SO2': 'SO2', 'SO4': 'SO4', 'HNO4': 'HNO4', 'C2H6': 'C2H6', 'RO2': 'RO2', 'MVK': 'MVK', 'PRPE': 'PRPE', 'OH': 'OH', 'ETO2': 'ETO2', 'MEK': 'MEK', 'MP': 'MP' , 'GMAO_TEMP':'GMAO_TEMP' },
#                    'GCFP_d2TRA_justTRA' : {'HIO3': 'TRA_80', 'OCPO': 'TRA_37', 'PPN': 'TRA_16', 'OCPI': 'TRA_35', 'O3': 'TRA_2', 'PAN': 'TRA_3', 'ACET': 'TRA_9', 'IEPOX': 'TRA_62', 'BrNO3': 'TRA_50', 'Br': 'TRA_45', 'HBr': 'TRA_48', 'HAC': 'TRA_58', 'ALD2': 'TRA_11', 'HNO3': 'TRA_7', 'HNO2': 'TRA_66', 'C2H5I': 'TRA_96', 'HNO4': 'TRA_23', 'OIO': 'TRA_70', 'MAP': 'TRA_63', 'PRPE': 'TRA_18', 'HI': 'TRA_71', 'CH2I2': 'TRA_76', 'IONO2': 'TRA_73', 'NIT': 'TRA_32', 'CH3Br': 'TRA_53', 'C3H7I': 'TRA_95', 'C3H8': 'TRA_19', 'DMS': 'TRA_25', 'CH2O': 'TRA_20', 'CH3IT': 'TRA_75','CH3I': 'TRA_75', 'NO2': 'TRA_64', 'NO3': 'TRA_65', 'N2O5': 'TRA_22', 'CHBr3': 'TRA_51', 'DST4': 'TRA_41', 'DST3': 'TRA_40', 'DST2': 'TRA_39', 'DST1': 'TRA_38', 'HOCl': 'TRA_89', 'NITs': 'TRA_33', 'RCHO': 'TRA_12', 'C2H6': 'TRA_21', 'MPN': 'TRA_54', 'INO': 'TRA_82', 'MP': 'TRA_24', 'CH2Br2': 'TRA_52', 'SALC': 'TRA_43', 'NH3': 'TRA_30', 'CH2ICl': 'TRA_93', 'RIP': 'TRA_61', 'ClO': 'TRA_90', 'NO': 'TRA_1', 'SALA': 'TRA_42', 'MOBA': 'TRA_56', 'R4N2': 'TRA_17', 'BrCl': 'TRA_92', 'OClO': 'TRA_91', 'PMN': 'TRA_15', 'CO': 'TRA_4', 'CH2IBr': 'TRA_94', 'ISOP': 'TRA_6', 'BCPO': 'TRA_36', 'MVK': 'TRA_13', 'BrNO2': 'TRA_49', 'IONO': 'TRA_72', 'Cl2': 'TRA_87', 'HOBr': 'TRA_47', 'PROPNN': 'TRA_57', 'Cl': 'TRA_88', 'I2O2': 'TRA_74', 'I2O3': 'TRA_83', 'I2O4': 'TRA_84', 'I2O5': 'TRA_85', 'MEK': 'TRA_10', 'MMN': 'TRA_60', 'ISOPN': 'TRA_55', 'SO4s': 'TRA_28', 'I2O': 'TRA_81', 'ALK4': 'TRA_5', 'MSA': 'TRA_29', 'I2': 'TRA_67', 'Br2': 'TRA_44', 'IBr': 'TRA_77', 'MACR': 'TRA_14', 'I': 'TRA_79', 'AERI': 'TRA_86', 'HOI': 'TRA_68', 'BrO': 'TRA_46', 'NH4': 'TRA_31', 'SO2': 'TRA_26', 'SO4': 'TRA_27', 'IO': 'TRA_69', 'H2O2': 'TRA_8', 'BCPI': 'TRA_34', 'ICl': 'TRA_78', 'GLYC': 'TRA_59','OH': 'OH','HO2': 'HO2',},
                    'GCFP_d2TRA_all_1.6' :{'HIO3': 'TRA_80', 'TRA_17': 'TRA_17', 'TRA_16': 'TRA_16', 'TRA_15': 'TRA_15', 'TRA_14': 'TRA_14', 'TRA_13': 'TRA_13', 'TRA_12': 'TRA_12', 'TRA_11': 'TRA_11', 'TRA_19': 'TRA_19', 'ACET': 'ACET', 'RIP': 'TRA_61', 'BrNO3': 'TRA_50', 'HAC': 'TRA_58', 'ALD2': 'ALD2', 'HNO3': 'HNO3', 'HNO2': 'HNO2', 'HNO4': 'HNO4', 'OIO': 'TRA_70', 'MAP': 'TRA_63', 'PRPE': 'PRPE', 'TRA_29': 'TRA_29', 'CH2I2': 'TRA_76', 'I2O2': 'TRA_74', 'NIT': 'TRA_32', 'CH3Br': 'TRA_53', 'C3H7I': 'TRA_95', 'MO2': 'MO2', 'C3H8': 'C3H8', 'I2O5': 'TRA_85', 'TRA_71': 'TRA_71', 'TRA_70': 'TRA_70', 'TRA_73': 'TRA_73', 'DMS': 'DMS', 'TRA_75': 'TRA_75', 'TRA_74': 'TRA_74', 'TRA_77': 'TRA_77', 'TRA_76': 'TRA_76', 'CH2O': 'CH2O', 'TRA_78': 'TRA_78', 'CH3IT': 'TRA_75', 'NO2': 'NO2', 'NO3': 'NO3', 'N2O5': 'N2O5', 'H2O2': 'H2O2', 'PAN': 'PAN', 'HOCl': 'TRA_89', 'TRA_18': 'TRA_18', 'GMAO_TEMP': 'GMAO_TEMP', 'RCHO': 'RCHO', 'C2H6': 'C2H6', 'INO': 'TRA_82', 'MP': 'MP', 'CH2Br2': 'TRA_52', 'CH2ICl': 'TRA_93', 'TRA_59': 'TRA_59', 'TRA_58': 'TRA_58', 'IEPOX': 'TRA_62', 'TRA_53': 'TRA_53', 'TRA_52': 'TRA_52', 'TRA_51': 'TRA_51', 'TRA_50': 'TRA_50', 'TRA_57': 'TRA_57', 'TRA_56': 'TRA_56', 'TRA_55': 'TRA_55', 'TRA_54': 'TRA_54', 'MOBA': 'TRA_56', 'CH3I': 'TRA_75', 'BrCl': 'TRA_92', 'OClO': 'TRA_91', 'CO': 'CO', 'BCPI': 'TRA_34', 'ISOP': 'ISOP', 'BCPO': 'TRA_36', 'MVK': 'MVK', 'TRA_28': 'TRA_28', 'Cl': 'TRA_88', 'TRA_26': 'TRA_26', 'TRA_27': 'TRA_27', 'TRA_24': 'TRA_24', 'I2O3': 'TRA_83', 'I2O4': 'TRA_84', 'TRA_23': 'TRA_23', 'TRA_20': 'TRA_20', 'TRA_21': 'TRA_21', 'MMN': 'TRA_60', 'I2O': 'TRA_81', 'HBr': 'TRA_48', 'ALK4': 'ALK4', 'I2': 'TRA_67', 'PPN': 'PPN', 'IBr': 'TRA_77', 'I': 'TRA_79', 'AERI': 'TRA_86', 'NH4': 'TRA_31', 'SO2': 'SO2', 'SO4': 'SO4', 'NH3': 'TRA_30', 'TRA_08': 'TRA_08', 'TRA_09': 'TRA_09', 'TRA_01': 'TRA_01', 'TRA_02': 'TRA_02', 'TRA_03': 'TRA_03', 'TRA_04': 'TRA_04', 'TRA_05': 'TRA_05', 'TRA_06': 'TRA_06', 'TRA_07': 'TRA_07', 'OCPI': 'TRA_35', 'OCPO': 'TRA_37', 'Br2': 'TRA_44', 'O3': 'O3', 'Br': 'TRA_45', 'TRA_96': 'TRA_96', 'TRA_95': 'TRA_95', 'TRA_94': 'TRA_94', 'TRA_93': 'TRA_93', 'TRA_92': 'TRA_92', 'TRA_91': 'TRA_91', 'TRA_90': 'TRA_90', 'TRA_62': 'TRA_62', 'TRA_63': 'TRA_63', 'TRA_60': 'TRA_60', 'TRA_61': 'TRA_61', 'TRA_66': 'TRA_66', 'TRA_67': 'TRA_67', 'C2H5I': 'TRA_96', 'TRA_65': 'TRA_65', 'TRA_68': 'TRA_68', 'TRA_69': 'TRA_69', 'OH': 'OH', 'IONO2': 'TRA_73', 'HI': 'TRA_71', 'CHBr3': 'TRA_51', 'TRA_46': 'TRA_46', 'DST4': 'TRA_41', 'DST3': 'TRA_40', 'DST2': 'TRA_39', 'DST1': 'TRA_38', 'NITs': 'TRA_33', 'TRA_48': 'TRA_48', 'TRA_49': 'TRA_49', 'TRA_44': 'TRA_44', 'TRA_45': 'TRA_45', 'RO2': 'RO2', 'TRA_47': 'TRA_47', 'TRA_40': 'TRA_40', 'TRA_41': 'TRA_41', 'TRA_42': 'TRA_42', 'TRA_43': 'TRA_43', 'MPN': 'TRA_54', 'ETO2': 'ETO2', 'IO': 'TRA_69', 'TRA_64': 'TRA_64', 'ClO': 'TRA_90', 'NO': 'NO', 'SALA': 'TRA_42', 'SALC': 'TRA_43', 'R4N2': 'R4N2', 'PMN': 'PMN', 'TRA_25': 'TRA_25', 'CH2IBr': 'TRA_94', 'TRA_22': 'TRA_22', 'BrNO2': 'TRA_49', 'IONO': 'TRA_72', 'Cl2': 'TRA_87', 'HOBr': 'TRA_47', 'PROPNN': 'TRA_57', 'MEK': 'MEK', 'TRA_72': 'TRA_72', 'ISOPN': 'TRA_55', 'SO4s': 'TRA_28', 'TRA_79': 'TRA_79', 'MSA': 'MSA', 'TRA_39': 'TRA_39', 'TRA_38': 'TRA_38', 'GLYC': 'TRA_59', 'TRA_35': 'TRA_35', 'TRA_34': 'TRA_34', 'TRA_37': 'TRA_37', 'TRA_36': 'TRA_36', 'TRA_31': 'TRA_31', 'TRA_30': 'TRA_30', 'TRA_33': 'TRA_33', 'TRA_32': 'TRA_32', 'HO2': 'HO2', 'MACR': 'TRA_14', 'HOI': 'TRA_68', 'BrO': 'TRA_46', 'ICl': 'TRA_78', 'TRA_80': 'TRA_80', 'TRA_81': 'TRA_81', 'TRA_82': 'TRA_82', 'TRA_83': 'TRA_83', 'TRA_84': 'TRA_84', 'TRA_85': 'TRA_85', 'TRA_86': 'TRA_86', 'TRA_87': 'TRA_87', 'TRA_88': 'TRA_88', 'TRA_89': 'TRA_89','GMAO_TEMP': 'GMAO_TEMP', 'GMAO_UWND': 'GMAO_UWND', 'GMAO_VWND': 'GMAO_VWND'},
#                    'GCFP_d2TRA_all_1.7' : {'TRA_74': 'ICl', 'TRA_25': 'DMS', 'TRA_68': 'CH2I2', 'TRA_44': 'Br2', 'TRA_70': 'CH2IBr', 'TRA_22': 'N2O5', 'TRA_76': 'IO', 'TRA_79': 'INO', 'TRA_23': 'HNO4', 'TRA_17': 'R4N2', 'TRA_16': 'PPN', 'TRA_15': 'PMN', 'TRA_14': 'MACR', 'TRA_13': 'MVK', 'TRA_12': 'RCHO', 'TRA_11': 'ALD2', 'TRA_10': 'MEK', 'TRA_53': 'CH3Br', 'TRA_52': 'CH2Br2', 'TRA_51': 'CHBr3', 'TRA_21': 'C2H6', 'TRA_57': 'PROPNN', 'TRA_56': 'MOBA', 'TRA_19': 'C3H8', 'TRA_18': 'PRPE', 'TRA_69': 'CH2ICl', 'TRA_50': 'BrNO3', 'TRA_39': 'DST2', 'TRA_38': 'DST1', 'TRA_73': 'IBr', 'TRA_35': 'OCPI', 'TRA_34': 'BCPI', 'TRA_37': 'OCPO', 'TRA_36': 'BCPO', 'TRA_31': 'NH4', 'TRA_30': 'NH3', 'TRA_33': 'NITs', 'TRA_32': 'NIT', 'TRA_77': 'HI', 'TRA_83': 'I2O3', 'TRA_55': 'ISOPN', 'TRA_54': 'MPN', 'TRA_72': 'I2', 'TRA_59': 'GLYC', 'TRA_62': 'IEPOX', 'TRA_63': 'MAP', 'TRA_60': 'MMN', 'TRA_61': 'RIP', 'TRA_48': 'HBr', 'TRA_49': 'BrNO2', 'TRA_64': 'NO2', 'TRA_65': 'NO3', 'TRA_20': 'CH2O', 'TRA_45': 'Br', 'TRA_46': 'BrO', 'TRA_47': 'HOBr', 'TRA_40': 'DST3', 'TRA_41': 'DST4', 'TRA_42': 'SALA', 'TRA_43': 'SALC', 'TRA_08': 'H2O2', 'TRA_09': 'ACET', 'TRA_75': 'I', 'TRA_28': 'SO4s', 'TRA_29': 'MSA', 'TRA_26': 'SO2', 'TRA_01': 'NO', 'TRA_02': 'O3', 'TRA_03': 'PAN', 'TRA_04': 'CO', 'TRA_05': 'ALK4', 'TRA_06': 'ISOP', 'TRA_07': 'HNO3', 'TRA_80': 'IONO', 'TRA_81': 'IONO2', 'TRA_82': 'I2O2', 'TRA_58': 'HAC', 'TRA_84': 'I2O4', 'TRA_85': 'AERI', 'TRA_27': 'SO4', 'TRA_78': 'OIO', 'TRA_66': 'HNO2', 'TRA_71': 'HOI', 'TRA_24': 'MP', 'TRA_67': 'CH3IT', 'TRA_9': 'ACET', 'TRA_8': 'H2O2', 'TRA_7': 'HNO3', 'TRA_6': 'ISOP', 'TRA_5': 'ALK4', 'TRA_4': 'CO', 'TRA_3': 'PAN', 'TRA_2': 'O3', 'TRA_1': 'NO'},           
                'GCFP_d2TRA_all_1.7' : {'TRA_25': 'DMS', 'TRA_77': 'HI', 'TRA_76': 'IO', 'TRA_23': 'HNO4', 'TRA_71': 'HOI', 'TRA_70': 'CH2IBr', 'TRA_15': 'PMN', 'TRA_14': 'MACR', 'TRA_13': 'MVK', 'TRA_12': 'RCHO', 'TRA_11': 'ALD2', 'TRA_10': 'MEK', 'TRA_79': 'INO', 'TRA_78': 'OIO', 'TRA_51': 'CHBr3', 'TRA_50': 'BrNO3', 'TRA_52': 'CH2Br2', 'TRA_46': 'BrO', 'TRA_19': 'C3H8', 'TRA_18': 'PRPE', 'TRA_47': 'HOBr', 'TRA_39': 'DST2', 'TRA_38': 'DST1', 'TRA_81': 'IONO2', 'TRA_35': 'OCPI', 'TRA_57': 'PROPNN', 'TRA_37': 'OCPO', 'TRA_36': 'BCPO', 'TRA_31': 'NH4', 'TRA_30': 'NH3', 'TRA_33': 'NITs', 'TRA_56': 'MOBA', 'TRA_83': 'I2O3', 'TRA_55': 'ISOPN', 'TRA_84': 'I2O4', 'TRA_54': 'MPN', 'TRA_5': 'ALK4', 'TRA_49': 'BrNO2', 'TRA_32': 'NIT', 'TRA_9': 'ACET', 'TRA_8': 'H2O2', 'TRA_7': 'HNO3', 'TRA_6': 'ISOP', 'TRA_59': 'GLYC', 'TRA_4': 'CO', 'TRA_3': 'PAN', 'TRA_2': 'O3', 'TRA_1': 'NO', 'TRA_62': 'IEPOX', 'TRA_63': 'MAP', 'TRA_60': 'MMN', 'TRA_61': 'RIP', 'TRA_66': 'HNO2', 'TRA_67': 'CH3IT', 'TRA_64': 'NO2', 'TRA_65': 'NO3', 'TRA_44': 'Br2', 'TRA_45': 'Br', 'TRA_68': 'CH2I2', 'TRA_69': 'CH2ICl', 'TRA_40': 'DST3', 'TRA_41': 'DST4', 'TRA_42': 'SALA', 'TRA_43': 'SALC', 'TRA_17': 'R4N2', 'TRA_28': 'SO4s', 'TRA_16': 'PPN', 'TRA_58': 'HAC', 'TRA_27': 'SO4', 'TRA_24': 'MP', 'TRA_29': 'MSA', 'TRA_22': 'N2O5', 'TRA_73': 'IBr', 'TRA_20': 'CH2O', 'TRA_21': 'C2H6', 'TRA_80': 'IONO', 'TRA_26': 'SO2', 'TRA_82': 'I2O2', 'TRA_72': 'I2', 'TRA_48': 'HBr', 'TRA_85': 'AERI', 'TRA_34': 'BCPI', 'TRA_75': 'I', 'TRA_53': 'CH3Br', 'TRA_74': 'ICl'}, 

                    'GCFP_d2TRA_all_1.7_EOH_actual_names' : {'HNO4': 'HNO4', 'PPN': 'PPN', 'TRA_17': 'R4N2', 'TRA_16': 'PPN', 'TRA_15': 'PMN', 'TRA_14': 'MACR', 'TRA_13': 'MVK', 'TRA_12': 'RCHO', 'TRA_11': 'ALD2', 'TRA_10': 'MEK', 'O3': 'O3', 'TRA_19': 'C3H8', 'TRA_18': 'PRPE', 'GMAO_UWND': 'GMAO_UWND', 'TRA_62': 'IEPOX', 'TRA_63': 'MAP', 'TRA_60': 'MMN', 'TRA_61': 'RIP', 'TRA_66': 'HNO2', 'TRA_67': 'CH3IT', 'TRA_65': 'NO3', 'TRA_68': 'CH2I2', 'TRA_69': 'CH2ICl', 'OH': 'OH', 'LAT': 'LAT', 'TRA_71': 'HOI', 'TRA_70': 'CH2IBr', 'TRA_73': 'IBr', 'TRA_72': 'I2', 'TRA_75': 'I', 'TRA_74': 'ICl', 'TRA_77': 'HI', 'TRA_76': 'IO', 'TRA_79': 'INO', 'TRA_78': 'OIO', 'NO2': 'NO2', 'NO3': 'NO3', 'N2O5': 'N2O5', 'H2O2': 'H2O2', 'GMAO_VWND': 'GMAO_VWND', 'PAN': 'PAN', 'GMAO_TEMP': 'GMAO_TEMP', 'TRA_48': 'HBr', 'TRA_49': 'BrNO2', 'TRA_44': 'Br2', 'TRA_45': 'Br', 'TRA_46': 'BrO', 'TRA_47': 'HOBr', 'TRA_40': 'DST3', 'TRA_41': 'DST4', 'TRA_42': 'SALA', 'TRA_43': 'SALC', 'TRA_59': 'GLYC', 'TRA_58': 'HAC', 'TRA_53': 'CH3Br', 'TRA_52': 'CH2Br2', 'TRA_51': 'CHBr3', 'TRA_50': 'BrNO3', 'TRA_57': 'PROPNN', 'TRA_56': 'MOBA', 'TRA_55': 'ISOPN', 'TRA_54': 'MPN', 'NO': 'NO', 'PMN': 'PMN', 'HNO3': 'HNO3', 'TRA_28': 'SO4s', 'TRA_29': 'MSA', 'TRA_26': 'SO2', 'TRA_27': 'SO4', 'TRA_24': 'MP', 'TRA_25': 'DMS', 'TRA_22': 'N2O5', 'TRA_23': 'HNO4', 'TRA_20': 'CH2O', 'TRA_21': 'C2H6', 'RO2': 'RO2', 'LON': 'LON', 'TRA_39': 'DST2', 'TRA_38': 'DST1', 'TRA_35': 'OCPI', 'TRA_34': 'BCPI', 'TRA_37': 'OCPO', 'TRA_36': 'BCPO', 'TRA_31': 'NH4', 'TRA_30': 'NH3', 'TRA_33': 'NITs', 'TRA_32': 'NIT', 'HO2': 'HO2', 'SO2': 'SO2', 'SO4': 'SO4', 'TRA_08': 'H2O2', 'TRA_09': 'ACET', 'HNO2': 'HNO2', 'TRA_03': 'PAN', 'TRA_04': 'CO', 'TRA_05': 'ALK4', 'TRA_06': 'ISOP', 'TRA_07': 'HNO3', 'TRA_80': 'IONO', 'TRA_81': 'IONO2', 'TRA_82': 'I2O2', 'TRA_83': 'I2O3', 'TRA_84': 'I2O4', 'TRA_85': 'AERI', 'TRA_86': 'EOH'},      
        'TRA_spec_met_all_1.7_EOH': {'MAO3': 'MAO3', 'DHMOB': 'DHMOB', 'ETP': 'ETP', 'RCO3': 'RCO3', 'MO2': 'MO2', 'EOH': 'EOH', 'MVKN': 'MVKN', 'R4P': 'R4P', 'ISNP': 'ISNP', 'RB3P': 'RB3P', 'MGLY': 'MGLY', 'MAOPO2': 'MAOPO2', 'RIO2': 'RIO2', 'PMNN': 'PMNN', 'PP': 'PP', 'VRP': 'VRP', 'RP': 'RP', 'MRO2': 'MRO2', 'HC5': 'HC5', 'ATO2': 'ATO2', 'PYAC': 'PYAC', 'R4N1': 'R4N1', 'DIBOO': 'DIBOO', 'LISOPOH': 'LISOPOH', 'HO2': 'HO2', 'ETHLN': 'ETHLN', 'ISNOOB': 'ISNOOB', 'ISNOOA': 'ISNOOA', 'ROH': 'ROH', 'MAN2': 'MAN2', 'B3O2': 'B3O2', 'INPN': 'INPN', 'MACRN': 'MACRN', 'PO2': 'PO2', 'VRO2': 'VRO2', 'MRP': 'MRP', 'PRN1': 'PRN1', 'ISNOHOO': 'ISNOHOO', 'MOBAOO': 'MOBAOO', 'MACRNO2': 'MACRNO2', 'ISOPND': 'ISOPND', 'HC5OO': 'HC5OO', 'ISOPNBO2': 'ISOPNBO2', 'RA3P': 'RA3P', 'ISOPNB': 'ISOPNB', 'ISOPNDO2': 'ISOPNDO2', 'PMNO2': 'PMNO2', 'IAP': 'IAP', 'MCO3': 'MCO3', 'IEPOXOO': 'IEPOXOO', 'MAOP': 'MAOP', 'INO2': 'INO2', 'OH': 'OH', 'PRPN': 'PRPN', 'GLYX': 'GLYX', 'A3O2': 'A3O2', 'ETO2': 'ETO2', 'R4O2': 'R4O2', 'ISN1': 'ISN1', 'KO2': 'KO2', 'ATOOH': 'ATOOH','GMAO_PSFC': 'GMAO_PSFC', 'GMAO_SURF': 'GMAO_SURF', 'GMAO_TEMP': 'GMAO_TEMP', 'GMAO_ABSH': 'GMAO_ABSH', 'GMAO_UWND': 'GMAO_UWND', 'GMAO_VWND': 'GMAO_VWND', 'TRA_9': 'ACET', 'TRA_8': 'H2O2', 'TRA_7': 'HNO3', 'TRA_6': 'ISOP', 'TRA_5': 'ALK4', 'TRA_4': 'CO', 'TRA_3': 'PAN', 'TRA_2': 'O3', 'TRA_1': 'NO', 'TRA_74': 'ICl', 'TRA_25': 'DMS', 'TRA_68': 'CH2I2', 'TRA_44': 'Br2', 'TRA_70': 'CH2IBr', 'TRA_22': 'N2O5', 'TRA_76': 'IO', 'TRA_79': 'INO', 'TRA_23': 'HNO4', 'TRA_17': 'R4N2', 'TRA_16': 'PPN', 'TRA_15': 'PMN', 'TRA_14': 'MACR', 'TRA_13': 'MVK', 'TRA_12': 'RCHO', 'TRA_11': 'ALD2', 'TRA_10': 'MEK', 'TRA_53': 'CH3Br', 'TRA_52': 'CH2Br2', 'TRA_51': 'CHBr3', 'TRA_21': 'C2H6', 'TRA_57': 'PROPNN', 'TRA_56': 'MOBA', 'TRA_19': 'C3H8', 'TRA_18': 'PRPE', 'TRA_69': 'CH2ICl', 'TRA_50': 'BrNO3', 'TRA_39': 'DST2', 'TRA_38': 'DST1', 'TRA_73': 'IBr', 'TRA_35': 'OCPI', 'TRA_34': 'BCPI', 'TRA_37': 'OCPO', 'TRA_36': 'BCPO', 'TRA_31': 'NH4', 'TRA_30': 'NH3', 'TRA_33': 'NITs', 'TRA_32': 'NIT', 'TRA_77': 'HI', 'TRA_83': 'I2O3', 'TRA_55': 'ISOPN', 'TRA_54': 'MPN', 'TRA_72': 'I2', 'TRA_59': 'GLYC', 'TRA_62': 'IEPOX', 'TRA_63': 'MAP', 'TRA_60': 'MMN', 'TRA_61': 'RIP', 'TRA_48': 'HBr', 'TRA_49': 'BrNO2', 'TRA_64': 'NO2', 'TRA_65': 'NO3', 'TRA_20': 'CH2O', 'TRA_45': 'Br', 'TRA_46': 'BrO', 'TRA_47': 'HOBr', 'TRA_40': 'DST3', 'TRA_41': 'DST4', 'TRA_42': 'SALA', 'TRA_43': 'SALC', 'TRA_08': 'H2O2', 'TRA_09': 'ACET', 'TRA_75': 'I', 'TRA_28': 'SO4s', 'TRA_29': 'MSA', 'TRA_26': 'SO2', 'TRA_01': 'NO', 'TRA_02': 'O3', 'TRA_03': 'PAN', 'TRA_04': 'CO', 'TRA_05': 'ALK4', 'TRA_06': 'ISOP', 'TRA_07': 'HNO3', 'TRA_80': 'IONO', 'TRA_81': 'IONO2', 'TRA_82': 'I2O2', 'TRA_58': 'HAC', 'TRA_84': 'I2O4', 'TRA_85': 'AERI', 'TRA_27': 'SO4', 'TRA_78': 'OIO', 'TRA_66': 'HNO2', 'TRA_71': 'HOI', 'TRA_24': 'MP', 'TRA_67': 'CH3IT' }, 
           'TRA_spec_met_all_1.7_EOH_no_trailing_zeroes':  {'EOH': 'EOH', 'TRA_17': 'R4N2', 'TRA_16': 'PPN', 'TRA_15': 'PMN', 'TRA_14': 'MACR', 'TRA_13': 'MVK', 'TRA_12': 'RCHO', 'TRA_11': 'ALD2', 'TRA_10': 'MEK', 'TRA_19': 'C3H8', 'TRA_18': 'PRPE', 'DHMOB': 'DHMOB', 'RP': 'RP', 'GMAO_UWND': 'GMAO_UWND', 'MAN2': 'MAN2', 'B3O2': 'B3O2', 'MRP': 'MRP', 'PRN1': 'PRN1', 'TRA_62': 'IEPOX', 'TRA_63': 'MAP', 'TRA_60': 'MMN', 'TRA_61': 'RIP', 'TRA_66': 'HNO2', 'TRA_67': 'CH3IT', 'TRA_64': 'NO2', 'TRA_65': 'NO3', 'TRA_68': 'CH2I2', 'TRA_69': 'CH2ICl', 'IAP': 'IAP', 'MCO3': 'MCO3', 'GMAO_SURF': 'GMAO_SURF', 'OH': 'OH', 'PRPN': 'PRPN', 'TRA7': 'HNO3', 'TRA6': 'ISOP', 'MAO3': 'MAO3', 'RCO3': 'RCO3', 'MO2': 'MO2', 'MACRNO2': 'MACRNO2', 'TRA_23': 'HNO4', 'TRA_71': 'HOI', 'TRA_70': 'CH2IBr', 'TRA_73': 'IBr', 'TRA_72': 'I2', 'TRA_75': 'I', 'TRA_74': 'ICl', 'ISNP': 'ISNP', 'TRA_76': 'IO', 'TRA_79': 'INO', 'RB3P': 'RB3P', 'TRA_51': 'CHBr3', 'ROH': 'ROH', 'PP': 'PP', 'ISOPNDO2': 'ISOPNDO2', 'HC5': 'HC5', 'TRA9': 'ACET', 'TRA_56': 'MOBA', 'MACRN': 'MACRN', 'DIBOO': 'DIBOO', 'MRO2': 'MRO2', 'INPN': 'INPN', 'GMAO_TEMP': 'GMAO_TEMP', 'PO2': 'PO2', 'ISOPND': 'ISOPND', 'TRA_48': 'HBr', 'TRA_1': 'NO', 'RA3P': 'RA3P', 'ISOPNB': 'ISOPNB', 'TRA_44': 'Br2', 'TRA_45': 'Br', 'TRA_46': 'BrO', 'TRA_47': 'HOBr', 'TRA_40': 'DST3', 'TRA_41': 'DST4', 'TRA_42': 'SALA', 'TRA_43': 'SALC', 'IEPOXOO': 'IEPOXOO', 'MAOP': 'MAOP', 'INO2': 'INO2', 'ETO2': 'ETO2', 'ISN1': 'ISN1', 'TRA_49': 'BrNO2', 'ETP': 'ETP', 'TRA5': 'ALK4', 'TRA4': 'CO', 'TRA_59': 'GLYC', 'TRA_58': 'HAC', 'TRA1': 'NO', 'R4P': 'R4P', 'TRA_3': 'PAN', 'TRA2': 'O3', 'TRA_53': 'CH3Br', 'TRA_52': 'CH2Br2', 'MAOPO2': 'MAOPO2', 'TRA_50': 'BrNO3', 'TRA_57': 'PROPNN', 'GMAO_VWND': 'GMAO_VWND', 'TRA_55': 'ISOPN', 'TRA_54': 'MPN', 'RIO2': 'RIO2', 'VRP': 'VRP', 'R4N1': 'R4N1', 'GMAO_PSFC': 'GMAO_PSFC', 'VRO2': 'VRO2', 'ISNOHOO': 'ISNOHOO', 'HC5OO': 'HC5OO', 'ETHLN': 'ETHLN', 'TRA_28': 'SO4s', 'TRA_29': 'MSA', 'TRA_26': 'SO2', 'TRA_27': 'SO4', 'TRA_24': 'MP', 'TRA_25': 'DMS', 'TRA_22': 'N2O5', 'R4O2': 'R4O2', 'TRA_20': 'CH2O', 'TRA_21': 'C2H6', 'TRA_77': 'HI', 'MVKN': 'MVKN', 'TRA_35': 'OCPI', 'TRA_78': 'OIO', 'MGLY': 'MGLY', 'PMNN': 'PMNN', 'TRA_39': 'DST2', 'TRA_38': 'DST1', 'ATO2': 'ATO2', 'TRA_34': 'BCPI', 'TRA_37': 'OCPO', 'TRA_36': 'BCPO', 'TRA_31': 'NH4', 'TRA_30': 'NH3', 'TRA_33': 'NITs', 'TRA_32': 'NIT', 'LISOPOH': 'LISOPOH', 'HO2': 'HO2', 'GMAO_ABSH': 'GMAO_ABSH', 'TRA8': 'H2O2', 'ISNOOB': 'ISNOOB', 'ISNOOA': 'ISNOOA', 'TRA_9': 'ACET', 'TRA_8': 'H2O2', 'TRA_7': 'HNO3', 'TRA_6': 'ISOP', 'TRA_5': 'ALK4', 'TRA_4': 'CO', 'TRA_3': 'PAN', 'TRA_2': 'O3', 'ISOPNBO2': 'ISOPNBO2', 'MOBAOO': 'MOBAOO', 'PMNO2': 'PMNO2', 'GLYX': 'GLYX', 'A3O2': 'A3O2', 'TRA_80': 'IONO', 'TRA_81': 'IONO2', 'TRA_82': 'I2O2', 'TRA_83': 'I2O3', 'TRA_84': 'I2O4', 'TRA_85': 'AERI', 'PYAC': 'PYAC', 'KO2': 'KO2', 'ATOOH': 'ATOOH', 'PRESS':'PRESS', 'U10M':'U10M'},
                    'red_specs_f_name': ['O3', 'NO2', 'NO', 'NO3', 'N2O5', 'HNO4', 'HNO3', 'HNO2', 'PAN', 'PPN', 'PMN', 'H2O2', 'HO2', 'OH', 'RO2', 'SO2', 'SO4', 'GMAO_TEMP', 'GMAO_UWND', 'GMAO_VWND', 'I2', 'HOI', 'IO', 'I', 'HI', 'OIO', 'INO', 'IONO', 'IONO2', 'I2O2', 'I2O4', 'I2O3', 'CH3IT', 'CH2I2', 'CH2ICl', 'CH2IBr'], 
                     
                    # Photolysis/Fast-J
                    'FastJ_lower' : [289.0, 298.25, 307.45, 312.45, 320.3, 345.0, 412.45],
                    'FastJ_upper' : [298.25, 307.45, 312.45, 320.3, 345.0, 412.45, 850.0],
                    'FastJ_mids' :  [294,303,310,316,333,380,574],
                    }  
    
    if rtn_dict:
        return GC_var_dict
    else:    
        return GC_var_dict[input_x]

# --------------
# 4.05 -  GEOS-Chem/ctm.bpch values
# --------------
def latex_spec_name(input_x, debug=False):
    spec_dict = {
            'OIO': 'OIO', 'C3H7I': 'C$_{3}$H$_{7}$I', 'IO': 'IO', 'I': 'I', 'I2': 'I$_{2}$', 'CH2ICl': 'CH$_{2}$ICl', 'HOI': 'HOI', 'CH2IBr': 'CH$_{2}$IBr', 'C2H5I': 'C$_{2}$H$_{5}$I', 'CH2I2': 'CH$_{2}$I$_{2}$', 'CH3IT': 'CH$_{3}$I', 'IONO': 'IONO','HIO3': 'HIO$_{3}$', 'ICl': 'ICl', 'I2O3': 'I$_{2}$O$_{3}$', 'I2O4': 'I$_{2}$O$_{4}$', 'I2O5': 'I$_{2}$O$_{5}$', 'INO': 'INO', 'I2O': 'I$_{2}$O', 'IBr': 'IBr','I2O2': 'I$_{2}$O$_{2}$', 'IONO2': 'IONO$_{2}$', 'HI':'HI', 'BrO':'BrO','Br':'Br','HOBr':'HOBr','Br2':'Br$_{2}$','CH3Br':'CH$_{3}$Br','CH2Br2':'CH$_{2}$Br$_{2}$', 'CHBr3':'CHBr$_{3}$','O3':'O$_{3}$', 'CO':'CO' , 'DMS':'DMS', 'NOx':'NOx', 'NO':'NO', 'NO2':'NO$_{2}$', 'NO3':'NO$_{3}$','HNO3':'HNO$_{3}$', 'HNO4':'HNO$_{4}$','PAN':'PAN', 'HNO2':'HNO$_{2}$', 'N2O5':'N$_{2}$O$_{5}$','ALK4':'>= C4 alkanes','ISOP':'Isoprene' ,'H2O2':'H$_{2}$O$_{2}$','ACET':'CH$_{3}$C(O)CH$_{3}$', 'MEK':'>C3 ketones', 'ALD2':'CH$_{3}$CHO', 'RCHO': 'CH$_{3}$CH$_{2}$CHO', 'MVK':'CH$_{2}$=CHC(O)CH$_{3}$', 'MACR':'Methacrolein', 'PMN':'CH$_{2}$=C(CH$_{3}$)C(O)OONO$_{2}$', 'PPN':'CH$_{3}$CH$_{2}$C(O)OONO$_{2}$', 'R4N2':'>= C4 alkylnitrates','PRPE':'>= C4 alkenes', 'C3H8':'C$_{3}$H$_{8}$','CH2O':'CH$_{2}$O', 'C2H6':'C$_{2}$H$_{6}$', 'MP':'CH$_{3}$OOH', 'SO2':'SO$_{2}$', 'SO4':'SO$_{4}$','SO4s':'SO$_{4}$ on SSA', 'MSA':'CH$_{4}$SO$_{3}$','NH3':'NH$_{3}$', 'NH4': 'NH$_{4}$', 'NIT': 'InOrg N', 'NITs': 'InOrg N on SSA', 'BCPI':'BCPI', 'OCPI':'OCPI', 'BCPO':'BCPO','OCPO':'OCPO', 'DST1':'DST1', 'DST2':'DST2','DST3':'DST3','DST4':'DST4','SALA':'SALA', 'SALC':'SALC',  'HBr':'HBr', 'BrNO2': 'BrNO$_{2}$', 'BrNO3': 'BrNO$_{3}$', 'MPN':'CH$_{3}$ON$_{2}$', 'ISOPN':'ISOPN', 'MOBA':'MOBA', 'PROPNN':'PROPNN', 'HAC':'HAC', 'GLYC':'GLYC', 'MMN':'MMN', 'RIP':'RIP', 'IEPOX':'IEPOX','MAP':'MAP', 'AERI':'Aerosol Iodine' , 'Cl2':'Cl$_{2}$', 'Cl':'Cl','HOCl':'HOCl','ClO':'ClO','OClO':'OClO','BrCl':'BrCl', 'HI+OIO+IONO+INO':'HI+OIO+IONO+INO','CH2IX':'CH$_{2}$IX', 'IxOy':'I$_{2}$O$_{X}$ ($_{X}$=2,3,4)', 'CH3I':'CH$_{3}$I', 'OH':'OH', 'HO2':'HO$_{2}$', 'MO2':'MO$_{2}$', 'RO2':'RO$_{2}$'

            ,'RD01':r'I + O$_{3}$ $\rightarrow$ IO + O$_{2}$'
            # Analysis names 
            ,'iodine_all':'All Iodine', 'Iy': 'I$_{Y}$', 'IOy': 'IO$_{Y}$', 'IyOx': 'I$_{Y}$O$_{X}$', 'IOx': 'IO$_{X}$','iodine_all_A':'All Iodine (Inc. AERI)', 'I2Ox': 'I$_{2}$O$_{X}$' , 'AERI/SO4': 'AERI/SO4', 'EOH':'Ethanol'
            , 'PSURF': 'Pressure at the bottom of level', 'GMAO_TEMP' : 'Temperature', 'TSKIN' : 'Temperature at 2m', 'GMAO_UWND':'Zonal Wind', 'GMAO_VWND':'Meridional Wind', 'U10M':'10m Meridional Wind', 'V10M': '10m Zonal Wind', 'CH2OO':'CH$_{2}$OO',
        # Family Names
           'N_specs' : 'NOy', 'NOy' :  'NO$_Y$', 
           'N_specs_no_I' : 'NOy exc. iodine',
        # typos
            'CH2BR2':'CH$_{2}$Br$_{2}$',
    

            }
    return spec_dict[input_x]

# --------------
# 4.06 - converts P/L tracer mulitpler to 1
# --------------
def p_l_unity(rxn, debug=False):
    p_l_dict = {
    'LR24': 1.0, 'LR25': 1.0, 'LR26': 1.0, 'LR27': 1.0, 'LR20': 1.0, 'LR21': 1.0, 'LR22': 1.0, 'LR30': 1.0, 'LR31': 1.0, 'LR23': 1.0, 'LR28': 1.0, 'LR29': 1.0, 'RD09': 1.0, 'PO3_46': 0.25, 'LR3': 1.0, 'LR2': 1.0, 'RD02': 1.0, 'PO3_03': 0.3, 'PO3_14': 1.0, 'PO3_02': 0.15, 'PO3_05': 0.15
    }
    return p_l_dict[rxn]

# --------------
#  4.07 - Returns tracers unit and scale (if requested)
# --------------
def tra_unit(x, scale=False, adjustment=False, adjust=True, \
            global_unit=False, ClearFlo_unit=False, debug=False ):
    tra_unit = {
    'OCPI': 'ppbv', 'OCPO': 'ppbv', 'PPN': 'ppbv', 'HIO3': 'pptv', 'O3': 'ppbv', 'PAN': 'ppbv', 'ACET': 'ppbC', 'RIP': 'ppbv', 'BrNO3': 'pptv', 'Br': 'pptv', 'HBr': 'pptv', 'HAC': 'ppbv', 'ALD2': 'ppbC', 'HNO3': 'ppbv', 'HNO2': 'ppbv', 'C2H5I': 'pptv', 'HNO4': 'ppbv', 'OIO': 'pptv', 'MAP': 'ppbv', 'PRPE': 'ppbC', 'HI': 'pptv', 'CH2I2': 'pptv', 'IONO2': 'pptv', 'NIT': 'ppbv', 'CH3Br': 'pptv', 'C3H7I': 'pptv', 'C3H8': 'ppbC', 'DMS': 'ppbv', 'CH2O': 'ppbv', 'CH3IT': 'pptv', 'NO2': 'ppbv', 'NO3': 'ppbv', 'N2O5': 'ppbv', 'CHBr3': 'pptv', 'DST4': 'ppbv', 'DST3': 'ppbv', 'DST2': 'ppbv', 'DST1': 'ppbv', 'HOCl': 'ppbv', 'NITs': 'ppbv', 'RCHO': 'ppbv', 'C2H6': 'ppbC', 'MPN': 'ppbv', 'INO': 'pptv', 'MP': 'ppbv', 'CH2Br2': 'pptv', 'SALC': 'ppbv', 'NH3': 'ppbv', 'CH2ICl': 'pptv', 'IEPOX': 'ppbv', 'ClO': 'ppbv', 'NO': 'pptv', 'SALA': 'ppbv', 'MOBA': 'ppbv', 'R4N2': 'ppbv', 'BrCl': 'pptv', 'OClO': 'ppbv', 'PMN': 'ppbv', 'CO': 'ppbv', 'CH2IBr': 'pptv', 'ISOP': 'ppbC', 'BCPO': 'ppbv', 'MVK': 'ppbv', 'BrNO2': 'pptv', 'IONO': 'pptv', 'Cl2': 'ppbv', 'HOBr': 'pptv', 'PROPNN': 'ppbv', 'Cl': 'ppbv', 'I2O2': 'pptv', 'I2O3': 'pptv', 'I2O4': 'pptv', 'I2O5': 'pptv', 'MEK': 'ppbC', 'MMN': 'ppbv', 'ISOPN': 'ppbv', 'SO4s': 'ppbv', 'I2O': 'pptv', 'ALK4': 'ppbC', 'MSA': 'ppbv', 'I2': 'pptv', 'Br2': 'pptv', 'IBr': 'pptv', 'MACR': 'ppbv', 'I': 'pptv', 'AERI': 'pptv', 'HOI': 'pptv', 'BrO': 'pptv', 'NH4': 'ppbv', 'SO2': 'ppbv', 'SO4': 'ppbv', 'IO': 'pptv', 'H2O2': 'ppbv', 'BCPI': 'ppbv', 'ICl': 'pptv', 'GLYC': 'ppbv',
    # Extra diagnostics to allow for simplified processing 
    'CH3I':'pptv', 'Iy':'pptv', 'PSURF': 'hPa', 'OH':'pptv', 'HO2':'pptv',
    'MO2': 'pptv', 'NOy':'pptbv','EOH': 'ppbv' , 'CO':'ppbv', 'CH4':'ppbv', 
    'TSKIN':'K', 'GMAO_TEMP': 'K', 'GMAO_VWND' :'m/s','GMAO_UWND': 'm/s', 'RO2': 'pptv', 'U10M':'m/s','V10M': 'm/s' , 'PRESS': 'hPa', 'CH2OO':'pptv',
    
    } 
    units = tra_unit[x]

    # Adjust to appropriate scale for pf analysis 
    if adjust:
        spec_2_pptv = GC_var('spec_2_pptv')
        spec_2_pptC = GC_var('spec_2_pptC') 
        if ( x in spec_2_pptv ):
            if (debug):
                print 'adjusting {} ({}) to {}'.format(x, units, 'pptv'   )
            units = 'pptv'
        if ( x in spec_2_pptC ):
            if (debug):
                print 'adjusting {} ({}) to {}'.format(x, units, 'pptC' )
            units = 'pptC'

    # Over ride adjustments for globally appro. units
    if global_unit:
        spec_2_ppbv = GC_var('spec_2_ppbv')                    
        spec_2_ppbC = GC_var('spec_2_ppbC')
        if ( x in spec_2_ppbv ):
            if (debug):
                print 'adjusting {} ({}) to {}'.format(x, units, 'ppbv'   )
            units = 'ppbv'
        if ( x in spec_2_ppbC ):
            if (debug):
                print 'adjusting {} ({}) to {}'.format(x, units, 'ppbC'   )
            units = 'ppbC'

    if ClearFlo_unit:
        if any( [ x == i for i in  [ 'NO', 'MACR', 'MVK' ] ] ):
            units = 'ppbv'
        if any( [ x == i for i in  [ 'PAN' ] ] ):
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

# --------------
# 4.09 -  Ox in species
# -------------
def Ox_in_species(in_=None, rxns=False, keys=False):
    species_Ox = {
    'HOIdf': 1.0, 'OIOdf': 2.0, 'BrNO3df': 2.0, 'HNO3df': 1.0, 'PPNdf': 1.0, 'IOdf': 1.0, 'N2O5df': 3.0, 'IONOdf': 1.0, 'PMNdf': 1.0, 'BrNO2df': 1.0, 'I2O4df': 4, 'MPNdf': 1.0, 'NO3df': 2.0, 'BrOdf': 1.0, 'HOBrdf': 1.0, 'HNO4df': 1.0, 'O3df': 1.0, 'I2O2df': 2.0, 'NO2df': 1.0, 'IONO2df': 2.0, 'PANdf': 1.0, 'OIO': 2.0, 'BrO': 1.0, 'HOBr': 1.0, 'N2O5': 3.0, 'IONO': 1.0, 'MPN': 1.0, 'BrNO2': 1.0, 'I2O2': 2.0, 'I2O4': 4, 'PPN': 1.0, 'HOI': 1.0, 'HNO3': 1.0, 'IONO2': 2.0, 'NO2': 1.0, 'IO': 1.0, 'HNO4': 1.0, 'PMN': 1.0, 'O3': 1.0, 'BrNO3': 2.0, 'PAN': 1.0, 'NO3': 2.0
    }
    rxn_Ox     = {
    'LO3_18': 1.0, 'LR25': 1.0, 'RD12': 2.0, 'LR21': 1.0, 'LO3_38': 1.0, 'LO3_10': 1.0, 'LO3_34': 1.0, 'LO3_35': 1.0, 'LO3_33': 1.0, 'LO3_30': 1.0, 'LR5': 2.0, 'LR6': 2.0, 'RD37': 2.0, 'LO3_05': 1.0, 'RD11': 2.0, 'LO3_06': 1.0, 'LO3_49': 1.0, 'LO3_04': 1.0, 'LO3_03': 1.0, 'LO3_02': 1.0, 'LO3_42': 1.0, 'LO3_41': 1.0, 'LO3_40': 1.0, 'LO3_47': 1.0, 'LO3_46': 1.0, 'LO3_09': 1.0, 'LO3_44': 1.0, 'LR30': 1.0, 'LO3_24': 1.0/2.0, 'LO3_21': 1.0, 'RD23': 2.0, 'LO3_54': 2.0, 'LO3_55': 1.0, 'LO3_08': 1.0, 'LO3_50': 1.0, 'LO3_51': 1.0, 'LO3_52': 1.0, 'LO3_53': 1.0, 'LR10': 1.0, 'LO3_36':1.0,
    # LO3_24 set to 1 (as 0.5*CoE) even though 2 Ox equivalents are lost, this allows for contribution to bromine and iodine loss to be inclued
    # LOX included for processing ease 
    'LOX':1.0, 'POX':1.0, 'PO3_14': 1.0, 'PO3_15':1.0 , 'RD98': 1.0, 'LO3_39':1.0 , 'RD63': 1.0, 
    # for prod analysis
    'PO3_69' : 1.0/2.0, 'PO3_35': 0.85, 'PO3_03':0.15/0.3, 'PO3_70': 0.4/1.4 , 'PO3_77': 1.0/2.0 , 'RD06':1.0, 'LR9':1.0
    }
    if (rxns):
        return rxn_Ox[ in_ ]
    if (keys):
        return species_Ox.keys()
    else:
        return species_Ox[ in_ ]

# ----
#  4.10 - Return dictionary of gaw sites
# ----
def gaw_2_name():
    wdf = get_dir('dwd') +'ozonesurface/' + 'gaw_site_list.h5'
    df= pd.read_hdf( wdf,  'wp', mode='r' )
    names = df.values[:,1]

    # alter long name for CVO
    ind=[ n for n, i in enumerate(names) if ( i =='Cape Verde Atmospheric Observatory' )]
    names[ind[0]] = 'Cape Verde'

    return dict( zip( df.index, names ))

# ----
#  4.11 - Returns list of gaw sites in HDF file of O3 surface data
# ----
def get_global_GAW_sites(f='gaw_site_list_global.h5'):
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
# 4.14 - Get scaling for a given unit
# --------
def get_unit_scaling( units ):

        misc = 'K', 'm/s', 'unitless', 'kg' ,'m', 'm2','kg/m2/s', \
            'molec/cm2/s', 'mol/cm3/s',  'kg/s', 'hPa'

        if any( [ (units ==  i) for i in 'pptv', 'pptC' ]):
            scaleby = 1E12
        elif any( [ (units ==  i) for i in 'ppbv', 'ppbC' ]):
            scaleby = 1E9
        elif any( [units ==i for i in misc ] ):
            scaleby = 1
        else:
            print 'WARNING: This unit is not in unit lists: ', units
        return scaleby


# ------------------------------------------- Section 5 -------------------------------------------
# --------------  Misc
#

# --------------
#  5.01 - Return MUTD runs  - no hal, just Br, just I, and I + Br
# --------------
def MUTD_runs( standard=True, sensitivity=False, titles=False, \
       IO_obs=False,no_I2Ox=False, respun=True,            \
       preindustrial=False, skip3=False, v10v92comp=False,              \
       nested_EU=False, just_bcase_no_hal=False, \
       just_bcase_std=False, ver='1.6', res='4x5',    \
       debug=False): 

    if debug:
        print standard, sensitivity, titles, IO_obs, preindustrial, skip3,\
                    nested_EU, just_bcase_std, ver
    pwd = get_dir( 'rwd' )
    l_dict= GC_var('latex_run_names')

    # Get version directories for versions
    d= { 
    '4x5': { 
    '1.6':'iGEOSChem_1.6_G5/',  
#    '1.6':'iGEOSChem_1.6.1_G5/',  
    '1.6.1':'iGEOSChem_1.6.1_G5/',  
    '1.5': 'iGEOSChem_1.5_G5/' ,
    '1.7': 'iGEOSChem_1.7_v10/'}, 
    '2x2.5': {
    '1.6':'iGEOSChem_1.6_G5_2x2.5/' 
    } }[res]
    d=d[ver]

    if sensitivity:
        l = [ 
        'no_hal', 'Just_Br',  'just_I', 'run', 'Just_I_org',  'I2Ox_double',\
        'I2Ox_half', 'het_double','het_half', 'no_het',  'Sulfate_up', \
        'MacDonald_iodide', 'I2Ox_phot_x2','I2Ox_phot_exp', 'no_I2Ox', \
        'BrO2pptv' ]                                     
        r = [ d + i for i in l ]

        # adjust names to inc. repsin
#        if respun:
#            r = [ i+'.respun' for i in r ]

    if standard and (not any( [preindustrial,sensitivity, v10v92comp]) ):
        if just_bcase_std:
            l = [ 'Just_Br',  'run' ]
#            l = [ 'Just_Br',  'no_I2Ox' ]
        elif just_bcase_no_hal:
            l = [ 'no_hal',  'run' ]
#            l = [ 'no_hal',  'no_I2Ox' ]
        else:
            if ver == '1.7':
#                l = ['no_hal', 'run' ]
                l = ['Just_Br', 'run' ]
            else:
                l = ['no_hal', 'Just_Br', 'just_I', 'run' ]
#        if any( [ (ver ==i) for i in  '1.5', '1.6' ] ) :
        r = [ d + i for i in  l ]

    if nested_EU:
        d= 'iGEOS_Chem_1.6_G5_NPOINTS/'
        if just_bcase_std:
            l = [ 'Just_Br',  'run' ]            
        else:
            l = [ 'no_hal',  'run' ]   
        r = [ d+i for  i in l ]

    # Setup latex titles list
    if titles and (not any( [preindustrial, v10v92comp]) ) and ( standard or \
        sensitivity or nested_EU):
        l = [ l_dict[i] for i in l ]       

    if v10v92comp:
        if just_bcase_std:
            l = [ 'Just_Br',  'run' ]            
        else:
            l = [ 'no_hal',  'run' ]        

        r= [ 'iGEOSChem_1.7_v10/'+ i for i in l ]+ \
             [ 'iGEOS_Chem_1.6_G5/'+ i for i in l ]
        l= [ 'GEOS-Chem v10 (no hal)', 'Iodine Sim. (v10)']  + \
            [ 'GEOS-Chem v9-2 (no hal)', 'Iodine Sim.'] 

    if IO_obs:
        l = 'run_CVO', 'run_GRO', 'run_MAL', 'run_HAL', 'run_TOR_updated.respun' 
        r = [d+ i for i in l]
        if no_I2Ox:
            r = [ i+'_no_I2Ox' for i in r ]

    if preindustrial:
            r = [ 
            'iGEOS_Chem_1.5_G5/no_hal', 'iGEOS_Chem_1.5_G5/run' ,   \
            'iGEOS_Chem_1.5_G5_preindustrial_no_hal/no_hal',    \
            'iGEOS_Chem_1.5_G5_preindustrial/run'  
            ]

            l =  [
            '(I-,Br-)', '(I+,Br+)' , '(I-,Br-) - 1750',  '(I+,Br+) - 1750'  
            ]

    if debug:                         
        print [pwd + i for i in r ], l

    if debug:
        [pwd + i for i in r ]
    if skip3:
        [ [ i.pop(0) for i in l, r ] for ii in range(3) ]
    if titles:
        return [pwd + i for i in r ], l
    else:
        return [pwd + i for i in r ]   

# --------------
# 5.02 - Store of  constants for use by funcs/progs
# --------------
def constants(input_x, debug=False):
    con_dict ={
    'RMM_air' : ( .78*(2.*14.)+.22*(2.*16.) )  ,
    'AVG' : 6.0221413E23, 
    'mol2DU': 2.69E20
    }
    return con_dict[input_x]


# ------------------------------------------- Section 6 -------------------------------------------
# -------------- Dynamic processing of p/l
#
    
# -------------
# 6.01 - Extract reactions to form a dictionary of active reactions (can be used as an external call)
# ------------- 
def rxn_dict_from_smvlog( wd, PHOTOPROCESS=None, ver='1.7', LaTeX=False, debug=False ):
    """ build a dictionary reaction of reaction details from smv.log"""
    
    if isinstance( PHOTOPROCESS, type(None) ):
        PHOTOPROCESS = {
        '1.6' : 457, '1.7' : 467
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
        for rxn in sorted(rdict.keys() ):
        
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
                    rxn_str = rxn_str.replace('++=', xarstr).replace('+', ' + ').replace('+ =', r' $\rightarrow$ ' )
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
                print '!'*100, 'ERROR HERE: {} {}'.format(  rxn, rxn_str )
        rdict = dict( zip(rxns, rxn_strs ) )
    return rdict
    
# -------------
# 6.02 - Extract reactions tracked by prod loss diag for a given p/l family
# ------------- 
def rxns_in_pl( wd, spec='LOX', debug=False ):
    fn =  'smv2.log'
    file_ =  open( wd+'/'+fn, 'rb' )
#    if debug:
    if True:
        print file_
    readrxn  = False
    for row in file_:
        row = row.split()
        if all( [ i in row for i in 'Family','coefficient' , 'rxns', spec] ):
            readrxn=True
        if len(row) < 1 :
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
    n = [int(rxn[1]) for rxn in rxns ]
    rxns = [rxn[2:] for rxn in rxns ]

    rdict = dict( zip(n, rxns) )
    return rdict

# ------------- 
# 6.03 - Extract reaction infomation for given tracer tags
# ------------- 
def rxn4pl( pls, wd='example/example', rdict=None, reduce_size=True, \
            debug=False ):
    # ---  Get Dict of reaction detail
    if debug:
        print 'rxn4pl called'
    if isinstance(rdict, type(None) ):
        #   Get Dict of all reactions, Keys = #s
        rdict = rxn_dict_from_smvlog(wd) 
        

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
         NOy_as_HOx=True, debug=False ):
    
    # assign family
#    famsn = [ 'Photolysis','HOx','NOy' ,'Bromine', 'Iodine' ]
    famsn = [ 'Photolysis','HOx' ,'Bromine', 'Iodine' ]
    fams = []
    for tag in tags:
        fams.append( get_tag_fam( tag) )
    # if rm NOy family (treat as NOx)
    if NOy_as_HOx:
        fams = [x if (x!='NOy') else 'HOx' for x in fams]

    fd = dict(zip(tags, fams) )

    ll =  []
    [ ll.append([]) for i in famsn ]
    for n, tag in enumerate( tags ) :
        for fn in range(len(famsn)):
            if fd[tag] == famsn[fn] :
                ll[fn].append( n)
    if fam:
        # Kludge - to allow for counting Ox loss via IO +BrO 50/50, 
        # ( add extra loss tag. )
        if IO_BrOx2:
            ll[famsn.index('Bromine')].append(  max([max(i) for i in ll])+1 )
            fams = fams +[ 'Bromine' ] 
            tags = tags +  [ 'LO3_24'] 
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

    # (PD??, RD??, LO3_??, PO3_??, LR??)
    for rxn in rxns:
        if debug:
            print [ i for i in rxn if  any( [x in i for x in 'PD', 'RD', 'PO3','LO3' , 'LR' ]) ]
        tags =  [i for i in rxn if any( [x in i for x in 'PD', 'RD', 'PO3','LO3' , 'LR' ]) ]
        try:
            tagsl.append( tags)
        except:
            tagsl = [tags]

    return tagsl

# -------------
# 6.06 - extract reactions tracked by prod loss diag in input.geos
# ------------- 
def p_l_species_input_geos( wd, ver='1.7', rm_multiple_tagged_rxs=False ):
    """ 
    Extract prod/loss species (input.geos) and reaction tags (globchem.dat) 
    """
    # find and open input.geos file
    fn = glob.glob(wd+'/*input.geos*')[0]
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

    # remove p/l with muliple values ( start from 12th input) - Kludge?
    if rm_multiple_tagged_rxs:
        PD, vars = [ i[11:] for i in PD, vars ]
        vars =  [ i[0] for i in  vars ]
        
    return PD, vars

# -------------
# 6.07 - extract all active tags from smv.log
# ------------- 
def tags_from_smvlog( wd ): #, spec='LOX' ):
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
    fn =  'smv2.log'
    file_ =  open( wd+'/'+fn, 'rb' )
    readrxn  = False
    leniency =  0
    for row in file_:
        row = row.split()
        if  all( [(i in row) for i in ['Families', 'for', 'prod', 'or', 'loss', 'output:'] ]):
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
    # --- get reaction dictionary 
    if rdict == None:
        rdict =  rxn_dict_from_smvlog( wd, ver=ver )
    
    # --- Caveats -  to adapt for long line errors in fortran written output
    errs = ['LO3_36']
    cerrs = ['RD95'] 
    if any([ (tag == i) for i in  errs ] ):
        tag  = cerrs[  errs.index( tag) ]
    
    # -- loop reactions, if tag in reaction return reaction
    rxns = []
    for n, rxn in enumerate( rdict.values() ):
#        if any( [tag in i for i in rxn]):
        if any( [(i.endswith(tag) ) for i in rxn]):
            rxns.append(  [rdict.keys()[n] ]+ rxn )
    return rxns

# -------------
# 6.10 - get details for a given tag
# ------------- 
def get_tag_details( wd, tag=None, PDs=None,  rdict=None, \
            PHOTOPROCESS=None, ver='1.7', LaTeX=False, print_details=False, \
            debug=False ):
    """ Retriveve prod/loss tag details from smv.log
        ( rxn number + reaction description)"""

    # what is the number of the first photolysis reaction?
    if isinstance( PHOTOPROCESS, type(None) ):
        PHOTOPROCESS = {
        '1.6' : 457, '1.7' : 467
        }[ver]

    # ---  get all reactions tags are active in smv.log    
    if rdict == None:
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
        print '!'*100, 'ERROR HERE: {} {}'.format(  tag, trxns  )
    if print_details:
        print dets
    # --- return a dictionary of all active tagged reactions details by tag : PD, number, rxn str, coeeffiecn
    else:
        return dets
        
# -------------
# 6.11 - Takes a reaction number add gives the Ox Coe
# ------------- 
def get_rxn_Coe(wd, num, tag, nums=None, rxns=None, tags=None, \
            Coe=None, spec='LOX', debug=False):

    # --- get dictionaries for reactions within
    if all( [ (i == None) for i in nums, rxns, tags, Coe ] ):
        nums, rxns, tags, Coe = prod_loss_4_spec( wd,  spec, all_clean=True )

    # Pull reaction coefficient 
    Coe_dict = dict(zip(nums, Coe) )             
    Coe = float(Coe_dict[ num ])
    if ('P' not in tag ): # consider all change positive - Kludge due to the assignment approach.
        Coe =  Coe*-1.0

    # --- over write Ox in species with prod_mod_tms values
    try:
        Coe = Ox_in_species(tag, rxns=True)
        if debug:
            print 'using values from Ox_in_species'
    except:
        if debug:
            print 'using values from smv.log @: {}'.format(wd)
    return Coe

# ------------------------------------------- Section 7 -------------------------------------------
# -------------- Observational Variables
#


# 5.02 - obs sites (e.g. deposition locs, lons, alts )
# 5.04 - Store site locations ( LAT,  LON,  ALT)

# --------------
#  7.01 - open ocean IO data  
# --------------
""" key = ref (e.g. name_year )     values =  alt, lon, lat, times, full_name  , IO (avg), BrO (avg), CAM-Chem IO, CAM-Chem BrO, group
"""
def IO_obs_data(just_IO=False):
    IO_obs_dict={'Read_2008': [0.0, -24.87, 16.85, 6.0, 'Read, 2008', '1', '2', '1', '2', 'Leeds'], 'Weddel_Sea': [0.03, -50.0, 75.0, 10.0, 'Weddel Sea', '-', '-', '-', '-', '-'], 'Oetjen_2009': [0.0, 73.5, 5.0, 2.0, 'Oetjen, 2009', '2.4', '-', '1', '-', 'Leeds'], 'Jones_2010_II': [0.0, -19.0, 30.0, 7.0, 'Jones, 2010, RHaMBLe II (26-36N)', '-', '-', '-', '-', '-'], 'Leser_2003': [0.0, -24.87, 16.85, 10.0, 'Leser, 2003', '-', '3.6', '-', '0.8', 'Heidelberg'], 'Jones_2010_I': [0.0, -19.0, 20.0, 7.0, 'Jones, 2010, RHaMBLe I (15-25N)', '-', '-', '-', '-', '-'], 'Halley_BAS': [0.03, -26.34, 75.35, 10.0, 'Halley, BAS', '-', '-', '-', '-', '-'], 'Theys_2007': [0.0, -20.9, -55.5, 8.0, 'Theys, 2007', '-', '<0.5', '-', '0.8', 'Belgium'], 'Dix_2013': [9.0, -160.0, 10.0, 1.0, 'Dix, 2013', '0.1', '-', '', '-', 'NCAR'], 'Allan_2000': [0.162, -16.6, 28.4, 6.0, 'Allan, 2000', '1.2', '-', '0.4', '-', 'Leeds'], 'Jones_2010_MAP': [0.0, -10.0, 55.0, 6.0, 'Jones, 2010, MAP', '-', '-', '-', '-', '-'], 'Grobmann_2013': [0.0, 150.0, 15.0, 1.0, 'Grobmann, 2013', '0.72-1.8', '-', '-', '-', 'Heidelberg'], 'Dix_2013_j_comp': [9.0, -160.0, 10.0, 5.0, 'Dix, 2013 (wrong month to allow with jones runs comparisons)', '-', '-', '-', '-', '-'], 'schonart_2009_II': [0.0, -80.0, -20.0, 10.0, 'Schonart, 2009 (satilitte)', '3.3', '-', '1', '-', 'Bremen'], 'Martin_2009': [0.0, -24.87, 16.85, 2.0, 'Martin, 2009', '-', '<3.0', '-', '1.2', 'Heidelberg'], 'schonart_2009_I': [0.0, -90.0, -5.0, 10.0, 'Schonart, 2009 (satilitte)', '3.3', '-', '1', '-', 'Bremen'], 'Yokouchi_2013': [0.0, 120.0, 15.0, 2.0, 'Yokouchi, 2013', '-', '-', '-', '-', '-'], 'Butz_2009': [0.045, -44.4, -2.4, 12.0, 'Butz, 2009', '0.1', '~1.0', '0.02', '0.5', 'Leeds']}
    IO_obs_dict_just = {'Read_2008': [0.0, -24.87, 16.85, 6.0, 'Read, 2008', '1', '2', '1', '2', 'Leeds'], 'schonart_2009_II': [0.0, -80.0, -20.0, 10.0, 'Schonart, 2009 (satilitte)', '3.3', '-', '1', '-', 'Bremen'], 'Oetjen_2009': [0.0, 73.5, 5.0, 2.0, 'Oetjen, 2009', '2.4', '-', '1', '-', 'Leeds'], 'Allan_2000': [0.162, -16.6, 28.4, 6.0, 'Allan, 2000', '1.2', '-', '0.4', '-', 'Leeds'], 'Leser_2003': [0.0, -24.87, 16.85, 10.0, 'Leser, 2003', '-', '3.6', '-', '0.8', 'Heidelberg'], 'Theys_2007': [0.0, -20.9, -55.5, 8.0, 'Theys, 2007', '-', '<0.5', '-', '0.8', 'Belgium'], 'Dix_2013': [9.0, -160.0, 10.0, 1.0, 'Dix, 2013', '0.1', '-', '', '-', 'NCAR'], 'Grobmann_2013': [0.0, 135.0, 15.0, 1.0, 'Grobmann, 2013', '1.095', '-', '-', '-', 'Heidelberg'], 'schonart_2009_I': [0.0, -90.0, -5.0, 10.0, 'Schonart, 2009 (satilitte)', '3.3', '-', '1', '-', 'Bremen'], 'Martin_2009': [0.0, -24.87, 16.85, 2.0, 'Martin, 2009', '-', '<3.0', '-', '1.2', 'Heidelberg'], 'Butz_2009': [0.045, -44.4, -2.4, 12.0, 'Butz, 2009', '0.1', '~1.0', '0.02', '0.5', 'Leeds'], 'Mahajan_2010': [0., -100, -10, 4.0, 'Mahajan, 2010', '0.58', '-', '-', '-', 'Leeds'] }

    if (just_IO):
        return IO_obs_dict_just
    else:
        return IO_obs_dict

# --------------
# 7.02 - BAE  flight ID to date
# -------------
def bae_flight_ID_2_date( f_ID, debug=False):
    if (debug):
        print 'bae_flight_ID_2_date called'
    f_ID_dict = {'847':'2014-02-18','846':'2014-02-17','845':'2014-02-17', '844':'2014-02-16','843':'2014-02-15','842':'2014-02-17','840':'2014-02-13','839':'2014-02-12','838':'2014-02-05','837':'2014-02-04','836':'2014-02-04','835':'2014-02-03','834':'2014-02-01','833':'2014-02-01','832':'2014-01-30','831':'2014-01-30','830':'2014-01-30','829':'2014-01-28','828':'2014-01-26','827':'2014-01-26','826':'2014-01-25','825':'2014-01-24','824':'2014-01-21','823':'2014-01-18'}
    return f_ID_dict[f_ID]

# --------------
# 7.03 - details on flights from FAAM's NETCDF core data files, 
# --------------
# removed 'CAST_flight, 'b822': '20140108''
def CAST_flight(all=True, CIMS=False, CIMSII=False):
    flight_dict={'b828': ['20140126', '0559', '20140126', '0930', 1390715993, 1390728623], 'b829': ['20140128', '1932', '20140129', '0224', 1390937523, 1390962280], 'b824': ['20140121', '1956', '20140122', '0413', 1390334165, 1390363990], 'b825': ['20140124', '1940', '20140125', '0140', 1390592404, 1390614058], 'b826': ['20140125', '0128', '20140125', '0637', 1390613301, 1390631839], 'b827': ['20140126', '0041', '20140126', '0445', 1390696903, 1390711511], 'b823': ['20140118', '1713', '20140119', '0147', 1390065184, 1390096069], 'b839': ['20140212', '0320', '20140212', '0824', 1392175209, 1392193471], 'b838': ['20140205', '2246', '20140206', '0545', 1391640363, 1391665549], 'b837': ['20140204', '1901', '20140205', '0318', 1391540463, 1391570305], 'b836': ['20140204', '0134', '20140204', '0627', 1391477699, 1391495264], 'b835': ['20140203', '2101', '20140204', '0104', 1391461263, 1391475894], 'b834': ['20140201', '0610', '20140201', '1217', 1391235033, 1391257044], 'b833': ['20140201', '0045', '20140201', '0606', 1391215547, 1391234764], 'b832': ['20140130', '0622', '20140130', '1131', 1391062943, 1391081485], 'b831': ['20140130', '0101', '20140130', '0614', 1391043666, 1391062492], 'b830': ['20140129', '0239', '20140129', '0830', 1390963178, 1390984251], 'b846': ['20140217', '2058', '20140218', '0510', 1392670683, 1392700201], 'b847': ['20140218', '0503', '20140218', '0857', 1392699809, 1392713860], 'b844': ['20140216', '1857', '20140217', '0324', 1392577026, 1392607467], 'b845': ['20140217', '0323', '20140217', '0753', 1392607430, 1392623590], 'b842': ['20140214', '0508', '20140214', '1035', 1392354530, 1392374130], 'b843': ['20140215', '1903', '20140216', '0400', 1392490986, 1392523256], 'b840': ['20140213', '0005', '20140213', '0658', 1392249904, 1392274686], 'b841': ['20140214', '0004', '20140214', '0458', 1392336265, 1392353901]}

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
# 7.05 - Stores locations for use by funcs/progs - LON, LAT,   ALT - double up? ( with 5.02 ?? )
# --------------
def get_loc(loc=None, rtn_dict=False, debug=False):
    loc_dict ={    
    'GUAM' :  ( 144.800, 13.500, 0  ),
    'CHUUK' : ( 151.7833, 7.4167, 0 ),
    'PILAU' : ( 134.4667,7.3500, 0 ),      
    'London': ( -0.1275, 51.5072, 0 ),
    'Weyborne' : (1.1380, 52.9420,  0 ), 
#    'Cape Verde': (16.848, -24.871, 0 ), 
#    'CVO': (16.848, -24.871, 0 ), 
    'Cape Verde': ( -24.871, 16.848, 0 ), 
    'CVO': (-24.871,16.848,  0 ), 
    'North Ken' :  (-0.214174, 51.520718, 0),
    'KEN' :  (-0.214174, 51.520718, 0),
    'BT tower': (-0.139055, 51.521556, 190),
    'BTT': (-0.139055, 51.521556, 190)
    }
    if rtn_dict:
        return loc_dict
    else:
        return loc_dict[loc]

# --------------
# 7.06 - Get Locations of observations (lats, lons, alts ) for given sites
# --------------
def get_obs_loc(loc, debug=False):
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
     7.8894] ],   
    'Weyborne' :[  [ 52.9420],  [1.1380 ]   ]    }
    return d[loc] 
    

    
# --------------
# 7.07 - sonde station variables (list of 432 sondes)
# -------------
def sonde_STNs():
    sonde_dict = {
    101: ['SYOWA', 101.0, -69.0, 39.58, 22.0, 'JPN', 'ANTARCTICA'], 104: ['BEDFORD', 104.0, 42.45, -71.267, 80.0, 'USA', 'IV'], 105: ['FAIRBANKS (COLLEGE)', 105.0, 64.817, -147.867, 138.0, 'USA', 'IV'], 107: ['WALLOPS ISLAND', 107.0, 37.898, -75.483, 13.0, 'USA', 'IV'], 108: ['CANTON ISLAND', 108.0, -2.76, -171.7, 3.0, 'USA', 'V'], 109: ['HILO', 109.0, 19.5735, -155.0485, 11.0, 'USA', 'V'], 111: ['AMUNDSEN-SCOTT (SOUTH POLE)', 111.0, -89.983, 0.0, 2820.0, 'ATA', 'ANTARCTICA'], 131: ['PUERTO MONTT', 131.0, -41.45, -72.833, 5.0, 'CHL', 'III'], 132: ['SOFIA', 132.0, 42.817, 23.383, 588.0, 'BGR', 'VI'], 137: ['TOPEKA', 137.0, 39.067, -95.633, 270.0, 'USA', 'IV'], 138: ['CHRISTCHURCH', 138.0, -43.483, 172.55, 34.0, 'NZL', 'V'], 149: ['OVEJUYO (LA PAZ)', 149.0, -16.517, -68.033, 3420.0, 'BOL', 'III'], 156: ['PAYERNE', 156.0, 46.49, 6.57, 491.0, 'CHE', 'VI'], 157: ['THALWIL', 157.0, 46.817, 8.455, 515.0, 'CHE', 'VI'], 163: ['WILKES', 163.0, -66.25, 110.517, 12.0, 'USA', 'ANTARCTICA'], 174: ['LINDENBERG', 174.0, 52.21, 14.12, 112.0, 'DEU', 'VI'], 175: ['NAIROBI', 175.0, -1.267, 36.8, 1745.0, 'KEN', 'I'], 181: ['BERLIN/TEMPLEHOF', 181.0, 52.467, 13.433, 50.0, 'DEU', 'VI'], 187: ['PUNE', 187.0, 18.553, 73.86, 559.0, 'IND', 'II'], 190: ['NAHA', 190.0, 26.2, 127.683, 27.0, 'JPN', 'II'], 191: ['SAMOA', 191.0, -14.25, -170.56, 82.0, 'ASM', 'V'], 194: ['YORKTON', 194.0, 51.263, -102.467, 504.0, 'CAN', 'IV'], 197: ['BISCARROSSE/SMS', 197.0, 44.367, -1.233, 18.0, 'FRA', 'VI'], 198: ['COLD LAKE', 198.0, 54.783, -110.05, 702.0, 'CAN', 'IV'], 199: ['BARROW', 199.0, 71.317, -156.635, 11.0, 'USA', 'IV'], 203: ['FT. SHERMAN', 203.0, 9.33, -79.983, 57.0, 'PAN', 'IV'], 205: ['THIRUVANANTHAPURAM', 205.0, 8.483, 76.97, 60.0, 'IND', 'II'], 206: ['BOMBAY', 206.0, 19.117, 72.85, 145.0, 'IND', 'II'], 210: ['PALESTINE', 210.0, 31.8, -95.717, 121.0, 'USA', 'IV'], 213: ['EL ARENOSILLO', 213.0, 37.1, -6.733, 41.0, 'ESP', 'VI'], 217: ['POKER FLAT', 217.0, 65.133, -147.45, 357.5, 'USA', 'IV'], 219: ['NATAL', 219.0, -5.71, -35.21, 30.5, 'BRA', 'III'], 221: ['LEGIONOWO', 221.0, 52.4, 20.967, 96.0, 'POL', 'VI'], 224: ['CHILCA', 224.0, -12.5, -76.8, -1.0, 'PER', 'III'], 225: ['KOUROU', 225.0, 5.333, -52.65, 4.0, 'GUF', 'III'], 227: ['MCDONALD OBSERVATORY', 227.0, 30.666, -90.933, 2081.0, 'USA', 'IV'], 228: ['GIMLI', 228.0, 50.633, -97.05, 228.0, 'CAN', 'IV'], 229: ['ALBROOK', 229.0, 8.983, -79.55, 66.0, 'PAN', 'IV'], 231: ['SPOKANE', 231.0, 47.667, -117.417, 576.0, 'USA', 'IV'], 233: ['MARAMBIO', 233.0, -64.233, -56.623, 196.0, 'ATA', 'ANTARCTICA'], 234: ['SAN JUAN', 234.0, 18.483, -66.133, 17.0, 'PRI', 'IV'], 235: ['LONG VIEW', 235.0, 32.5, -94.75, 103.0, 'USA', 'IV'], 236: ['COOLIDGE FIELD', 236.0, 17.283, -61.783, 10.0, 'ATG', 'IV'], 237: ['GREAT FALLS', 237.0, 47.483, -111.35, 1118.0, 'USA', 'IV'], 238: ['DENVER', 238.0, 39.767, -104.883, 1611.0, 'USA', 'IV'], 239: ['SAN DIEGO', 239.0, 32.76, -117.19, 72.5, 'USA', 'IV'], 242: ['PRAHA', 242.0, 50.02, 14.45, 304.0, 'CZE', 'VI'], 254: ['LAVERTON', 254.0, -37.867, 144.75, 21.0, 'AUS', 'V'], 255: ['AINSWORTH (AIRPORT)', 255.0, 42.583, -100.0, 789.0, 'USA', 'IV'], 256: ['LAUDER', 256.0, -45.03, 169.683, 370.0, 'NZL', 'V'], 257: ['VANSCOY', 257.0, 52.115, -107.165, 510.0, 'CAN', 'IV'], 260: ['TABLE MOUNTAIN (CA)', 260.0, 34.4, -117.7, 2286.0, 'USA', 'IV'], 262: ['SODANKYLA', 262.0, 67.335, 26.505, 179.0, 'FIN', 'VI'], 265: ['IRENE', 265.0, -25.91, 28.211, 1524.0, 'ZAF', 'I'], 280: ['NOVOLASAREVSKAYA / FORSTER', 280.0, -70.767, 11.867, 110.0, 'ATA', 'ANTARCTICA'], 297: ['S.PIETRO CAPOFIUME', 297.0, 44.65, 11.617, 11.0, 'ITA', 'VI'], 303: ['IQALUIT', 303.0, 63.75, -68.55, 20.0, 'CAN', 'IV'], 308: ['MADRID / BARAJAS', 308.0, 40.46, -3.65, 650.0, 'ESP', 'VI'], 315: ['EUREKA / EUREKA LAB', 315.0, 80.04, -86.175, 310.0, 'CAN', 'IV'], 316: ['DE BILT', 316.0, 52.1, 5.18, 4.0, 'NLD', 'VI'], 318: ['VALENTIA OBSERVATORY', 318.0, 51.93, -10.25, 14.0, 'IRL', 'VI'], 323: ['NEUMAYER', 323.0, -70.65, -8.25, 42.0, 'ATA', 'ANTARCTICA'], 328: ['ASCENSION ISLAND', 328.0, -7.98, -14.42, 91.0, 'SHN', 'I'], 329: ['BRAZZAVILLE', 329.0, -4.28, 15.25, 314.0, 'COG', 'I'], 330: ['HANOI', 330.0, 21.033, 105.84, 5.0, 'VNM', 'II'], 333: ['PORTO NACIONAL', 333.0, -10.8, -48.4, 240.0, 'BRA', 'III'], 334: ['CUIABA', 334.0, -15.6, -56.1, 990.0, 'BRA', 'III'], 335: ['ETOSHA PAN', 335.0, -19.2, 15.9, 1100.0, 'NAM', 'I'], 336: ['ISFAHAN', 336.0, 32.477, 51.425, 1550.0, 'IRN', 'II'], 338: ['BRATTS LAKE (REGINA)', 338.0, 50.205, -104.705, 592.0, 'CAN', 'IV'], 339: ['USHUAIA', 339.0, -54.85, -68.308, 15.0, 'ARG', 'III'], 344: ['HONG KONG OBSERVATORY', 344.0, 22.31, 114.17, 66.0, 'HKG', 'II'], 348: ['ANKARA', 348.0, 39.95, 32.883, 896.0, 'TUR', 'VI'], 360: ['PELLSTON (MI)', 360.0, 45.56, -84.67, 238.0, 'USA', 'IV'], 361: ['HOLTVILLE (CA)', 361.0, 32.81, -115.42, -18.0, 'USA', 'IV'], 394: ['BROADMEADOWS', 394.0, -37.6914, 144.9467, 108.0, 'AUS', 'V'], 400: ['MAITRI', 400.0, -70.46, 11.45, 223.5, 'ATA', 'ANTARCTICA'], 401: ['SANTA CRUZ', 401.0, 28.42, -16.26, 36.0, 'ESP', 'I'], 404: ['JOKIOINEN', 404.0, 60.81, 23.5, 103.0, 'FIN', 'VI'], 406: ['SCORESBYSUND', 406.0, 70.49, -21.98, 50.0, 'GRL', 'VI'], 418: ['HUNTSVILLE', 418.0, 34.72, -86.64, 196.0, 'USA', 'IV'], 420: ['BELTSVILLE (MD)', 420.0, 39.02, -76.74, 64.0, 'USA', 'IV'], 432: ['PAPEETE (TAHITI)', 432.0, -18.0, -149.0, 2.0, 'PYF', 'V'], 434: ['SAN CRISTOBAL', 434.0, -0.92, -89.6, 8.0, 'ECU', 'III'], 435: ['PARAMARIBO', 435.0, 5.81, -55.21, 22.5, 'SUR', 'III'], 436: ['LA REUNION ISLAND', 436.0, -20.99, 55.48, 61.5, 'REU', 'I'], 437: ['WATUKOSEK (JAVA)', 437.0, -7.57, 112.65, 50.0, 'IDN', 'V'], 438: ['SUVA (FIJI)', 438.0, -18.13, 178.315, 6.0, 'FJI', 'V'], 439: ['KAASHIDHOO', 439.0, 5.0, 73.5, 1.0, 'MDV', 'V'], 441: ['EASTER ISLAND', 441.0, -27.17, -109.42, 62.0, 'CHL', 'III'], 443: ['SEPANG AIRPORT', 443.0, 2.73, 101.7, 17.0, 'MYS', 'V'], 444: ['CHEJU', 444.0, 33.5, 126.5, 300.0, 'KOR', 'II'], 445: ['TRINIDAD HEAD', 445.0, 40.8, -124.16, 55.0, 'USA', 'IV'], 448: ['MALINDI', 448.0, -2.99, 40.19, -6.0, 'KEN', 'I'], 450: ['DAVIS', 450.0, -68.577, 77.973, 16.0, 'ATA', 'ANTARCTICA'], 456: ['EGBERT', 456.0, 44.23, -79.78, 253.0, 'CAN', 'IV'], 457: ['KELOWNA', 457.0, 49.93, -119.4, 456.0, 'CAN', 'IV'], 458: ['YARMOUTH', 458.0, 43.87, -66.1, 9.0, 'CAN', 'IV'], 459: ['TBD', 459.0, 0.0, 0.0, 0.0, '', 'VI'], 460: ['THULE', 460.0, 76.53, -68.74, 57.0, 'GRL', 'VI'], 466: ['MAXARANGUAPE (SHADOZ-NATAL)', 466.0, -5.445, -35.33, 32.0, 'BRA', 'III'], 472: ['COTONOU', 472.0, 6.21, 2.23, 10.0, 'BEN', 'I'], 477: ['HEREDIA', 477.0, 10.0, -84.11, 1176.0, 'CRI', 'IV'], 480: ['SABLE ISLAND', 480.0, 43.93, -60.02, 4.0, 'CAN', 'IV'], 482: ['WALSINGHAM', 482.0, 42.6, -80.6, 200.0, 'CAN', 'IV'], 483: ['BARBADOS', 483.0, 13.16, -59.43, 32.0, 'BRB', 'III'], 484: ['HOUSTON (TX)', 484.0, 29.72, -95.4, 19.0, 'USA', 'IV'], 485: ['TECAMEC (UNAM)', 485.0, 19.33, -99.18, 2272.0, 'MEX', 'IV'], 487: ['NARRAGANSETT', 487.0, 41.49, -71.42, 21.0, 'USA', 'IV'], 488: ['PARADOX', 488.0, 43.92, -73.64, 284.0, 'USA', 'IV'], 489: ['RICHLAND', 489.0, 46.2, -119.16, 123.0, 'USA', 'IV'], 490: ['VALPARAISO (IN)', 490.0, 41.5, -87.0, 240.0, 'USA', 'IV'], 494: ['ALAJUELA', 494.0, 9.98, -84.21, 899.0, 'CRI', 'IV']
    }
    return sonde_dict 

# ----
#  7.08 - returns  (lat, lon, alt (press), timezone (UTC) ) for a given site
# ----
def gaw_2_loc(site,  f =  'GLOBAL_SURFACE_O3_2006_2012.nc'):#, f 
    """
    
    - Another file is availible with just GAW sites:
          'GAW_SURFACE_O3_2006_2012.nc'  
    """
#    from AC_tools.funcs4generic import hPa_to_Km
    try:
        gaw_sites= {
        'SMO': (-14.247,  -170.565,1002.7885270480558, -11), 'MNM':(24.285, 153.981, 1011.9342452324959, 9), 'BMW':(32.27, -64.88, 1008.6109830510485, -4 ) ,'CVO': (16.848, -24.871, 1011.6679817831093, -1), 'RPB':(13.17000, -59.43000, 1007.0196960034474, -4 ), 'ogasawara': (26.38, 142.10,996.08181619552602, 9 ), 'OGA': (26.38, 142.10,996.08181619552602, 9 ) ,

        # Add extras for ease of analysis (e.g. Roscoff ... )
        'ROS': (48.433, -3.5904, 1011.6679817831093, +1)       
        }
        return gaw_sites[ site ] 
    except:
        wd= get_dir('dwd') +'ozonesurface/' 
        with  Dataset(wd+f, 'r', format='NETCDF4') as f:
            lon= f.groups[site].longitude
            alt =  f.groups[site].altitude /1E3
            lat =  f.groups[site].latitude
            print [ (i, type(i) ) for i in lat, lon, alt ]
        return (lat, lon, float( hPa_to_Km([alt], reverse=True)[0] ), -9999 )
