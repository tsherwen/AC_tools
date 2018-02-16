#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Automatically tag reactions in a GEOS-Chem (KPP) mechansim for a given family


Notes
-------
 - This script using AC_tools with contains functions to extract data from KPP mechanism files. This utilizes the fact that:
     - KPP mechanisms are contructed and information can be extracted from both the produced files (e.g. gckpp_Monitor.F90) and those the mechanism is constructed from (e.g. *.eqn)
     - When using P/L tagging, the output P/L infomation (e.g. rxns in family and their stiochmetry in gckpp_Monitor.F90) can be used to add more specific tags/for post-processing

"""
# Compatibility with both python 2 and 3
from __future__ import print_function
# Import modules used within the script
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
import AC_tools as AC
import matplotlib.pyplot as plt
import datetime as datetime
import glob


def main( folder=None, print_formatted_KPP_file=True, GC_version=None,
        verbose=True, mechanism = 'Standard', debug=False):
    """
    Parse combined KPP file (>= v11-2d) to dictionary of pd.DataFrames, then
    print these into files than can be pasted into a GEOS-Chem KPP mechanism

    Parameters
    -------
    print_formatted_KPP_file (boolean): Save the uniformly formated .eqn file
    folder (str): folder of GEOS-Chem code directory
    mechanism (str): mechanism to create tag files for ?
    GC_version (str): version of GEOS-Chem (e.g. v11-2)
    debug (boolean): print out extra infomation for debugging
    verbose (boolean): print out extra information during processing

    Returns
    -------
    (None)

    Notes
    -----
    The script below assumes this workflow:
     (1a) Make sure formatting is uniforming .eqn file (if not do step 1b)
     ( check this by seeing if diff between EXISTING_MECH_???.eqn current .eqn)
     (1b) Over write the .eqn in the mechanism directory with the new .eqn file
     (2) Compile the via build_mechanism.sh script (or "kpp gckpp.kpp")
     (3) If no compile issues, then tag the mechansim
     (4) Update the gckpp.kpp file to include the names of the new P/L tags
     ( lines outputted in 'gckpp.kpp_extra_lines_for_tagged_mech_?' file )
     (5) Run "make realclean" in the code directory and compile
     (7) Run the model with the compilie executable inc. P/L tags
     (6) P/L ouput can then be analysed by standard approaches
     ( AC_tools contains some functions that can automate this too )
    """
    # ---- Local settings
    # GEOS-chem version?
    if isinstance( GC_version, type(None) ):
        # TODO - get this online from log files
        # (already implemented in AC_Tools)
        GC_version = 'v11-2'

    # Add mechanism name to string (and KPP folder)
    folder += 'KPP/{}/'.format( mechanism )

    # ---- Process input files
    # Create a dictionary to store various (e.g. gas-phase, Heterogeneous)
    KPP_dicts = {}
    # Get the mechanism KPP file
    filename = glob.glob( folder+'/*.eqn')[0].split('/')[-1]
    # Get header lines from *.eqn file
    headers = AC.KPP_eqn_file_headers(folder=folder, filename=filename)
    # Get species and details on species as a DataFrame
    species_df = AC.KPP_eqn_file_species(folder=folder, filename=filename )
    # Get dictionaries of all reactions
    rxn_dicts = AC.get_dicts_of_KPP_eqn_file_reactions(folder=folder, \
        filename=filename )
    # Process rxns to be in dictionaries of DataFrames
    # (with extra diagnostic columns, inc. reactants, products, metadata,...)
    rxn_dicts = AC.process_KPP_rxn_dicts2DataFrames(rxn_dicts=rxn_dicts)

    # Update the numbering of indexes...
    Gas_dict = rxn_dicts['Gas-phase']
    Het_dict = rxn_dicts['Heterogeneous']
    Hv_dict = rxn_dicts['Photolysis']
    # Update index start points
    Het_dict.index = Het_dict.index + Gas_dict.shape[0]
    Hv_dict.index = Hv_dict.index + Gas_dict.shape[0] + Het_dict.shape[0]
    rxn_dicts['Heterogeneous'] = Het_dict
    rxn_dicts['Photolysis'] =  Hv_dict

    # --- Print out input KPP files with updated formatting (prior to tagging)
    # (Uniform formatting required for parsing - this step may not be required)
    if print_formatted_KPP_file:
        AC.print_out_dfs2KPP_eqn_file( headers=headers, species_df=species_df, \
            rxn_dicts=rxn_dicts, extr_str='EXISTING_MECH_{}'.format(mechanism) )

    # ---- Get outputted KPP files and process these...
    # Get outputted KPP mechanism
    KPP_output_mech = AC.get_dict_of_KPP_mech( wd=folder, GC_version=GC_version)

    # ---------------------- Tagging of Mechanism
    # Initialise dictionary to store tags used for reactions
    dict_of_tags = {}
    current_tag = 'T000'
    tag_prefix ='T'

    # --- Tag LOx reactions
    # Get tagged reactions (as a dictionary)
    fam = 'LOx'
    df_fam = AC.get_reactants_and_products4tagged_fam(folder=folder,
        KPP_output_mech=KPP_output_mech, fam=fam )

    # Loop reaction indexes for LOx family
    for key_ in rxn_dicts.keys():
        df_tmp = rxn_dicts[key_].copy()
        # Get indices of rxnes in tagged family
        rxns_in_mech = sorted(df_tmp[df_tmp.index.isin(df_fam.index)].index)
        print( len(rxns_in_mech) )
        for ix in rxns_in_mech:
            rxn_str =  df_tmp.loc[ix]['rxn_str']
            # Get a new tag and add to the reaction string
            current_tag = AC.get_KPP_PL_tag(current_tag, tag_prefix=tag_prefix)
            rxn_str += ' + '+ current_tag
            df_tmp.loc[ix, 'rxn_str'] = rxn_str
            # Save the tag for later processing
            dict_of_tags[current_tag] = {'rxn_str':rxn_str, 'fam': fam}
        # Now update the DataFrame in the rxn_dicts dictionary
        rxn_dicts[key_] = df_tmp

    # --- Add tags for Other families too?
    # ( Or just all reactions that contain a species of interest )
#     fams = 'ClNO2', 'POx'
#     for fam in fams:
#         for key_ in rxn_dicts.keys():
#             df_tmp = rxn_dicts[key_]
#             for ix in df_tmp.index:
#                 rxn_str =  df_tmp.loc[ix]['rxn_str']
#                 if (fam in rxn_str):
#                     # Get a new tag and add to the reaction string
#                     current_tag = AC.get_KPP_PL_tag(current_tag)
#                     rxn_str += ' + '+ current_tag
#                     df_tmp.loc[ix,'rxn_str'] = rxn_str
#                     # save the tag for later processing
#                     dict_of_tags[current_tag] = {'rxn_str':rxn_str, 'fam': fam}
#             # Now update the DataFrame in the rxn_dicts dictionary
#             rxn_dicts[key_] = df_tmp

    # --- Add the species to the species_df
    tags = sorted( dict_of_tags.keys() )
    fams = [ r'{'+dict_of_tags[i]['fam']+'}' for i in tags]
    df_spec_tmp = pd.DataFrame( [[False]]*len(tags), columns=['inactive'] )
    df_spec_tmp['Description'] = fams
    df_spec_tmp.index = tags
    df_spec_tmp.sort_values('Description')
    # Add to existing DataFrame
    species_df = pd.concat( [species_df, df_spec_tmp] )

    # --- Print out updated KPP .eqn file (with tags)
    AC.print_out_dfs2KPP_eqn_file( headers=headers, species_df=species_df, \
        rxn_dicts=rxn_dicts, extr_str='TAGGED_MECH_{}'.format(mechanism) )

    # --- Save out the tags and the reactions tagged
    # (to use for post-processing of tagged output)
    df = pd.DataFrame( dict_of_tags )
    savetitle = 'Tagged_reactions_in_{}'.format( mechanism )
    df.to_csv( savetitle )

    # --- Save out lines that need to be added to the gckpp.kpp file
    AC.print_out_lines_for_gckpp_file( dict_of_tags=dict_of_tags,
        extr_str=mechanism)


if __name__ == "__main__":
    main()

