from __future__ import print_function
#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Functions for KPP input/output file Processing

Notes
-------
These functions are specifically for GEOS-Chem versions v11-1g and later. They have been tested on v11-2d-R1,v11-2d-R2, v11-2d-R3, v12.0, v12.1, v12.2., v12.3 ...

"""
import os
import sys
import glob
import pandas as pd
import xarray as xr
import re
from netCDF4 import Dataset
import numpy as np
import math
import time
import calendar
import datetime as datetime
from datetime import datetime as datetime_
import logging

# AC_tools imports
from . generic import myround
from .GEOSChem_bpch import *
from .GEOSChem_nc import *


def get_fam_prod_loss4tagged_mech(wd=None, fam='LOx', ref_spec='O3',
                                  tags=None, RR_dict=None, Data_rc=None,
                                  Var_rc=None, tags2_rxn_num=None,
                                  RR_dict_fam_stioch=None, region=None,
                                  rm_strat=False,
                                  weight_by_molecs=False, verbose=True,
                                  debug=False):
    """
    Extract prod/loss for family from wd and code directory

    Parameters
    -------
    fam (str): tagged family to track (already compiled in KPP mechanism)
    ref_spec (str): reference species to normalise to
    region (str): region to consider (by masking all but this location)
    wd (str): working directory ("wd") of model output
    debug, verbose (bool): switches to turn on/set verbosity of output to screen
    RR_dict (dict): dictionary of KPP rxn. mechanism (from get_dict_of_KPP_mech)
    RR_dict_fam_stioch (dict): dictionary of stoichiometries for
        get_stioch4family_rxns
    tags2_rxn_num (dict): dicionary mapping tags to reaction numbers
    tags (list): list of prod/loss tags to extract
    Var_rc (dict): dictionary containing variables for working directory ('wd')
        ( made by AC_tools' get_default_variable_dict function)
    Data_rc (dict): dictionary containing model data
        (made by AC_tools' get_shared_data_as_dict function )
    rm_strat (bool): (fractionally) replace values in statosphere with zeros
    weight_by_molecs (bool): weight grid boxes by number of molecules

    Returns
    -------
    (None)
    """
    # --- Get tags
    # Get dictionaries (reaction, tags, stoichiometry ) if not provided
    if isinstance(RR_dict, type(None)):
        RR_dict = get_dict_of_KPP_mech(wd=CODE_wd,
                                       GC_version=Data_rc['GC_version'],
                                       Mechanism=Mechanism)
        if debug:
            print(len(RR_dict), [RR_dict[i] for i in RR_dict.keys()[:5]])
    if isinstance(tags2_rxn_num, type(None)):
        tags2_rxn_num = {v: k for k, v in tags_dict.items()}
    if isinstance(RR_dict_fam_stioch, type(None)):
        RR_dict_fam_stioch = get_stioch4family_rxns(
            fam=fam, RR_dict=RR_dict, Mechanism=Mechanism)
    # --- Get data
    # Get prod/loss arrays
    ars = get_GC_output(wd=Var_rc['wd'], r_list=True,
                        vars=['PORL_L_S__'+i for i in tags],
                        trop_limit=Var_rc['trop_limit'])
    # Limit prod/loss vertical dimension?
    limit_Prod_loss_dim_to = Var_rc['limit_Prod_loss_dim_to']
    # Covert units based on whether model output is monthly
    if Data_rc['output_freq'] == 'Monthly':
        month_eq = True  # use conversion in convert_molec_cm3_s_2_g_X_s
    else:
        month_eq = False
    # Now convert the units (to G/s)
    ars = convert_molec_cm3_s_2_g_X_s(ars=ars, ref_spec=ref_spec,
                                      # Shared settings...
                                      months=Data_rc['months'],
                                      years=Data_rc['years'],
                                      vol=Data_rc['vol'], t_ps=Data_rc['t_ps'],
                                      trop_limit=Var_rc['trop_limit'],
                                      rm_strat=Var_rc['rm_strat'],
                                      # There are 59 levels of computation for P/l in
                                      # v11-1+ (so limit to 59)
                                      limit_Prod_loss_dim_to=limit_Prod_loss_dim_to,
                                      # ... and function specific settings...
                                      month_eq=month_eq,
                                      conbine_ars=False)
    # Add stoichiometric scaling (# of Ox losses per tagged rxn. )
    ars = [i*RR_dict_fam_stioch[tags2_rxn_num[tags[n]]]
           for n, i in enumerate(ars)]
    # Scale to annual
    if Data_rc['output_freq'] == 'Monthly':
        # Should this be summated then divided adjusted to time points.
        # Sum over time
        ars = [i.sum(axis=-1) for i in ars]
        # Adjust to equivalent months.
        ars = [i/len(Data_rc['months'])*12 for i in ars]
    else:
        # Average over time?
        ars = [i.mean(axis=-1) for i in ars]
        # Scale to annual
        ars = [i*60.*60.*24.*365. for i in ars]
    # Get time in troposphere diagnostic
    print(Data_rc['t_ps'].shape)
    print(Data_rc['t_ps'].mean(axis=-1).shape)
    t_ps = Data_rc['t_ps'].mean(axis=-1)[..., :limit_Prod_loss_dim_to]
    # Check tropospheric LOx total
    arr = np.ma.array(ars)
    print(arr.shape, t_ps.shape, limit_Prod_loss_dim_to)
    LOx_trop = (arr * t_ps[None, ...]).sum() / 1E12
    if verbose:
        print('Annual tropospheric Ox loss (Tg O3): ', LOx_trop)
    # Remove the stratosphere by multiplication through by "time in troposphere"
    if rm_strat:
        ars = [i*t_ps for i in ars]
    # Select data by location or average globally?
    if not isinstance(region, type(None)):
        # Also allow for applying masks here...
        print('NOT SETUP!!!')
        sys.exit()
    else:
        if weight_by_molecs:
            ars = [molec_weighted_avg(i, weight_lon=True, res=Data_rc['res'],
                                      weight_lat=True, wd=Var_rc['wd'],
                                      trop_limit=Var_rc['trop_limit'],
                                      rm_strat=Var_rc['rm_strat'], \
                                      # provide shared data arrays averaged over time...
                                      molecs=Data_rc['molecs'].mean(axis=-1),
                                      t_p=t_ps) \
                   for i in ars]
    if debug:
        print([i.shape for i in ars])
    return ars


# -------------- Functions to produce KPP ***input*** file(s)


def print_out_dfs2KPP_eqn_file(species_df=None, headers=None,
                               extr_str='OUTPUT', rxn_dicts=None):
    """
    Print out DataFrames to *.eqn file
    """
    # ---- Create KPP Mechanism file to save lines to
    a = open('{}.eqn'.format(extr_str), 'w')
    # ---  Print headers to file (list of strings)
    for line in headers:
        print(line.replace('\n', ''), file=a)
    # ----  Print species to file
    print('', file=a)
    print('#include atoms', file=a)
    print('', file=a)
    print('#DEFVAR', file=a)
    print('', file=a)
    # Now print active species
    ptr_str = '{:<11}= IGNORE; '
    for spec in species_df[species_df['inactive'] == False].index:
        df_tmp = species_df[species_df.index == spec]
#        print( spec)
        dsp_str = df_tmp['Description'].values[0]
        if ('{' not in dsp_str):
            dsp_str = '{'+dsp_str+'}'
#        print( ptr_str.format(spec)+dsp_str)
        print(ptr_str.format(spec)+dsp_str, file=a)
#        print( df_tmp['Description'].values[0] )
    # Now print inactive species
    print('', file=a)
    print('#DEFFIX', file=a)
    print('', file=a)
    for spec in species_df[species_df['inactive'] == True].index:
        df_tmp = species_df[species_df.index == spec]
        dsp_str = df_tmp['Description'].values[0]
        if ('{' not in dsp_str):
            dsp_str = '{'+dsp_str+'}'
        print(ptr_str.format(spec)+dsp_str, file=a)
    # ----  Print equations to file
    print('', file=a)
    print('#EQUATIONS', file=a)
#    rxn_types = sorted(rxn_dicts.keys())
    rxn_types = ['Gas-phase', 'Heterogeneous', 'Photolysis']
    for rxn_type in rxn_types:
        # Setup print strings
        KPP_line_max = 43
        if (rxn_type == 'Photolysis'):
            func_max_len = 15
            ptr_str = '{:<44} {:<15} {}'
        else:  # Heterogeneous and Gas-phase
            func_max_len = 36.
            ptr_str = '{:<44} {:<36} {}'
        # print out headers for mech
        print('//', file=a)
        print('// {}'.format(rxn_type) + ' reactions', file=a)
        print('//', file=a)
        # Now loop by reaction
        rxn_dict = rxn_dicts[rxn_type]
        for ind in rxn_dict.index:
            print
            df_tmp = rxn_dict[rxn_dict.index == ind]
            rxn = df_tmp['rxn_str'].values[0]
            eqn = df_tmp['eqn'].values[0]
            metadata = df_tmp['metadata'].values[0].strip()
            # Check if rxn is longer than KPP rxn strin max...
            if (len(rxn) < KPP_line_max):
                ptr_vars = rxn+' :', eqn+';', metadata
                print(ptr_str.format(*ptr_vars), file=a)
            else:
                print(rxn, KPP_line_max)
                # Just split substrings ? (NO, must split by reactant)
#                 frac_of_max = len( rxn ) / KPP_line_max
#                 nchunks = myround( frac_of_max, base=1, round_up=True )
#                 sub_strs = chunks(l=rxn, n=int(KPP_line_max)-2 )
                sub_strs = split_KPP_rxn_str_into_chunks(rxn)
                # Loop sub strings
                for n, sub_str in enumerate(sub_strs):
                    # If not the final sub string, then just print sub string
                    if (n+1 != len(sub_strs)):
                        print(sub_str, file=a)
                    # Now print substring
                    else:
                        ptr_vars = sub_str+' :', eqn+';', metadata
                        print(ptr_str.format(*ptr_vars), file=a)
    # ---  Now close the file...
    a.close()


# -------------- Extract rxn data from KPP ***output*** files


def KPP_eqn_file_headers(folder=None, filename=None):
    """
    Get headers from KPP *.eqn file
    """
    # Open file and loop by line
    header_lines = []
    with open(folder+filename, 'r') as file_:
        for line in file_:
            header_lines += [line]
            if (line.strip() == '}'):
                break
    return header_lines


def get_reactants_and_products4tagged_fam(fam='LOx', KPP_output_mech=None,
                                          folder=None):
    """ Return a list of reactants and products of reactions making fam tag """
    import pandas as pd
    # Get the outputted KPP mechanism if not provided.
    if isinstance(KPP_output_mech, type(None)):
        KPP_output_mech = get_dict_of_KPP_mech(wd=folder)
    # Make a DataFrame from the dictionary
    s = pd.Series(KPP_output_mech)
    df = pd.DataFrame()
    df['rxn str'] = s
#    df.index = df.index - 1 # Adjust Fortran numbering in output to Pythonic #
    # Split apart reaction str to give products and reactants
    def extract_products(input):
        return str(input).split(' --> ')[-1].strip()

    def extract_reactants(input):
        return str(input).split(' --> ')[0].strip()
    df['prod'] = df['rxn str'].map(extract_products)
    df['react'] = df['rxn str'].map(extract_reactants)
    # Make the formating the same as in input files
    # - Check if the file mapping can work between KPP input & output
    df['KPP input react'] = df['react'].map(update_KPP_rxn_str)
#    KPP_Input_df['KPP input react'] = KPP_Input_df['react'].map(update_rxn_str)
    find_breakpoint = False
    if find_breakpoint:  # Seems to be C3H8 + OH = A3O2 :
        #        for s,e in [ [0, 50], [] ]:
        for s, e in zip(range(0, 700, 50), range(50, 750, 50)):
            bool = df['react'][s:e] == KPP_Input_df['KPP input react'][s:e]
            a_ = df['react'][s:e][~bool]
            b_ = KPP_Input_df['KPP input react'][s:e][~bool]
#            print dict( zip( a_, b_ )  )
            # Turn into sorted lists
            a_ = [list(sorted([ii.strip()for ii in i.split('+')])) for i in a_]
            b_ = [list(sorted([ii.strip()for ii in i.split('+')])) for i in b_]
            # Compare these and only print if they differ
            for n, a__ in enumerate(a_):
                if a__ != b_[n]:
                    print(a__, b_[n])
    # Only consider reaction that include family in the products
    def fam_in_rxn(input):
        return (fam in input)
    rtn_vars = ['react', 'prod', 'KPP input react']
    df = df[df['prod'].map(fam_in_rxn)][rtn_vars]
    return df


def get_KPP_tagged_rxns(fam='LOx', filename='gckpp_Monitor.F90',
                        Mechanism='Halogens', RR_dict=None, wd=None):
    """
    Search compiled KPP mechanism for reactions with a certain tag ("fam") in
    their mechanism (e.g. LOx for loss of Ox)

    Parameters
    -------
    fam (str): prod/loss family in KPP file (e.g. Ox loss or "LOx" )
    wd (str): the working (code) directory to search for files in
    Mechanism (str): name of mechanism (e.g. dir in KPP folder)
    filename (str): name of KPP monitor file

    Returns
    -------
    (list)
    """
    # Get dictionary of reactions in Mechanism
    if isinstance(RR_dict, type(None)):
        RR_dict = get_dict_of_KPP_mech(Mechanism=Mechanism, filename=filename,
                                       wd=wd)
    # Loop dictionary of reactions and save those that contain tag
    tagged_RR_dummies = []
    for key_ in list(RR_dict.keys()):
        # Get reaction string
        rxn_str = RR_dict[key_]
        # Collection reactions with tag
        if fam in rxn_str:
            tagged_RR_dummies += [key_]
    return tagged_RR_dummies


# -------------- Extract rxn data from KPP ***input*** file(s)


def KPP_eqn_file_species(folder=None, filename=None, debug=False):
    """
    Get species from *.eqn file
    """
    import pandas as pd
    import numpy as np
    # ----open file and loop by line
    with open(folder+filename, 'r') as file_:
        specs = []
        num2read_line_from = 999999
        species_is_not_active = False
        inactive = []
        for n_line, line_ in enumerate(file_):
            # Are we at the species section?
            if ('#DEFVAR' in line_):
                num2read_line_from = n_line+2
            # Now save the species detail
            if (n_line >= num2read_line_from) and (len(line_) > 1):
                specs += [line_]
            # Also include inactive species? - Yes
            # (but track that these are inactive)
            if ("#DEFFIX" in line_):
                species_is_not_active = True
            if ('IGNORE' in line_):
                inactive += [species_is_not_active]
            # Stop when got to lines of equations
            if ("#EQUATIONS" in line_):
                break
            if debug:
                print(num2read_line_from)
    # Process extracted lines to dict of names and descriptions...
    def strip_line(line_):
        try:
            name, descrip = line_.strip().split('= IGNORE;')
            line_ = [name.strip(), descrip.strip()]
        except:
            print('ERROR for: {}'.format(line_))
        return line_
    specs = [strip_line(i) for i in specs if (len(i) > 12)]
    # Return as a DataFrame
    maxlen = 400
    df = pd.DataFrame([i[1] for i in specs], dtype=(np.str, maxlen))
    df.index = [i[0] for i in specs]
    df.columns = ['Description']
    # Add activity to DataFrame
    df['inactive'] = inactive
    return df


def get_dicts_of_KPP_eqn_file_reactions(folder=None, filename=None,
                                        debug=False):
    """
    Get reactions (Heterogeneous, Photolysis, Gas-phase) from *.eqn file
    """
    import string
    # Local vars
    KPP_rxn_funcs = (
    # Main KPP functions
    'HET', 'PHOTOL', 'GCARR','GCJPLPR',
#    'GC_',
    # Specialist KPP functions for mechanism
    'GC_HO2HO2', 'GC_OHCO',
    'GC_RO2NO', 'GC_OHHNO3', 'GC_RO2HO2', 'GC_HACOHA',
    'GC_RO2HO2', 'GC_HACOHB', 'GC_TBRANCH', 'GCJPLEQ',
    'GC_GLYCOHA', 'GC_GLYCOHB', 'GC_DMSOH', 'GC_GLYXNO3',
    # Include ISOP reaction functions (inc. GC)
    'GC_ALK', 'GC_NIT', 'GC_PAN', 'GC_EPO', 'GC_ISO1', 'GC_ISO2',
    # KPP function without GC prefix
#   'NIT', 'PAN', 'ALK', 'EPO',
    'ARRPLUS', 'TUNPLUS',
    '1.33E-13+3.82E-11*exp', # Why is this function not in gckpp.kpp?
    )
    # Loop lines in file
    with open(folder+filename, 'r') as file_:
        rxns_dict = {}
        for rxns in ('Gas-phase', 'Heterogeneous', 'Photolysis',):
            #        for rxns in ('Gas-phase',):
            num2read_line_from = 999999
            eqns = []
            # Tmp variables to catch KPP reactions over more than one line
            tmp_line_ = ''
            partial_rxn_str = False
            # Now loop
            for n_line, line_ in enumerate(file_):
                # Are we at the species section?
                if ('//' in line_) and (rxns in line_):
                    num2read_line_from = n_line+2
                # Now save the species detail
                if (n_line >= num2read_line_from) and (len(line_) > 1):
                    # Is there a rxn function in  the line?
                    if not any([(i in line_) for i in KPP_rxn_funcs]):
                        # If not then str/line only contains a partial strin
                        partial_rxn_str = True
                        rxn_func_in_str = False
                    else:
                        rxn_func_in_str = True
                    # If only part of rxn str, then save and add to next line
                    if partial_rxn_str:
                        tmp_line_ += line_.strip()
                        if debug:
                            print('added to tmp_line_:', tmp_line_, line_)
                    # If the tmp str is empty, just fill it with the line_
                    if (tmp_line_ == ''):
                        tmp_line_ = line_.strip()
                    # (extra) check: is there is rxn func in the tmp_line_ str?
                    if any([(i in tmp_line_) for i in KPP_rxn_funcs]):
                        rxn_func_in_str = True
                    if rxn_func_in_str:
                        eqns += [tmp_line_]
                        # Reset the tmp str
                        tmp_line_ = ''
                        partial_rxn_str = False
                # Stop reading lines at end of section
                if ("//" in line_) and len(eqns) > 3:
                    break
                print(n_line, line_)
            # Remove spacing errors in KPP eqn entries

            def remove_KPP_spacing_errors(input):
                """ Remove differences in spacing in KPP """
                # Numbers
                for num in [str(i) for i in range(0, 10)]:
                    input = input.replace(' +'+num, ' + '+num)
                # Lower case letters
                for let in string.ascii_lowercase:
                    input = input.replace(' +'+let, ' + '+let)
                # Upper case letters
                for let in string.ascii_uppercase:
                    input = input.replace(' +'+let, ' + '+let)
                return input
            # Remove KPP spacing errors
            eqns = [remove_KPP_spacing_errors(i) for i in eqns]
            # Check for "eqns" with more than one eqn in
            eqns = split_combined_KPP_eqns(eqns)
            # Save out
            rxns_dict[rxns] = eqns
    return rxns_dict


# -------------- Misc. helper functions for KPP file processing/analysis


def split_KPP_rxn_str_into_chunks(rxn, KPP_line_max=43, debug=False):
    """
    Split a rxn str so that it is not longer than KPP max line length
    """
    print(rxn)
    # Sub-function to cut strings to last " +"

    def return_string_in_parts_ending_with_plus(input):
        """ return string upto last ' +'  in string """
        # Find the last ' + ' in string
        ind = input[::-1].find(' +')
        return input[:-(ind)-1]
    # Now loop the sub_Strs
    rxn_len = len(rxn)
    remainder = rxn_len
    sub_str_lens = [0]
    sub_strs = []
    # Loops to do ()
    times2loop = rxn_len / float(KPP_line_max)
    times2loop = myround(times2loop, round_up=True, base=1) + 1
    for n in range(times2loop):
        if debug:
            print(n, remainder, sub_str_lens)
        if (remainder > 0):  # And len(sub_str):
            #        while (remainder > 0):
            # Cut reaction string from length taken
            new_str = rxn[sum(sub_str_lens):][:int(KPP_line_max)]
            if debug:
                print(new_str)
            # If the new string is the end of the rxn string
            if len(new_str) < int(KPP_line_max):
                sub_strs += [new_str]
                sub_str_lens += [len(new_str)]
                remainder = 0
                break
            # If not chop str to las "+"
            else:
                new_str = return_string_in_parts_ending_with_plus(new_str)
                # Now save to sub str list
                sub_strs += [new_str]
                sub_str_lens += [len(new_str)]
                remainder = rxn_len - sum(sub_str_lens)
            if debug:
                print(new_str)
        else:
            break
    print(sub_strs)
    return sub_strs


def get_KPP_PL_tag(last_tag, tag_prefix='T'):
    """
    Get the next P/L tag in a format T???
    """
    assert (len(last_tag) == 4), "Tag must be 4 characers long! (e.g. T???)"
    last_tag_num = int(last_tag[1:])
    return '{}{:0>3}'.format(tag_prefix, last_tag_num+1)


def print_out_lines_for_gckpp_file(tags=None, extr_str=''):
    """
    Print lines to gckpp.kpp for added for tags in now in the .eqn file
    """
    # Create a *.kpp file for lines to be saved to
    a = open('gckpp.kpp_extra_lines_for_tagged_mech_{}'.format(extr_str), 'w')
    print('#FAMILIES', file=a)
    for tag in tags:
        print('{} : {};'.format('P'+tag, tag), file=a)
    # Now close the file...
    a.close()


def update_KPP_rxn_str(input, rtn_as_list=False):
    """
    Update reaction string in KPP reaction string
    """
    # Remove +hv from phtolysis reactions
    input = input.replace('+hv', '')
    input = input.replace('+ hv', '')
    # Remove {M} from 3 body reactions
    input = input.replace('{+M}', '')
    # Remove any trailing/leading  white space
    input = input.strip()
    # If a self reaction, then  write  Y
    if '+' in input:
        tmp_ = input.split('+')
        if len(tmp_) == 2:
            if tmp_[0].strip() == tmp_[1].strip():
                input = '{} {}'.format(len(tmp_), tmp_[0])
            elif rtn_as_list:
                input = sorted(tmp_)
        else:
            print('ERROR!', tmp_, input)
    return input.strip()


def process_KPP_rxn_dicts2dfs(rxn_dicts=None, debug=False,
                                     Use_FORTRAN_KPP_numbering=True):
    """
    Process lists of strings of KPP mechanism into organised DataFrame
    """
    import pandas as pd
    # Loop dictionary of reaction mechanisms
    for key_ in rxn_dicts.keys():
        rxn_list = rxn_dicts[key_]
        rxn_list_processed = []
        # Loop reactions and extract and clean meta data
        for rxn in rxn_list:
            rxn_str = rxn.split(':')[0].strip()
            eqn = rxn.split(':')[1].split(';')[0].strip()
            metadata = rxn[rxn.find(';')+1:]
#            if len(metadata) > 0:
#                print( metadata )
#                metadata_start = metadata.find('{')
#                metadata_end = metadata[::-1.find('}')
#                metadata = metadata[metadata_start:metadata_end+1]
            # Split off reactants and products...
            react = rxn_str.split('=')[0].strip()
            prod = rxn_str.split('=')[-1].strip()
            rxn_list_processed += [[rxn_str, react, prod, eqn, metadata]]
        # Create a dataframe
        df = pd.DataFrame(rxn_list_processed)
        df.columns = ['rxn_str', 'react', 'prod', 'eqn', 'metadata']
        # Force use of index numbers from FORTRAN / KPP
        if Use_FORTRAN_KPP_numbering:
            df.index = df.index + 1
            print('WARNING - updated DataFrame dictionary to start @1 (not 0)')
        # Overwrite dictionary with dataframe
        rxn_dicts[key_] = df
    return rxn_dicts


def split_combined_KPP_eqns(list_in):
    """ Split combined KPP eqn strings  """
    # Indices of KPP eqn strings with more than one "="
    inds = []
    for n, str in enumerate(list_in):
        if str.count('=') > 1:
            inds += [n]
    strs = [list_in[i] for i in inds]
    print('WARNING: found {} combined strings:'.format(len(inds)), strs)
    # If there are any combined eqn strings, then split them
    if len(inds) > 0:
        new_eqns = []
        for ind in inds:
            tmp_str = list_in[ind]
            print(tmp_str, inds)
            ind_1st_equals = tmp_str.find('=')
            ind_2nd_equals = tmp_str[::-1].find('=')
            tmp_str_short = tmp_str[ind_1st_equals:len(tmp_str)-ind_2nd_equals]
            try:
                last_bracket_of_1st_eqn = tmp_str_short[::-1].find('}')
                # Calculate the split point for the reactions
                ind_2nd_eqn = ind_1st_equals+len(tmp_str_short)
                ind_2nd_eqn = ind_2nd_eqn - last_bracket_of_1st_eqn
                #
                new_eqns += [[tmp_str[:ind_2nd_eqn], tmp_str[ind_2nd_eqn:]]]
            except:
                print('Issue with inputted file (meta data... )')
                sys.ext()
        # Drop old ones
        [list_in.pop(i) for i in sorted(inds)[::-1]]
        # Add news ones and return
#        if len(inds) == 1:
        for n, ind in enumerate(inds):
            list_in = list_in[:ind] + new_eqns[n] + list_in[ind:]
            print(list_in[:ind][-5:])
            print(list_in[ind+1:][:5])
    return list_in


# -------------- Funcs specifically for v11-1g and later ends here
# (funcs still do use funcs from elsewhere)


def get_stioch4family_rxns(fam='LOx', filename='gckpp_Monitor.F90',
                           Mechanism='Halogens', RR_dict=None, wd=None,
                           debug=False):
    """
    Get the stiochmetery for each reaction in family from compiled KPP file

    Parameters
    -------
    fam (str): prod/loss family in KPP file (e.g. Ox loss or "LOx" )
    wd (str): the working (code) directory to search for files in
    Mechanism (str): name of mechanism (e.g. dir in KPP folder)
    filename (str): name of KPP monitor file

    Returns
    -------
    (dict)
    """
    # Get dictionary of reactions in Mechanism
    if isinstance(RR_dict, type(None)):
        RR_dict = get_dict_of_KPP_mech(Mechanism=Mechanism, filename=filename,
                                       wd=wd)
    # Assume unity if no coeffeicent
    unity = 1.0
    # Loop dictionary of reactions and save those that contain tag
    tagged_rxns = []
    tagged_rxn_stioch = []
    for key_ in list(RR_dict.keys()):
        # Get reaction string
        rxn_str = RR_dict[key_]
        # Collection reactions with tag
        if fam in rxn_str:
            tagged_rxns += [key_]
            # Split reaction str by '+' (excluding reactants part)
            rxn_str = rxn_str[18:].split('+')
            if debug:
                print(rxn_str)
            # Get product
            product_str = [i for i in rxn_str if (fam in i)]
            if debug:
                print(product_str)
            # Find index of reaction
#            ind_of_rxn = [ n for n,i in product_str if (fam in i) ]
            # Split coefficient ("Coe") from reaction
            product_str = product_str[0].strip().split()
            if len(product_str) > 1:
                tagged_rxn_stioch += [float(product_str[0])]
            else:
                tagged_rxn_stioch += [unity]
    return dict(list(zip(tagged_rxns, tagged_rxn_stioch)))


def get_tags4family(fam='LOx', filename='gckpp_Monitor.F90',
                    Mechanism='Halogens', tag_prefix='PT', RR_dict=None,
                    wd=None,
                    get_one_tag_per_fam=True, debug=False):
    """
    For a P/L family tag (e.g. LOx), if there are individual tags then these
    are extracted

    Parameters
    -------
    fam (str): prod/loss family in KPP file (e.g. Ox loss or "LOx" )
    wd (str): the working (code) directory to search for files in
    Mechanism (str): name of mechanism (e.g. dir in KPP folder)
    filename (str): name of KPP monitor file

    Returns
    -------
    (dict)
    """
    # Get dictionary of reactions in Mechanism
    if isinstance(RR_dict, type(None)):
        RR_dict = get_dict_of_KPP_mech(Mechanism=Mechanism, filename=filename,
                                       wd=wd)
    # Loop dictionary of reactions and save those that contain tag
    tagged_rxns = []
    tagged_rxn_tags = []
    for key_ in list(RR_dict.keys()):
        # Get reaction string
        rxn_str = RR_dict[key_]
        # Collection reactions with tag
        if fam in rxn_str:
            #            if debug: print( key_, fam, rxn_str )
            tagged_rxns += [key_]
            # Split reaction str by '+' (excluding reactants part)
            rxn_str = rxn_str[18:].split('+')
            # Look for tag(s)... - should only be one per reaction!
            tags = [i.strip() for i in rxn_str if (tag_prefix in i)]
            if debug:
                print(rxn_str, tags)
            if len(tags) >= 1:
                if len(tags) == 1:
                    tagged_rxn_tags += tags
                else:
                    prt_str = 'WARNING: {} tags for rxn! - {} - {}'
                    if debug:
                        print(prt_str.format(len(tags), tags, rxn_str))
                    if get_one_tag_per_fam:
                        tagged_rxn_tags += [tags[0]]
                    else:
                        tagged_rxn_tags += tags
            else:
                tagged_rxn_tags += [
                    'WARNING: RXN. NOT TAGGED! ({})'.format(key_)]
                if debug:
                    print((key_, fam, 'ERROR!', tags, rxn_str))
    assert len(tagged_rxns) == len(tagged_rxn_tags), "# tags doesn't = # rxns!"
    return dict(list(zip(tagged_rxns, tagged_rxn_tags)))


def get_dict_of_KPP_mech(filename='gckpp_Monitor.F90',
                         Mechanism='Halogens', GC_version='v11-01', wd=None):
    """
    Get a dictionary of KPP mechansim from compile mechanism

    Parameters
    -------
    wd (str): the working (code) directory to search for files in
    Mechanism (str): name of mechanism (e.g. dir in KPP folder)
    filename (str): name of KPP monitor file
    GC_version (str): name of GEOS-Chem version

    Returns
    -------
    (dict)
    """
    log_str = 'get_dict_of_KPP_mech called for Mech:{} (ver={}, wd={})'
    logging.info(log_str.format(Mechanism, GC_version, wd))
    # Get base working code directory
    if isinstance(wd, type(None)):
        #        wd = sys.agrv[1]
        # Hardwire for testing
        #         wd = '/work/home/ts551/data/all_model_simulations/iodine_runs/'
        #         wd += 'iGEOSChem_5.0/code_TMS_new/'
        #         wd += '/KPP/{}/'.format(Mechanism)
        print('wd with code must be provided')
    MECH_start_str = 'INTEGER, DIMENSION(1) :: MONITOR'
#    rxn_line_ind = '! index'
    strs_in_non_rxn_lines = 'N_NAMES_', 'CHARACTER', 'DIMENSION', 'PARAMETER'
    RR_dict = {}
    with open(wd+filename, 'r') as file:
        # Loop by line in file
        start_extracting_mech_line = 1E99
        for n, line in enumerate(file):
            # Check for mechanism identifier
            if (MECH_start_str in line):
                start_extracting_mech_line = n+3
            # Extract reaction mechanism info
            if n >= start_extracting_mech_line:
                # break out of loop after empyty line (len<3)
                if len(line) < 3:
                    break
                # Check if the line contains a reaction str
#                if (rxn_line_ind in line) or start_extracting_mech_line:
                # Check if line contrains strings not in reaction lines
                if not any([(i in line) for i in strs_in_non_rxn_lines]):
                    rxn_str = line[6:106]
                    # RR??? dummy tags only available on v11-01 patches!
                    if GC_version == 'v11-01':
                        RR_dummy = 'RR{}'.format(line[118:].strip())
                        # Add to dictionary
                        RR_dict[RR_dummy] = rxn_str
                    # Use just reaction number for other versions...
                    else:
                        try:
                            RR_num = int(line[118:].strip())
                            # Add to dictionary
                            RR_dict[RR_num] = rxn_str
                        except ValueError:
                            # Add gotcha for lines printed without an index in
                            # KPP
                            # e.g. final line is missing index in KPP output.
                            # CHARACTER(LEN=???), PARAMETER, DIMENSION(??) ::
                            #  EQN_NAMES_?? = (/ &
                            #  ' ... ', & ! index ???
                            #  ' ... ' /)
                            # (no index given for the final output. so assume
                            # N+1)
                            logging.debug('rxn # assumed as {} for {}'.format(
                                RR_num+1, rxn_str))
                            # use the value of the line previous and add to
                            # dictionary
                            RR_dict[RR_num+1] = rxn_str
    return RR_dict


def prt_families4rxns_to_input_to_PROD_LOSS(fam='LOx', rxns=None, wd=None,
                                            filename='gckpp_Monitor.F90',
                                            Mechanism='Halogens'):
    """
    Takes a list of reaction numbers or (fam) and prints out prod/loss in
    a form that can be copy/pasted into KPP (gckpp.kpp)

    Parameters
    -------
    rnxs (list): lisst of reactions to use for
    fam (str): prod/loss family in KPP file (e.g. Ox loss or "LOx" )
    wd (str): the working (code) directory to search for files in
    Mechanism (str): name of mechanism (e.g. dir in KPP folder)
    filename (str): name of KPP monitor file

    Returns
    -------
    (None)
    """
    # Get list of tagged reactions for family
    if isinstance(rxns, type(None)):
        rxns = get_KPP_tagged_rxns(fam=fam, filename=filename,
                                   Mechanism=Mechanism, wd=wd)
    # --- Print to screen
    header_Str = '#FAMILIES'
    pstr = 'P{} : {};'
    print('>>>> Copy and paste the below text into gckpp.kpp <<<')
    print(header_Str)
    for rxn in rxns:
        print((pstr.format(rxn, rxn)))


def prt_families4rxns_to_input_to_PROD_LOSS_globchem_spec(fam='LOx',
                                                          rxns=None, wd=None,
                                                          filename='gckpp_Monitor.F90',
                                                          Mechanism='Halogens'):
    """
    Takes a list of reaction numbers or (fam) and prints out prod/loss in
    a form that can be copy/pasted into KPP (globchem.spc)

    Parameters
    -------
    rnxs (list): lisst of reactions to use for
    fam (str): prod/loss family in KPP file (e.g. Ox loss or "LOx" )
    wd (str): the working (code) directory to search for files in
    Mechanism (str): name of mechanism (e.g. dir in KPP folder)
    filename (str): name of KPP monitor file

    Returns
    -------
    (None)

    Notes
    -------
     - From v11-2d the KPP mechanism is in a single *.eqn file (globchem.spc was
        retired).
    """
    # Get list of tagged reactions for family
    if isinstance(rxns, type(None)):
        rxns = get_KPP_tagged_rxns(fam=fam, filename=filename,
                                   Mechanism=Mechanism, wd=wd)
    # --- Print to screen
    header_Str = '#DEFVAR'
    pstr = '     {:<10}      =          IGNORE;'
    print('>>>> Copy and paste the below text into globchem.spc <<<')
    print(header_Str)
    for rxn in rxns:
        print((pstr.format(rxn, rxn)))


def get_KKP_mech_from_eqn_file_as_dicts(folder=None, Mechanism='Tropchem',
                                        filename=None):
    """
    Get KPP mechanism as a dictionary of dictionaries
    """
    # Set the filename if not provided
    if isinstance(filename, type(None)):
        filename = '{}.eqn'.format(Mechanism)
    # Get dictionaries of all reactions
    print('get_KKP_mech_from_eqn_file_as_dicts', folder, filename)
    rxn_dicts = get_dicts_of_KPP_eqn_file_reactions(
        folder=folder, filename=filename)
    # Process rxns to be in dictionaries of DataFrames
    # (with extra diagnostic columns, inc. reactants, products, metadata,...)
    rxn_dicts = process_KPP_rxn_dicts2dfs(rxn_dicts=rxn_dicts)
    # Update the numbering of DataFrame indexes...
    Gas_dict = rxn_dicts['Gas-phase']
    Het_dict = rxn_dicts['Heterogeneous']
    Hv_dict = rxn_dicts['Photolysis']
    # Update index start points
    Het_dict.index = Het_dict.index + Gas_dict.shape[0]
    Hv_dict.index = Hv_dict.index + Gas_dict.shape[0] + Het_dict.shape[0]
    rxn_dicts['Heterogeneous'] = Het_dict
    rxn_dicts['Photolysis'] = Hv_dict
    # Return  dictionary of dictionaries
    return rxn_dicts


def get_KKP_mech_from_eqn_file_as_df(folder=None, Mechanism='Tropchem',
                                     filename=None):
    """
    Get KPP mechanism as a pandas DataFrame
    """
    # Set the filename if not provided
    if isinstance(filename, type(None)):
        filename = '{}.eqn'.format(Mechanism)
    print('get_KKP_mech_from_eqn_file_as_df', folder, filename)
    # Get dictionary of dictionaries
    dicts = get_KKP_mech_from_eqn_file_as_dicts(folder=folder,
                                                Mechanism=Mechanism,
                                                filename=filename)
    # Add a flag for reaction type
    for key in dicts.keys():
        df = dicts[key]
        df['Type'] = key
        dicts[key] = df
    # Combine into a single data frame
    df = pd.concat([dicts[i] for i in dicts.keys()], join='outer')
    return df


def get_dictionary_of_tagged_reactions(filename='globchem.eqn',
                                       Mechanism='Halogens', wd=None):
    """
    Construct a dictionary of reaction strings for tagged reactions in
    globchem.eqn

    Notes
    -------
     - From v11-2d the KPP mechanism is in a single *.eqn file (globchem.eqn was
        retired and replaced with <mechanism name>.eqn).
    """
    # Get base working code directory
    if isinstance(wd, type(None)):
        #        wd = sys.agrv[1]
        # Hardwire for testing
        wd = '/work/home/ts551/data/all_model_simulations/iodine_runs/iGEOSChem_5.0/code_TMS_new/'
        wd += '/KPP/{}/'.format(Mechanism)
    # Local variable
    tracer_names_like_tag = ['TRO2']
    # Initialise lists to store tags
    tags_l = []
    rxn_l = []
    # Now open file and look for tags
    with open(wd+filename) as file_:
        # Loop lines and store reaction strs (1st line) that contain tags
        for line in file_:
            try:
                # What does the reaction string begin?
                start_str_index = line.index('}')+2
                # Select reaction string and just products.
                reaction_str = line[start_str_index:51]
                product_str = reaction_str.split(' =')[-1]
                # Assume Tag is 1st and in form TG??? (where ??? are numerical )
                # update this to use regular expressions?
                part_of_str_that_would_contain_tag = product_str[:6]
                if ' T' in part_of_str_that_would_contain_tag:
                    # Strip tag from string.
                    tag = part_of_str_that_would_contain_tag.strip()
                    # Add if actually a tag.
                    if tag not in tracer_names_like_tag:
                        tags_l += [tag]
                        # Select tag
                        rxn_l += [reaction_str.strip()]
            # If the line doesn't contain a reaction # ("{?}"), then pass
            except ValueError:
                pass
    return dict(list(zip(tags_l, rxn_l)))


def get_Ox_fam_based_on_reactants(filename='gckpp_Monitor.F90',
                                  Mechanism='Halogens', wd=None,  fam='LOx',
                                  tag_prefix='PT',
                                  tags=None, input_KPP_mech=None, RR_dict=None,
                                  GC_version='v11-01',
                                  debug=False):
    """
    Get Ox famiy of reaction tag using KPP reaction string.

    Parameters
    -------
    rxns (list): lisst of reactions to use for
    fam (str): prod/loss family in KPP file (e.g. Ox loss or "LOx" )
    wd (str): the working (code) directory to search for files in
    Mechanism (str): name of mechanism (e.g. dir in KPP folder)
    filename (str): name of KPP monitor file

    Returns
    -------
    (None)
    """
    # --- Local variables
    # Get dictionary of reactions in Mechanism
    if isinstance(RR_dict, type(None)):
        RR_dict = get_dict_of_KPP_mech(Mechanism=Mechanism,
                                       filename=filename, wd=wd)
    if isinstance(input_KPP_mech, type(None)):
        input_KPP_mech = get_KKP_mech_from_eqn_file_as_df(
            Mechanism=Mechanism, folder=wd)
    # Get tags for reaction unless already provided
    if isinstance(tags, type(None)):
        tags = get_tags4family(wd=wd, fam=fam, filename=filename,
                               Mechanism=Mechanism, tag_prefix=tag_prefix,
                               RR_dict=RR_dict)
    # --- Extra local variables
    HOx = ['OH', 'HO2', 'H2O2']
    # Temporarily add to definition of HOx
    HOx += ['O1D', 'O', 'O3', 'NO3']
    ClOx = ['CFC', 'Cl', 'ClO']
    NOx = ['NO', 'NO2', 'NO3', 'N2O5']
    non_I_specs_with_I_char = [
        'INO2', 'ISN1', 'ISNOOA', 'ISNOHOO', 'ISOPNB', 'ISOPND', 'ISOP', 'ISNP',
        # Add species in v11-2d
        'IONITA',
        # Add species in v11-2d (benchmarks)
        'IMAO3', 'IPMN', 'MONITS', 'MONITU', 'ISNP', 'LIMO', 'ISNOOB', 'PIO2',
        'LIMO2'
    ]
    # Get list of tagged hv reactions
#    all_rxns = sorted(RR_dict.keys())
#    all_rxns = sorted(RR_hv_dict.keys()) # This also give the same list as above line
#     hv_rxns_tags = []
#     for rxn_ in all_rxns:
#         rxn_str = input_KPP_mech.loc[input_KPP_mech.index == rxn_,:]['rxn_str'].values[0]
#         if ('hv' in rxn_str):
#             hv_rxns_tags += [rxn_]
#         else:
#             pass
    hv_rxns = input_KPP_mech.loc[input_KPP_mech['Type'] == 'Photolysis', :]
    hv_rxns = hv_rxns.index.astype(list)

    # The below line will not work as "hv" not in rxn strings in gckpp_Monitor.F90
#    hv_rxns_tags = [i for i in all_rxns if ('hv' in RR_dict[i])]
    # Loop tagged reactions and assign family
    tagged_rxns = sorted(tags.keys())
    fam_l = []
    for n_rxn, rxn_ in enumerate(tagged_rxns):
        _debug = False
        # Initialise booleans
        Br_is_reactant = False
        Cl_is_reactant = False
        I_is_reactant = False
        NOx_is_reactant = False
        HOx_is_reactant = False
        hv_is_reactant = False
        IONITA_is_formed = False
        # --- Get reaction string  and check within
        rxn_str = RR_dict[rxn_]
        reactant_str = rxn_str.split('-->')[0]
        product_str = rxn_str.split('-->')[1]
        # - Check if species are in reactions
        # BrOx reaction?
        if 'Br' in reactant_str:
            Br_is_reactant = True
        # IOx reaction?
        if 'I' in reactant_str:
            # Account for non-iodine species with I Character in...
            if not any([(i in reactant_str) for i in non_I_specs_with_I_char]):
                I_is_reactant = True
        # ClOx reaction?
        if any([(i in reactant_str) for i in ClOx]):
            Cl_is_reactant = True
        # Gotcha for heterogenous N2O5 breakdown by Cl-
#        if (tags[rxn_] == 'T172') and (ver='v11-01')::
#        if (tags[rxn_] == 'T164') and (GC_version=='v11-01'):
#            Cl_is_reactant=True
        if ('N2O5' in reactant_str) and ('ClNO2' in reactant_str):
            Cl_is_reactant = True
        # NOx reaction?
        if any([(i in reactant_str) for i in NOx]):
            NOx_is_reactant = True
#        print (tags[rxn_][1:] in hv_rxns_tags), tags[rxn_][1:], hv_rxns_tags
#        print RR_dict[rxn_]
        # hv reaction?
        if (rxn_ in hv_rxns):
            hv_is_reactant = True
        # HOx reaction?
        if any([(i in reactant_str) for i in HOx]):
            HOx_is_reactant = True
        # Gotcha for IONITA in v11-01
#         if (GC_version=='v11-01'):
#             IONITA_rxns = ['T040', 'T039', 'T037']
#             if any( [(tags[rxn_] == i) for i in IONITA_rxns] ):
#                 IONITA_is_formed = True
        # IONITA in  reaction?
        # (IONITA = Aer-phase organic nitrate from isoprene precursors)
        if ('IONITA' in product_str):
            IONITA_is_formed = True
        # --- Assign familes
        # 1st check if halogen crossover reaction...
        if Cl_is_reactant and Br_is_reactant:
            fam_l += ['ClOx-BrOx']
        elif Cl_is_reactant and I_is_reactant:
            fam_l += ['ClOx-IOx']
        elif Br_is_reactant and I_is_reactant:
            fam_l += ['BrOx-IOx']
        # Now just check if single halogen family...
        elif Cl_is_reactant:
            fam_l += ['ClOx']
        elif Br_is_reactant:
            fam_l += ['BrOx']
        elif I_is_reactant:
            fam_l += ['IOx']
        # Now check if photolysis reaction...
        elif hv_is_reactant:
            fam_l += ['hv']
        # HOx?
        elif NOx_is_reactant:
            fam_l += ['NOx']
        # HOx?
        elif HOx_is_reactant:
            fam_l += ['HOx']
        # If IONITA reaction
        elif IONITA_is_formed:
            fam_l += ['NOx']
        # Not assigned?!
        else:
            fam_l += ['NOT ASSIGNED!!!']
            _debug = True
        # Special cases?
        # O1D + H2O --> 2 OH
        if ('O1D' in reactant_str) and ('H2O' in reactant_str):
            fam_l[n_rxn] = 'hv'
        # debug print ?
        if _debug:
            specs = 'Br', 'I', 'Cl', 'NOx', 'hv', 'HOx'
            booleans = Br_is_reactant, I_is_reactant, Cl_is_reactant,  \
                NOx_is_reactant, hv_is_reactant, HOx_is_reactant
            print((list(zip(specs, booleans))))
        if debug:
            print((rxn_str, fam_l[n_rxn]))
    # Nicer family names?
#    if mk_family_names_readable:
    name_dict = {
        'ClOx-BrOx': 'Cl+Br', 'BrOx-IOx': 'Br+I', 'BrOx': 'Bromine',
        'ClOx-IOx': 'Cl+I', 'IOx': 'Iodine', 'ClOx': 'Chlorine', 'HOx': 'HO$_{\\rm x}$',
        'hv': 'Photolysis', 'NOx': 'NO$_{\\rm x}$',
        # Allow unassigned reactions to pass (to be caught later)
        'NOT ASSIGNED!!!': 'NOT ASSIGNED!!!'
    }
    fam_l = [name_dict[i] for i in fam_l]
    tag_l = [tags[i] for i in tagged_rxns]
    return dict(list(zip(tag_l, fam_l)))


def get_tags_in_rxn_numbers(rxn_nums=[], RR_dict=None,
                            filename='gckpp_Monitor.F90', Mechanism='Halogens',
                            wd=None,
                            debug=False):
    """
    Get tags in given list of reactions

    Parameters
    -------
    rxns (list): lisst of reactions to use for
    Mechanism (str): name of mechanism (e.g. dir in KPP folder)
    filename (str): name of KPP monitor file

    Returns
    -------
    (dict)
    """
    # Get dictionary of reactions in Mechanism
    if isinstance(RR_dict, type(None)):
        RR_dict = get_dict_of_KPP_mech(Mechanism=Mechanism, filename=filename,
                                       wd=wd)
    # Loop dictionary of reactions and save those that contain tag
    tagged_rxns = []
    tags_for_rxns = []
    # Create a sub dictionary for reaction numbers provided
    sub_dict = dict([(i, RR_dict[i]) for i in rxn_nums])
    print(sub_dict)
    for key_ in list(sub_dict.keys()):
        # Get reaction string
        rxn_str = RR_dict[key_]
        # Collection reactions with tag
        tag_strs = [' + T', '> T']
        # Replace this with reg ex?
        if debug:
            print((key_, rxn_str, any([i in rxn_str for i in tag_strs])))
        if any([i in rxn_str for i in tag_strs]):
            # Split reaction str by '+' (excluding reactants part)
            product_str_list = rxn_str[18:].split('+')
            # Select tag(s)
            tag_ = [i for i in product_str_list if (' T' in i)]
            if len(tag_) > 1:
                print(('WARNING - more than one tag for reaction? :', tag_))
                sys.exit()
            else:
                tags_for_rxns += [tag_[0].strip()]
                tagged_rxns += [key_]
    return dict(list(zip(tagged_rxns, tags_for_rxns)))


def get_oxidative_release4specs(filename='gckpp_Monitor.F90',
                                Mechanism='Halogens', wd=None):
    """
    Certain species have a fixed concentration in the model, for these speices
    the source value is calculated by the summ of the oxidative release
    ( OH, Br, Cl )

    Parameters
    -------
    wd (str): the working (code) directory to search for files in
    Mechanism (str): name of mechanism (e.g. dir in KPP folder)
    filename (str): name of KPP monitor file

    Returns
    -------
    (list)
    """
    # Species ?
    specs = ['CHBr3', 'CH3Cl', 'CH2Cl2', 'CHCl3']
    # Get dictionary of reactions in Mechanism
    RR_dict = get_dict_of_KPP_mech(Mechanism=Mechanism, filename=filename,
                                   wd=wd)
    # Loop species
    RR_rxn_dummies = []
    for spec in specs:
        # Loop reactions
        for key_ in RR_dict:
            rxn_str = RR_dict[key_]
            # Species in reaction?
            if spec in rxn_str:
                print(rxn_str)
                RR_rxn_dummies += [key_]
    return RR_rxn_dummies


def get_Ox_loss_dicts(fam='LOx', ref_spec='O3',
                      wd=None, Mechanism='Halogens', rm_strat=False,
                      weight_by_molecs=True, CODE_wd=None,
                      dpi=320, suffix='', full_vertical_grid=False,
                      save_plot=True, show_plot=False,
                      verbose=True, debug=False):
    """
    Get Ox loss data and variables as a dictionary (rxn. by rxn.)

    Parameters
    -------
    fam (str): tagged family to track (already compiled in KPP mechanism)
    ref_spec (str): reference species to normalise to
    wd (str): working directory ("wd") of model output
    CODE_wd (str): root of code directory containing the tagged KPP mechanism
    Mechanism (str): name of the KPP mechanism (and folder) of model output
    weight_by_molecs (bool): weight grid boxes by number of molecules
    rm_strat (bool): (fractionally) replace values in statosphere with zeros
    debug, verbose (bool): switches to turn on/set verbosity of output to screen
    full_vertical_grid (bool): use the full vertical grid for analysis

    Returns
    -------
    (None)

    Notes
    -----

    """
    # - Get key model variables, model settings, etc
    # Setup dictionary for shared variables
    Var_rc = get_default_variable_dict(full_vertical_grid=full_vertical_grid)
    # Get locations of model output/core
    assert os.path.exists(wd), 'working directory not found @: {}'.format(wd)
    CODE_wd = '/{}/KPP/{}/'.format(CODE_wd, Mechanism)
    assert os.path.exists(CODE_wd), 'code directory not found @: ' + CODE_wd
    # Use provided working directory
    if not isinstance(wd, type(None)):
        Var_rc['wd'] = wd
    # Setup dictionary for shared data
    var_list = [
        'generic_4x5_wd', 'months', 'years', 'datetimes', 'output_freq',
        'output_vertical_grid', 's_area', 'vol', 't_ps', 'n_air', 'molecs',
        'alt'
    ]
    Data_rc = get_shared_data_as_dict(Var_rc=Var_rc, var_list=var_list)
    # Get reaction dictionary
    RR_dict = get_dict_of_KPP_mech(wd=CODE_wd,
                                   GC_version=Data_rc['GC_version'],
                                   Mechanism=Mechanism)
    if debug:
        print((len(RR_dict), [RR_dict[i] for i in list(RR_dict.keys())[:5]]))
    # Get tags for family
    tags_dict = get_tags4family(fam=fam, wd=CODE_wd, RR_dict=RR_dict)
    tags = list(sorted(tags_dict.values()))
    tags2_rxn_num = {v: k for k, v in list(tags_dict.items())}
    if debug:
        print(tags)
    # Get stiochiometry of reactions for family
    RR_dict_fam_stioch = get_stioch4family_rxns(fam=fam,
                                                RR_dict=RR_dict,
                                                Mechanism=Mechanism)
    # - Extract data for Ox loss for family from model
    ars = get_fam_prod_loss4tagged_mech(RR_dict=RR_dict,
                                        tags=tags, rm_strat=rm_strat,
                                        Var_rc=Var_rc, Data_rc=Data_rc, wd=wd,
                                        fam=fam, ref_spec=ref_spec,
                                        tags2_rxn_num=tags2_rxn_num,
                                        RR_dict_fam_stioch=RR_dict_fam_stioch,
                                        weight_by_molecs=weight_by_molecs)
    # Get dictionary of families that a tag belongs to
    fam_dict = get_Ox_fam_based_on_reactants(fam=fam, tags=tags_dict,
                                             RR_dict=RR_dict, wd=CODE_wd,
                                             Mechanism=Mechanism,
                                             GC_version=Data_rc['GC_version'])
    # Get indices of array for family.
    sorted_fam_names = ['Photolysis', 'HO$_{\\rm x}$', 'NO$_{\\rm x}$']
    halogen_fams = ['Chlorine', 'Cl+Br', 'Bromine', 'Br+I', 'Cl+I', 'Iodine', ]
    sorted_fam_names += halogen_fams
    # - Place all variables/data of share us into a dictionary and return this
    d = {
        'sorted_fam_names': sorted_fam_names,
        'fam_dict': fam_dict,
        'ars': ars,
        'RR_dict_fam_stioch': RR_dict_fam_stioch,
        'RR_dict': RR_dict,
        'tags2_rxn_num': tags2_rxn_num,
        'tags': tags,
        'tags_dict': tags_dict,
        'Data_rc': Data_rc,
        'Var_rc': Var_rc,
        'halogen_fams': halogen_fams,
    }
    return d


def prt_lines4species_database_yml(tags, extr_str=''):
    """
    Print lines for tags to paste into the species_database.yml file
    """
    # Create a *.kpp file for lines to be saved to
    FileStr = 'species_database_extra_lines_for_tagged_mech_{}'
    a = open(FileStr.format(extr_str), 'w')
    pstr1 = '{}:'
    pstr2 = '  Is_Gas: true'
    for tag in tags:
        print(pstr1.format(tag), file=a)
        print(pstr2, file=a)
