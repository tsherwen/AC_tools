#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Functions for SMVGEAR input/output file Processing

Notes
-------
 - These functions are specifically for GEOS-Chem versions prior to v11-01.
 - They are no longer maintained, as KPP is default ODE solve for GEOS-Chem now
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

# AC_tools imports
from .. GEOSChem_bpch import *
from .. GEOSChem_nc import *


# -------------- Smvgear input/output file Processing
# NOTE: this is now redundant as GEOS-Chem uses KPP going forward. procssing/parsing
#       scripts are  also present in AC_tools for KPP.
def prt_Prod_Loss_list4input_geos(spec_list=None, prefix_list=None,
                                  start_num=1):
    """
    Print to screen lines to be inserted into input.geos to track prod/loss of species

    Parameters
    -------
    spec_list (list): list of species to print lines for
    start_num (int): number to start labeling reaction index for input.geos
    prefix_list (list): list of prefixs to dictacte if prod or loss calc ("P" or "L")

    Returns
    -------
    (None)

    Notes
    -------
    """
    # Formatted string to be filled by variables
    line_str = '{:<6}chemical family   : {}: {}'
    # Activity for species?
    if isinstance(prefix_list, type(None)):
        prefix_list = ['P'] * len(spec_list)
    # Loop species (specs) provided
    for n, spec in enumerate(spec_list):
        index = start_num+n
        print((line_str.format(get_suffix(index), prefix_list[n]+spec, spec)))


def prt_Species_List_lines4globchem_dat(spec_list=None, activty_list=None):
    """
    Print to screen lines to be inserted into globchem.dat for given
    list of species

    Parameters
    -------
    spec_list (list): list of species to print lines for
    acvtivty_list (list): list of activitys ("A"=Active, "L"=live, "D"=Dead)

    NOTES
    -------
    - Globchem.dat contains the chemistry for the smvgear solver.
    - These lines must be present at the top to initialise each species
    """
    line1_str = '{} {:<10}         1.00 1.000E-20 1.000E-20 1.000E-20 1.000E-20'
    line2_str = '                     0   0   0   0   0   0   0   0   0     '
    # Activity for species?
    if isinstance(activty_list, type(None)):
        activty_list = ['A'] * len(spec_list)
    # Print lines to insert into globchem.dat
    for n, spec in enumerate(spec_list):
        print((line1_str.format(acvtivty_list[n], spec)))
        print(line2_str)


# -------------- Dynamic processing of p/l
#
# NOTE - below functions for automated processing/parsing of smvgear
# mechanisms for tagging. It is redundant and far more coherent version
# for processing KPP mechanisms and tags is included in funcs4GEOS.py


def rxn_dict_from_smvlog(wd, PHOTOPROCESS=None, ver='1.7',
                         LaTeX=False, debug=False):
    """
    Build a dictionary reaction of reaction details from smv.log.
    This can be used as an external call to analyse other prod/loss
    reactions through smvgear

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    debug (bool): legacy debug option, replaced by python logging
    PHOTOPROCESS (int): smvgear index of 1st photochemical reaction
    LaTeX (bool): convert species to LaTeX style formatting
    ver (str): The GEOS-Chem halogen version that is being used

    Returns
    -------
    (dict) dictionary of active reactions in smvgear output

    Notes
    -----
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """

    if isinstance(PHOTOPROCESS, type(None)):
        PHOTOPROCESS = {
            '1.6': 457, '1.6.2': 452, '1.6.3': 461, '1.7': 467, '2.0': 555,
            '3.0': 547,
            # placeholder for v4.0
            '4.0': 547
        }[ver]

    fn = 'smv2.log'
    if debug:
        print((wd+'/'+fn))
    file_ = open(wd+'/'+fn, 'r')
    readrxn = False
    for row in file_:
        row = row.split()
        if 'NMBR' in row:
            readrxn = True
        if len(row) < 1:
            readrxn = False
        if readrxn:
            try:
                rxns.append(row)
            except:
                rxns = [row]

    # -- remove 'NMBR'
    rxns = [i for i in rxns if ('NMBR' not in i)]
    n = [int(rxn[0]) for rxn in rxns]
    rxns = [rxn[1:] for rxn in rxns]
    rdict = dict(list(zip(n, rxns)))

    # --- Process to Latex
    if LaTeX:
        for rxn in sorted(rdict.keys()):

            # -- Use Latex formatting?
            if rxn > PHOTOPROCESS-1:
                #                xarstr = r' $\xrightarrow{hv}$ '
                xarstr = r' + hv $\rightarrow$ '
            else:
                xarstr = r' + M $\rightarrow$ '
            # --- get all print on a per tag basis the coe, rxn str
            try:
                rxn_str = ''.join(rdict[rxn][4:])
                if LaTeX:
                    rxn_str = rxn_str.replace('++=', xarstr)
                    rxn_str = rxn_str.replace('+', ' + ')
                    rxn_str = rxn_str.replace('+ =', r' $\rightarrow$ ')
                else:
                    pass
                if debug:
                    print(rxn_str)
                try:
                    rxn_strs += [rxn_str]
                    rxns += [rxn]
                except:
                    rxn_strs = [rxn_str]
                    rxns = [rxn]
            except:
                print(('!'*100, 'ERROR HERE: >{}<  >{}<'.format(rxn, rxn_str)))
        rdict = dict(list(zip(rxns, rxn_strs)))
    return rdict


def rxns_in_pl(wd, spec='LOX', debug=False):
    """
    Extract reactions tracked by p/l family in smvgear

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    debug (bool): legacy debug option, replaced by python logging
    spec (str): prod-loss variable name (from input.geos)

    Returns
    -------
    (list) reaction in production/loss family

    Notes
    -----
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """

    fn = 'smv2.log'
    file_ = open(wd+'/'+fn, 'r')
    if debug:
        print(file_)
    readrxn = False
    # required strings in file line
    conditions = 'Family', 'coefficient', 'rxns', spec
    for row in file_:
        row = row.split()
        if debug:
            print((row, spec, all([i in row for i in conditions])))
        if all([i in row for i in conditions]):
            readrxn = True
        if (len(row) < 1) or ('REACTANTS:' in row):
            readrxn = False
        if readrxn:
            try:
                rxns.append(row)
            except:
                rxns = [row]

    # -- Check that rxns ahave been found?
    if len(rxns) < 1:
        print(('ERROR: No rxns. found for >{}<, correct family?'.format(spec)))
        sys.exit(0)

    # -- remove 'Family'
    rxns = [i for i in rxns if ('Family' not in i)]
    if debug:
        print(('number (len of list) of reacitons: ', len(rxns)))
    n = [int(rxn[1]) for rxn in rxns]
    rxns = [rxn[2:] for rxn in rxns]

    rdict = dict(list(zip(n, rxns)))
    return rdict


def rxn4pl(pls, wd='example/example', rdict=None, reduce_size=True,
           ver='1.7', debug=False):
    """
    Get information on reaction in smvgear from a provide reaction tag

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    ver (str): The GEOS-Chem halogen version that is being used
    debug (bool): legacy debug option, replaced by python logging
    spec (str): prod-loss variable name (from input.geos)
    rdict (dict): dictionary of prod-loss reaction (made by rxn_dict_from_smvlog)

    Returns
    -------
    (list) reaction in production/loss family

    Notes
    -----
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """

    # ---  Get Dict of reaction detail
    if debug:
        print('rxn4pl called')
    if isinstance(rdict, type(None)):
        #   Get Dict of all reactions, Keys = #s
        rdict = rxn_dict_from_smvlog(wd, ver=ver)

    if debug:
        for i in list(rdict.keys()):
            if any([(s_ in ''.join(rdict[i])) for s_ in pls]):
                print((i, 'yes'))

    # --- Indices for
    # reduce dict size
    if reduce_size:
        keys = [i for i in list(rdict.keys()) if
                any([(s_ in ''.join(rdict[i])) for s_ in pls])]

        # re-make dictionary
        rdict = dict(list(zip(keys, [rdict[i] for i in keys])))

    # loop via pl
    keys = np.array([[pl, [i for i in list(rdict.keys())
                           if any([(pl in ''.join(rdict[i]))])][0]] for pl in pls])

    # --- Return as reactions referenced by tag
    return dict(list(zip(keys[:, 0],  [rdict[int(i)] for i in keys[:, 1]])))


def get_indicies_4_fam(tags, fam=False, IO_BrOx2=False, rtnspecs=False,
                       NOy_as_HOx=True, Include_Chlorine=False, debug=False):
    """
    Return indicies (in list form) for in a given family

    Parameters
    ----------
    tags (list): list of tags(str) from globchem.dat
    fam (bool): also return which family the tag belongs to
    IO_BrOx2 (bool): include an addtion tracer for XO+XO reactions?
    Include_Chlorine (bool): include an extra tracer for ClO+XO
    NOy_as_HOx (Boolean): include NOy losses in HOx family.
    debug (bool): legacy debug option, replaced by python logging
    rtnspecs(bool): return tags as part of function output

    Returns
    -------
    (list) or indices for globche.dat/prod and loss tags

    Notes
    -----
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """
    # assign family
#    famsn = [ 'Photolysis','HOx','NOy' ,'Bromine', 'Iodine' ]
    if Include_Chlorine:
        famsn = ['Photolysis', 'HOx', 'Chlorine', 'Bromine', 'Iodine']
    else:
        famsn = ['Photolysis', 'HOx', 'Bromine', 'Iodine']
    fams = []
    for tag in tags:
        fams.append(get_tag_fam(tag))
    # if rm NOy family (treat as NOx)
    if NOy_as_HOx:
        fams = [x if (x != 'NOy') else 'HOx' for x in fams]

    # Create dictionary from tags and fam assignment
    fd = dict(list(zip(tags, fams)))

    # Select tags with assigned family in list ("famsn")
    ll = []
    [ll.append([]) for i in famsn]
    for n, tag in enumerate(tags):
        for fn in range(len(famsn)):
            if fd[tag] == famsn[fn]:
                ll[fn].append(n)

    # Consider Ox loss for
    if fam:
        #    if False: # Kludge for test.
        # Kludge - to allow for counting Ox loss via XO +XO 50/50 between fams,
        # ( add extra loss tag. )
        if IO_BrOx2:
            #        if False:
            # Add extra tag for reaction ( IO + BrO )
            ll[famsn.index('Bromine')].append(max([max(i) for i in ll])+1)
            fams = fams + ['Bromine']
            tags = tags + ['LO3_24']
        if Include_Chlorine:
            #        if False:
            # Add extra tag for reaction ( ClO + BrO )
            ll[famsn.index('Chlorine')].append(max([max(i) for i in ll])+1)
            fams = fams + ['Chlorine']
            tags = tags + ['LO3_82']
            # Add extra tag for reaction  ( ClO + IO )
            ll[famsn.index('Chlorine')].append(max([max(i) for i in ll])+1)
            fams = fams + ['Chlorine']
            tags = tags + ['LO3_87']

        if rtnspecs:
            return ll, fams, tags
        else:
            return ll, fams
    else:
        return ll


def get_p_l_tags(rxns, debug=False):
    """
    Get p/l tags for a given smvgear reaction

    Parameters
    ----------
    rxns (list): list of rxn numbers(int) from smvgear
    debug (bool): legacy debug option, replaced by python logging

    Returns
    -------
    (list) of tags for globchem.dat/prod and loss tags

    Notes
    -----
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.

    """

    # (PD??, RD??, LO3_??, PO3_??, LR??)
    prefixs = 'PD', 'RD', 'PO3', 'LO3', 'LR'

    for rxn in rxns:
        if debug:
            print([i for i in rxn if any([(x in i) for x in prefixs])])
        tags = [i for i in rxn if any([(x in i) for x in prefixs])]

        try:
            tagsl.append(tags)
        except:
            tagsl = [tags]

    return tagsl


def p_l_species_input_geos(wd, ver='1.7', rm_multiple_tagged_rxs=False,
                           debug=False):
    """
    Extract prod/loss species (input.geos) and reaction tags (globchem.dat)

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    debug (bool): legacy debug option, replaced by python logging
    ver (str): The GEOS-Chem halogen version that is being used
    rm_multiple_tagged_rxs(bool): only return one tag per rxn.

    Returns
    -------
    (list) globchem.dat tags and prod/loss ("PD") vars from input.geos

    Notes
    -----
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """
    # find and open input.geos file
    fn = glob.glob(wd+'/*input.geos*')[0]
    if any([(i in fn) for i in ('~', '#')]):
        print(('Trying next "input.geos" file - as FAIL for :', fn))
        fn = glob.glob(wd+'/*input.geos*')[1]

    if debug:
        print(('p_l_species_input_geos called using : ', wd, fn))
    file_ = open(fn, 'r')

    # Read in just the prod loss section
    strs_in_1st_line = 'Number', 'of', 'P/L', 'families'
    section_line_divider = '------------------------+----------' + \
        '--------------------------------------------'
    readrxn = False
    for row in file_:
        row = row.split()

        # once at prod/loss section, start added to list
        if all([i in row for i in strs_in_1st_line]):
            readrxn = True

        # if not at end of prod/loss section, add to list
        if section_line_divider in row:
            readrxn = False
        if readrxn:
            try:
                rxns.append(row)
            except:
                rxns = [row]

    # -- Only consider 'Family' ( no headers e.g. 'families' )
    rxns = [i for i in rxns if ('families' not in i)]
    rxns = [[i.replace(':', '') for i in r] for r in rxns]

    # Kludge, adjust for extra space 12-99
    # ( This is no longer required for 1.7 + )
    if ver == '1.6':
        [i.pop(0) for i in rxns if ('th' not in i[0])]

    # Extract just PD (input.geos) and vars (globchem.dat vars )
    PD = [rxn[4] for rxn in rxns]
    vars = [rxn[5:] for rxn in rxns]
    if debug:
        print((rxns, PD, vars, ver))

    # remove p/l with muliple values ( start from 12th input) - Kludge?
    if rm_multiple_tagged_rxs:
        PD, vars = [i[11:] for i in (PD, vars)]
        vars = [i[0] for i in vars]

    return PD, vars


def tags_from_smvlog(wd):  # , spec='LOX' ):
    """
    Get all active p/l tags in smvgear ( from smv2.log )

    Notes
    -----
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """
    fn = 'smv2.log'
    file_ = open(wd+'/'+fn, 'r')
    readrxn = False
    for row in file_:
        row = row.split()
        if all([(i in row) for i in ['NBR', 'NAME', 'MW', 'BKGAS(VMRAT)']]):
            readrxn = True
        if len(row) < 1:
            readrxn = False
        if readrxn:
            try:
                rxns.append(row)
            except:
                rxns = [row]

    # Remove 'NMBR'
    rxns = [i for i in rxns if ('NBR' not in i)]
    rxns = [rxn[1] for rxn in rxns]
    # Only consider tags and return
    inclusions = ('PD', 'RD', 'PO3', 'LO3', 'LR')
    return [i for i in rxns if any([x in i for x in inclusions])]


def PDs_from_smvlog(wd, spec='LOX'):
    """
    Get all active PDs tags in smvgear ( from smv2.log )

    Notes
    -----
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """
    fn = 'smv2.log'
    file_ = open(wd+'/'+fn, 'r')
    readrxn = False
    leniency = 0
    entries_in_title = ['Families', 'for', 'prod', 'or', 'loss', 'output:']
    for row in file_:
        row = row.split()
        if all([(i in row) for i in entries_in_title]):
            readrxn = True
            leniency = 1
        if len(row) < 1:
            if leniency < 0:
                readrxn = False
            leniency -= 1
        if readrxn:
            try:
                rxns.append(row)
            except:
                rxns = [row]

    # Remove 'NMBR'
    exceptions = ['SPECIES', '='*79, 'Families']
    rxns = [i for i in rxns if all([(ii not in i) for ii in exceptions])]
    rxns = [j for k in rxns for j in k]
    return rxns


def rxns4tag(tag, rdict=None, ver='1.7', wd=None):
    """
    Get a list of all reactions with a given p/l tag

    Notes
    -----
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """
    # --- get reaction dictionary
    if isinstance(rdict, type(None)):
        rdict = rxn_dict_from_smvlog(wd, ver=ver)

    # --- Caveats -
    # to adapt for long line errors in fortran written output
    errs = ['LO3_36']  # + ['LO3_87']
    cerrs = ['RD95']  # + ['LR48']
    # To account for reaction where not all channels result in Ox loss
    errs += ['RD48']
    cerrs += ['LO3_87']
    if any([(tag == i) for i in errs]):
        tag = cerrs[errs.index(tag)]

    # -- loop reactions, if tag in reaction return reaction
    rxns = []
    for n, rxn in enumerate(rdict.values()):

        expanded_rxn_str = [i.split('+') for i in rxn]
        expanded_rxn_str = [
            item for sublist in expanded_rxn_str for item in sublist]

        # ( Issue) Why endswith? Restore to use if contains any tag
#        if any( [ (i.endswith(tag) ) for i in rxn]):
        # This is because otherwise 'LR10' would be read as 'LR100'
#        if any( [tag in i for i in rxn]): # <= This will lead to false +ve
        # However, fortran print statment err for (  LO3_87 )
        if any([i.endswith(tag) for i in expanded_rxn_str]):
            rxns.append([list(rdict.keys())[n]] + rxn)

    return rxns


def get_tag_details(wd, tag=None, PDs=None,  rdict=None, PHOTOPROCESS=None,
                    ver='1.7',
                    LaTeX=False, print_details=False, debug=False):
    """
    Retriveve prod/loss tag details from smv.log
        ( rxn number + reaction description)

    Notes
    -----
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """

    # what is the number of the first photolysis reaction?
    if isinstance(PHOTOPROCESS, type(None)):
        PHOTOPROCESS = {
            '1.6': 457,  '1.6.2': 452, '1.7': 467, '2.0': 555, '3.0': 547, '4.0': 547
        }[ver]

    # ---  get all reactions tags are active in smv.log
    if isinstance(rdict, type(None)):
        rdict = rxn_dict_from_smvlog(wd, ver=ver)
    trxns = rxns4tag(tag, wd=wd, rdict=rdict)

    # --- get all print on a per tag basis the coe, rxn str
    try:
        rxn_str = ''.join(trxns[0][5:9])

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
            print(rxn_str)
        dets = [tag, trxns[0][0], rxn_str]
    except:
        print(('!'*100, 'ERROR HERE: >{}< >{}<'.format(tag, trxns)))
    if print_details:
        print(dets)
    # --- return a dictionary of all active tagged reactions details by tag : PD, number, rxn str, coeeffiecn
    else:
        return dets


def get_rxn_Coe(wd, num, tag, nums=None, rxns=None, tags=None, Coe=None, spec='LOX',
                ver='1.6', debug=False):
    """
    Retrieve given reaction coefficient for smvgear (from smv2.log)

    Notes
    -----
     - This is no longer the case. However, previously if using dev.
    Iy scheme, then the listed values from fuction Ox_in_species() will be used.
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """

    # --- get dictionaries for reactions within
    if all([(i == None) for i in (nums, rxns, tags, Coe)]):
        nums, rxns, tags, Coe = prod_loss_4_spec(wd,  spec, all_clean=True,
                                                 ver=ver)
    if debug:
        print((nums, Coe))

    # Pull reaction coefficient  from dictionary
    Coe_dict = dict(list(zip(nums, Coe)))
    Coe = float(Coe_dict[num])
    # Consider all change positive - Kludge
    # ( This is due to the assignment approach, where P=prod, L=loss )
    if ('P' not in tag):
        Coe = Coe*-1.0

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


def get_pldict_reactants(pl_dict=None, only_rtn_tracers=True, rm_OH=True,
                         rm_Cl=True,
                         tags=None, debug=False):
    """
    Get reactant from smv2.log dictionary

    Notes
    -----
     - some reactants are not tracers ( see list: non_TRAs )
    ( to remove these set only_rtn_tracers=True )
     - to remove OH from reactant list set rm_OH=True
     - This function is useful, but update to GEOS-Chem flexchem ( in >v11)
    will make it redundant and therefore this is not being maintained.
    """
#    non_TRAs = ['CH4', '', 'ETHLN', 'ISOPND', 'E', 'M', 'HCO', 'MVKN', 'ACTA']

    # Get reaction strings
    if isinstance(tags, type(None)):
        tags = list(pl_dict.keys())
    strs = [pl_dict[i][1] for i in tags]
    if debug:
        print(('rxn strs: ',  strs, len(strs)))
    # remove arrows from reactions
    strs = [i.replace('+M=', '+=').replace('+O2=', '+=').replace('+N2=', '+=')
            for i in strs]
    if debug:
        print((strs, len(strs)))
    # select reactants
    strs = [i.split('+=')[0] for i in strs]
    if debug:
        print((strs, len(strs)))
    if rm_OH:     # remove OH from reaction strings
        strs = [i.replace('+OH', '+') for i in strs]
        for n, str in enumerate(strs):
            if str.startswith('OH+'):
                strs[n] = str.replace('OH', '+')
    if rm_Cl:     # remove Cl from reaction strings
        strs = [i.replace('+Cl', '+') for i in strs]
        for n, str in enumerate(strs):
            if str.startswith('Cl+'):
                strs[n] = str.replace('Cl+', '+')

    # remove "+" punctuation
    strs = [i.replace('+', '').strip() for i in strs]
    if debug:
        print((strs, len(strs)))

    return strs


def PLO3_to_PD(PL, fp=True, wd=None, ver='1.6', res='4x5', verbose=False,
               debug=False):
    """
    Converts globchem.dat tracer to PD/LD from prod/loss diag in input.geos

    Notes
    -----
     - 'fp' option is now obselete.
     - UPDATED NEEDED: MORE DETAILED DESCRIPT.
    """

    if verbose:
        print(('PLO3_to_PD called for wd = ', wd))

    versions = [
        '1.3', '1.4', '1.5', '1.6', '1.6.1', '1.6.2', '1.6.3', '1.7', '2.0', '3.0', '4.0']
    if any([(ver == i) for i in versions]):

        if isinstance(wd, type(None)):
            print('WARNING: Using MUTD wd')
            wd = MUTD_runs(ver=ver, res=res, debug=debug)[0]

        # Get list of assigned PDs for vars
        PDs, vars = p_l_species_input_geos(wd, ver=ver,
                                           rm_multiple_tagged_rxs=True,
                                           debug=debug)

        # Add other (non 'PD') vars for ease of processing
        non_PDs = [
            'PIOx', 'iLOX', 'LIOx', 'iPOX', 'POX', 'LOX', 'LOx', 'L_Iy', 'LOH',
            'LCl', 'POH', 'PCl', 'P_Iy', 'L_Bry', 'P_Bry', 'L_Cly', 'P_Cly',
            'LSBrA', 'LSBrC', 'PSBrA', 'PSBrC'
        ]
        vars += non_PDs
        PDs += non_PDs

        if debug:
            print((vars, PDs))

        return dict(list(zip(vars, PDs)))[PL]
    else:
        print('update programme - manual PLdict now obsolete. ')


def get_pl_dict(wd, spec='LOX', rmx2=False, ver='1.7',
                rm_redundant_ClBrI_tags=False,
                debug=False):
    """
    Get reaction IDs for each rxn. in spec (p/l, e.g. LOX). This is the driver for
    the prod/loss programmes

    Notes
    -----
     - UPDATED NEEDED: MORE DETAILED DESCRIPT.
    """
    if debug:
        print(('get_pl_dict called for ', ver, spec, wd))

    # Extract details on reactio in p/l family
    nums, rxns, tags, Coe = prod_loss_4_spec(wd,  spec, all_clean=True,
                                             ver=ver, debug=debug)

    # Make a dictionary of coeffiecnts of reaction
    Coe_dict = dict(list(zip(nums, Coe)))

    # unpack for mulutple tags of same reactions, then get details
    unpacked_tags = [j for k in tags for j in k]

    # Kludge - remove 'PO3_10' temporarily from dictionary as these contain
    # partial tag names (  Fortran print statment cut off )
    unpacked_tags = [i for i in unpacked_tags if ('PO3_10' not in i)]

    if debug:
        print(unpacked_tags)

    details = [get_tag_details(wd, tag, ver=ver) for tag in unpacked_tags]

    # Kludge - 1 rxn missing from POx tracking? - 999  + "'ISOPND+OH', '+', '=1.0ISOPND'"
    [details.pop(n) for n, i in enumerate(details) if i[1] == 364]
    ind = [n for n, i in enumerate(nums) if i == 354]

    # Get coefficients for change in reaction family
    # NOTE: This does not exclude adjustment for non unity globchem.dat tags
    Coes = [get_rxn_Coe(wd, d[1], unpacked_tags[n], nums=nums,
                        rxns=rxns, tags=tags, Coe=Coe, spec=spec, debug=debug)
            for n, d in enumerate(details)]

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
            ['RD95', 'LO3_36'],   \
            # OIO +hv =>  => AERI (LOSS) ( either RD67 or LO3_35 fine. )
            ['RD67', 'LO3_35'],
            # HOBr + hv=>  ( either LR25 or LO3_84 fine. )
            ['LR25', 'LO3_84']
        ]
        # Use Justin's 'LR??'/Johan's JTO1 tags in preference to 'LO3??' tags
        # This is an Kludge that only is necessary for "NOHAL" runs
        if rm_redundant_ClBrI_tags:
            d += [
                ['LR6', 'LO3_73'], ['LR5',  'LO3_73'], ['LR10', 'LO3_74'],
                ['LR25', 'LO3_84'], ['LO3_82', 'LO3_82'], ['LO3_76', 'LO3_76'],
                ['LO3_75', 'LO3_75']
            ]
            # make sure 'LO3_73' is dropped regardless of index selection
            d += [['LO3_73', 'LO3_73'],  ['LO3_73', 'LO3_73'],
                  ['LO3_74', 'LO3_74'], ['LO3_84', 'LO3_84']
                  ]

        # Either one can be removed as currently two will be present and both #
        # are  equally weighted by use of the Ox_in_species diagnostic ...
        # ( However, this  only if all are equiv. )
        if rm_redundant_ClBrI_tags:
            # Drop the 2nd element in "d" list #
            d = [i[1] for i in d]
        else:
            # Drop the 1st element in "d" list # <= MUST USE for Halogen sims.
            d = [i[0] for i in d]

        ind = [n for n, i in enumerate(details) if any([i[0] == ii
                                                        for ii in d])]
        # make sure indices ('ind') are only removed once
        ind = list(sorted(set(ind)))

        if debug:
            print((d, ind, [len(i) for i in (details, Coes)],
                   [[i[0] for i in details][ii] for ii in ind][::-1]))
        # If cases have been found, remove these
        if len(ind) > 0:
            [l.pop(i) for i in ind[::-1] for l in (details, Coes)]
        if debug:
            print([len(i) for i in (details, Coes)])

    # return a dictionary indexed by p/l tracer, with rxn #,
    # reaction str and Coe of rxn.
    return dict(list(zip([i[0] for i in details], [i[1:] + [Coes[n]]
                                                   for n, i in enumerate(details)])))


def prod_loss_4_spec(wd, fam, all_clean=True, ver='1.7', debug=False):
    """
    Retrieve reaction numbers for family of tags

    Notes
    -----
     - coefficecents ("Coe") returned are for the family (e.g LOX)
    within a reaciton. ( aka not for the tag )
     -  UPDATED NEEDED: MORE DETAILED DESCRIPT.
    """

    # ---  Get Dict of all reactions, Keys = #s
    rdict = rxn_dict_from_smvlog(wd, ver=ver)

    # ---  Get reaction # tracked by p/l diag for spec and coefficient.
    rxns = rxns_in_pl(wd, fam)
    nums = list(rxns.keys())
    Coe = [rxn[-1] for rxn in list(rxns.values())]

    # --- get all details from full reaction dictionary
    rxns = [rdict[i] for i in nums]

    # --- get tags for tracked reactions, state where reactions are un tracked
    tags = get_p_l_tags(rxns)

    # --- cleaned tags
    if all_clean:
        tags = [[re.sub('\+\d.\d', '',  i) for i in u] for u in tags]
        tags = [[re.sub('\=\d.\d', '',  i) for i in u] for u in tags]

        # -- remove erroneous read/  Kludge on val
        # ---  Fortran write error leads to combination of species at the
        #  end of long line of a chemical reaction in globchem.dat
        if debug:
            print([i[:3] for i in (nums, rxns, tags, Coe)])
            print([len(i) for i in (nums, rxns, tags, Coe)])

        # LO3_36RD95 is present twice as this reaction is present 2 times in the code
        # update 16 01 11: LO3_36RD95 is now present x3 in version 3.0
        # ( due split uptake in iodine to aerosol )
        errs = [
            'LO3_36RD95', 'LO3_36RD95', 'LO3_36RD95',
            'ISOPNDPO3_50',
            'ISOPNDLR40',
            'LO3_30LR42',
            'LO3_30LR43',
            'LO3_39LR46',
            'LO3_39LR47',
            'LO3_87LR48', 'LO3_87LR48',
            'LISOPOLR86',
            'LO3_50LR113',
            'LO3_50LR114',
            # final v3.0 runs
            'LO3_65LR126'
        ]
        cerrs = [
            ['LO3_36', 'RD95'], ['LO3_36', 'RD95'], ['LO3_36', 'RD95'],
            ['PO3_50'],
            ['LR40'],
            ['LO3_30', 'LR42'],
            ['LO3_30', 'LR43'],
            ['LO3_39', 'LR46'],
            ['LO3_39', 'LR47'],
            ['LO3_87'], ['LO3_87'],
            ['LR86'], \
            #        ['LO3_50', 'LR113' ], \
            #        ['LO3_50', 'LR114' ]
            # revserse LR tags also ( LO3_50 now redundant? ) - NO
            ['LR113', 'LO3_50'], \
            ['LR114', 'LO3_50'],
            #        ['LR113'],
            #        ['LR114'],
            # final v3.0 runs
            ['LR126', 'LO3_65'],

        ]
#        errs = ['LO3_36RD95' , 'ISOPNDPO3_50', 'ISOPNDLR40']
#        cerrs = [ ['RD95'], ['PO3_50'], ['LR40'] ]
        for n, e in enumerate(errs):
            try:
                # Get index of erroneous rxn tag
                ind = [nn for nn, i in enumerate(tags) if
                       any([(e in ii) for ii in i])][0]
                # Extract vars for a given index
                vars = [i[ind] for i in (nums, rxns, tags, Coe)]
                if debug:
                    print((3, [i[-1] for i in (nums, rxns, tags, Coe)], vars,
                           [len(i) for i in (nums, rxns, tags, Coe)]))
                # remove index ( "ind" ) value from nums, rxns, tags, and Coe
                [i.pop(ind) for i in (nums, rxns, tags, Coe)]

                # Add the cerrs values on the end
                if debug:
                    print((4, [i[-1] for i in (nums, rxns, tags, Coe)],
                           [len(i) for i in (nums, rxns, tags, Coe)]))
                nums += [vars[0]]
                rxns += [vars[1]]
                tags += [cerrs[n]]
                Coe += [vars[-1]]

                if debug:
                    print((6, [i[-1] for i in (nums, rxns, tags, Coe)],
                           [len(i) for i in (nums, rxns, tags, Coe)]))
                    print(('->'*30,  'SUCCESS >{}<  >{}<'.format(n, e)))
            except:
                print(('>'*50, 'FAIL (NOT REPLACED) >{}< >{}<'.format(n, e)))

    # KLUDGE! - rm empty list values of ones that contain errs
#    ind = [ n for n,i in enumerate(tags) if ( (len(i)==0) or (i[0] in errs) ) ]
#    [ [ l.pop(i) for i in sorted(ind)[::-1] ] for  l in nums, rxns, tags, Coe ]
    if debug:
        print(tags)

#    print '1'*300, tags

    return nums, rxns, tags, Coe
