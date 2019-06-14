#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Plot Vertical loss of Ox by family. Function can also print budget to csv.

This is an example script to use AC_tools KPP mechanism parsing/tagging functions. The

python AC_tools/Scripts/analyse_vertical_Ox_loss_by_route_in_KPP_mech.py <working directory with NetCDF of GEOS-Chem output>

"""
import AC_tools as AC
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from netCDF4 import Dataset
import os


def main( wd=None, CODE_wd=None ):
    """
    Driver for analysis of LOx via KPP in GEOS-Chem

    Notes
    -----
     - comment/uncommet functions as required
    """
    # Manually let locations of Ox loss here
    root = '/users/ts551/scratch/GC/'
    CODE_wd = root+'/Code/Code.v11-02_Cl_v3_0/'
    wd = root+'rundirs/GC_v11_2d_plus_Clv3/geosfp_4x5_tropchem_Cl.v3_0.1year.2016.tagged/'
    Mechanism = 'Tropchem'
    # Get all the necessary data as as a dictionary object
    Ox_loss_dict = AC.get_Ox_loss_dicts(wd=wd, CODE_wd=CODE_wd, Mechanism=Mechanism)
    # Plot vertical odd oxygen (Ox) loss via route (chemical family)
    plot_vertical_fam_loss_by_route(Ox_loss_dict=Ox_loss_dict, Mechanism=Mechanism)
    # Analyse odd oxygen (Ox) loss budget via route (chemical family)
    calc_fam_loss_by_route(Ox_loss_dict=Ox_loss_dict, Mechanism=Mechanism)


def plot_vertical_fam_loss_by_route(fam='LOx', ref_spec='O3',
                                    wd=None, Mechanism='Halogens', rm_strat=False,
                                    weight_by_molecs=True, CODE_wd=None,
                                    full_vertical_grid=True, dpi=320, suffix='',
                                    save_plot=True, show_plot=False,
                                    limit_plotted_alititude=True, lw=16,
                                    Ox_loss_dict=None, fontsize=10, cmap=plt.cm.jet,
                                    verbose=True, debug=False):
    """
    Plot vertical odd oxygen (Ox) loss via route (chemical family)

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
    limit_plotted_alititude (bool): limit the plotted vertical extend to troposphere
    suffix (str): suffix in filename for saved plot
    dpi (int): resolution to use for saved image (dots per square inch)
    Ox_loss_dict (dict), dictionary of Ox loss variables/data (from get_Ox_loss_dicts)

    Returns
    -------
    (None)

    Notes
    -----
     - AC_tools includes equivlent functions for smvgear mechanisms
    """
    # - Local variables/ Plot extraction / Settings
    if isinstance(Ox_loss_dict, type(None)):
        Ox_loss_dict = AC.get_Ox_loss_dicts(wd=wd, CODE_wd=CODE_wd, fam=fam,
                                           ref_spec=ref_spec,
                                           Mechanism=Mechanism, rm_strat=rm_strat,
                                           weight_by_molecs=weight_by_molecs,
                                           full_vertical_grid=full_vertical_grid,
                                           )
    # extract variables from data/variable dictionary
    sorted_fam_names = Ox_loss_dict['sorted_fam_names']
    fam_dict = Ox_loss_dict['fam_dict']
    ars = Ox_loss_dict['ars']
    RR_dict_fam_stioch = Ox_loss_dict['RR_dict_fam_stioch']
    RR_dict = Ox_loss_dict['RR_dict']
    tags2_rxn_num = Ox_loss_dict['tags2_rxn_num']
    tags = Ox_loss_dict['tags']
    tags_dict = Ox_loss_dict['tags_dict']
    Data_rc = Ox_loss_dict['Data_rc']
    # Combine to a single array
    arr = np.array(ars)
    if debug:
        print((arr.shape))
    # - Process data for plotting
    fam_tag = [fam_dict[i] for i in tags]
    fam_ars = []
    for fam_ in sorted_fam_names:
        # Get indices for routes of family
        fam_ind = [n for n, i in enumerate(fam_tag) if (i == fam_)]
        if debug:
            print((fam_ind, len(fam_ind)))
        # Select these ...
        fam_ars += [arr[fam_ind, ...]]
    # Recombine and sum by family...
    if debug:
        print(([i.shape for i in fam_ars], len(fam_ars)))
    arr = np.array([i.sum(axis=0) for i in fam_ars])
    if debug:
        print((arr.shape))
    # - Plot up as a stack-plot...
    # Normalise to total and conver to % (*100)
    arr = (arr / arr.sum(axis=0)) * 100
    # Add zeros array to beginning (for stack/area plot )
    arr_ = np.vstack((np.zeros((1, arr.shape[-1])), arr))
    # Setup figure
    fig, ax = plt.subplots(figsize=(9, 6), dpi=dpi, \
                           facecolor='w', edgecolor='w')
    # Plot by family
    for n, label in enumerate(sorted_fam_names):
        # Print out some summary stats
        if verbose:
            print(n, label, arr[:n, 0].sum(axis=0), arr[:n+1, 0].sum(axis=0), end=' ')
            print(arr[:n, :].sum(), arr[:n+1, :].sum())
            print([i.shape for i in (Data_rc['alt'], arr)])
        # Fill between X
        plt.fill_betweenx(Data_rc['alt'], arr[:n, :].sum(axis=0),
                          arr[:n+1, :].sum(axis=0),
                          color=cmap(1.*n/len(sorted_fam_names)))
        # Plot the line too
        plt.plot(arr[:n, :].sum(axis=0), Data_rc['alt'], label=label,
                 color=cmap(1.*n/len(sorted_fam_names)), alpha=0,
                 lw=lw,)
    # Beautify the plot
    plt.xlim(0, 100)
    xlabel = '% of total O$_{\\rm x}$ loss'
    plt.xlabel(xlabel, fontsize=fontsize*.75)
    plt.yticks(fontsize=fontsize*.75)
    plt.xticks(fontsize=fontsize*.75)
    plt.ylabel('Altitude (km)', fontsize=fontsize*.75)
    leg = plt.legend(loc='upper center', fontsize=fontsize)
    # Update lengnd line sizes ( + update line sizes)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(lw/2)
        legobj.set_alpha(1)
    plt.ylim(Data_rc['alt'][0], Data_rc['alt'][-1])
    # Limit plot y axis to 12km?
    if limit_plotted_alititude:
        plt.ylim(Data_rc['alt'][0], 12)
    # Show plot or save?
    if save_plot:
        filename = 'Ox_loss_plot_by_vertical_{}_{}'.format(Mechanism, suffix)
        plt.savefig(filename, dpi=dpi)
    if show_plot:
        plt.show()


def calc_fam_loss_by_route(wd=None, fam='LOx', ref_spec='O3',
                           rm_strat=True, Mechanism='Halogens', Ox_loss_dict=None,
                           weight_by_molecs=False, full_vertical_grid=False,
                           CODE_wd=None, verbose=True, debug=False):
    """
    Build an Ox budget table like table 4 in Sherwen et al 2016b

    Parameters
    -------
    fam (str): tagged family to track (already compiled in KPP mechanism)
    ref_spec (str): reference species to normalise to
    wd (str): working directory ("wd") of model output
    CODE_wd (str): root of code directory containing the tagged KPP mechanism
    rm_strat (bool): (fractionally) replace values in statosphere with zeros
    Ox_loss_dict (dict), dictionary of Ox loss variables/data (from get_Ox_loss_dicts)
    Mechanism (str): name of the KPP mechanism (and folder) of model output
    weight_by_molecs (bool): weight grid boxes by number of molecules
    full_vertical_grid (bool): use the full vertical grid for analysis
    debug, verbose (bool): switches to turn on/set verbosity of output to screen

    Returns
    -------
    (None)

    Notes
    -----
     - AC_tools includes equivlent functions for smvgear mechanisms
    """
    # - Local variables/ Plot extraction / Settings
    if isinstance(Ox_loss_dict, type(None)):
        Ox_loss_dict = AC.get_Ox_loss_dicts(wd=wd, CODE_wd=CODE_wd, fam=fam,
                                            ref_spec=ref_spec,
                                            Mechanism=Mechanism, rm_strat=rm_strat,
                                            weight_by_molecs=weight_by_molecs,
                                            full_vertical_grid=full_vertical_grid,
                                            )
    # Extract variables from data/variable dictionary
    fam_dict = Ox_loss_dict['fam_dict']
    ars = Ox_loss_dict['ars']
    RR_dict_fam_stioch = Ox_loss_dict['RR_dict_fam_stioch']
    RR_dict = Ox_loss_dict['RR_dict']
    tags2_rxn_num = Ox_loss_dict['tags2_rxn_num']
    tags = Ox_loss_dict['tags']
    tags_dict = Ox_loss_dict['tags_dict']
    halogen_fams = Ox_loss_dict['halogen_fams']
    # --- Do analysis on model output
    # Sum the total mass fluxes for each reaction
    ars = [i.sum() for i in ars]
    # Sum all the Ox loss routes
    total = np.array(ars).sum()
    # Create a dictionary of values of interest
    dict_ = {
        'Total flux': [i.sum(axis=0) for i in ars],
        'Total of flux (%)': [i.sum(axis=0)/total*100 for i in ars],
        'Family': [fam_dict[i] for i in tags],
        'rxn #': [tags2_rxn_num[i] for i in tags],
        'rxn str': [RR_dict[tags2_rxn_num[i]] for i in tags],
        'stoich': [RR_dict_fam_stioch[tags2_rxn_num[i]] for i in tags],
        'tags': tags,
    }
    # Create pandas dataframe
    df = pd.DataFrame(dict_)
    # Sort the data and have a look...
    df = df.sort_values('Total flux', ascending=False)
    if debug:
        print(df.head())
    # Sort values again and save...
    df = df.sort_values(['Family', 'Total flux'], ascending=False)
    if debug:
        print(df.head())
    df.to_csv('Ox_loss_budget_by_rxn_for_{}_mechanism.csv'.format(Mechanism))
    # Now select the most important routes
    grp = df[['Family', 'Total flux']].groupby('Family')
    total = grp.sum().sum()
    # Print the contribution by family to screen
    print((grp.sum() / total * 100))
    # Print the contribution of all the halogen routes
    hal_LOx = (grp.sum().T[halogen_fams].sum().sum() / total * 100).values[0]
    if verbose:
        print(('Total contribution of halogens is: {:.2f} %'.format(hal_LOx)))
    # Add Halogen total and general total to DataFrame
    dfFam = grp.sum().T
    dfFam['Total'] =dfFam.sum().sum()
    dfFam['Halogens'] = dfFam[halogen_fams].sum().sum()
    # Update units to Tg O3
    dfFam = dfFam.T /1E12
    # return dictionaries of LOx by reaction or by family (in Tg O3)
    if rtn_by_rxn:
        return df /1E12
    if rtn_by_fam:
        return dfFam


if __name__ == "__main__":
    main()
