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


def main(wd=None, CODE_wd=None, verbose=False, debug=False):
    """
    Driver for Ox loss analysis via KPP in GEOS-Chem from NetCDF output

    Notes
    -----
     - comment/uncomment functions as required
    """
    # - Local variables
    fam = 'LOx'
    ref_spec = 'O3'
    # Manually set locations of model output with the Ox loss diagnostic
    root = '/users/ts551/scratch/GC/'
    CODE_wd = root + '/Code/Code.BleedingEdge/'
    wd = root + '/rundirs/'
    wd += 'merra2_4x5_standard.v12.9.1.BASE.Oi.MacDonald2014.tagged/'
    wd += '/OutputDir/'
#    Mechanism = 'Tropchem'
    Mechanism = 'Standard'
    # Get locations of model output/core
    assert os.path.exists(wd), 'working directory not found @: {}'.format(wd)
    CODE_wd = '/{}/KPP/{}/'.format(CODE_wd, Mechanism)
    assert os.path.exists(CODE_wd), 'code directory not found @: ' + CODE_wd

    # - Model output from tagged mechanism run

    # Get all the necessary data as as a dictionary object?
    # NOTE: this is a previous approach using dictionaries

    # Get a dictionary of information about the model run
    # NOTE: Should this approach be retired? Yes.
#    full_vert_grid = True
#    VarDict = AC.get_default_variable_dict(full_vert_grid=full_vert_grid,
#                                           wd=wd)
    # Now Get the StateMet object... for time in troposphere diagnostic
    StateMet = AC.get_StateMet_ds(wd=wd)

    # - KPP mechanism
    # - Analyse the Ox loss budget's numerical terms
    # Get the dictionary of the KPP mechanism.
    Ox_fam_dict = AC.get_Ox_fam_dicts(fam=fam, ref_spec=ref_spec,
                                      Mechanism=Mechanism,
#                                      tag_prefix=tag_prefix,
                                      wd=wd, CODE_wd=CODE_wd,
                                      StateMet=StateMet,
                                      rm_strat=True,
                                      weight_by_molecs=True,
                                      )

    LatLonAlt_dict = AC.gchemgrid(rtn_dict=True)
    alt_array = LatLonAlt_dict['c_km_geos5']
    # Plot vertical odd oxygen (Ox) loss via route (chemical family)
    suffix = 'v12.9.1_.png'
    AC.plot_vertical_fam_loss_by_route(Ox_fam_dict=Ox_fam_dict,
                                       alt_array=alt_array,
                                       Mechanism=Mechanism,
                                       suffix=suffix)

    # - Analyse the Ox loss budget's numerical terms
    # Get the dictionary of the KPP mechanism.
    Ox_fam_dict = AC.get_Ox_fam_dicts(fam=fam, ref_spec=ref_spec,
                                      Mechanism=Mechanism,
#                                      tag_prefix=tag_prefix,
                                      wd=wd, CODE_wd=CODE_wd,
                                      StateMet=StateMet,
                                      rm_strat=True,
                                      weight_by_molecs=False,
                                      )


    # Analyse odd oxygen (Ox) loss budget via route (chemical family)
    suffix = 'v12.9.1'
    df = AC.calc_fam_loss_by_route(Ox_fam_dict=Ox_fam_dict,
                                   Mechanism=Mechanism,
                                   suffix=suffix)







