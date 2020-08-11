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


def main(wd=None, CODE_wd=None):
    """
    Driver for Ox loss analysis via KPP in GEOS-Chem from NetCDF output

    Notes
    -----
     - comment/uncomment functions as required
    """
    # Manually set locations of model output with the Ox loss diagnostic
    root = '/users/ts551/scratch/GC/'
    CODE_wd = root+'/Code/Code.v11-02_Cl_v3_0/'
    wd = root + '/rundirs/'
    wd += 'merra2_4x5_standard.v12.9.1.BASE.Oi.MacDonald2014.tagged/'
#    Mechanism = 'Tropchem'
    Mechanism = 'Standard'
    # Get a dictionary of information about the model run
    # NOTE: Should this approach be retired? Yes.
    full_vert_grid = True
    VarDict = AC.get_default_variable_dict(full_vert_grid=full_vert_grid,
                                           wd=wd)
    # Now Get the StateMet object...



    # Get all the necessary data as as a dictionary object
    Ox_loss_dict = AC.get_Ox_loss_dicts(
        wd=wd, CODE_wd=CODE_wd, Mechanism=Mechanism)
    # Plot vertical odd oxygen (Ox) loss via route (chemical family)
    plot_vertical_fam_loss_by_route(
        Ox_loss_dict=Ox_loss_dict, Mechanism=Mechanism)
    # Analyse odd oxygen (Ox) loss budget via route (chemical family)
    calc_fam_loss_by_route(Ox_loss_dict=Ox_loss_dict, Mechanism=Mechanism)
