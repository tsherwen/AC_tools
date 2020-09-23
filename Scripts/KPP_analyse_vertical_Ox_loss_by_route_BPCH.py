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
    Driver for Ox loss analysis via KPP in GEOS-Chem from bpch output

    Notes
    -----
     - comment/uncomment functions as required
    """
    # Manually let locations of Ox loss here
    root = '/users/ts551/scratch/GC/'
    CODE_wd = root+'/Code/Code.v11-02_Cl_v3_0/'
    wd = root+'rundirs/GC_v11_2d_plus_Clv3/geosfp_4x5_tropchem_Cl.v3_0.1year.2016.tagged/'
    Mechanism = 'Tropchem'
    # Get all the necessary data as as a dictionary object
    Ox_fam_dict = AC.get_Ox_fam_dicts_BPCH(wd=wd, CODE_wd=CODE_wd,
                                           weight_by_molecs=True,
                                           Mechanism=Mechanism)
    # Plot vertical odd oxygen (Ox) loss via route (chemical family)
    Data_rc = Ox_fam_dict['Data_rc']
    alt_array = Data_rc['alt']
    AC.plot_vertical_fam_loss_by_route(Ox_fam_dict=Ox_fam_dict,
                                       Mechanism=Mechanism,
                                       alt_array=alt_array)

    # Get all the necessary data as as a dictionary object
    # (Not weighted by molecules)
    Ox_fam_dict = AC.get_Ox_fam_dicts_BPCH(wd=wd, CODE_wd=CODE_wd,
                                           weight_by_molecs=False,
                                           rm_strat=True,
                                           Mechanism=Mechanism)
    # Analyse odd oxygen (Ox) loss budget via route (chemical family)
    AC.calc_fam_loss_by_route(Ox_fam_dict=Ox_fam_dict, Mechanism=Mechanism)


if __name__ == "__main__":
    main()
