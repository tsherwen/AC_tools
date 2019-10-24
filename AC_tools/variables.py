#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Variable store/dictionarys for use in AC_tools

Use help(<name of function>) to get details on a particular function.

Notes
 - This code will be updated to use a user configuration approach (*.rc file) shortly
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
from . core import *


def get_unit_scaling(units, scaleby=1):
    """
    Get scaling for a given unit string

    Parameters
    -------
    units (str): units
    scaleby (float): scaling factor for unit

    Returns
    -------
    (float)
    """
    logging.debug("Getting unit scaling for {units}".format(units=units))
    misc = (
        'K', 'm/s', 'unitless', 'kg', 'm', 'm2', 'kg/m2/s', 'molec/cm2/s',
        'mol/cm3/s',  'kg/s', 'hPa', 'atoms C/cm2/s', 'kg S', 'mb', 'atoms C/cm2/s',
        'molec/cm3', 'v/v', 'cm/s', 's-1', 'molec/m3', 'W/m2', 'unitless',
        '$\mu$g m$^{-3}$',
    )
    # parts per trillion
    if any([(units == i) for i in ('pptv', 'pptC', 'ppt')]):
        scaleby = 1E12
    # parts per billion
    elif any([(units == i) for i in ('ppbv', 'ppbC', 'ppb')]):
        scaleby = 1E9
    elif any([units == i for i in misc]):
        scaleby = 1.
        logging.debug('WARNING: {} is not in unit lists: '.format(units))
    return scaleby


class species:
    """
    Class for holding infomation about chemical specie

    Notes
    -----
     - the class is built from a csv file in the data folder of the repository
     - This file was built using the table of Species in GEOS-Chem on the wiki linked
     below. It was updated on
    http://wiki.seas.harvard.edu/geos-chem/index.php/Species_in_GEOS-Chem
    """
    def __repr__(self):
        rtn_str = "This is a class to hold chemical species information for {}"
        return rtn_str.format(self.name)

    def __str__(self):
        for key in self.__dict__.keys():
            print( "{:<20}:".format(key, self.key) )

    def __init__(self, name):
        self.name = name
        self.help = ("""This is a class to get information on a species from a CSV file
   It might contain the following information:
   self.RMM        = Molecular weight of the species in g mol-1
   self.latex      = The latex name of the species
   self.Smiles     = The smiles string of the species
   self.InChI      = The InChI string of the species
   self.Phase      = Denotes whether the species is in the gas or aerosol phase.
   self.Formula    = Chemical formula of the species
   self.long_name  = A longer, descriptive name for the species
   self.Chem       = Is the species contained in one of the pre-built KPP chemistry mechanisms used by GEOS-Chem?
   self.Advect     = Is the species subject to advection, cloud convection, + BL mixing?
   self.Drydep     = Is the species subject to dry deposition?
   self.Wetdep     = Is the species soluble and subject to wet deposition?
   self.Phot       = Is the species included in the FAST-JX photolysis mechanism?
   self.Mechanisms = the GEOS-Chem chemistry mechanisms to which the species belongs
   self.Ox         = Shows the number of molecules of the species that will be included in the computation of the P(Ox) and L(Ox) diagnostic families
   self.Version    = The version this species was added to GEOS-Chem or the latest version this species was updated (if available)
    self.Carbons   = number of carbon atoms in species
   """)
        # Set folder and filename to use, then check they are present
        # NOTE: "Python AC_tools/Scripts/get_data_files.py" retrieves data files
        folder = os.path.dirname(__file__) + '/../data/'
        filename = 'GEOS_ChemSpecies_fullchem_v0.1.0.csv'
        assert os.path.exists(folder+filename), "Error: Species csv not found!"
        # Open species csv file as dataframe
        dfM = pd.read_csv(folder+filename)
        # see if species in df
        df = dfM.loc[dfM['Species'] == str(self.name), : ]
        # Add properties from csv file
        if df.shape[0] != 0:
            # - Now attributes for species
            self.formula = str(df['Formula'].values[0])
            self.long_name = str(df['Full name'].values[0])
            self.RMM = str(df['Molec wt\n(g/mol)'].values[0])
            self.Phase = str(df['Gas or Aer'].values[0])
            self.Chem = str(df['Chem'].values[0])
            self.Advect = str(df['Advect'].values[0])
            self.Drydep = str(df['Drydep'].values[0])
            self.Wetdep = str(df['Wetdep'].values[0])
            self.Phot = str(df['Phot'].values[0])
            self.Mechanisms = str(df['Mechanisms'].values[0])
            self.Ox = str(df['Ox?'].values[0])
            self.Version = str(df['Version\nadded/\nupdated'].values[0])
            # inc smiles, InChI, and latex add offline from MChem_tools file.
            # NOTE: smiles/InChI/latex need updating
            self.InChI = str(df['InChI'].values[0])
            self.smiles = str(df['smiles'].values[0])
            self.LaTeX = str(df['LaTeX'].values[0])

            # - Update the formating of columns of DataFrame
            # TODO: update to read a pre-processed file for speed.
            # Add the number of carbons in a species
            def add_carbon_column(x):
                try:
                    return float(x.split('(12, ')[-1][:-2])
                except:
                    return np.NaN
#            df['Carbons'] = df['Molec wt\n(g/mol)'].map(add_carbon_column)
#            val = df['Molec wt\n(g/mol)'].values[0]
#            df['Carbons'] = np.NaN
#            df.loc[:,'Carbons'] = add_carbon_column(val)
            self.Carbons = add_carbon_column(self.RMM)
            # Make sure mass is shown as RMM
            def mk_RMM_a_float(x):
                try:
                    return float(x.split('(12, ')[0].strip())
                except:
                    return float(x)
#            val = df['Molec wt\n(g/mol)'].values[0]
#            df.loc[:,'Molec wt\n(g/mol)'] = mk_RMM_a_float(val)
            self.RMM = mk_RMM_a_float(self.RMM)
            # Convert booleans from crosses to True or False
            def update_X_to_bool(x):
                if x == 'X':
                    return True
                else:
                    return False
            booleans = 'Chem', 'Advect', 'Drydep', 'Wetdep', 'Phot'
#            for col in booleans:
#                val = df[col].values[0]
#                df.loc[:,col] = update_X_to_bool(val)
            self.Chem = update_X_to_bool(self.Chem)
            self.Advect = update_X_to_bool(self.Advect)
            self.Drydep = update_X_to_bool(self.Drydep)
            self.Wetdep = update_X_to_bool(self.Wetdep)
            self.Phot = update_X_to_bool(self.Phot)

        else:
            print("Species not found in CSV file ({})".format(filename))
        # Remove the DataFrame from memory
        del dfM

    def help(self):
        '''
        Another way to get help for the class.
        '''
        help(self)
        return


def constants(input_x, rtn_dict=False):
    """
    Dictionary storing commonly used constants

    Parameters
    ----------
    input_x (str): name of contant of interest
    rtn_dict (bool): retrun complete dictionary instead of single value

    Returns
    -------
    value of constant (float) or dictionary on constants and values (dict)

    Notes
    -----
    """
    con_dict = {
        # Relative atomic mass of air (Just considering N2+O2)
        'RMM_air': (.78*(2.*14.)+.22*(2.*16.)),
        # Avogadro constant (Mol^-1)
        'AVG': 6.0221413E23,
        # Dobson unit - (molecules per meter squared)
        'mol2DU': 2.69E20,
        # Specific gas constant for dry air (J/(kg·K))
        'Rdry': 287.058,
    }
    if rtn_dict:
        return con_dict
    else:
        return con_dict[input_x]


def get_loc(loc=None, rtn_dict=False, debug=False):
    """
    Dictionary to store locations (lon., lat., alt.)

    Returns
    -------
    (tuple)

    Notes
    -----
     - Data arranged: LON, LAT, ALT
    ( LON in deg E, LAT in deg N, ALT in metres a.s.l. )
     - double up? ( with 5.02 ?? )
     - Now use Class of GEO_site in preference to this func?
     - UPDATE NEEDED: move this to func_vars4obs
    """
    loc_dict = {
        # - CAST/CONTRAST
        'GUAM':  (144.800, 13.500, 0),
        'CHUUK': (151.7833, 7.4167, 0),
        'PILAU': (134.4667, 7.3500, 0),
        'London': (-0.1275, 51.5072, 0),
        'Weybourne': (1.1380, 52.9420,  0),
        'WEY': (1.1380, 52.9420,  0),  # Weyboure ID
        'Cape Verde': (-24.871, 16.848, 0),
        'CVO': (-24.871, 16.848,  0),  # Cape Verde ID
        'CVO1': (-24.871, 16.848,  0),  # Cape Verde ID
        'CVO (N)': (-24.871, 16.848+4,  0),  # Cape Verde (N)
        'CVO (N) 2x2.5': (-24.871, 16.848+2,  0),  # Cape Verde (N)
        'CVO2': (-24.871, 16.848+4,  0),  # Cape Verde (N)
        'CVO (NNW)': (-24.871, 16.848+8.,  0),  # Cape Verde (NNW)
        'CVO3': (-24.871, 16.848+8.,  0),  # Cape Verde (NNW)
        'CVO (NW)': (-24.871-4, 16.848+4.,  0),  # Cape Verde (NW)
        'CVO (NW) 2x2.5': (-24.871-2, 16.848+2.,  0),  # Cape Verde (NW)
        'CVO4': (-24.871-4, 16.848+4.,  0),  # Cape Verde (NW)
        'CVO (W)': (-24.871-4, 16.848,  0),  # Cape Verde (W)
        'CVO (W) 2x2.5': (-24.871-2, 16.848,  0),  # Cape Verde (W)
        'CVO5': (-24.871-4, 16.848,  0),  # Cape Verde (W)
        'CVO (S)': (-24.871, 16.848-4,  0),  # Cape Verde (S)
        'CVO6': (-24.871, 16.848-4,  0),  # Cape Verde (S)
        'CVO (SW)': (-24.871-4, 16.848-4,  0),  # Cape Verde (SW)
        'CVO7': (-24.871-4, 16.848-4,  0),  # Cape Verde (SW)
        # - ClearFlo
        'North Ken':  (-0.214174, 51.520718, 0),
        'KEN':  (-0.214174, 51.520718, 0),
        'BT tower': (-0.139055, 51.521556, 190),
        'BTT': (-0.139055, 51.521556, 190),
        # - ClNO2 sites
        'HOU': (-95.22, 29.45, 0),
        'BOL': (-105.27, 40.0, 1655 + 150),
        'LAC': (-118.23, 34.05, 	0),
        'HES': (8.45, 50.22, 	825),
        'SCH': (114.25, 22.22, 60),
        'TEX': (-95.425000, 30.350278,  60),
        'CAL': (-114.12950, 51.07933,  1100),
        'PAS':  (-118.20, 34.23, 246),
        # - ClNO2 (UK) sites
        'PEN':  (-4.1858, 50.3214, 0.),
        'LEI_AUG':  (-1.127311, 52.619823, 0.),
        'LEI_MAR':  (-1.127311, 52.619823, 0.),
        'LEI':  (-1.127311, 52.619823, 0.),
        'Leicester':  (-1.127311, 52.619823, 0.),
        'Mace_head_M3': (-10.846408, 53.209003, 0.),
        'Penlee':  (-4.1858, 50.3214, 0.),
        'Penlee_M2': (-2.0229414, 49.7795272, 0.),
        'Penlee_M3': (-5.3652425, 49.8370764, 0.),
        'Penlee_M4': (-4.15,  50.25, 0.),
        'Penlee_M5': (-0.85, 50.25, 0.),
        'Penlee_M6': (-7.05, 50.25, 0.),
        'Penlee_M7': (-4.1858, 50.1, 0.),
        # - Europe sites
        'DZK':  (4.5000, 52.299999237, 4),
        # - sites with preindustrial ozone observations
        'MON':  (2.338333, 48.822222,  75+5),
        #    'MON' : (2.3, 48.8, 80), # Monsoursis
        # Pavelin  et al. (1999) / Mickely  et al. (2001)
        # 0 m. a. l. assumed following ( "<500m") in Mickely  et al. (2001)
        'ADE': (138.0, -35.0, 0),  # Adelaide
        'COI': (-8.0, 40.0, 0),  # Coimbra
        'HIR': (132.0, 34.0, 0),  # Hiroshima
        'HOB': (147.0, -43.0, 0),  # Hobart
        'HOK': (114.0, 22.0, 0),  # Hong Kong
        'LUA': (14.0, -9.0, 0),  # Luanda
        'MAU': (57.0, -20.0, 0),  # Mauritius
        'MOV': (-56.0, -35.0, 0),  # Montevideo
        'MVT': (4.0, 44.0, 1900),  # Mont-Ventoux
        'NEM': (145.0, 43.0, 0),  # Nemuro
        'TOK': (139.0, 35.0, 0),  # Tokyo
        'VIE': (16.0, 48.0, 0),  # Vienna
        'PDM': (0.0, 43.0, 1000),  # Pic du midi
        # - Miscellaneous
        #    'MAC' : ( -10.846408, 53.209003, 0 ) # Mace Head.
        'MAC': (-9.9039169999999999, 53.326443999999995, 0),  # Mace Head.
        'Mace Head': (-9.9039169999999999, 53.326443999999995, 0),  # .
        'Brittany': (-4.0, 48.7, 0),  # Brittany, France
        'Ria de Arousa': (-8.87, 42.50, 0),  # Ria de Arousa, Spain
        'Mweenish Bay': (-9.83, 53.31, 0),  # Ireland
        'Harestua': (10.7098608, 60.2008617, 0),  # Norway
        'Cartagena': (-1.0060599, 37.6174104, 0),  # Spain
        'Malasapina - final day': (-8.338, 35.179, 0),  # final day of cruise
        'Dagebull': (8.69, 54.73, 0),
        'Lilia': (-4.55, 48.62, 0),
        'Heraklion': (25.1, 35.3, 0),  # Heraklion, Crete
        'Sylt': (8.1033406, 54.8988164, 0),
        'Sicily': (14.2371407,  38.5519809, 0),  # Sicily
        #    'Frankfurt' : ( 8.45,50.22, )
        # - Global GAW sites (from GAWSIS)
        'Barrow': (-156.6114654541, 71.3230133057,  11),
        'Ascension Island': (-14.3999996185, -7.9699997902, 91),
        'Neumayer': (-8.265999794, -70.6660003662, 42),
        'Hilo': (-155.0700073242, 19.5799999237,  11),
        'Samoa': (-170.5645141602, -14.2474746704, 77),
        'Assekrem': (5.6333332062, 23.2666664124,  2710),
        # - Misc
        'UoM_Chem': (-2.2302418, 53.4659844, 38),
        'CDD': (6.83333333, 45.8333333, 4250),
        'NEEM': (-51.12, 77.75, 2484),
        'Welgegund': (26.939311,  -26.570146, 1480),
        'WEL': (26.939311,  -26.570146, 1480),  # abrev. Welgegund
        'WEL-W': (26.939311-1.5,  -26.570146, 1480),  # abrev. Welgegund
        'WEL-SW': (26.939311-1.5,  -26.570146-1.5, 1480),  # abrev. Welgegund
        'Botsalano': (25.75, -25.54, 1420),
        'BOT': (25.75, -25.54, 1420),  # abrev. Botsalano
        'Marikana': (27.48, -25.70, 1170),
        'MAR': (27.48, -25.70, 1170),  # abrev. Marikana
        'Elandsfontein': (29.42, -26.25, 1750),
        'ELA': (29.42, -26.25, 1750),  # abrev. Elandsfontein
        # - Global GAW sites
        'ASK': (5.63, 23.27, 2710.0000000000005),
        'BRW': (-156.6, 71.32, 10.999999999999746),
        'CGO': (144.68, -40.68, 93.99999999999973),
        'CMN': (10.7, 44.18, 2165.0),
        'CPT': (18.48, -34.35, 229.99999999999997),
#        'CVO': (-24.871, 16.848, 10.000000000000103), # Already present.
        'JFJ': (7.987, 46.548, 3580.0),
        'LAU': (169.67, -45.03, 369.99999999999983),
        'MHD': (-9.9, 53.33, 4.999999999999905),
        'MLO': (-155.578, 19.539, 3397.0),
        'MNM': (153.981, 24.285, 7.999999999999767),
        'NMY': (-8.25, -70.65, 41.99999999999969),
        'SMO': (-170.565, -14.247, 77.00000000000001),
        'SPO': (-24.8, -89.98, 2810.0),
        'THD': (-124.15, 41.05, 119.99999999999997),
        # - NOAA sites
        # https://www.esrl.noaa.gov/gmd/grad/antuv/Palmer.jsp
        'Palmer Station' : (64.05, -64.767, 21.),
        'PSA' : (64.05, -64.767, 21.),
        # https://www.esrl.noaa.gov/gmd/dv/site/LEF.html
        'Park Falls Wisconsin' : (-90.2732, 45.9451, 472.00),
        'LEF' : (-90.2732, 45.9451, 472.00),
        # https://www.esrl.noaa.gov/gmd/obop/mlo/aboutus/siteInformation/kumukahi.html
        'Cape Kumukahi' : (-154.82, 19.54, 15.0),
        'KUM' : (-154.82, 19.54, 15.0),
        # https://www.esrl.noaa.gov/gmd/dv/site/NWR.html
        # NOTE: altitude in the GAW NetCDF is differrent (680.7339159020077m)
        'Niwot Ridge' : (-105.5864, 40.0531, 3523.00),
        'NWR' : (-105.5864, 40.0531, 3523.00),
        # https://gawsis.meteoswiss.ch/GAWSIS/#/search/station/stationReportDetails/487
        'Alert' : (-62.3415260315, 82.4991455078, 210.),
        'ALT' : (-62.3415260315, 82.4991455078, 210.),
        # https://gawsis.meteoswiss.ch/GAWSIS/#/search/station/stationReportDetails/312
        'Summit' : (-38.4799995422, 72.5800018311, 3238.),
        'SUM' : (-38.4799995422, 72.5800018311, 3238.),
        # https://www.esrl.noaa.gov/gmd/hats/stations/hfm.html
        # https://gawsis.meteoswiss.ch/GAWSIS/#/search/station/stationReportDetails/173
        # https://www.esrl.noaa.gov/gmd/dv/site/HFM.html
        'Havard Forest' :  ( -72.3000030518, 42.9000015259, 340.),
        'HFM' : ( -72.3000030518, 42.9000015259, 340.),
    }
    if rtn_dict:
        return loc_dict
    else:
        return loc_dict[loc]


def site_code2name(code):
    """
    Get the full site name from a given code (e.g. GAW ID)
    """
    d = {
    'CMN' : 'Monte Cimone',
    'CGO' : 'Cape Grim',
    'BRW' : 'Barrow',
    'HFM' : 'Havard Forest',
    'SUM' : 'Summit',
    'NWR' : 'Niwot Ridge',
    'KUM' : 'Cape Kumukahi',
    'LAU' : 'Lauder',
    'ASK' : 'Assekrem',
    'JFJ' : 'Jungfraujoch',
    'MHD' : 'Mace Head',
    'MLO' : 'Mauna Loa',
    'MNM' : 'Minamitorishima',
    'NMY' : 'Neumayer',
    'SMO' : 'Samoa',
    'SPO' : 'South Pole',
    'THD' : 'Trinidad Head',
    'ALT' : 'Alert',
    'CPT' : 'Cape Point',
    'LEF' : 'Park Falls Wisconsin',
    'PSA' : 'Palmer Station',
    }
    return d[code]


def sort_locs_by_lat(sites):
    """
    Order given list of sties by their latitudes
    """
    # Get info
    vars = [get_loc(s) for s in sites]  #  lon, lat, alt,
    # Sort by lat, index orginal sites list and return
    lats = [i[1] for i in vars]
    slats = sorted(lats)[::-1]
    return [sites[i] for i in [lats.index(ii) for ii in slats]]



