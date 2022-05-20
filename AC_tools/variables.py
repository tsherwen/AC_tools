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
import inspect
import yaml
import pandas as pd
from netCDF4 import Dataset
import sys
import glob
# Math/Analysis
import numpy as np
# the below import needs to be updated,
# imports should be specific and in individual functions
# import tms modules with shared functions
from . core import *
#from . utils import read_yaml_file


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
        'mol/cm3/s',  'kg/s', 'hPa', 'atoms C/cm2/s', 'kg S', 'mb',
        'atoms C/cm2/s',
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
    Class for holding infomation about chemical species

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
            print("{:<20}:".format(key, self.key))

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
        df = dfM.loc[dfM['Species'] == str(self.name), :]
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
        'BMW':  (64.75, 32.32, 53),  #  Tudor Hill Atmospheric observatory
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
        'Palmer Station': (64.05, -64.767, 21.),
        'PSA': (64.05, -64.767, 21.),
        # https://www.esrl.noaa.gov/gmd/dv/site/LEF.html
        'Park Falls Wisconsin': (-90.2732, 45.9451, 472.00),
        'LEF': (-90.2732, 45.9451, 472.00),
        # https://www.esrl.noaa.gov/gmd/obop/mlo/aboutus/siteInformation/kumukahi.html
        'Cape Kumukahi': (-154.82, 19.54, 15.0),
        'KUM': (-154.82, 19.54, 15.0),
        # https://www.esrl.noaa.gov/gmd/dv/site/NWR.html
        # NOTE: altitude in the GAW NetCDF is differrent (680.7339159020077m)
        'Niwot Ridge': (-105.5864, 40.0531, 3523.00),
        'NWR': (-105.5864, 40.0531, 3523.00),
        # https://gawsis.meteoswiss.ch/GAWSIS/#/search/station/stationReportDetails/487
        'Alert': (-62.3415260315, 82.4991455078, 210.),
        'ALT': (-62.3415260315, 82.4991455078, 210.),
        # https://gawsis.meteoswiss.ch/GAWSIS/#/search/station/stationReportDetails/312
        'Summit': (-38.4799995422, 72.5800018311, 3238.),
        'SUM': (-38.4799995422, 72.5800018311, 3238.),
        # https://www.esrl.noaa.gov/gmd/hats/stations/hfm.html
        # https://gawsis.meteoswiss.ch/GAWSIS/#/search/station/stationReportDetails/173
        # https://www.esrl.noaa.gov/gmd/dv/site/HFM.html
        'Havard Forest':  (-72.3000030518, 42.9000015259, 340.),
        'HFM': (-72.3000030518, 42.9000015259, 340.),
        # - ARNA locations
        'Dakar': (-17.467686, 14.716677, 22),
        'DSS': (-17.467686, 14.716677, 22),  # Dakar airport code (as above)
        'Sao Vicente Airport': (-25.0569, 16.8331, 20),
        'VXE': (-25.0569, 16.8331, 20),  # Sao Vincite code (as above)
        'Praia Airport': (-23.4939, 14.9242, 70),
        'RAI': (-23.4939, 14.9242, 70),  # Praia airport code (as above)
        # Other "nearby" airports
        'Gran Canaria Airport': (-15.386667, 27.931944, 24),
        # Gran Canaria airport code (as above)
        'LPA': (-15.386667, 27.931944, 24),
        'Lisbon Airport': (-9.134167, 38.774167, 114),
        'LIS': (-9.134167, 38.774167, 114),  # Lisbon airport code (as above)
        'Paris (Charles de Gaulle) Airport': (-2.547778, 49.009722, 119),
        'CDG': (-2.547778, 49.009722, 119),  # Paris airport code (as above)
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
        'CMN': 'Monte Cimone',
        'CGO': 'Cape Grim',
        'BRW': 'Barrow',
        'HFM': 'Havard Forest',
        'SUM': 'Summit',
        'NWR': 'Niwot Ridge',
        'KUM': 'Cape Kumukahi',
        'LAU': 'Lauder',
        'ASK': 'Assekrem',
        'JFJ': 'Jungfraujoch',
        'MHD': 'Mace Head',
        'MLO': 'Mauna Loa',
        'MNM': 'Minamitorishima',
        'NMY': 'Neumayer',
        'SMO': 'Samoa',
        'SPO': 'South Pole',
        'THD': 'Trinidad Head',
        'ALT': 'Alert',
        'CPT': 'Cape Point',
        'LEF': 'Park Falls Wisconsin',
        'PSA': 'Palmer Station',
    }
    return d[code]


def sort_locs_by_lat(sites):
    """
    Order given list of sties by their latitudes
    """
    # Get info
    vars = [get_loc(s) for s in sites]  # lon, lat, alt,
    # Sort by lat, index orginal sites list and return
    lats = [i[1] for i in vars]
    slats = sorted(lats)[::-1]
    return [sites[i] for i in [lats.index(ii) for ii in slats]]


def latex_spec_name(input_x, debug=False):
    """
    Formatted (Latex) strings for species and analysis names


    Returns
    -------
    (tuple)

    Notes
    -----
     - Redundant? Can use class structure ( see species instance )
    """
    spec_dict = {
        'OIO': 'OIO', 'C3H7I': 'C$_{3}$H$_{7}$I', 'IO': 'IO', 'I': 'I',
        'I2': 'I$_{2}$', 'CH2ICl': 'CH$_{2}$ICl', 'HOI': 'HOI',
        'CH2IBr': 'CH$_{2}$IBr', 'C2H5I': 'C$_{2}$H$_{5}$I',
        'CH2I2': 'CH$_{2}$I$_{2}$', 'CH3IT': 'CH$_{3}$I',
        'IONO': 'INO$_{2}$', 'HIO3': 'HIO$_{3}$', 'ICl': 'ICl',
        'I2O3': 'I$_{2}$O$_{3}$', 'I2O4': 'I$_{2}$O$_{4}$',
        'I2O5': 'I$_{2}$O$_{5}$', 'INO': 'INO', 'I2O': 'I$_{2}$O',
        'IBr': 'IBr', 'I2O2': 'I$_{2}$O$_{2}$', 'IONO2': 'INO$_{3}$',
        'HI': 'HI',
        'BrO': 'BrO', 'Br': 'Br', 'HOBr': 'HOBr', 'Br2': 'Br$_{2}$',
        'CH3Br': 'CH$_{3}$Br', 'CH2Br2': 'CH$_{2}$Br$_{2}$',
        'CHBr3': 'CHBr$_{3}$', 'O3': 'O$_{3}$', 'CO': 'CO', 'DMS': 'DMS',
        'NO': 'NO', 'NO2': 'NO$_{2}$',
        'NO3': 'NO$_{3}$', 'HNO3': 'HNO$_{3}$', 'HNO4': 'HNO$_{4}$',
        'PAN': 'PAN', 'HNO2': 'HONO', 'N2O5': 'N$_{2}$O$_{5}$',
        'ALK4': '$\geq$C4 alkanes', 'ISOP': 'Isoprene',
        'H2O2': 'H$_{2}$O$_{2}$',
        'ACET': 'CH$_{3}$C(O)CH$_{3}$', 'MEK': '>C3 ketones',
        'RCHO': 'CH$_{3}$CH$_{2}$CHO',
        'MVK': 'CH$_{2}$=CHC(O)CH$_{3}$', 'MACR': 'Methacrolein',
        'PMN': 'PMN', 'PPN': 'PPN',
        'R4N2': '$\geq$C4 alkylnitrates', 'PRPE': '$\geq$C3 alkenes',
        'C3H8': 'C$_{3}$H$_{8}$', 'CH2O': 'CH$_{2}$O',
        'C2H6': 'C$_{2}$H$_{6}$', 'MP': 'CH$_{3}$OOH', 'SO2': 'SO$_{2}$',
        'SO4': 'SO${_4}{^{2-}}$', 'SO4s': 'SO${_4}{^{2-}}$ on SSA',
        'MSA': 'CH$_{4}$SO$_{3}$', 'NH3': 'NH$_{3}$', 'NH4': 'NH${_4}{^+}$',
        'NIT': 'NO${_3}{^-}$', 'NITs': 'NO${_3}{^-}$ on SSA', 'BCPI': 'BCPI',
        'OCPI': 'OCPI', 'BCPO': 'BCPO', 'OCPO': 'OCPO', 'DST1': 'DST1',
        'DST2': 'DST2', 'DST3': 'DST3', 'DST4': 'DST4', 'SALA': 'SALA',
        'SALC': 'SALC',  'HBr': 'HBr', 'BrNO2': 'BrNO$_{2}$',
        'BrNO3': 'BrNO$_{3}$', 'MPN': 'CH$_{3}$O$_{2}$NO$_{2}$',
        'ISOPN': 'ISOPN', 'MOBA': 'MOBA', 'PROPNN': 'PROPNN',
        'HAC': 'HAC', 'GLYC': 'GLYC', 'MMN': 'MMN', 'RIP': 'RIP',
        'IEPOX': 'IEPOX', 'MAP': 'MAP', 'AERI': 'Aerosol Iodine',
        'Cl2': 'Cl$_{2}$',
        'Cl': 'Cl', 'HOCl': 'HOCl', 'ClO': 'ClO', 'OClO': 'OClO',
        'BrCl': 'BrCl',
        'HI+OIO+IONO+INO': 'HI+OIO+INO$_{2}$+INO',
        'CH2IX': 'CH$_{2}$IX (X=Cl, Br, I)',
        'IxOy': 'I$_{2}$O$_{X}$ ($_{X}$=2,3,4)',
        'CH3I': 'CH$_{3}$I', 'OH': 'OH', 'HO2': 'HO$_{2}$', 'MO2': 'MO$_{2}$',
        'RO2': 'RO$_{2}$', 'ISALA': 'Iodine on SALA',
        'ISALC': 'Iodine on SALC', 'CH4': 'CH$_{4}$', 'MOH': 'Methanol',
        'RD01': r'I + O$_{3}$ $\rightarrow$ IO + O$_{2}$',
        # Adjusted names
        'ALD2': 'Acetaldehyde',
        # Analysis names
        'iodine_all': 'All Iodine', 'Iy': 'I$_{\\rm y}$',\
        'IOy': 'IO$_{\\rm y}$', \
        'IyOx': 'I$_{y}$O$_{x}$',
        'IOx': 'IO$_{\\rm x}$', \
        'iodine_all_A': 'All Iodine (Inc. AERI)',  \
        'I2Ox': 'I$_{2}$O$_{\\rm X}$', 'AERI/SO4': 'AERI/SO4', \
        'EOH': 'Ethanol', 'OH reactivity / s-1': 'OH reactivity / s$^{-1}$', \
        'PSURF': 'Pressure at the bottom of level', \
        'GMAO_TEMP': 'Temperature', 'TSKIN': 'Temperature at 2m', \
        'GMAO_UWND': 'Zonal Wind', 'GMAO_VWND': 'Meridional Wind', \
        'U10M': '10m Meridional Wind', 'V10M': '10m Zonal Wind', \
        'CH2OO': 'CH$_{2}$OO', 'Sulfate': 'Sulfate', 'VOCs': 'VOCs', \
        'GMAO_ABSH': 'Absolute humidity', 'GMAO_SURF': 'Aerosol surface area', \
        'GMAO_PSFC': 'Surface pressure',
        # Family/group species/tracer Names
        'N_specs': 'NO$_{\\rm y}$', 'NOy': 'NO$_{\\rm y}$',
        'Bry': 'Br$_{\\rm y}$', 'Cly': 'Cl$_{\\rm y}$', 'NIT+NITs': 'NIT+NITs',  \
        'N_specs_no_I': 'NO$_{\\rm y}$ exc. iodine', 'TSO4': 'TSO$_4$',
        'NOx': 'NO$_{\\rm x}$', 'HOx': 'HO$_{\\rm x}$', 'TNO3': 'TNO$_3$', \
        'SOx': 'SO$_{\\rm x}$', 'PM2.5': 'PM$_{2.5}$', 'VOC': 'VOC', \
        'NIT+NH4+SO4': 'NO${_3}{^-}$+NH${_4}{^+}$+SO${_4}{^{2-}}$',
        'NIT_ALL': 'NIT (all)', 'PM10': 'PM$_{10}$', 'BC': 'BC',
        # typos
        'CH2BR2': 'CH$_{2}$Br$_{2}$',\
        # Cly names
        'ClOO': 'ClOO', 'Cl2': 'Cl$_{2}$', \
        'BrCl': 'BrCl', 'ICl': 'ICl', 'HOCl': 'HOCl', 'ClO': 'ClO',
        'ClOO': 'ClOO', \
        'OClO': 'OClO', 'Cl2O2': 'Cl$_{2}$O$_{2}$', 'HCl': 'HCl', \
        'ClNO2': 'ClNO$_{2}$', 'ClNO3': 'ClNO$_{3}$', 'Cl': 'Cl',\
        'CH3Cl': 'CH$_{3}$Cl',  'CH2Cl2': 'CH$_{2}$Cl$_{2}$', \
        'CHCl3': 'CHCl$_{3}$',
        # Bry names
        'BrSALC': 'Br- on SALC', 'BrSALA': 'Br- on SALA',
        # add planeflight variabels  for ease of processing
        'LON': 'Lon.', 'LAT': 'Lat.', 'PRESS': 'Press.',
        # add combined species for easy of processing
        'HOCl+Cl2': 'HOCl+Cl$_2$', 'HOBr+Br2': 'HOBr+Br$_2$',
        'HNO3/NOx': 'HNO$_3$/NO$_{\\rm x}$',
        'HNO3+NIT': 'HNO$_3$+NIT', 'HNO3+NO3': 'HNO$_3$+NO$_3$',
        'NIT/NOx': 'NIT/NO$_{\\rm x}$', 'HNO3/NIT': 'HNO$_3$/NIT',
        # Add pseudo species
        'Cl-': 'Cl$^{-}$',
        # Extra sub species from GEOS-CF
        'PM2.5(dust)': 'PM$_{2.5}$(dust)',
        'PM2.5(SO4)': 'PM$_{2.5}$(SO${_4}{^{2-}}$)',
        'PM2.5(NIT)': 'PM$_{2.5}$(NO${_3}{^-}$)',
        'PM2.5(SOA)': 'PM$_{2.5}$(SOA)',
        'PM2.5(SSA)': 'PM$_{2.5}$(SSA)',
        'PM2.5(BC)': 'PM$_{2.5}$(BC)',
        'PM2.5(OC))': 'PM$_{2.5}$(OC)',
        # Extra GEOS-chem advected tracers - in standard as of v12.9.1
        # TODO - complete expanding species
        # http://wiki.seas.harvard.edu/geos-chem/index.php/Species_in_GEOS-Chem
        'ACTA': 'ACTA', 'ATOOH': 'ATOOH', 'BENZ': 'Benzene',
        'CCl4': 'CCl4', 'CFC11': 'CFC11', 'CFC113': 'CFC113',
        'CFC114': 'CFC114', 'CFC115': 'CFC115', 'CFC12': 'CFC12',
        'CH3CCl3': 'CH3CCl3', 'ETHLN': 'ETHLN', 'ETNO3': 'ETNO3',
        'ETP': 'ETP', 'GLYX': 'GLYX', 'H1211': 'H1211',
        'H1301': 'H1301', 'H2402': 'H2402', 'H2O': 'H2O',
        'HC5A': 'HC5A', 'HCFC123': 'HCFC123', 'HCFC141b': 'HCFC141b',
        'HCFC142b': 'HCFC142b', 'HCFC22': 'HCFC22', 'HCOOH': 'HCOOH',
        'HMHP': 'HMHP', 'HMML': 'HMML', 'HONIT': 'HONIT',
        'HPALD1': 'HPALD1', 'HPALD2': 'HPALD2', 'HPALD3': 'HPALD3',
        'HPALD4': 'HPALD4', 'HPETHNL': 'HPETHNL', 'ICHE': 'ICHE',
        'ICN': 'ICN', 'ICPDH': 'ICPDH', 'IDC': 'IDC',
        'IDCHP': 'IDCHP', 'IDHDP': 'IDHDP', 'IDHPE': 'IDHPE',
        'IDN': 'IDN', 'IEPOXA': 'IEPOXA', 'IEPOXB': 'IEPOXB',
        'IEPOXD': 'IEPOXD', 'IHN1': 'IHN1', 'IHN2': 'IHN2',
        'IHN3': 'IHN3', 'IHN4': 'IHN4', 'INDIOL': 'INDIOL',
        'INPB': 'INPB', 'INPD': 'INPD', 'IONITA': 'IONITA',
        'IPRNO3': 'IPRNO3', 'ITCN': 'ITCN', 'ITHN': 'ITHN',
        'LIMO': 'LIMO', 'LVOC': 'LVOC', 'LVOCOA': 'LVOCOA',
        'MACR1OOH': 'MACR1OOH', 'MCRDH': 'MCRDH', 'MCRENOL': 'MCRENOL',
        'MCRHN': 'MCRHN', 'MCRHNB': 'MCRHNB', 'MCRHP': 'MCRHP',
        'MENO3': 'MENO3', 'MGLY': 'MGLY', 'MONITA': 'MONITA',
        'MONITS': 'MONITS', 'MONITU': 'MONITU', 'MPAN': 'MPAN',
        'MTPA': 'MTPA', 'MTPO': 'MTPO', 'MVKDH': 'MVKDH',
        'MVKHC': 'MVKHC', 'MVKHCB': 'MVKHCB', 'MVKHP': 'MVKHP',
        'MVKN': 'MVKN', 'MVKPC': 'MVKPC', 'N2O': 'N2O',
        'NPRNO3': 'NPRNO3', 'OCS': 'OCS', 'pFe': 'pFe',
        'PIP': 'PIP', 'PP': 'PP', 'PRPN': 'PRPN',
        'PYAC': 'PYAC', 'R4P': 'R4P', 'RA3P': 'RA3P',
        'RB3P': 'RB3P', 'RIPA': 'RIPA', 'RIPB': 'RIPB',
        'RIPC': 'RIPC', 'RIPD': 'RIPD', 'RP': 'RP',
        'SALAAL': 'SALAAL', 'SALACL': 'SALACL', 'SALCAL': 'SALCAL',
        'SALCCL': 'SALCCL', 'SOAGX': 'SOAGX', 'SOAIE': 'SOAIE',
        'SOAP': 'SOAP', 'SOAS': 'SOAS', 'TOLU': 'TOLU',
        'XYLE': 'Xylene',
        # Extra GEOS-chem species - in standard as of v12.9.1
        'O1D': 'O($^{1}$D)', 'O': 'O', 'hv': '$h\\nu$',
    }
    return spec_dict[input_x]


def species_mass(spec):
    """
    Function to get species relative molecular mass (RMM) in g/mol

    Parameters
    ----------
    spec (str): species/tracer/variable name

    Returns
    -------
    (float)

    Notes
    -----
     - Use legacy AC_tools yml as default offline yml file for now.
     - aim switch to GC JSON file outputed by GEOS-Chem
     - C3H5I == C2H5I (this is a vestigle typo, left in to allow for
    use of older model run data  )
    """
    try:
        d = read_yaml_file('species_mass.yml')
        return d[spec]
    except KeyError:
        d = read_GC_species_database()
        SpecMassStr = 'MW_g'
        return d[spec][SpecMassStr]


def spec_stoich(spec, IO=False, I=False, NO=False, OH=False, N=False,
                C=False, Br=False, Cl=False, S=False, ref_spec=None,
                debug=False):
    """
    Returns unit equivalent of X ( e.g. I ) for a given species

    Parameters
    ----------
    spec (str): species/tracer/variable name
    ref_spec (str): species which number of spec equiv. in is being sought
    wd (str): Specify the wd to get the results from a run.
    res (str): the resolution if wd not given (e.g. '4x5' )
    debug (bool): legacy debug option, replaced by python logging
    IO, I, NO, OH, N, C, Br, Cl, S (bool): reference species to use (default=I)

    Returns
    -------
    (float) number equivlent units of ref_spec in spec.

    Notes
    -----
     - This can be automatically set by providing a reference species or by
     setting boolean input parametiers.
     - Update Needed: re-write to take stioch species
     (e.g. OH, I instead of booleans )
     - asssume I == True as default
     - C3H5I == C2H5I
     (this is a vestigle typo, left in to allow for use of older model runs)
     - aerosol cycling specs
    # 'LO3_36' : (2.0/3.0) , 'LO3_37' : (2.0/4.0),  # aersol loss rxns... 'LO3_37' isn't true loss, as I2O4 is regen. temp
     - Aerosol loss rxns ( corrected stoich. for Ox, adjsutment need for I )
    """
    # If reference species provided automatically select family
    if not isinstance(ref_spec, type(None)):
        if any([(ref_spec == i) for i in ('I', 'Iy', 'Iodine')]):
            I = True
        if any([(ref_spec == i) for i in ('Br', 'Bry', 'Bromine')]):
            Br = True
        if any([(ref_spec == i) for i in ('Cl', 'Cly', 'Chlorine')]):
            Cl = True
        if any([(ref_spec == i) for i in ('C', 'VOC')]):
            C = True
        if any([(ref_spec == i) for i in ('N', 'NOy', 'NOx')]):
            N = True
        if any([(ref_spec == i) for i in ('OH', 'HO2')]):
            OH = True
        if ref_spec == 'IO':
            IO = True
        if ref_spec == 'NO':
            NO = True
        if any([(ref_spec == i) for i in ('S', 'SOx', 'Sulfate')]):
            S = True

    if debug:
        vars = ref_spec, IO, I, NO, OH, N, C, Br, Cl
        varsn = 'ref_spec', 'IO', 'I', 'N', 'OH', 'N', 'C', 'Br', 'Cl'
        print(("'spec_stoich'  called for: ", list(zip(varsn, vars))))

    # Select dictionary ( I=True is the default... )
    if IO:
        d = {
            'RD11': 2.0, 'RD10': 1.0, 'RD12': 2.0, 'LO3_36': 1./3.,
            'RD09': 1.0,
            'RD66': 1.0, 'RD23': 1.0, 'RD37': 1.0, 'LO3_24': 1.0/2.0,
            'RD56': 1.0,
            'RD01': 1.0, 'RD08': 1.0, 'RD46': 2.0, 'RD30': 1.0, 'RD25': 1.0,
            'RD27': 1.0, 'RD97': 1.0
        }
    elif NO:
        d = {
            'NO2': 1.0, 'NO3': 1.0, 'N2O5': 2.0, 'NO': 1.0, 'PPN': 1.0,
            'R4N2': 1.0,
            'BrNO3': 1.0, 'INO': 1.0, 'PAN': 1.0, 'PMN': 1.0, 'HNO3': 1.0,
            'HNO2': 1.0, 'NH3': 1.0, 'HNO4': 1.0, 'BrNO2': 1.0,
            'IONO': 1.0, 'PROPNN': 1.0, 'NH4': 1.0, 'MPN': 1.0, 'MMN': 1.0,
            'ISOPN': 1.0, 'IONO2': 1.0
        }
    elif OH:
        d = {
            'LO3_18': 2.0, 'LO3_03': 1.0,  'PO3_14': 1.0, 'RD65': 1.0,
            'LR25': 1.0,
            'LOH': 1.0, 'POH': 1.0, 'LO3_86': 1.0, 'RD98': 1.0, \
            # Redundant: 'RD95': 1.0,
            # also include HO2 and OH for HOx calculations
            'OH': 1.0, 'HO2': 1.0
        }
    elif S:
        d = {
            'S': 1.0, 'SO4': 1.0, 'SO4s': 1.0, 'SO4S': 1.0, 'SO2': 1.0,
            'DMS': 1.0,
            'SO4D1': 1.0, 'SO4D2': 1.0, 'SO4D3': 1.0, 'SO4D4': 1.0,
        }
    elif N:
        d = {
            'RD10': 1.0, 'LR26': 1.0, 'LR27': 1.0, 'LR20': 1.0, 'RD17': 1.0,
            'RD16': 1.0, 'RD19': 1.0, 'RD18': 2.0, 'LR28': 1.0, 'LO3_30': 1.0,
            'RD75': 1.0, 'LR7': 1.0, 'LR8': 1.0, 'RD56': 1.0, 'RD24': 1.0,
            'LO3_39': 1.0, 'RD25': 1.0, 'RD81': 1.0, 'LR35': 1.0, 'LR18': 1.0,
            'LR17': 1.0, 'LR11': 1.0, 'LR39': 1.0, 'RD20': 1.0, 'RD21': 2.0,
            'RD22': 1.0, 'RD23': 1.0, 'RD68': 1.0, 'RD69': 1.0, \
            # NOy ( N in 'NOy')
            'NO2': 1.0, 'NO3': 1.0, 'N2O5': 2.0, 'NO': 1.0, 'PPN': 1.0, \
            'R4N2': 2.0, 'BrNO3': 1.0, 'INO': 1.0, 'PAN': 1.0, 'PMN': 1.0, \
            'HNO3': 1.0, 'HNO2': 1.0, 'NH3': 1.0, 'HNO4': 1.0, 'BrNO2': 1.0, \
            'IONO': 1.0, 'PROPNN': 1.0, 'NH4': 1.0, 'MPN': 1.0, 'MMN': 1.0, \
            'ISOPN': 1.0, 'IONO2': 1.0, 'ClNO2': 1.0, 'ClNO3': 1.0,
            'NIT': 1.0, 'NITs': 1.0, 'NITS': 1.0, \
            #
            'NITD1': 1.0, 'NITD2': 1.0, 'NITD3': 1.0, 'NITD4': 1.0,
        }
    elif C:
        d = {
            'ACET': 3.0, 'ALD2': 2.0, 'C2H6': 2.0, 'C3H8': 3.0, 'ISOP': 5.0,
            'PRPE': 3.0, 'ALK4': 4.0, 'MEK': 4.0,
            'APINE': 10.0, 'BPINE': 10.0, 'LIMON': 10.0, 'SABIN': 10.0,
            'MYRCN': 10.0,
            'CAREN': 10.0, 'OCIMN': 10.0, 'XYLE': 8.0,
        }
    elif Br:
        d = {
            'CH3Br': 1.0, 'HOBr': 1.0, 'BrO': 1.0, 'CHBr3': 3.0, 'Br2': 2.0,
            'BrSALC': 1.0, 'CH2IBr': 1.0, 'BrCl': 1.0, 'Br': 1.0,
            'CH2Br2': 2.0,
            'IBr': 1.0, 'BrSALA': 1.0, 'BrNO2': 1.0, 'BrNO3': 1.0, 'HBr': 1.0,
            # for ease of processing also include Seasalt Br2
            'SSBr2': 2.0,
            # Also have reaction tracers
            'LR73': 1.0,
            # Note: stoichometry is for **GAS** phase Br (aka not SSA )
            # ( Aka JT03s == Br2 ( ==2 ), but one is BrSALA/BrSALC therefore =1)
            'JT03s': 1.0, 'JT04s': 1.0, 'JT05s': 1.0,
            # BrCl from HOBr or hv
            'JT02s': 1.0, 'JT08': 1.0,
            # v11 KPP Tags
            'T149': 3.0, 'T127': 0.680+1.360, 'T126': 0.440+0.560, 'T071': 3.0,
            'T198': 0.150, 'T199': 0.150, 'T200': 0.150, 'T082': 1.0
        }
    elif Cl:
        d = {
            'ClO': 1.0, 'Cl': 1.0, 'ClOO': 1.0, 'ClNO3': 1.0, 'ClNO2': 1.0,
            'Cl2': 2.0, 'OClO': 1.0, 'HOCl': 1.0, 'HCl': 1.0, 'Cl2O2': 2.0,
            'BrCl': 1.0, 'ICl': 1.0, 'CH2Cl2': 2.0, 'CHCl3': 3.0,
            'CH2ICl': 1.0,
            'CH3Cl': 1.0,
            # Also have reaction tracers
            'LR62': 3.0, 'LR107': 3.0,
            'LR74': 1.0, 'LR106': 1.0, 'LR103': 1.0,
            'LR75': 2.0, 'LR105': 2.0, 'LR104': 2.0,
            # BrCl from HOBr or hv
            'JT02s': 1.0, 'JT08': 1.0,
            # ICl  (assuming 0.85:0.15 )
            'RD59': 0.15, 'RD92': 0.15, 'RD63': 0.15,
            # N2O5+SSA=>ClNO2
            'LR114': 1.0,
            # v11 KPP Tags
            'T174': 3.0, 'T203': 3.0,
            'T173': 2.0, 'T201': 2.0, 'T202': 2.0,
            'T172': 1.0, 'T171': 1.0, 'T143': 1.0,
            'T155': 1.0, 'T135': 1.0, 'T212': 1.0,
            'T198': 0.850, 'T199': 0.850, 'T200': 0.850,
            'PT213': 1.0, 'PT214': 1.0,
        }
    else:  # ( I=True is the default... )
        d = {
            'RD11': 1.0, 'RD10': 1.0, 'HIO3': 1.0, 'RD15': 1.0, 'RD62': 2.0,
            'RD17': 1.0, 'RD16': 1.0, 'RD19': 1.0, 'LO3_37': 0.5, 'CH2I2': 2.0,
            'AERII': 1.0, 'CH2ICl': 1.0, 'PIOx': 1.0, 'C3H7I': 1.0,
            'RD73': 1.0,
            'RD72': 2.0, 'RD71': 1.0, 'RD70': 1.0, 'C3H5I': 1.0, 'RD57': 1.0,
            'CH3IT': 1.0, 'IO': 1.0, 'LO3_38': 1.0, 'RD61': 1.0, 'RD68': 1.0,
            'I2': 2.0, 'IONO': 1.0, 'LO3_36': 0.6666666666666666, 'INO': 1.0,
            'RD88': 1.0, 'RD89': 1.0, 'LOx': 1.0, 'RD06': 1.0, 'RD07': 1.0,
            'RD02': 1.0, 'RD01': 1.0, 'I': 1.0,  'LO3_24': 0.5, 'AERI': 1.0,
            'HOI': 1.0, 'RD64': 2.0, 'RD65': 1.0, 'RD66': 1.0, 'RD67': 1.0,
            'RD60': 1.0, 'RD47': 1.0, 'C2H5I': 1.0, 'RD63': 1.0, 'RD20': 1.0,
            'RD22': 1.0, 'RD24': 1.0, 'RD69': 1.0, 'RD27': 1.0, 'OIO': 1.0,
            'CH2IBr': 1.0, 'LIOx': 1.0, 'L_Iy': 1.0, 'ICl': 1.0, 'IBr': 1.0,
            'RD95': 2.0, 'I2O2': 2.0, 'I2O3': 2.0, 'I2O4': 2.0, 'I2O5': 2.0,
            'HI': 1.0, 'I2O': 2.0, 'RD59': 1.0, 'RD93': 2.0, 'RD92': 1.0,
            'IONO2': 1.0, 'RD58': 1.0, 'ISALA': 1.0, 'ISALC': 1.0, 'CH3I': 1.0, \
            # p/l for: IO, I
            'RD15': 1.0, 'RD17': 1.0, 'RD75': 1.0, 'RD72': 2.0, 'RD71': 1.0, \
            'RD70': 1.0, 'RD56': 1.0, 'RD69': 1.0, 'RD88': 1.0, 'RD89': 1.0, \
            'RD06': 1.0, 'RD07': 1.0, 'RD08': 1.0, 'RD64': 2.0, 'RD65': 1.0, \
            'RD67': 1.0, 'RD46': 2.0, 'RD47': 1.0, 'RD20': 1.0, 'RD22': 1.0, \
            'RD68': 1.0, 'RD25': 1.0, 'RD96': 1.0, 'RD11': 1.0, 'RD12': 2.0, \
            'RD02': 1.0, 'RD16': 1.0, 'RD19': 1.0, 'RD24': 1.0, 'RD09': 1.0, \
            'RD23': 1.0, 'RD37': 1.0, 'RD97': 1.0, \
            # kludge for test analysis (HEMCO emissions )
            'ACET': 1.0, 'ISOP': 1.0, 'CH2Br2': 1.0, 'CHBr3': 1.0,
            'CH3Br': 1.0, \
            # Iodine in het loss/cycling reactions
            # loss to SSA/other aerosols
            # HOI
            'LR44': 1.0, 'LR45': 1.0, 'LR32': 1.0,   \
            # HI other
            'LR34': 1.0, \
            # IONO2
            'LR42': 1.0, 'LR43': 1.0, 'LR35': 1.0, \
            # IONO
            'LR46': 1.0, 'LR47': 1.0, 'LR39': 1.0,
            # --- KPP tags
            # Iy cycling sinks...
            'T217': 1.0, 'T216': 1.0, 'T198': 1.0, 'T199': 1.0, 'T196': 1.0,
            'T183': 1.0, 'T195': 1.0, 'T184': 1.0,  'T215': 1.0, 'T197': 1.0,
            # I2Oy
            'T190': 2.0, 'T193': 2.0, 'T187': 2.0,
            'T186': 2.0, 'T189': 2.0, 'T192': 2.0,
            'T185': 2.0, 'T188': 2.0, 'T191': 2.0,
        }

    # Kludge for testing. Allow values to equal 1.0 if not defined.
    try:
        if debug:
            prt_str = '{} (ref_spec: {}) stoichiometry : {}'
            print((prt_str.format(spec, ref_spec,  d[spec])))
        return d[spec]

    except:
        prt_str = '!!!!!!! WARNING - Kludge assumming stoichiometry = 1.0, for'
        prt_str += ' {} (ref_spec given as: {})'
        print((prt_str.format(spec, ref_spec)))
        return 1.0


def tra_unit(x, scale=False, adjustment=False, adjust=True, global_unit=False,
             ClearFlo_unit=False, IUPAC_unit=False, use_pf_species_units=False,
             debug=False):
    """
    Get appropirate unit for Tracer (and scaling if requested)

    Parameters
    ----------
    x (str): species/tracer/variable name
    adjustment (float): adjustment to give unit (+,- value etc )
    scale (float): scaling factor for unit
    adjust (bool): set==True to adjust input.geos unit values to adjusted values
    global_unit (bool): set units to globally relevent ones
    IUPAC_unit (bool): set units to IUPAC uints
    ClearFlo_unit (bool): set units to those used in the ClearFlo campaign
    debug (bool): legacy debug option, replaced by python logging

    Returns
    -------
    units (str), adjustment to give units (float), and scaling for units (float)

    Notes
    -----
     - Is this redundant now with the species class?
     - "Appropirate" unit is taken from GEOS-Chem input.geos
     - Option to use IUPAC unit. ( set IUPAC_unit==True )
    """
    d = read_yaml_file('species_units.yml')
    try:
        units = d[x]
    except KeyError:
        LogStr = 'provided species/tracer ({}) not in unit dict (assuming {})'
        units = 'pptv'
        logging.info(LogStr.format(x, units))

    # Adjust to appropriate scale for pf analysis
    if adjust:
        spec_2_pptv = GC_var('spec_2_pptv')
        spec_2_pptC = GC_var('spec_2_pptC')
        if (x in spec_2_pptv):
            if debug:
                print(('adjusting {} ({}) to {}'.format(x, units, 'pptv')))
            units = 'pptv'
        if (x in spec_2_pptC):
            if debug:
                print(('adjusting {} ({}) to {}'.format(x, units, 'pptC')))
            units = 'pptC'

    # Over ride adjustments for globally appro. units
    if global_unit:
        spec_2_ppbv = GC_var('spec_2_ppbv')
        spec_2_ppbC = GC_var('spec_2_ppbC')
        if (x in spec_2_ppbv):
            if debug:
                print(('adjusting {} ({}) to {}'.format(x, units, 'ppbv')))
            units = 'ppbv'
        if (x in spec_2_ppbC):
            if debug:
                print(('adjusting {} ({}) to {}'.format(x, units, 'ppbC')))
            units = 'ppbC'

    if ClearFlo_unit:
        units2ppbv = ['NO', 'MACR', 'MVK', 'PRPE', 'ALK4', 'ALD2']
        if any([x == i for i in units2ppbv]):
            units = 'ppbv'
        if any([x == i for i in ['PAN', ]]):
            units = 'pptv'
        if any([x == i for i in ['ISOP']]):
            units = 'ppbC'

    if use_pf_species_units:
        TRA_v10_species = [
            'A3O2', 'ATO2', 'B3O2', 'EOH', 'ETO2', 'ETP', 'GLYX', 'HO2', 'IAP', 'INO2', 'INPN', 'ISN1', 'ISNOOA', 'ISNOOB', 'ISNOHOO', 'ISNP', 'KO2', 'MAN2', 'MAO3', 'MAOP', 'MAOPO2', 'MCO3', 'MGLY', 'MO2', 'MRO2', 'MRP', 'OH', 'PO2', 'PP', 'PRN1', 'PRPN', 'R4N1', 'R4O2', 'R4P', 'RA3P', 'RB3P', 'RCO3', 'RIO2', 'ROH', 'RP', 'VRO2', 'VRP', 'LISOPOH', 'ISOPND', 'ISOPNB', 'HC5', 'DIBOO', 'HC5OO', 'DHMOB', 'MOBAOO', 'ISOPNBO2', 'ISOPNDO2', 'ETHLN', 'MACRN', 'MVKN', 'PYAC', 'IEPOXOO', 'ATOOH', 'PMNN', 'MACRNO2', 'PMNO2'
        ]
        if (x in TRA_v10_species):
            units = 'molec/cm3'

    if scale:
        scaleby = get_unit_scaling(units)

    if adjustment:
        if units == 'K':
            units = 'Deg. Celcuis'
            adjustby = -273.15
        else:
            adjustby = 0

    if IUPAC_unit:
        if units == 'ppbv':
            units = 'nmol mol$^{-1}$'
        if units == 'ppbC':
            units = 'nmol (C) mol$^{-1}$'
        if units == 'pptv':
            units = 'pmol mol$^{-1}$'
        if units == 'pptC':
            units = 'pmol (C) mol$^{-1}$'

    if scale and (not adjustment):
        return units, scaleby
    elif (scale and adjustment):
        return units, scaleby, adjustby
    else:
        return units


def get_ref_spec(spec='LIOx'):
    """
    Store of reference species for families

    Parameters
    ----------
    spec (str): species/tracer/variable name

    Returns
    -------
    ref_spec (str) reference species for a given family

    Notes
    -----
    This is for use in conbination  with functions that calculate relative values
    (e.g. in units of Ox, I, etc)
    """
    d = read_yaml_file('reference_families_for_species.yml')
    try:
        return d[spec]
    except KeyError:
        pstr = "WARNING: Just returning provided species ('{}') as ref_spec"
        print(pstr.format(spec))
        return spec


def GC_var(input_x=None, rtn_dict=False, debug=False):
    """
    return common variables variables lists (e.g. species in a family like NOx)

    Parameters
    ----------
    input_x (str): species/tracer/variable/etc name

    Returns
    -------
    values (e.g. list) from dictionary store

    Notes
    -----
     - Some older less self explanatory variables include:
    f_var = GC flux (EW, NS , UP) variables
    Ox = 'Ox', 'POX', 'LOX' + list of drydep species
    ( # not inc. 'NO3df', 'HNO4df', 'BrOdf' , 'BrNO2', 'IO', 'IONO', 'OIO', )
    Ox_p = Ox prod list
    Ox-l = Ox loss list
    d_dep  = dry dep (category, name = species)
    w_dep  = wet dep ( 3x categories
    ( WETDCV = rain out loss in convective updrafts (kg/s),
    WETDLS = rainout in large scale precip (kg/s),
    CV-FLX = Mass change due to cloud convection (kg/s); name = species)
    BL_m = UPWARD MASS FLUX FROM BOUNDARY-LAYER MIXING,
    (category, name = species)
    f_strat  = strat flux (to tropsosphere) (category, name = species)
    """
    if debug:
        print('GC_var called for {}'.format(input_x))
    d = read_yaml_file('family_variables.yml')
    if rtn_dict:
        return d
    else:
        return d[input_x]


def get_conversion_factor_kgX2kgREF(spec=None, ref_spec=None, stioch=None,
                                    debug=False):
    """
    Return conversion factor for mass (e.g. kg) X to mass ref_spec
    """
    # Get the reference species if not provided and calculate the mass ratio
    if isinstance(ref_spec, type(None)):
        ref_spec = get_ref_spec(spec)
        if debug:
            print("ref_spec for '{}': {}".format(spec, ref_spec))
    factor = 1/species_mass(spec)*species_mass(ref_spec)
    # Consider if there is a stoichiometry between ref_spec and species
    if isinstance(stioch, type(None)):
        stioch = spec_stoich(spec, ref_spec=ref_spec)
        if debug:
            print("stoich for '{}': {}".format(spec, stioch))
    if stioch != 1.0:
        factor = factor * stioch
    return factor


def get_GC_aerosol_species(YAML_filename='species_database_GCv12_9.yml',
                           path=None, AerosolVar='Is_Aerosol',
                           debug=False):
    """
    Retrieve a list of GEOS-Chem aerosol species
    """
    d = read_GC_species_database(path=path)
    AerosolSpecies = []
    for key in d.keys():
        try:
            Is_Aerosol = d[key][AerosolVar]
            if Is_Aerosol:
                AerosolSpecies += [key]
        except KeyError:
            PrtStr = "WARNING: Aerosol variable '{}' not present for '{}'"
            if debug:
                print(PrtStr.format(AerosolVar, key))
    # Include aerosol families in this list too
    AerosolSpecies += ['NIT-all', 'SO4-all', 'DST-all', 'DSTAL-all', 'SAL-all']
    return AerosolSpecies


def read_GC_species_database(YAML_filename='species_database_GCv12_9.yml',
                             path=None, debug=False
                             ):
    """
    Read in the GEOS-Chem YAML file
    """
    d = read_yaml_file(path=path, YAML_filename=YAML_filename, debug=debug)
    return d


def get_spec_properties():
    """
    Get the species properties using the GEOS-Chem json files in GCPy

    Notes
    ----
     - TODO: remove this function and other duplicates?
    """
    return read_GC_species_database()


def read_yaml_file(YAML_filename, path=None, debug=False):
    """
    Helper function to read yaml files
    """
    assert type(YAML_filename) == str, 'Passed filename argument must be str!'
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    if isinstance(path, type(None)):
        path = os.path.dirname(os.path.abspath(filename))
    with open(os.path.join(path, YAML_filename), "r") as stream:
        try:
            d = yaml.safe_load(stream)
            if debug:
                print(d)
        except yaml.YAMLError as exc:
            print(exc)
    return d
