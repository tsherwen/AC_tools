#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Functions for use with the GEOS-Chem chemical transport model (CTM).

Use help(<name of function>) to get details on a particular function.

"""
# - Required modules:
# I/O / Low level
import os
import sys
import csv
import glob
import pandas as pd
import xarray as xr
import re
from netCDF4 import Dataset
#try:
#    import iris
#except ImportError:
#    print('WARNING iris not imported')
import logging
# Math/Analysis
import numpy as np
from time import mktime
import scipy.stats as stats
import math
# Time
import time
import calendar
import datetime as datetime
from datetime import datetime as datetime_
# The below imports need to be updated,
# imports should be specific and in individual functions
# import tms modules with shared functions
from .funcs4core import *
from .funcs4generic import *
from .funcs4time import *
from .funcs4pf import *
from .funcs_vars import *
#from .Scripts.bpch2netCDF import convert_to_netCDF



def GetGEOSChemFilesAsDataset( FileStr='GEOSChem.SpeciesConc.*.nc4', wd=None,
                               debug=False ):
    """
    Extract GEOS-Chem NetCDF files that match file string format to a xr.dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    StateMet (dataset): Dataset object containing time in troposphere

    Returns
    -------
    (dataset)

    Notes
    """
    import glob
    # Check input
    assert type(wd) == str, 'Working directory (wd) provided must be a string!'
    # Get files
    files = glob.glob( wd+FileStr )
    assert len(files) >= 1, 'No files found matching-{}'.format( wd+FileStr )
    # open all of these files as single Dataset
    ds = xr.open_mfdataset( files )
    return ds


def GetTropBurdenInGg( ds=None, SpecVar=None, StateMet=None, wd=None, TimeInTropVar='N/A',
        TropLevelVar='Met_TropLev', AirMassVar='Met_AD', AvgOverTime=False, TropMask=None,
        SumSpatially=True, RmTroposphere=True, use_time_in_trop=False,
        SpeciesConcPrefix='SpeciesConc_'  ):
    """
    Get Tropospheric burden for a/all species in dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    StateMet (dataset): Dataset object containing time in troposphere

    Returns
    -------
    (dataset or pandas.DataFrame)

    Notes
    -----
     - A pandas dataframe is returned
    """
    # Only setup to take xarray datasets etc currently...
    assert type(StateMet) != None, 'Func. just setup to take StateMet currently'
    assert type(ds) == xr.Dataset, 'Func. just setup to take a xr.dataset currently'
    # Setup a dataset to process and return
    dsL = ds.copy()
    # Extract local variables
    AirMass = StateMet[AirMassVar]
    TropLevel = StateMet[TropLevelVar]
    # Create mask for stratosphere if not provided
    if isinstance( TropMask, type(None) ):
        TropMask = Create4DMask4TropLevel( StateMet=StateMet,
                                      use_time_in_trop=use_time_in_trop,
                                      )
    # only allow "SpeciesConc" species
    Specs2Convert = [i for i in dsL.data_vars if 'SpeciesConc' in i]
    dsL = dsL[Specs2Convert]
    # Covert all species into burdens (Gg)
    if isinstance( SpecVar, type(None) ):
        # Loop by spec
        for SpecVar in Specs2Convert:
            # Remove diagnostic prefix from chemical species
            spec = SpecVar.replace(SpeciesConcPrefix, '')
            # Check units
            SpecUnits = dsL[SpecVar].units
            MixingRatioUnits = SpecUnits ==  'mol mol-1 dry'
            assert MixingRatioUnits, 'Units must be in v/v terms! ({})'.format(SpecUnits)
            # v/v * (mass total of air (kg)/ 1E3 (converted kg to g)) = moles of tracer
            dsL[SpecVar] = dsL[SpecVar] * (AirMass*1E3 / constants('RMM_air'))
            # Convert moles to mass (* RMM) , then to Gg
            dsL[SpecVar]  = dsL[SpecVar] * float(species_mass(spec)) / 1E9
        # Return values averaged over time if requested
        if AvgOverTime:
            dsL = dsL.mean(dim='time')
        # Remove the tropospheric values?
        if RmTroposphere:
            # Loop by spec
            for SpecVar in Specs2Convert:
                 dsL[SpecVar] = dsL[SpecVar].where(TropMask)

    else:
        # Just consifer the species of interest
        dsL = dsL[[SpecVar]]
        # Remove diagnostic prefix from chemical species
        spec = SpecVar.replace(SpeciesConcPrefix, '')
        # Check units
        SpecUnits = dsL[SpecVar].units
        MixingRatioUnits = SpecUnits == 'mol mol-1 dry'
        assert MixingRatioUnits, 'Units must be in v/v terms! ({})'.format(SpecUnits)
        # v/v * (mass total of air (kg)/ 1E3 (converted kg to g)) = moles of tracer
        dsL[SpecVar] = dsL[SpecVar] * (AirMass*1E3 / constants('RMM_air'))
        # Convert moles to mass (* RMM) , then to Gg
        dsL[SpecVar]  = dsL[SpecVar] * float(species_mass(spec)) / 1E9
        # Return values averaged over time if requested
        if AvgOverTime:
            dsL = dsL.mean(dim='time')
        # remove the tropospheric values?
        if RmTroposphere:
            dsL[SpecVar] = dsL[SpecVar].where(TropMask)

    # Sum the values spatially?
    if SumSpatially:
        dsL= dsL.sum()
        return dsL.to_array().to_pandas()
    else:
        return dsL


def Create4DMask4TropLevel( StateMet=None,
                            TropLevelVar='Met_TropLev',
                            DynTropPressVar='Met_TropP',
                            PmidPress = 'Met_PMID',
                            use_time_in_trop=False,
                            ):
    """
    Create a mask to remove the stratosphere from GEOSChem output
    """
    # Extract local variables
#    TropLevel = StateMet[TropLevelVar]
    DynTropPress = StateMet[DynTropPressVar]
    # PMID: Pressure at average pressure level
    PmidPress = StateMet[PmidPress]
    # Just use an integer value for this for now
#    TropLevel = TropLevel.astype(int)
#    MASK = (StateMet['Met_PMID'] > )
    # Just mask as middle pressures values above dynamic troposphere for now
    # this can then be used like ds[VarName].where( MASK )
    # and summed via np.nansum( ds[VarName].where(MASK).values )
    MASK = PmidPress > DynTropPress
    return MASK


def GetSpeciesConcDataset( FileStr='GEOSChem.SpeciesConc.*.nc4', wd=None ):
    """
    Wrapper to retrive GEOSChem SpeciesConc NetCDFs as a xr.dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    FileStr (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return GetGEOSChemFilesAsDataset( FileStr=FileStr, wd=wd )


def GetStateMetDataset( FileStr='GEOSChem.StateMet.*', wd=None ):
    """
    Wrapper to get NetCDF StateMet output as a Dataset
    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    FileStr (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return GetGEOSChemFilesAsDataset( FileStr=FileStr, wd=wd )


def GetProdLossDataset( FileStr='GEOSChem.ProdLoss..*', wd=None ):
    """
    Wrapper to get NetCDF StateMet output as a Dataset
    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    FileStr (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return GetGEOSChemFilesAsDataset( FileStr=FileStr, wd=wd )

