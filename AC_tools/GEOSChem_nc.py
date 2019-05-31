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
#import csv
import glob
import pandas as pd
import xarray as xr
import re
from netCDF4 import Dataset
# try:
#    import iris
# except ImportError:
#    print('WARNING iris not imported')
import logging
# Math/Analysis
import numpy as np
#from time import mktime
#import scipy.stats as stats
#import math
# Time
import time
#import calendar
#import datetime as datetime
#from datetime import datetime as datetime_
# The below imports need to be updated,
# imports should be specific and in individual functions
# import tms modules with shared functions
from .core import *
from .generic import *
from .AC_time import *
from .planeflight import *
from .variables import *
#from .Scripts.bpch2netCDF import convert_to_netCDF


def GetGEOSChemFilesAsDataset(FileStr='GEOSChem.SpeciesConc.*.nc4', wd=None,
                              debug=False):
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
    files = glob.glob(wd+FileStr)
    assert len(files) >= 1, 'No files found matching-{}'.format(wd+FileStr)
    # open all of these files as single Dataset
    ds = xr.open_mfdataset(files)
    return ds


def GetTropBurdenInGg(ds=None, spec=None, SpecVar=None, StateMet=None, wd=None,
                      TropLevelVar='Met_TropLev', AirMassVar='Met_AD', AvgOverTime=False,
                      SumSpatially=True, RmTroposphere=True, use_time_in_trop=False,
                      TropMask=None, SpeciesConcPrefix='SpeciesConc_',
                      TimeInTropVar='N/A',
                      ):
    """
    Get Tropospheric burden for a/all species in dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    StateMet (dataset): Dataset object containing time in troposphere
    TropMask (nd.array): 3D or 4D boolean array where stratosphere is False
    spec (str): Name of the species (optional)
    SpecVar (str):  Name of the species inc. Diagnostic prefix (optional)
    SpeciesConcPrefix (str): the diagnostic prefix for concentration

    Returns
    -------
    (xr.dataset or pandas.DataFrame)

    Notes
    -----
     - A pandas dataframe is returned if values are requested to be summed spatially
     (e.g. SumSpatially=True), otherwise a dataset xr.dataset is returned.
    """
    # Only setup to take xarray datasets etc currently...
    assert type(StateMet) != None, 'Func. just setup to take StateMet currently'
    assert type(
        ds) == xr.Dataset, 'Func. just setup to take a xr.dataset currently'
    # Setup a dataset to process and return
    dsL = ds.copy()
    # Extract local variables
    AirMass = StateMet[AirMassVar]
#    TropLevel = StateMet[TropLevelVar]
    # Define SpecVar as the diga. prefix + "spec" if "spec" provided
    if not isinstance(spec, type(None)):
        if isinstance(SpecVar, type(None)):
            SpecVar = SpeciesConcPrefix+spec
    # Create mask for stratosphere if not provided
    if isinstance(TropMask, type(None)):
        TropMask = Create4DMask4TropLevel(StateMet=StateMet)
    # only allow "SpeciesConc" species
    Specs2Convert = [i for i in dsL.data_vars if 'SpeciesConc' in i]
    dsL = dsL[Specs2Convert]
    MXUnits = 'mol mol-1 dry'
    # Covert all species into burdens (Gg)
    if isinstance(SpecVar, type(None)) and isinstance(spec, type(None)):
        # Loop by spec
        for SpecVar in Specs2Convert:
            # Remove diagnostic prefix from chemical species
            spec = SpecVar.replace(SpeciesConcPrefix, '')
            # Check units
            SpecUnits = dsL[SpecVar].units
            MixingRatioUnits = MXUnits == SpecUnits
            assert_str = "Units must be in '{}' terms! (They are: '{}')"
            assert MixingRatioUnits, assert_str.format(MXUnits, SpecUnits)
            # v/v * (mass total of air (kg)/ 1E3 (converted kg to g)) = moles of tracer
            dsL[SpecVar] = dsL[SpecVar] * (AirMass*1E3 / constants('RMM_air'))
            # Convert moles to mass (* RMM) , then to Gg
            dsL[SpecVar] = dsL[SpecVar] * float(species_mass(spec)) / 1E9
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
        MixingRatioUnits = MXUnits == SpecUnits
        assert_str = "Units must be in '{}' terms! (They are: '{}')"
        assert MixingRatioUnits, assert_str.format(MXUnits, SpecUnits)
        # v/v * (mass total of air (kg)/ 1E3 (converted kg to g)) = moles of tracer
        dsL[SpecVar] = dsL[SpecVar] * (AirMass*1E3 / constants('RMM_air'))
        # Convert moles to mass (* RMM) , then to Gg
        dsL[SpecVar] = dsL[SpecVar] * float(species_mass(spec)) / 1E9
        # Return values averaged over time if requested
        if AvgOverTime:
            dsL = dsL.mean(dim='time')
        # remove the tropospheric values?
        if RmTroposphere:
            dsL[SpecVar] = dsL[SpecVar].where(TropMask)
    # Sum the values spatially?
    if SumSpatially:
        dsL = dsL.sum()
        return dsL.to_array().to_pandas()
    else:
        return dsL

def plot_up_surface_changes_between2runs( ds_dict=None, levs=[1], specs=[],
        BASE='', NEW='', prefix='IJ_AVG_S__', update_PyGChem_format2COARDS=False ):
    """
    Compare BASE and NEW datasets for given species using GCPy

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    BASE (str): name of the dataset (key) in ds_dict to have as the reference for changes
    NEW (str): name of the dataset (key) in ds_dict to compare against BASE
    specs (list): list of the names of the species to plot
    ds_dict (dict): dictionary of xr.datasets objects
    update_PyGChem_format2COARDS (bool): update the dataset names to be COARDS?
    prefix (str): category string that proceeds all variable names of specs
    levs (list): levels to plot spatial changes for

    Returns
    -------
    (xr.dataset or pandas.DataFrame)

    Notes
    -----

    """
    import gcpy
    # Species to plot
    vars2use = [prefix+i for i in specs]
    unit=None
    PDFfilenameStr = 'Oi_surface_change_{}_vs_{}_lev_{:0>2}'
    # Set datasets to use and  Just include the variables to plot in the dataset
    title1 = BASE
    title2 = NEW
    ds1 = ds_dict[ BASE ][ vars2use ].copy()
    ds2 = ds_dict[ NEW ][ vars2use ].copy()
    # Average over time
    print( ds1, ds2 )
    ds1 = ds1.mean(dim='time')
    ds2 = ds2.mean(dim='time')
    # Remove vestigial coordinates.
    # (e.g. the time_0 coord... what is this?)
    vars2drop = ['time_0']
    dsL = [ds1,ds2]
    for var2drop in vars2drop:
        for n, ds in enumerate(dsL):
            CoordVars = [ i for i in ds.coords ]
            if var2drop in CoordVars:
                ds = ds.drop(var2drop)
                dsL[n] = ds
    ds1, ds2 = dsL
    # Update dimension names
    if update_PyGChem_format2COARDS:
        ds1 = Convert_PyGChem_Iris_DataSet2COARDS_NetCDF(ds=ds1)
        ds2 = Convert_PyGChem_Iris_DataSet2COARDS_NetCDF(ds=ds2)
    # Now plot this using the compare_single_level script
    for lev in levs:
        # Just select surface (default) or lev in list provided
        ds1 = ds1.isel(lev=lev)
        ds2 = ds2.isel(lev=lev)
        # Make sure the units are present (xarray loses these after some actions)
        for var_ in vars2use:
            ds1[var_].attrs = ds_dict[title1][var_].attrs
            ds2[var_].attrs = ds_dict[title2][var_].attrs
        # Plot and save through GCPy
        PDFfilename = PDFfilenameStr.format( BASE, NEW, lev )
        gcpy.compare_single_level( ds1, title1, ds2, title2, varlist=vars2use,
                 ilev=0, pdfname=PDFfilename+'.pdf',)


def Create4DMask4TropLevel(StateMet=None,
                           TropLevelVar='Met_TropLev',
                           DynTropPressVar='Met_TropP',
                           PmidPress='Met_PMID',
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


def ReadInInstfilesSaveOnlySurfaceValues(wd=None, FileStr='GEOSChem.inst1hr.*',
                                         FileExtension='.nc4', SaveNewNetCDF=True,
                                         DeleteExistingNetCDF=True):
    """
    Extract just surface values and save as NetCDF (& DELETE old NetCDF)
    """
    import glob
    # Check input
    assert type(wd) == str, 'Working directory (wd) provided must be a string!'
    # Get files
    files = glob.glob(wd+FileStr)
    assert len(files) >= 1, 'No files found matching-{}'.format(wd+FileStr)
    for FullFileRoot in sorted(files):
        print(FullFileRoot)
        ds = xr.open_dataset(FullFileRoot)
        ds = ds.sel(lev=ds['lev'][0].values)
        # Create new file name
        Suffix = '_Just_surface'+FileExtension
        FileStrNew = FullFileRoot.replace(FileExtension, Suffix)
        print(FileStrNew)
        # Save and close file
        if SaveNewNetCDF:
            ds.to_netcdf(FileStrNew, engine='scipy')
        ds.close()
        # Delele old file?
        if DeleteExistingNetCDF:
            os.remove(FullFileRoot)


def GetSpeciesConcDataset(FileStr='GEOSChem.SpeciesConc.*.nc4', wd=None):
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
    return GetGEOSChemFilesAsDataset(FileStr=FileStr, wd=wd)


def GetInst1hrDataset(FileStr='GEOSChem.inst1hr.*', wd=None):
    """
    Wrapper to get NetCDF 1hr instantaneous (Inst1hr) output as a Dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    FileStr (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return GetGEOSChemFilesAsDataset(FileStr=FileStr, wd=wd)


def GetStateMetDataset(FileStr='GEOSChem.StateMet.*', wd=None):
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
    return GetGEOSChemFilesAsDataset(FileStr=FileStr, wd=wd)


def GetProdLossDataset(FileStr='GEOSChem.ProdLoss.*', wd=None):
    """
    Wrapper to get NetCDF ProdLoss output as a Dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    FileStr (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return GetGEOSChemFilesAsDataset(FileStr=FileStr, wd=wd)


def GetJValuesDataset(FileStr='GEOSChem.JValues.*', wd=None):
    """
    Wrapper to get NetCDF photolysis rates (Jvalues) output as a Dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    FileStr (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return GetGEOSChemFilesAsDataset(FileStr=FileStr, wd=wd)


def GetHEMCODiagnostics_AsDataset(FileStr='HEMCO_diagnostics.*', wd=None):
    """
    Wrapper to get HEMCO diagnostics NetCDF output as a Dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    FileStr (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return GetGEOSChemFilesAsDataset(FileStr=FileStr, wd=wd)


def Convert_PyGChem_Iris_DataSet2COARDS_NetCDF(ds=None, transpose_dims=True):
    """
    Convert a PyChem/Iris dataset into a COARDS compliant xr.dataset/NetCDF

    Parameters
    ----------
    ds (dataset): input Dataset object
    transpose_dims (bool): transpose the dimension order?

    Returns
    -------
    (dataset)
    """
    # PyGChem NetCDF (via iris backend ordering)
    PyGChem_Iris_order = ('time', 'longitude', 'latitude', 'dim3')
    # Make sure the Datasets are using the correct names
    rename_dims = {
    'dim3': 'lev', 'latitude': 'lat', 'longitude': 'lon',
    'model_level_number':'lev'
    }
    for CoordVar in ds.coords:
        try:
            ds = ds.rename(name_dict={CoordVar: rename_dims[CoordVar]})
        except KeyError:
            print('Not renamed {} as not in Dataset'.format(CoordVar))
    # Update the ordering to follow COARDS
    if transpose_dims:
        CoordVars = [ i for i in ds.coords ]
        if 'time' in CoordVars:
            COARDS_order = ('time', 'lev', 'lat', 'lon')
        else:
            COARDS_order = ('lev', 'lat', 'lon')
        # Transpose to ordering
        ds = ds.transpose(*COARDS_order)
    return ds


def Convert_HEMCO_ds2Gg_per_yr( ds, vars2convert=None, var_species_dict=None,
                          Output_freq='End', verbose=False, debug=False):
    """
    Convert emissions in HEMCO dataset to mass/unit time

    vars2convert (list), NetCDF vairable names to convert
    var_species_dict (dict), dictionary to map variables names to chemical species
    Output_freq (str), output frequency dataset made from HEMCO NetCDF file output

    """
    # Get chemical species for each variable name
    var_species = {}
    for var in vars2convert:
        try:
            var_species[var] = var_species_dict[var]
        except:
#            if verbose:
            PrtStr = "WARNING - using variable name '{}' as chemical species!"
            print( PrtStr.format( var) )
            var_species[var] = var
    # Print assumption about end units.
    if Output_freq == 'End':
        print( "WARNING - Assuming Output frequnecy ('End') is monthly")

    # Get equivalent unit for chemical species (e.g. I, Br, Cl, N, et c)
    ref_specs = {}
    for var in vars2convert:
        ref_specs[var] = get_ref_spec( var_species[var] )
    # Loop dataset by variable
    for var_n, var_ in enumerate(vars2convert):
        # extract var arr
        arr = ds[var_].values
        # --- Adjust units to be in kg/gridbox
        # remove area units
        if ds[var_].units == 'kg/m2/':
            arr = arr * ds['AREA']
        elif ds[var_].units == 'kg/m2/s':
            arr = arr * ds['AREA']
            # now remove seconds
            if Output_freq == 'Hourly':
                arr = arr*60.*60.
            elif Output_freq == 'Daily':
                arr = arr*60.*60.*24.
            elif Output_freq == 'Weekly':
                arr = arr*60.*60.*24.*(365./52.)
            elif (Output_freq == 'Monthly') or (Output_freq == 'End'):
                arr = arr*60.*60.*24.*(365./12.)
            else:
                print('WARNING: ({}) output convert. unknown'.format(
                    Output_freq))
                sys.exit()
        elif ds[var_].units  == 'kg':
            pass # units are already in kg .
        else:
            print('WARNING: unit convert. ({}) unknown'.format(unit_dict[var_]))
            sys.exit()
        # --- convert to Gg species
        # get spec name for output variable
        spec = var_species[var_]
        # Get equivalent unit for species (e.g. I, Br, Cl, N, et c)
        ref_spec = ref_specs[var_]
        # get stiochiometry of ref_spec in species
        stioch = spec_stoich(spec, ref_spec=ref_spec)
        RMM_spec = species_mass(spec)
        RMM_ref_spec = species_mass(ref_spec)
        # update values in array
        arr = arr / RMM_spec * RMM_ref_spec * stioch
        # (from kg=>g (*1E3) to g=>Gg (/1E9))
        arr = arr*1E3 / 1E9
        if set(ref_specs) == 1:
            units = '(Gg {})'.format(ref_spec)
        else:
            units = '(Gg X)'
        if debug:
            print(arr.shape)
        # reassign arrary
        ds[var_].values = arr
        # Update units too
        attrs = ds[var_].attrs
        attrs['units'] = units
        ds[var_].attrs = attrs
    return ds


def GetSummaryStatsOnHEMCOds_in_Gg_per_yr( ds, vars2use=None ):
    """
    Get summary statistics on dataframe of data
    """
    # master list to hold values
    master_l = []

    # loop by species
    for var_n, var_ in enumerate(vars2use):

        # sub list to hold calculated values
        sub_l = []
        headers = []

        # --- process data to summary statistics
        sum_over_lat_lon = arr.sum(axis=-1).sum(axis=-1).copy()

        # If monthly... process useful summary stats...
        Monthly_Output_freqs = 'Monthly', 'End'
        if Output_freq in Monthly_Output_freqs:
            if Output_freq == 'End':
                print(('!'*100, 'WARNING: End output assumed to monthly!'))

            # Add a monthly total
            sub_l += [sum_over_lat_lon.mean(axis=0)]
            headers += ['Mon. avg {}'.format(units)]

            # Add a monthly max
            sub_l += [sum_over_lat_lon.max(axis=0)]
            headers += ['Mon. max {}'.format(units)]

            # Add a monthly min
            sub_l += [sum_over_lat_lon.min(axis=0)]
            headers += ['Mon. max {}'.format(units)]

            # Annual equi.
            sub_l += [arr.sum() / len(arr[:, 0, 0])*12]
            headers += ['Ann. equiv. {}'.format(units)]

            # Annual equi.
            sub_l += [arr.sum() / len(arr[:, 0, 0])*12/1E3]
            header_ = 'Ann. equiv. (Tg {})'.format(ref_spec)
            if len(set(ref_specs)) > 1:
                header_ = 'Ann. equiv. (Tg X)'.format(ref_spec)
            headers += [header_]

        # If daily?!
        elif Output_freq == 'Daily':

            # Add a monthly total
            sub_l += [sum_over_lat_lon.mean(axis=0)]
            headers += ['Daily avg {}'.format(units)]

            # Annual equi.
            sub_l += [arr.sum() / len(arr[:, 0, 0])*365.]
            headers += ['Ann. equiv. {}'.format(units)]

            # Annual equi.
            sub_l += [arr.sum() / len(arr[:, 0, 0])*365./1E3]
            header_ = ['Ann. equiv. (Tg {})'.format(ref_spec)]
            if len(set(ref_specs)) > 1:
                header_ = 'Ann. equiv. (Tg X)'.format(ref_spec)
            headers += [header_]

        else:
            prt_str = 'WARNING: no processing setup for {} output'
            print(prt_str.format(Output_freq))
            sys.exit()

        # save to master list
        master_l += [sub_l]

    # Make into a DataFrame
    df = pd.DataFrame(master_l)
    df.columns = headers
    df.index = variables
    if verbose:
        print(df)
    return df

