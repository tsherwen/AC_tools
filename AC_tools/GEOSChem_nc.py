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
import datetime as datetime
from datetime import datetime as datetime_
# The below imports need to be updated,
# imports should be specific and in individual functions
# import tms modules with shared functions
from .core import *
from .generic import *
from .AC_time import *
from .planeflight import *
from .variables import *
#from .Scripts.bpch2netCDF import convert_to_netCDF
import gc


def get_GEOSChem_files_as_ds(file_str='GEOSChem.SpeciesConc.*.nc4', wd=None,
                             collection=None, dates2use=None,
                             parallel=True, data_vars="minimal",
                             coords="minimal", compat="override",
                             combine='by_coords',
                             open_with_coords_dropped=False,
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
    # Get files (and check if the a HEMCO collection is being requested)
    if isinstance(collection, str):
        glob_pattern = '{}/*.{}.*'.format(wd, collection)
        is_HEMCO_collection = ('hemco' in collection.lower())
    else:
        glob_pattern = '{}/{}'.format(wd, file_str)
        is_HEMCO_collection = ('hemco' in file_str.lower())
    files = glob.glob(glob_pattern)
    assert len(files) >= 1, 'No files found matching-{}'.format(glob_pattern)
    # Sort the files based on their name (which contains a regular datastring)
    files = list(sorted(files))
    # Only open dates for certain dates?
    if not isinstance(dates2use, type(None)):
        FileRootsVar = 'FileRoots'
        df = pd.DataFrame(files)
        df = pd.DataFrame({FileRootsVar: files})
        # Which date format to look for in filenames?
        if is_HEMCO_collection:
            format = '%Y%m%d%H%M'
        else:
            format = '%Y%m%d_%H%Mz'
        # Setup a helper function to extract dates from file strings

        def get_date_from_filename(x, format=format):
            """
            Extract Dates from filenames

            Notes
            -------
             - It is assumed that the date ends the file string before the
             format identifier
            """
            date_str = x.split('.')[-2]
            dt = datetime_.strptime(date_str, format)
            return dt
        dtVar = 'datetime'
        df[dtVar] = df[FileRootsVar].map(get_date_from_filename)
        bool = df[dtVar].isin(dates2use)
        files = list(df.loc[bool, FileRootsVar].values)
    # Open all of these files as single Dataset
    if open_with_coords_dropped:
        def drop_all_coords(ds):
            return ds.reset_coords(drop=True)
        ds = xr.open_mfdataset(files, combine='by_coords',
                               preprocess=drop_all_coords)
    else:
        # NOTE: Updated to use faster opening settings for files sharing the same coords
        # https://github.com/pydata/xarray/issues/1823
        ds = xr.open_mfdataset(files,
                               #                           concat_dim='time',
                               combine=combine,
                               data_vars=data_vars, coords=coords,
                               compat=compat, parallel=parallel)
        #

    return ds


def get_Gg_trop_burden(ds=None, spec=None, spec_var=None, StateMet=None,
                       wd=None,
                       trop_level_var='Met_TropLev', air_mass_var='Met_AD',
                       avg_over_time=False,
                       sum_patially=True, rm_strat=True, use_time_in_trop=True,
                       trop_mask=None, spec_conc_prefix='SpeciesConc_',
                       time_in_trop_var='N/A', vars2use=None,
                       sum_spatially=True,
                       debug=False,
                       ):
    """
    Get Tropospheric burden for given/or all species in dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    StateMet (dataset): Dataset object containing time in troposphere
    trop_mask (nd.array): 3D or 4D boolean array where stratosphere is False
    rm_strat (bool): remove the stratospheric values
    spec (str): Name of the species (optional)
    spec_var (str):  Name of the species inc. Diagnostic prefix (optional)
    spec_conc_prefix (str): the diagnostic prefix for concentration
    vars2use (list): list of variables to calculate burden for

    Returns
    -------
    (xr.dataset or pandas.DataFrame)

    Notes
    -----
     - A pandas dataframe is returned if values are requested to be summed spatially
     (e.g. sum_patially=True), otherwise a dataset xr.dataset is returned.
    """
    # Only setup to take xarray datasets etc currently...
    ass_Str = 'WATNING: Func. just setup to take StateMet currently'
    assert type(StateMet) != None, ass_Str
    ass_Str = 'WATNING: Func. just setup to take a xr.dataset currently'
    assert type(ds) == xr.Dataset, ass_Str
    # Setup a dataset to process and return
    dsL = ds.copy()
    # Extract local variables
    AirMass = StateMet[air_mass_var]
#    TropLevel = StateMet[trop_level_var]
    # Define spec_var as the diga. prefix + "spec" if "spec" provided
    if not isinstance(spec, type(None)):
        if isinstance(spec_var, type(None)):
            spec_var = spec_conc_prefix+spec
    # Create mask for stratosphere if not provided
    if isinstance(trop_mask, type(None)) and not use_time_in_trop:
        trop_mask = create4Dmask4trop_level(StateMet=StateMet)
    # Only allow "SpeciesConc" species
    ass_Str = 'WARNING: Duplicates found in vars2use list!'
    assert len(set(vars2use)) == len(vars2use), ass_Str
    if isinstance(vars2use, type(None)):
        vars2use = [i for i in dsL.data_vars if 'SpeciesConc' in i]
    # Only try to extract requested species if they are in the dataset
    NotInDataset = [i for i in vars2use if (i not in dsL.data_vars)]
    if len(NotInDataset) > 1:
        print('WARNING: not extracing variables not in ds: ', NotInDataset)
        vars2use = [i for i in vars2use if (i not in NotInDataset)]
    dsL = dsL[vars2use]
    MXUnits = 'mol mol-1 dry'
    # Covert all species into burdens (Gg)
    if isinstance(spec_var, type(None)) and isinstance(spec, type(None)):
        # Loop by spec
        for spec_var in vars2use:
            # Remove diagnostic prefix from chemical species
            spec = spec_var.replace(spec_conc_prefix, '')
            if debug:
                PStr = 'Attempting in ds conversion of {} ({}, type: {})'
                print(PStr.format(spec, spec_var, type(dsL[spec_var])))
                print(dsL[spec_var].attrs)
            # Check units
            SpecUnits = dsL[spec_var].units
            MixingRatioUnits = MXUnits == SpecUnits
            assert_str = "Units must be in '{}' terms! (They are: '{}')"
            assert MixingRatioUnits, assert_str.format(MXUnits, SpecUnits)
            # v/v * (mass total of air (kg)/ 1E3 (converted kg to g)) = moles of tracer
            conversion_factor = (AirMass*1E3 / constants('RMM_air'))
            dsL[spec_var] = dsL[spec_var] * conversion_factor
            # Convert moles to mass (* RMM) , then to Gg
            dsL[spec_var] = dsL[spec_var] * float(species_mass(spec)) / 1E9
        # Remove the non-tropospheric values?
        if rm_strat:
            # Loop by spec
            if use_time_in_trop:
                dsL = rm_fractional_troposphere(dsL, vars2use=vars2use,
                                                StateMet=StateMet)
            else:
                for spec_var in vars2use:
                    dsL[spec_var] = dsL[spec_var].where(trop_mask)
    else:
        # Just consifer the species of interest
        dsL = dsL[[spec_var]]
        # Remove diagnostic prefix from chemical species
        spec = spec_var.replace(spec_conc_prefix, '')
        # Check units
        SpecUnits = dsL[spec_var].units
        MixingRatioUnits = MXUnits == SpecUnits
        assert_str = "Units must be in '{}' terms! (They are: '{}')"
        assert MixingRatioUnits, assert_str.format(MXUnits, SpecUnits)
        # v/v * (mass total of air (kg)/ 1E3 (converted kg to g)) = moles of tracer
        dsL[spec_var] = dsL[spec_var] * (AirMass*1E3 / constants('RMM_air'))
        # Convert moles to mass (* RMM) , then to Gg
        dsL[spec_var] = dsL[spec_var] * float(species_mass(spec)) / 1E9
        # Remove the non-tropospheric values?
        if rm_strat:
            # Loop by spec
            if use_time_in_trop:
                dsL = rm_fractional_troposphere(dsL, vars2use=vars2use,
                                                StateMet=StateMet)
            else:
                for spec_var in vars2use:
                    dsL[spec_var] = dsL[spec_var].where(trop_mask)
    # Return values averaged over time if requested
    if avg_over_time:
        dsL = dsL.mean(dim='time')
    # Sum the values spatially?
    if sum_spatially:
        dsL = dsL.sum()
        return dsL.to_array().to_pandas()
    else:
        return dsL


def plot_up_surface_changes_between2runs(ds_dict=None, levs=[1], specs=[],
                                         BASE='', NEW='', prefix='IJ_AVG_S__',
                                         update_PyGChem_format2COARDS=False):
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
    unit = None
    PDFfilenameStr = 'Oi_surface_change_{}_vs_{}_lev_{:0>2}'
    # Set datasets to use and  Just include the variables to plot in the dataset
    title1 = BASE
    title2 = NEW
    ds1 = ds_dict[BASE][vars2use].copy()
    ds2 = ds_dict[NEW][vars2use].copy()
    # Average over time
    print(ds1, ds2)
    ds1 = ds1.mean(dim='time')
    ds2 = ds2.mean(dim='time')
    # Remove vestigial coordinates.
    # (e.g. the time_0 coord... what is this?)
    vars2drop = ['time_0']
    dsL = [ds1, ds2]
    for var2drop in vars2drop:
        for n, ds in enumerate(dsL):
            CoordVars = [i for i in ds.coords]
            if var2drop in CoordVars:
                ds = ds.drop(var2drop)
                dsL[n] = ds
    ds1, ds2 = dsL
    # Update dimension names
    if update_PyGChem_format2COARDS:
        ds1 = convert_pyGChem_iris_ds2COARDS_ds(ds=ds1)
        ds2 = convert_pyGChem_iris_ds2COARDS_ds(ds=ds2)
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
        PDFfilename = PDFfilenameStr.format(BASE, NEW, lev)
        gcpy.compare_single_level(ds1, title1, ds2, title2, varlist=vars2use,
                                  ilev=0, pdfname=PDFfilename+'.pdf',)


def create4Dmask4trop_level(StateMet=None,
                            trop_level_var='Met_TropLev',
                            dyn_trop_press_var='Met_TropP',
                            pmid_press='Met_PMID',
                            use_time_in_trop=False,
                            ):
    """
    Create a mask to remove the stratosphere from GEOSChem output
    """
    # Extract local variables
#    TropLevel = StateMet[trop_level_var]
    DynTropPress = StateMet[dyn_trop_press_var]
    # PMID: Pressure at average pressure level
    pmid_press = StateMet[pmid_press]
    # Just use an integer value for this for now
#    TropLevel = TropLevel.astype(int)
#    MASK = (StateMet['Met_PMID'] > )
    # Just mask as middle pressures values above dynamic troposphere for now
    # This can then be used like ds[VarName].where( MASK )
    # and summed via np.nansum( ds[VarName].where(MASK).values )
    MASK = pmid_press > DynTropPress
    return MASK


def rm_fractional_troposphere(ds, vars2use=None, StateMet=None,
                              TropFracDA=None,
                              TropFracVar='FracOfTimeInTrop'):
    """
    rm stratospheric contribution by multiplying by grid box's time in troposphere
    """
    # Get the fraction of boxes in the troposphere as a data array
    if isinstance(TropFracDA, type(None)):
        TropFracDA = StateMet[TropFracVar]
    # Set the variables to use as all data_vars if list not provided
    if isinstance(vars2use, type(None)):
        vars2use = list(ds.data_vars)
    # Multiply through by time in troposphere to remove stratosphere.
    for var in vars2use:
        ds[var] = ds[var] * TropFracDA
    return ds


def read_inst_files_save_only_surface(wd=None, file_str='GEOSChem.inst1hr.*',
                                      file_extension='.nc4',
                                      save_new_NetCDF=True,
                                      delete_existing_NetCDF=True):
    """
    Extract just surface values and save as NetCDF (& DELETE old NetCDF)
    """
    import glob
    # Check input
    assert type(wd) == str, 'Working directory (wd) provided must be a string!'
    # Get files
    files = glob.glob(wd+file_str)
    assert len(files) >= 1, 'No files found matching-{}'.format(wd+file_str)
    for FullFileRoot in sorted(files):
        print(FullFileRoot)
        ds = xr.open_dataset(FullFileRoot)
        ds = ds.sel(lev=ds['lev'][0].values)
        # Create new file name
        Suffix = '_Just_surface'+file_extension
        file_strNew = FullFileRoot.replace(file_extension, Suffix)
        print(file_strNew)
        # Save and close file
        if save_new_NetCDF:
            ds.to_netcdf(file_strNew, engine='scipy')
        ds.close()
        # Delele old file?
        if delete_existing_NetCDF:
            os.remove(FullFileRoot)


def GetSpeciesConcDataset(file_str='GEOSChem.SpeciesConc.*.nc4', wd=None,
                          dates2use=None):
    """
    Wrapper to retrive GEOSChem SpeciesConc NetCDFs as a xr.dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    file_str (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return get_GEOSChem_files_as_ds(file_str=file_str, wd=wd,
                                    dates2use=dates2use)


def get_Inst1hr_ds(file_str='GEOSChem.inst1hr.*', wd=None,
                   dates2use=None):
    """
    Wrapper to get NetCDF 1hr instantaneous (Inst1hr) output as a Dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    file_str (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return get_GEOSChem_files_as_ds(file_str=file_str, wd=wd,
                                    dates2use=dates2use)


def get_StateMet_ds(file_str='GEOSChem.StateMet.*', wd=None,
                    dates2use=None):
    """
    Wrapper to get NetCDF StateMet output as a Dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    file_str (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return get_GEOSChem_files_as_ds(file_str=file_str, wd=wd,
                                    dates2use=dates2use)


def get_DryDep_ds(file_str='GEOSChem.DryDep.*', wd=None,
                  dates2use=None):
    """
    Wrapper to get NetCDF dry deposition output as a dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    file_str (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return get_GEOSChem_files_as_ds(file_str=file_str, wd=wd,
                                    dates2use=dates2use)


def get_WetLossConv_ds(file_str='GEOSChem.WetLossConv.*', wd=None,
                       dates2use=None):
    """
    Wrapper to get NetCDF Wet Loss via Convection output as a dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    file_str (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return get_GEOSChem_files_as_ds(file_str=file_str, wd=wd,
                                    dates2use=dates2use)


def get_WetLossLS_ds(file_str='GEOSChem.WetLossLS.*', wd=None,
                     dates2use=None):
    """
    Wrapper to get Wet Loss via large-scale Convection NetCDF output as a ds

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    file_str (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return get_GEOSChem_files_as_ds(file_str=file_str, wd=wd,
                                    dates2use=dates2use)


def get_ProdLoss_ds(file_str='GEOSChem.ProdLoss.*', wd=None,
                    dates2use=None):
    """
    Wrapper to get NetCDF ProdLoss output as a Dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    file_str (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return get_GEOSChem_files_as_ds(file_str=file_str, wd=wd,
                                    dates2use=dates2use)


def GetJValuesDataset(file_str='GEOSChem.JValues.*', wd=None,
                      dates2use=None):
    """
    Wrapper to get NetCDF photolysis rates (Jvalues) output as a Dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    file_str (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return get_GEOSChem_files_as_ds(file_str=file_str, wd=wd,
                                    dates2use=dates2use)


def get_HEMCO_diags_as_ds(file_str='HEMCO_diagnostics.*', wd=None,
                          dates2use=None):
    """
    Wrapper to get HEMCO diagnostics NetCDF output as a Dataset

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    file_str (str): a str for file format with wildcards (?, *)

    Returns
    -------
    (dataset)
    """
    return get_GEOSChem_files_as_ds(file_str=file_str, wd=wd,
                                    dates2use=dates2use)


def convert_pyGChem_iris_ds2COARDS_ds(ds=None, transpose_dims=True):
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
        'model_level_number': 'lev'
    }
    for CoordVar in ds.coords:
        try:
            ds = ds.rename(name_dict={CoordVar: rename_dims[CoordVar]})
        except KeyError:
            print('Not renamed {} as not in Dataset'.format(CoordVar))
    # Update the ordering to follow COARDS
    if transpose_dims:
        CoordVars = [i for i in ds.coords]
        if 'time' in CoordVars:
            COARDS_order = ('time', 'lev', 'lat', 'lon')
        else:
            COARDS_order = ('lev', 'lat', 'lon')
        # Transpose to ordering
        ds = ds.transpose(*COARDS_order)
    return ds


def convert_HEMCO_ds2Gg_per_yr(ds, vars2convert=None, var_species_dict=None,
                               output_freq='End',
                               convert_unaveraged_time=False,
                               verbose=False, debug=False):
    """
    Convert emissions in HEMCO dataset to mass/unit time

    vars2convert (list), NetCDF vairable names to convert
    var_species_dict (dict), dictionary to map variables names to chemical species
    output_freq (str), output frequency dataset made from HEMCO NetCDF file output

    """
    # Get chemical species for each variable name
    var_species = {}
    for var in vars2convert:
        try:
            var_species[var] = var_species_dict[var]
        except:
            #            if verbose:
            PrtStr = "WARNING - using variable name '{}' as chemical species!"
            print(PrtStr.format(var))
            var_species[var] = var
    # Print assumption about end units.
    if output_freq == 'End':
        print("WARNING - Assuming Output frequnecy ('End') is monthly")

    # Get equivalent unit for chemical species (e.g. I, Br, Cl, N, et c)
    ref_specs = {}
    for var in vars2convert:
        try:
            ref_specs[var] = get_ref_spec(var_species[var])
        except KeyError:
            PStr = "WARNING: Using '{}' as reference species for '{}'"
            print(PStr.format(var, var))
    # Loop dataset by variable
    for var_n, var_ in enumerate(vars2convert):
        if debug:
            print('{:<2} {} '.format(var_n, var_))
        # Extract the variable array
        try:
            arr = ds[var_].values
        except KeyError:
            PStr = "WARNING: skipping variable '({})' as not in dataset"
            print(PStr.format(var_))
            continue

        # --- Adjust units to be in kg/gridbox
        # Remove area units
        if ds[var_].units == 'kg/m2/':
            arr = arr * ds['AREA']
        elif ds[var_].units == 'kg/m2/s':
            arr = arr * ds['AREA']
            # now remove seconds
            if convert_unaveraged_time:
                if output_freq == 'Hourly':
                    arr = arr*60.*60.
                elif output_freq == 'Daily':
                    arr = arr*60.*60.*24.
                elif output_freq == 'Weekly':
                    arr = arr*60.*60.*24.*(365./52.)
                elif (output_freq == 'Monthly') or (output_freq == 'End'):
                    arr = arr*60.*60.*24.*(365./12.)
                else:
                    print('WARNING: ({}) output convert. unknown'.format(
                        output_freq))
                    sys.exit()
            else:
                arr = arr*60.*60.*24.*365.
        elif ds[var_].units == 'kg':
            pass  # Units are already in kg .
        else:
            print('WARNING: unit convert. ({}) unknown'.format(ds[var_].units))
            sys.exit()
        # - Convert to Gg species
        # Get spec name for output variable
        spec = var_species[var_]
        # Get equivalent unit for species (e.g. I, Br, Cl, N, et c)
        ref_spec = ref_specs[var_]
        # Get stoichiometry of ref_spec in species
        stioch = spec_stoich(spec, ref_spec=ref_spec)
        RMM_spec = species_mass(spec)
        RMM_ref_spec = species_mass(ref_spec)
        # Update values in array
        arr = arr / RMM_spec * RMM_ref_spec * stioch
        # (from kg=>g (*1E3) to g=>Gg (/1E9))
        arr = arr*1E3 / 1E9
        if set(ref_specs) == 1:
            units = '(Gg {})'.format(ref_spec)
        else:
            units = '(Gg X)'
        if debug:
            print(arr.shape)
        # Reassign arrary
        ds[var_].values = arr
        # Update units too
        attrs = ds[var_].attrs
        attrs['units'] = units
        ds[var_].attrs = attrs
    return ds


def get_HEMCO_ds_summary_stats_Gg_yr(ds, vars2use=None, ref_spec=None):
    """
    Get summary statistics on dataframe of data
    """
    # Master list to hold values
    master_l = []

    # Loop by species
    for var_n, var_ in enumerate(vars2use):

        # Sub list to hold calculated values
        sub_l = []
        headers = []
        # Get reference species for variable

        # Process data to summary statistics
        sum_over_lat_lon = arr.sum(axis=-1).sum(axis=-1).copy()

        # If monthly... process useful summary stats...
        Monthly_output_freqs = 'Monthly', 'End'
        if output_freq in Monthly_output_freqs:
            if output_freq == 'End':
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
        elif output_freq == 'Daily':

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
            print(prt_str.format(output_freq))
            sys.exit()

        # Save to master list
        master_l += [sub_l]

    # Make into a DataFrame
    df = pd.DataFrame(master_l)
    df.columns = headers
    df.index = variables
    if verbose:
        print(df)
    return df


def AddChemicalFamily2Dataset(ds, fam='NOy', prefix='SpeciesConc_'):
    """
    Add a variable to dataset for a chemical family (e.g. NOy, NOx...)
    """
    if fam == 'NOx':
        fam2use = 'NOx'
        CopyVar = 'NO'
        ds[prefix+fam2use] = ds[prefix+CopyVar].copy()
        ds[prefix+fam2use] = ds[prefix+fam2use] + ds[prefix+'NO2']
        # Copy and update attributes
        attrs = ds[prefix+CopyVar].attrs
        attrs['long_name'] = 'Dry mixing ratio of species {}'.format(fam2use)
        ds[prefix+fam2use].attrs = attrs

    elif fam == 'NOy':
        # Add NOy as defined in GEOS-CF, NOy =
        # no_no2_hno3_hno4_hono_2xn2o5_pan_organicnitrates_aerosolnitrates
        vars2use = GC_var('NOy-all')
        fam2use = 'NOy'
        CopyVar = 'N2O5'
        # 2 N2O5 in NOy, so 2x via template
        ds[prefix+fam2use] = ds[prefix+CopyVar].copy()
        for var in vars2use:
            try:
                ds[prefix+fam2use] = ds[prefix+fam2use] + ds[prefix+var]
            except KeyError:
                pass
        # Copy and update attributes
        attrs = ds[prefix+CopyVar].attrs
        attrs['long_name'] = 'Dry mixing ratio of species {}'.format(fam2use)
        ds[prefix+fam2use].attrs = attrs

    elif fam == 'NOy-gas':
        # Add a variable for gas-phase NOy
        fam2use = 'NOy-gas'
        CopyVar = 'N2O5'
        vars2use = GC_var(fam2use)
        ds[prefix+fam2use] = ds[prefix+CopyVar].copy()
        for var in vars2use:
            try:
                ds[prefix+fam2use] = ds[prefix+fam2use] + ds[prefix+var]
            except KeyError:
                pass
        # Copy and update attributes
        attrs = ds[prefix+CopyVar].attrs
        attrs['long_name'] = 'Dry mixing ratio of species {}'.format(fam2use)
        ds[prefix+fam2use].attrs = attrs

    elif fam == 'NIT-all':
        fam2use = 'NIT-all'
        CopyVar = 'NIT'
        vars2use = GC_var(fam2use)
        ds[prefix+fam2use] = ds[prefix+CopyVar].copy()
        vars2use.pop(vars2use.index(CopyVar))
        # Add dust nitrates if present.
        for var in vars2use:
            try:
                ds[prefix+fam2use] = ds[prefix+fam2use] + ds[prefix+var]
            except KeyError:
                pass
        # Copy and update attributes
        attrs = ds[prefix+CopyVar].attrs
        attrs['long_name'] = 'Dry mixing ratio of species {}'.format(fam2use)
        ds[prefix+fam2use].attrs = attrs

    elif fam == 'SO4-all':
        fam2use = 'SO4-all'
        CopyVar = 'SO4'
        vars2use = GC_var(fam2use)
        ds[prefix+fam2use] = ds[prefix+CopyVar].copy()
        vars2use.pop(vars2use.index(CopyVar))
        # Add dust sulfates if present.
        for var in vars2use:
            try:
                ds[prefix+fam2use] = ds[prefix+fam2use] + ds[prefix+var]
            except KeyError:
                pass
        # Copy and update attributes
        attrs = ds[prefix+CopyVar].attrs
        attrs['long_name'] = 'Dry mixing ratio of species {}'.format(fam2use)
        ds[prefix+fam2use].attrs = attrs

    elif fam == 'SOx':
        fam2use = 'SOx'
        CopyVar = 'SO2'
        vars2use = GC_var(fam2use)
        ds[prefix+fam2use] = ds[prefix+CopyVar].copy()
        vars2use.pop(vars2use.index(CopyVar))
        # Add extra variables if present (check family definition for SOx)
        for var in vars2use:
            try:
                ds[prefix+fam2use] = ds[prefix+fam2use] + ds[prefix+var]
            except KeyError:
                pass
        # Copy and update attributes
        attrs = ds[prefix+CopyVar].attrs
        attrs['long_name'] = 'Dry mixing ratio of species {}'.format(fam2use)
        ds[prefix+fam2use].attrs = attrs

    elif fam == 'Cly':
        fam2use = 'Cly'
        ref_spec = 'Cl'
        vars2use = GC_var(fam2use)
        CopyVar = vars2use[0]
        stoich = spec_stoich(CopyVar, ref_spec=ref_spec)
        vars2use = GC_var(fam2use)
        ds[prefix+fam2use] = ds[prefix+CopyVar].copy() * stoich
        # Add dust nitrates if present.
        for var in vars2use[1:]:
            stoich = spec_stoich(var, ref_spec=ref_spec)
            try:
                ds[prefix+fam2use] = ds[prefix+fam2use] + ds[prefix+var]*stoich
            except KeyError:
                pass
        # Copy and update attributes
        attrs = ds[prefix+CopyVar].attrs
        attrs['long_name'] = 'Dry mixing ratio of species {}'.format(fam2use)
        ds[prefix+fam2use].attrs = attrs

    elif fam == 'Bry':
        fam2use = 'Bry'
        ref_spec = 'Br'
        vars2use = GC_var(fam2use)
        CopyVar = vars2use[0]
        stoich = spec_stoich(CopyVar, ref_spec=ref_spec)
        vars2use = GC_var(fam2use)
        ds[prefix+fam2use] = ds[prefix+CopyVar].copy() * stoich
        # Add dust nitrates if present.
        for var in vars2use[1:]:
            stoich = spec_stoich(var, ref_spec=ref_spec)
            try:
                ds[prefix+fam2use] = ds[prefix+fam2use] + ds[prefix+var]*stoich
            except KeyError:
                pass
        # Copy and update attributes
        attrs = ds[prefix+CopyVar].attrs
        attrs['long_name'] = 'Dry mixing ratio of species {}'.format(fam2use)
        ds[prefix+fam2use].attrs = attrs

    elif fam == 'Iy':
        fam2use = 'Iy'
        ref_spec = 'I'
        vars2use = GC_var(fam2use)
        CopyVar = vars2use[0]
        stoich = spec_stoich(CopyVar, ref_spec=ref_spec)
        vars2use = GC_var(fam2use)
        ds[prefix+fam2use] = ds[prefix+CopyVar].copy() * stoich
        # Add dust nitrates if present.
        for var in vars2use[1:]:
            stoich = spec_stoich(var, ref_spec=ref_spec)
            try:
                ds[prefix+fam2use] = ds[prefix+fam2use] + ds[prefix+var]*stoich
            except KeyError:
                pass
        # Copy and update attributes
        attrs = ds[prefix+CopyVar].attrs
        attrs['long_name'] = 'Dry mixing ratio of species {}'.format(fam2use)
        ds[prefix+fam2use].attrs = attrs

    else:
        print('TODO - setup family and stoich conversion for {}'.format(fam))

    gc.collect()
    return ds


def get_general_stats4run_dict_as_df(run_dict=None, extra_str='', REF1=None,
                                     REF2=None, REF_wd=None, res='4x5',
                                     trop_limit=True,
                                     save2csv=True, prefix='GC_',
                                     extra_burden_specs=[],
                                     extra_surface_specs=[],
                                     GC_version='v12.6.0',
                                     use_time_in_trop=True, rm_strat=True,
                                     dates2use=None, round=3,
                                     use_REF_wd4Met=False,
                                     verbose=False, debug=False):
    """
    Get various stats on a set of runs in a dictionary ({name: location})

    Parameters
    ----------
    run_dict (dict): dicionary of run names and locations
    use_REF_wd4Met (bool): use a reference working directory for shared values?
    REF_wd (str): directory to use to extract shared variables
    REF1 (str): name of (1st) run in dictionary to to % change calculations from
    REF2 (str): name of (2nd) run in dictionary to to % change calculations from
    prefix (str):  string to include as a prefix in saved csv's filename
    extra_str (str):  string to include as a suffx in saved csv's filename
    save2csv (bool): save dataframe as a csv file
    trop_limit (bool): limit analysis to the troposphere?
    extra_burden_specs (list): list of extra species to give trop. burden stats on
    extra_surface_specs (list): list of extra species to give surface conc. stats on
    res (str): resolution of the modul output (e.g. 4x5, 2x2.5, 0.125x0.125)
    round (int): number of decimal places to round dataframe too

    Returns
    -------
    (pd.DataFrame)
    """
    # - Define local variables
    # Mass unit scaling
    mass_scale = 1E3
    mass_unit = 'Tg'
    # Mixing ratio (v/v) scaling?
    ppbv_unit = 'ppbv'
    ppbv_scale = 1E9
    pptv_unit = 'pptv'
    pptv_scale = 1E12
    # Setup lists and check for families (NOy, NIT-all, Iy, Cly, Bry, ... )
    core_specs = ['O3',  'NO', 'NO2']
    prefix = 'SpeciesConc_'
    ALLSp = list(set(core_specs+extra_burden_specs+extra_surface_specs))
    FamilyNames = GC_var('FamilyNames')
    families2use = []
    for FamilyName in FamilyNames:
        if FamilyName in ALLSp:
            families2use += [FamilyName]
    if verbose or debug:
        print(families2use)
    # Core dataframe for storing calculated stats on runs
    df = pd.DataFrame()
    # - Get core data required
    # Get all of the speciesConcs for runs as list of datasets
    dsD = {}
    for key in run_dict.keys():
        ds = GetSpeciesConcDataset(wd=run_dict[key], dates2use=dates2use)
        # Add families to Dataset
        if len(families2use) >= 1:
            for fam in families2use:
                ds = AddChemicalFamily2Dataset(ds, fam=fam, prefix=prefix)
        dsD[key] = ds
    # - Get burdens for core species
    avg_over_time = True  # Note: burdens area averaged overtime
    specs2use = list(set(core_specs+['CO']+extra_burden_specs))
    vars2use = [prefix+i for i in specs2use]
    for key in run_dict.keys():
        # Get StateMet object for 1st of the runs and use this for all runs
        if use_REF_wd4Met:
            # Set working directory for shared variables
            if isinstance(REF_wd, type(None)):
                REF_wd = run_dict[list(run_dict.keys())[0]]
            StateMet = get_StateMet_ds(wd=REF_wd, dates2use=dates2use)
        else:
            StateMet = get_StateMet_ds(wd=run_dict[key], dates2use=dates2use)
        # Average burden over time
        ds = dsD[key]  # .mean(dim='time', keep_attrs=True)
        S = get_Gg_trop_burden(ds, vars2use=vars2use, StateMet=StateMet,
                               use_time_in_trop=use_time_in_trop,
                               avg_over_time=avg_over_time,
                               rm_strat=rm_strat,
                               debug=debug)
        # convert to ref spec equivalent (e.g. N for NO2, C for ACET)
        for spec in specs2use:
            ref_spec = get_ref_spec(spec)
            val = S[prefix+spec]
            S[prefix+spec] = val/species_mass(spec)*species_mass(ref_spec)
        # Upate varnames
        varnames = ['{} burden ({})'.format(i, mass_unit) for i in specs2use]
        S = S.rename(index=dict(zip(list(S.index.values), varnames)))
        # Save the values for run to central DataFrame
        df[key] = S

    # Transpose dataframe
    df = df.T

    # Scale units
    for col_ in df.columns:
        if 'Tg' in col_:
            df.loc[:, col_] = df.loc[:, col_].values/mass_scale
    # Transpose back to variables as index
    df = df.T

    # - Add Ozone production and loss...

    # - Add lightning source if in HEMCO output
#     var2use = 'EmisNO_Lightning'
#     varName = 'Lightning (Tg N/yr)'
#     try:
#     # TODO -  add check for HEMCO NetCDF files in the output folder
# #    if True:
#         dsH = {}
#         for key in run_dict.keys():
#             dsH[key] = get_HEMCO_diags_as_ds(wd=run_dict[key])
#             ds = ds[key]
#             val = (ds[var2use].mean(dim='time').sum(dim='lev') * ds['AREA'] )
#             val2 = val.values.sum() * 60 * 60 * 24 * 365 # => /yr
#             df.loc[key,varName] = val2*1E3/1E12
#     except KeyError:
#         pass
#         df.loc[key,varName] = np.nan

    # - Surface concentrations
    specs2use = list(set(core_specs+['N2O5']+extra_surface_specs))
    prefix = 'SpeciesConc_'
    # Loop by run and get stats
    for key in run_dict.keys():
        if debug:
            print(key)
        ds = dsD[key].copy()
        # Select surface and average over time
        ds = ds.mean(dim='time')
        ds = ds.isel(lev=ds.lev == ds.lev[0])
        for spec in specs2use:
            # Get units and scaling
            units, scale = tra_unit(spec, scale=True)
            # Surface ozone
            varname = '{} surface ({})'.format(spec, units)
            var = prefix+spec
            # Save values on a per species basis to series
            val = get_avg_2D_conc_of_X_weighted_by_Y(ds, Xvar=var, Yvar='AREA')
            # Save calculated values to dataframe
            df.loc[varname, key] = val

    # Transpose dataframe
    df = df.T
    # Scale units
    for col_ in df.columns:
        if 'ppb' in col_:
            df.loc[:, col_] = df.loc[:, col_].values*ppbv_scale
        if 'ppt' in col_:
            df.loc[:, col_] = df.loc[:, col_].values*pptv_scale

    # - OH concentrations

    # - Add methane lifetime calculation

    # - Processing and save?
    # Calculate % change from base case for each variable
    if not isinstance(REF1, type(None)):
        for col_ in df.columns:
            pcent_var = col_+' (% vs. {})'.format(REF1)
            df[pcent_var] = (df[col_]-df[col_][REF1]) / df[col_][REF1] * 100
    if not isinstance(REF2, type(None)):
        for col_ in df.columns:
            pcent_var = col_+' (% vs. {})'.format(REF2)
            df[pcent_var] = (df[col_]-df[col_][REF2]) / df[col_][REF2] * 100

    # Transpose back to variables as index
    df = df.T
    # Re-order columns
    df = df.reindex(sorted(df.columns), axis=1)
    # Reorder index
    df = df.T.reindex(sorted(df.T.columns), axis=1).T
    # Now round the numbers
    df = df.round(round)
    # Save csv to disk
    if save2csv:
        csv_filename = '{}_summary_statistics{}.csv'.format(prefix, extra_str)
        df.to_csv(csv_filename)
    # Return the DataFrame too
    return df


def add_molec_den2ds(ds, MolecVar='Met_MOLCES', AirDenVar='Met_AIRDEN'):
    """
    Add molecules/cm3 to xr.dataset (must contain AirDenVar)
    """
    # Calculate number of molecules
    try:
        ds[MolecVar]
    except KeyError:
        RMM_air = constants('RMM_air') / 1E3  # g/mol => kg/mol
        ds[MolecVar] = ds[AirDenVar].copy() / 1E6  # kg/m3 => kg/cm3
        values = ds[MolecVar].values
        values = values / RMM_air  # kg/cm3 / kg/mol
        values = values * constants('AVG')  # mol/cm3 => molecules/cm3
        ds[MolecVar].values = values
    return ds


def add_Cly_Bry_Iy_2ds(ds=None, prefix='SpeciesConc_',
                       add_ind_specs2ds=False, verbose=False):
    """
    Add all Xy (X=Cl, Br, I) variables to xr.dataset (e.g. SpeciesConc* )
    """
    fams = 'Cly', 'Bry', 'Iy'
    for fam in fams:
        ds = add_Xy_2ds(ds, var2add=fam, prefix=prefix,
                        add_ind_specs2ds=add_ind_specs2ds,
                        verbose=verbose)
    return ds


def add_Xy_2ds(ds=None, var2add='Cly', prefix='SpeciesConc_',
               add_ind_specs2ds=False,  verbose=False):
    """
    Add an Xy (X=Cl, Br, I) to xr.dataset (e.g. SpeciesConc* )
    """
    # Get the reference species
    ref_spec = get_ref_spec(var2add)
    # Get Xy species
    specs2use = GC_var(var2add)
    # Remove HCl from list of species (if present) as it is a reservoir species
    if var2add == 'Cly':
        try:
            idx = specs2use.index('HCl')
            specs2use.pop(idx)
            print('WARNING: removed HCl from Cly definition')
        except ValueError:
            pass
    if verbose:
        Pstr = "Using species for '{}' family (ref_spec: {}): {}"
        print(Pstr.format(var2add, ref_spec, specs2use))
    # Setup Xy variable as template of 1st Xy species
    spec2use = specs2use[0]
    var2use = '{}{}'.format(prefix, spec2use)
    stioch = spec_stoich(spec2use, ref_spec=ref_spec)
    ds[var2add] = ds[var2use].copy() * stioch
    # Also save the individual species in reference species terms?
    if add_ind_specs2ds:
        Var2Save = '{}-in-{}-units'.format(spec2use, ref_spec)
        ds[Var2Save] = ds[var2use].copy() * stioch
    # Add rest of Xy species and scale to stoichiometry
    for spec2use in specs2use[1:]:
        var2use = '{}{}'.format(prefix, spec2use)
        values = ds[var2use].values * stioch
        ds[var2add].values = ds[var2add].values + values
        # Also save the individual species in reference species terms?
        if add_ind_specs2ds:
            Var2Save = '{}-in-{}-units'.format(spec2use, ref_spec)
            ds[Var2Save] = ds[var2use].copy() * stioch
    return ds


def get_specieslist_from_input_geos(folder=None, filename='input.geos'):
    """
    Extract the species list from the input.geos file
    """
    line2start_read = '%%% ADVECTED SPECIES MENU %%%'
    line2end_read = '------------------------+--------------------------------'
    species_lines = []
    with open(folder+filename, 'r') as file:
        save_line = False
        for line in file:
            if line2end_read in line:
                save_line = False
            if save_line:
                species_lines += [line]
            if line2start_read in line:
                save_line = True

    species = [i.split(':')[-1].strip() for i in species_lines]
    return species
