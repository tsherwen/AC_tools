#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Functions for use with the NASA models in the GEOS family

Use help(<name of function>) to get details on a particular function.

"""
# - Required modules:
# compatibility with both python 2 and 3
from __future__ import print_function
# I/O / Low level
import os
import sys
import re
import glob
import logging
import wget
import requests
import pandas as pd
import xarray as xr
import xesmf as xe
from netCDF4 import Dataset
from bs4 import BeautifulSoup
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
from time import gmtime, strftime
# The below imports need to be updated,
# imports should be specific and in individual functions
# import tms modules with shared functions
from .core import *
from .generic import *
from .AC_time import *
from .planeflight import prt_PlaneFlight_files_v12_plus
#from .variables import *
#from .bpch2netCDF import convert_to_netCDF


def get_GEOSCF_as_ds_via_OPeNDAP(collection='chm_inst_1hr_g1440x721_p23',
                                 mode='fcast', date=None):
    """
    Get the GEOS Composition Forecast (GEOS-CF) as a xr.Dataset (using OPeNDAP)

    Parameters
    ----------
    mode (str): retrieve the forecast (fcast) or assimilated fields (assim)
    date (datetime.datetime): date to retrieve forecast from or assimilation for
    collection (str): data collection to access (e.g. chm_inst_1hr_g1440x721_p23)

    Returns
    -------
    (xr.dataset)

    NOTES
    ---
     - default is to get the latest forecast for chemistry
     - See documentation for details:
    https://gmao.gsfc.nasa.gov/weather_prediction/GEOS-CF/
     - Collections include:
    chm_inst_1hr_g1440x721_p23
    chm_tavg_1hr_g1440x721_v1
    htf_inst_15mn_g1440x721_x1
    met_inst_1hr_g1440x721_p23
    met_tavg_1hr_g1440x721_x1
    xgc_tavg_1hr_g1440x721_x1
    """
    # Root OPeNDAP directory
    root_url = 'https://opendap.nccs.nasa.gov/dods/gmao/geos-cf/{}/'.format(
        mode)
    # Make up the complete URL for a forecast or assimilation field
    if mode == 'fcast':
        # Which date to use?
        if isinstance(date, type(None)):
            # Use the lastest file (default)
            URL = '{}/{}.latest'.format(root_url, collection)
        else:
            # Use a file specified in arguments
            correct_type = type(date) == datetime.datetime
            AssStr = "'date' variable must be a datetime.datetime object"
            assert correct_type, AssStr
            # Use the lastest file (default)
            dstr = date.strftime(format='%Y%m%d')
            URL = '{}/{}/{}.{}_12z'.format(root_url,
                                           collection, collection, dstr)
    elif mode == 'assim':
        # Just retrieve an OPeNDAP pointer to the entire dataset for now
        URL = '{}/{}'.format(root_url, collection)
    else:
        print("WARNING: GEOS-CF mode provided ('{}') not known".format(mode))
        sys.exit()
    # Open the dataset via OPeNDAP and return
    ds = xr.open_dataset(URL)
    return ds


def get_GEOS5_as_ds_via_OPeNDAP(collection='inst3_3d_aer_Nv',
                                fcast_start_hour=12,
                                mode='seamless', dt=None):
    """
    Get the GEOS-5 model product (GEOS-5) as a xr.Dataset (using OPeNDAP)

    Parameters
    ----------
    mode (str): retrieve the forecast (fcast) or assimilated fields (assim) or
                both (seemless)
    dt (datetime.datetime): date to retrieve forecast from or assimilation for
    collection (str): data collection to access
                      (e.g. chm_inst_1hr_g1440x721_p23)
    fcast_start_hour (int): hour the forcast started on a given day

    Returns
    -------
    (xr.dataset)

    NOTES
    ---
     - default is to get the latest forecast for chemistry (via seamless route)
     - See documentation for details:    https://geos5.org/wiki/index.php?title=GEOS-5_Earth_System_Modeling_and_Data_Assimilation
     - Collections include:
     - The forecast set going at different times are for different length.
     00 - ~10 days
     06 - ~1.5 days
     12 - ~5 days
     18 - ~1.5 days
     - the 5 day forecast for a given day is selected as default
       (fcast_start_hour)

    """
    # Root OPeNDAP directory
    root_url = 'https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/{}/'
    root_url = root_url.format(mode)
    # Make up the complete URL for a forecast or assimilation field
    if (mode == 'fcast') or (mode == 'seamless'):
        # Which date to use?
        if isinstance(dt, type(None)):
            # Use the lastest file (default)
            URL = '{}/{}.latest'.format(root_url, collection)
        else:
            # Use a file specified in arguments
            correct_type = type(dt) == datetime.datetime
            AssStr = "'date' variable must be a datetime.datetime object"
            assert correct_type, AssStr
            # Use the 'lastest' file (default)
            # NOTE: lastest 12 file is used to match up with GEOS-CF
            # TODO: update this. this will not give enough data
            dstr = dt.strftime(format='%Y%m%d')
            URL = '{}/{}/{}.{}_{:0>2}'
            URL = URL.format(root_url, collection, collection, dstr,
                             fcast_start_hour)
    elif mode == 'assim':
        # Just retrieve an OPeNDAP pointer to the entire dataset for now
        URL = '{}/{}'.format(root_url, collection)
    else:
        print("WARNING: GEOS-5 mode provided ('{}') not known".format(mode))
        sys.exit()
    # Open the dataset via OPeNDAP and return
    ds = xr.open_dataset(URL)
    return ds


def get_GEOS5_online_diagnostic_plots(dt=None, ptype='wxmaps',
                                      region='atlantic', field='precip',
                                      fcst=None, stream='G5FPFC',
                                      taus=[0], level='0',
                                      folder='./', prefix='ARNA',
                                      verbose=False):
    """
    Get the static meteorological plots for GEOS5 from NASA's NCCS/Fluid

    Parameters
    -------
    dt (datetime): datetime for start forecast to access
    ptype (str): type of plot (e.g. wxmaps or chem2d)
    field (str): data field to access
    region (str): Plotted region (e.g. atlantic)
    fcst (str): forecast string to use (else evaluated from 'dt')
    taus (list of int): list of integer timesteps to save from
    level (list): level to extract (only = 0 currently setup)
    stream (str): which model data stream to use
    folder (str): where to save plots
    prefix (Str): prefix to append to beginning of all filenames

    Returns
    -------
    (None)

    Notes
    -----
     - Please refer to NASA's fluid website for further infomation
    https://fluid.nccs.nasa.gov/
    """
    # if no day given, use previous day
    if isinstance(dt, type(None)):
        # Just use yesterday for Now
        TNow = time2datetime( [gmtime()] )[0]
        dt = datetime.datetime(TNow.year, TNow.month, TNow.day-1)

    # Which forecast to use
    if isinstance(fcst, type(None)):
        fcst_str = '{}{:0>2}{:0>2}T{:0>2}0000'
        fcst = fcst_str.format( dt.year, dt.month, dt.day, dt.hour )
    # What is the website location for the data?
    site = 'https://fluid.nccs.nasa.gov/'
    if ptype == 'wxmaps':
        type_root = '{}/{}'.format( site, ptype)
    else:
        type_root = '{}/wxmaps/{}'.format( site, ptype)

    # What is the outline file structure?
    urlstr = '{}/?one_click=1&tau={:0>3}&stream={}&level={}&region={}&fcst={}&field={}'
    # Loop by requested time from forecast start
    for tau in taus:
        # Compile the URL
        URL = urlstr.format(type_root, tau, stream, level, region, fcst, field)
        # Request URL and then get images from location
        r = requests.get(URL)
        html = r.text
        soup = BeautifulSoup(html, "html.parser")
        images = soup.findAll('img')
        # check that the URL was correct
        AssStr = 'WARNING: No images found @ URL ({})'
        assert len(images) > 0, AssStr.format(URL)
        # Get the one image
        for img in images:
            src = img.get('src')
            title = img.get('title')
            if field in src:
                image_URL = site+src
                if verbose:
                    print(img, src, title, image_URL)
                # Download using Python wget
                f = '{}_{}_{}_fcast_{}_{}_{}_{:0>2}_{:0>3}.png'
                name = f.format(prefix, ptype, stream, fcst, field, region,
                                level, tau )
                wget.download(image_URL, folder+name)


def get_GEOS5_datagram_plots( dt=None, stream='G5FPFC', folder=None,
                              prefix='ARNA', plts2get=['du_16.7_-23.0'],
                              verbose=False, debug=False ):
    """
    Get GEOS5 datagram plots and save these locally

    Parameters
    -------
    dt (datetime): datetime for start forecast to access
    stream (str): which model data stream to use
    folder (str): where to save plots
    prefix (Str): prefix to append to beginning of all filenames
    plts2get (list): list of plots to retrieve (example is dust at Cape Verde)

    Returns
    -------
    (None)

    Notes
    -----
     - Please refer to NASA's fluid website for further infomation
    https://fluid.nccs.nasa.gov/
    """
    # if no day given, use previous day
    if isinstance(dt, type(None)):
        # Just use yesterday for Now
        TNow = time2datetime( [gmtime()] )[0]
        dt = datetime.datetime(TNow.year, TNow.month, TNow.day-1)
    date_str = '{}_{:0>2}_{:0>2}'.format(dt.year, dt.month, dt.day)
    # get plots in list
    site = 'https://fluid.nccs.nasa.gov/'
    gram_root = site+'/gram/static/plots/'
    # loop and retrieve the files
    for plt2get in plts2get:
        url = '{}{}.png'.format( gram_root, plt2get )
        # Download using wget through python
        fstr = '{}_{}_{}_datagram_{}.png'
        fstr = fstr.format( prefix, stream, date_str, plt2get )
        filename = fstr.format(date_str, stream, )
        if debug:
            pstr = 'Getting {} and saving here: {}'
            print(pstr.format(url, folder+filename))
        wget.download(url, folder+filename)


def get_GEOSCF_vertical_levels(print_equivalents=False, native_levels=False):
    """
    get a dictionary of GEOS Composition Forecast (GEOS-CF) vertical levels
    """
    # Get a list of the pressure levels in GEOS-CF
    if native_levels:
        HPa_l  = [
        0.01, 0.02, 0.0327, 0.0476, 0.066, 0.0893, 0.1197, 0.1595, 0.2113, 0.2785, 0.365, 0.4758, 0.6168, 0.7951, 1.0194, 1.3005, 1.6508, 2.085, 2.6202, 3.2764, 4.0766, 5.0468, 6.2168, 7.6198, 9.2929, 11.2769, 13.6434, 16.4571, 19.7916, 23.7304, 28.3678, 33.81, 40.1754, 47.6439, 56.3879, 66.6034, 78.5123, 92.3657, 108.663, 127.837, 150.393, 176.93, 208.152, 244.875, 288.083, 337.5, 375.0, 412.5, 450.0, 487.5, 525.0, 562.5, 600.0, 637.5, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 820.0, 835.0, 850.0, 865.0, 880.0, 895.0, 910.0, 925.0, 940.0, 955.0, 970.0, 985.0
        ]
    else:
        HPa_l = [
            1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400,
            350, 300, 250, 200, 150, 100, 50
        ]

    # Get altitudes in km, then convert to metres
    Alt = hPa_to_Km(HPa_l)
    Alt = np.array(Alt) / 1E3
    # Print out a summary
    if print_equivalents:
        # For just HPa and km
        pstr = 'A pressure {:>4}HPa of is equiv to {:>4,.3f}km'
        Alt_dict = dict(zip(HPa_l, Alt*1E3))
        for HPa in HPa_l:
            print(pstr.format(HPa, Alt_dict[HPa]))
        # Also for kft
        pstr = 'A press. {:>4} HPa of is equiv to {:>4,.3f} km ({:>4,.3f} kft)'
        for HPa in HPa_l:
            print(pstr.format(HPa, Alt_dict[HPa], Alt_dict[HPa]/304.8))
    return HPa_l


def extract_GEOSCF4FAAM_flight(folder=None, flight_ID='C216', folder4csv=None,
                               PressVar="PS_RVSM",
                               LonVar='LON_GIN',
                               LatVar='LAT_GIN', TimeVar='Time',
                               testing_mode=True, csv_suffix='',
                               inc_ds_vars_in_csv=False):
    """
    Extract the GEOS-CF model for a given FAAM BAe146 flight
    """
    # Retrieve FAAM BAe146 Core NetCDF files
    filename = 'core_faam_*_{}_1hz.nc'.format(flight_ID.lower())
    file2use = glob.glob(folder+filename)
    if len(file2use) > 1:
        print( 'WARNING: more that one file found! (so using latest file)' )
        print(file2use)
    ds = xr.open_dataset( file2use[0] )
    # Only select the variable of intereest and drop where these are NaNs
    df = ds[ [PressVar, LatVar, LonVar, TimeVar] ].to_dataframe()
    df = df.dropna()
    # If doing a test, then just extract the first 150 points of flight track
    if testing_mode:
        df = df.head(150)
    # Get variables from chemistry collection
    collection = 'chm_inst_1hr_g1440x721_p23'
    df1 = extract_GEOSCF_assim4flight(df, PressVar=PressVar, LonVar=LonVar,
                                      LatVar=LatVar, TimeVar=TimeVar,
                                      collection=collection)
    # Get variables from meterology collection
    collection = 'met_inst_1hr_g1440x721_p23'
    df2 = extract_GEOSCF_assim4flight(df, PressVar=PressVar, LonVar=LonVar,
                                      LatVar=LatVar, TimeVar=TimeVar,
                                      collection=collection)
    # Combine dataframes and remove duplicate columns
    dfs = [df1, df2]
    df = pd.concat(dfs, axis=1)
    if inc_ds_vars_in_csv:
        FAAM_df = ds.to_dataframe()
        df = pd.concat([df, FAAM_df], axis=1)
        duplicates = [i for i in FAAM_df.columns if i in df.columns ]
        if len(duplicates)>1:
            print('WARNING: duplicates in FAAM and GEOS-CF dataframe headers!')
    df = df[ list(set(df.columns)) ]
    # Save dateframe to a csv file
    if isinstance(folder4csv, type(None)):
        folder4csv = './'
    filename = 'FAAM_flightpath_extracted_from_GEOSCF_4_{}{}'
    filename = filename.format(flight_ID, csv_suffix)
    filename = rm_spaces_and_chars_from_str(filename)
    df.to_csv(folder4csv+filename+'.csv')



def extract_GEOSCF_assim4flight(df=None, ds=None,
                                LatVar='lat', vars2extract=None,
                                collection='met_inst_1hr_g1440x721_p23',
                                LonVar='lon', PressVar='hPa',
                                mode='assim', TimeVar='time',
                                dsAltVar='lev', dsTimeVar='time',
                                dsLonVar='lon', dsLatVar='lat',
                                testing_mode=False, inc_attrs=True,
                                TEMP_nc_name=None,
                                debug=False):
    """
    Extract GEOS-CF collection for flightpath locations in dataframe
    """
    # Get the start and end date of flight (with a 1/4 day buffer)
    sdate = add_days(df.index.min(), -0.25)
    edate = add_days(df.index.max(), 0.25)
    # Retrieve the 3D fields from OPenDAP for given collection
    ds = get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode)
    # Extract all of the data variables unless a specific list is provided
    if isinstance(vars2extract, type(None)):
        vars2extract = list(ds.data_vars)
    # Restrict the dataset to the day(s) of the flight
    time_bool = [((i>=sdate) & (i<=edate)) for i in ds.time.values]
    ds = ds.isel(time=time_bool)
    # Reduce the dataset size to the spatial locations of the flight (+ buffer)
    spatial_buffer = 2 # degrees lat / lon
    lat_min = df[LatVar].values.min() - spatial_buffer
    lat_max = df[LatVar].values.max() + spatial_buffer
    lat_bool = [((i>=lat_min) & (i<=lat_max)) for i in ds[dsLatVar].values]
    ds = ds.isel(lat=lat_bool)
    lon_min = df[LonVar].values.min() - spatial_buffer
    lon_max = df[LonVar].values.max() + spatial_buffer
    lon_bool = [((i>=lon_min) & (i<=lon_max)) for i in ds[dsLonVar].values]
    ds = ds.isel(lon=lon_bool)
    # Get a list of the levels that the data is interpolated onto
    HPa_l = get_GEOSCF_vertical_levels(native_levels=False)
    # Get nearest indexes in 4D data from locations in dataframe
    idx_dict = calc_4D_idx_in_ds(ds_hPa=HPa_l, ds=ds, df=df,
                                 TimeVar=TimeVar,
                                 AltVar=PressVar,
                                 LatVar=LatVar, LonVar=LonVar,
                                 dsLonVar=dsLonVar, dsLatVar=dsLatVar,
                                 dsTimeVar=dsTimeVar,
                                 )
    # Make a dictionaries to convert between ds and df variable names
    df2ds_dict = {
    LonVar:dsLonVar, LatVar:dsLatVar, TimeVar:dsTimeVar, PressVar:dsAltVar,
    }
    df2ds_dict_r = {v: k for k, v in list(df2ds_dict.items())}
    # Save subset of dataset locally and then reload
    if isinstance(TEMP_nc_name, type(None)):
         TEMP_nc_name = 'TEMP_NetCDF_{}.nc'.format(collection)
    ds = save_ds2disk_then_reload(ds, savename=TEMP_nc_name)
    # Create a data frame for values
    dfN = pd.DataFrame()
    # Extraction of data points in a bulk manner
    for nval, var in enumerate( vars2extract ):
        print(var)
        # Now extract values
        dims2use = list(ds[var].coords)
        idx_list = [idx_dict[df2ds_dict_r[i]] for i in dims2use]
        vals = ds[var].values[tuple(idx_list)]
        dfN[vars2extract[nval]] = vals
    # Also save model time variable to dataframe
    time_idx = idx_dict[df2ds_dict_r[dsTimeVar]]
    dfN['model-time'] = ds[dsTimeVar].values[time_idx]
    # Add model lat, lon and pressure to dataframe
    lon_idx = idx_dict[df2ds_dict_r[dsLonVar]]
    dfN['model-lon'] = ds[dsLonVar].values[lon_idx]
    lat_idx = idx_dict[df2ds_dict_r[dsLatVar]]
    dfN['model-lat'] = ds[dsLatVar].values[lat_idx]
    alt_idx = idx_dict[df2ds_dict_r[dsAltVar]]
    dfN['model-alt'] = ds[dsAltVar].values[alt_idx]
    # Include variable attributes from original dataset
    if inc_attrs:
        for col in dfN.columns:
            try:
                attrs = ds[col].attrs.copy()
                dfN[col].attrs = attrs
            except KeyError:
                pass
    # Save the datetime as a column too
    dfN['Datetime'] = df.index.values
    # Update the model datetime to be in datetime units
    dfN['model-time'] = pd.to_datetime(dfN['model-time'].values)
    # Remove the temporary NetCDF file (of OPenDap dataset subset) from disk
    rm_file(folder='./', filename=TEMP_nc_name)
    # Return the extracted dataframe of flighttrack points
    return dfN


def mk_planeflight_input4FAAM_flight(folder=None, flight_ID='C216',
                                     folder4csv=None,
                                     PressVar="PS_RVSM",
                                     LonVar='LON_GIN',
                                     LatVar='LAT_GIN', TimeVar='Time',
                                     LocVar='TYPE', LocName='B-146',
                                     DateVar='datetime',
                                     testing_mode=True, csv_suffix='',
                                     num_tracers=203,
                                     Username='Tomas Sherwen',
                                     slist=None,
                                     Extra_spacings=False
                                     ):
    """
    Extract the GEOS-CF model for a given FAAM BAe146 flight
    """
    # Retrieve FAAM BAe146 Core NetCDF files
    filename = 'core_faam_*_{}_1hz.nc'.format(flight_ID.lower())
    file2use = glob.glob(folder+filename)
    if len(file2use) > 1:
        print( 'WARNING: more that one file found! (so using latest file)' )
        print(file2use)
    ds = xr.open_dataset( file2use[0] )
    # Only select the variable of intereest and drop where these are NaNs
    df = ds[ [PressVar, LatVar, LonVar, TimeVar] ].to_dataframe()
    df = df.dropna()
    # Add a location name (as Type)
    df[LocVar] = LocName
    # remove the index name and add index values to a column
    df.index.name = None
    try:
        df[DateVar]
    except KeyError:
        df['datetime'] = df.index.values
    # If doing a test, then just extract the first 150 points of flight track
    if testing_mode:
        df = df.head(150)
    # Call planeflight maker...
    prt_PlaneFlight_files_v12_plus(df=df, slist=slist,
                                   Extra_spacings=Extra_spacings,
                                   LON_var=LonVar, LAT_var=LatVar,
                                   PRESS_var=PressVar, loc_var=LocVar,
                                   num_tracers=num_tracers,
                                   Username=Username,)


def regrid_restart_file4flexgrid(dsA, OutFile=None, lons=None,
                                lats=None, res='1x1', folder='./',
                                vars2regrid=None, rm_regridder=True,
                                save2netcdf=True,
                                debug=False):
    """
    Regrid a GEOS-Chem restart file to a given resolution using xEMSF

    Parameters
    -------
    dsA (xr.dataset): dataset to regrid
    OutFile (str): name to save regretted netCDF to
    folder (str): directory address to save file to
    res (str): resolution to regrid to (TODO: add this functionality)
    vars2regrid (list): list of variables to regrid (or all data_vars used)
    rm_regridder (bool): remove the regridded file afterwards?
    save2netcdf (Boolean): Save the regridded file as a NetCDF?
    debug (bool): print debug statements to screen

    Returns
    -------
    (xr.dataset)

    Notes
    -----
     - Important note: Extra dimensions must be on the left, i.e. (time, lev, lat, lon) is correct but (lat, lon, time, lev) would not work. Most data sets should have (lat, lon) on the right (being the fastest changing dimension in the memory). If not, use DataArray.transpose or numpy.transpose to preprocess the data.
    """
    # TODO - setup using AC_Tools to import standard resolutions?
    # Create a dataset to re-grid into
    ds_out = xr.Dataset(dsA.coords)
    # Replace the lat and lon coords
    del ds_out['lat']
    ds_out = ds_out.assign_coords({'lat':lats})
    del ds_out['lon']
    ds_out = ds_out.assign_coords({'lon':lons})
    # Create a regidder (to be reused )
    regridder = xe.Regridder(dsA, ds_out, 'bilinear', reuse_weights=True)
    # Loop and regrid variables
    ds_l = []
    if isinstance(vars2regrid, type(None)):
        vars2regrid = dsA.data_vars
    # Only regrid variables wit lon and lat in them
    vars2leave = [i for i in vars2regrid if 'lat' not in dsA[i].coords.keys()]
    vars2regrid = [i for i in vars2regrid if 'lat' in dsA[i].coords.keys() ]
    for var2use in vars2regrid:
        if debug:
            print(var2use)
        # Create a dataset to re-grid into
        ds_out = xr.Dataset(dsA.coords)
        # Replace the lat and lon coords with the ones to regrid to
        del ds_out['lat']
        ds_out = ds_out.assign_coords({'lat':lats})
        del ds_out['lon']
        ds_out = ds_out.assign_coords({'lon':lons})
        # Get a DataArray
        dr = dsA[var2use]
        # Build regridder
        dr_out = regridder(dr)
        # Exactly the same as input?
        if ('time' in dsA[var2use].coords):
            xr.testing.assert_identical(dr_out['time'], dsA['time'])
        else:
            pstr = 'WARNING: time variable not coord. for var. - not checked'
            if debug:
                print(pstr)
        # Save variable
        ds_l += [dr_out]
    # Combine variables
    ds = xr.Dataset()
    for n, var2use in enumerate(vars2regrid):
        ds[var2use] = ds_l[n]
    # Transfer attributes
    for coord in ds.coords:
        ds[coord].attrs = dsA[coord].attrs.copy()
    for var2use in ds.data_vars:
        ds[var2use].attrs = dsA[var2use].attrs.copy()
    # Add non re-regridded variables into returned dataset
    for var2use in vars2leave:
        ds[var2use] = dsA[var2use]
    # Clean up
    if rm_regridder:
        regridder.clean_weight_file()
    # Save the file (and return)
    if save2netcdf:
        if isinstance(OutFile, type(None)):
            OutFile = '{}.out.nc'.format(InFile)
        ds.to_netcdf(os.path.join(folder, OutFile))
    return ds
