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
import glob
import pandas as pd
import xarray as xr
import re
from netCDF4 import Dataset
import logging
import wget
from bs4 import BeautifulSoup
import requests
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
#from .planeflight import *
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
            assert correct_type, "'date' variable must be a datetime.datetime object"
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
    fcst (str): forecast string to use (if not provided, then evaluated from 'dt')
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
        fcst = '{}{:0>2}{:0>2}T{:0>2}0000'.format( dt.year, dt.month, dt.day, dt.hour )
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
        URL = urlstr.format(type_root, tau, stream, level, region, fcst, field )
        # Request URL and then get images from location
        r = requests.get(URL)
        html = r.text
        soup = BeautifulSoup(html, 'lxml')
        images = soup.findAll('img')
        # check that the URL was correct
        assert len(images) > 0, 'WARNING: No images found @ URL ({})'.format(URL)
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
                name = f.format(prefix, ptype, stream, fcst, field, region, level, tau )
                wget.download(image_URL, folder+name)


def get_GEOS5_datagram_plots( dt=None, stream='G5FPFC', folder=None,
                              prefix='ARNA', plts2get=['du_16.7_-23.0'],
                              verbose=False ):
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
        fstr = '{}_{}_{}_datagram_{}.png'.format( prefix, stream, date_str, plt2get )
        filename = fstr.format(date_str, stream, )
        wget.download(url, folder+filename)

