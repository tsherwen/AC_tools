#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Core functions used to be used by all function levels of
GEOS-Chem/Data analysis in AC_Tools

NOTES:
 - user specific functions need to be migrated to a seperate fuction
"""
# - Required modules:
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import platform
import sys
import logging
import os
import inspect
from math import log10, floor
import math


def get_dir(input, loc='earth0'):
    """
    Retrieves directories within structure on a given platform
        ( e.g. York computer (earth0), Mac, Old cluster (atmosviz1... etc)  )

    Notes
    -----
     - This function is not general enough to be transferrable.
     - Update to use $USER flag.
    """
    import getpass
    import platform
    host = platform.node()
    user = getpass.getuser()
    tms_users = ['Tomas', 'tomassherwen', 'ts551']

    # Mac setup
    if (host == 'tomasmbp13.york.ac.uk') or ('Tomas-13-MBP'):
        home = '/Users/{}/'.format(user)
        if user in tms_users:
            d = {
                'rwd': home+'PhD/Data/MUTD_iGEOS-Chem_output/',
                'dwd':  home+'PhD/Data/',
                'npwd': home+'PhD/Data/np_arrs/',
                'tpwd': home+'GITHub/PhD_progs/',
                'ppwd': home+'Pictures/'
            }
            d = d[input]
        else:
            d = home

    # Earth0 setup
    if host == 'earth0':
        home = '/work/home/{}/'.format(user)
        if user in tms_users:
            d = {
                'rwd': home + 'data/all_model_simulations/iodine_runs/',
                'dwd': home + 'data/',
                'fwd': home + 'labbook/Python_progs/d_fast-J_JX/data/',
                'lwd': home + 'labbook/',
                'npwd': home + 'data/np_arrs/',
                'tpwd': home + 'labbook/Python_progs/',
                'ppwd': home + 'labbook/plots_images/'
            }
            d = d[input]
        else:
            d = home

    # Viking setup
    if 'viking' in host:
        home = '/users/{}/scratch/'.format(user)
        if user in tms_users:
            d = {
                'rwd': home + '/GC/rundirs/',
                'dwd': home + 'data/',
                'fwd': home + 'labbook/Python_progs/d_fast-J_JX/data/',
                'lwd': home + 'labbook/',
                'npwd': home + 'data/np_arrs/',
                'tpwd': home + 'labbook/Python_progs/',
                'ppwd': home + 'labbook/plots_images/'
            }
            d = d[input]
        else:
            d = home

    # Atmosviz1 setup
    if host == 'atmosviz1':
        home = '/home/{}/'.format(user)
        if user in tms_users:
            d = {
                'rwd': home + 'data/model/',
                'dwd':  home + 'data/',
                'lwd':  home + 'labbook/',
                'npwd': home + 'data/np_arrs/',
                'tpwd': home + 'labbook/PhD_progs/'
            }
            d = d[input]
        else:
            d = home

    # YARCC setup
    YARCC_login_nodes = ['login{}.york.ac.uk'.format(i) for i in range(1, 4)]
    if host in YARCC_login_nodes:
        home = '/shared/earth_home/{}/'.format(user)  # just use mounted drive
        if user in tms_users:
            d = {
                'rwd': home + 'data/all_model_simulations/iodine_runs/',
                'dwd': home + 'data/',
                'fwd': home + 'labbook/Python_progs/d_fast-J_JX/data/',
                'lwd': home + 'labbook/',
                'npwd': home + 'data/np_arrs/',
                'tpwd': home + 'labbook/Python_progs/',
                'ppwd': home + 'labbook/plots_images/'
            }
            d = d[input]
        else:
            d = home

    return d


def get_gc_lat(lat, res='4x5', wd=None, filename='ctm.nc', debug=False):
    """
    Get index of lat for given resolution

    Parameters
    ----------
    lat (float): latitude to convert
    wd (str): the directory to search for file in
    filename (Str): name of NetCDF file (e.g. ctm.nc or ts_ctm.nc)
    res (str): the resolution if wd not given (e.g. '4x5' )
    debug (bool): legacy debug option, replaced by python logging

    Returns
    -------
    (float)
    """
    NIU, lat_c, NIU = get_latlonalt4res(res=res, wd=wd, filename=filename)
    del NIU
    return find_nearest_value(lat_c, lat)


def get_gc_lon(lon, res='4x5', wd=None, filename='ctm.nc', debug=False):
    """
    Get index of lon for given resolution

    Parameters
    ----------
    lon (float): longitude to convert
    wd (str): the directory to search for file in
    filename (Str): name of NetCDF file (e.g. ctm.nc or ts_ctm.nc)
    res (str): the resolution if wd not given (e.g. '4x5' )
    debug (bool): legacy debug option, replaced by python logging

    Returns
    -------
    (float)
    """
    lon_c, NIU, NIU = get_latlonalt4res(res=res, wd=wd, filename=filename)
    del NIU
    return find_nearest_value(lon_c, lon)


def get_dims4res(res=None, r_dims=False, invert=True, trop_limit=False,
                 just2D=False, full_vert_grid=False, add_n_time_dims=None,
                 debug=False):
    """
    Get dimension of GEOS-Chem output for given resolution

    Parameters
    ----------
    invert (bool): invert dictionary keys and values
    trop_limit (bool): limit 4D arrays to troposphere
    r_dims (bool): return dicionary of dimensions
    just2D (bool): just return horizontal dimensions
    debug (bool): legacy debug option, replaced by python logging
    add_n_time_dims (int): number of time dimensions to add

    Returns
    -------
    (tuple)
    """
    # Dictionary of max dimensions of standard GEOS-Chem output
    dims = {
        '4x5': (72, 46, 47),
        '2x2.5': (144, 91, 47),
        '1x1': (360, 181, 47),
        '0.5x0.5': (720, 361, 47),
        '0.5x0.666': (121, 81, 47),
        # tms - update to be '0.25x0.3125_EU' for consistancy?
        '0.25x0.3125': (177, 115, 47),
        '0.25x0.3125_CH': (225, 161, 47),
        '0.25x0.3125_WA': (145, 89, 47),
        '0.5x0.625': (145, 133, 47),
        '0.083x0.083': (4320, 2160, 72),  # 9km resolution?
        '0.125x0.125': (2880, 1441, 72),  # nature run (~12km globally)
    }
    if debug:
        print(dims)

    # If full using full vertical
    if full_vert_grid:
        vals = []
        for i in list(dims.values()):
            vals += [(i[0], i[1], 72)]
        dims = dict(list(zip(list(dims.keys()), vals)))
        if debug:
            print(dims)

    # Only consider the GEOS-Chem chemical troposphere
    if trop_limit:
        vals = []
        for i in list(dims.values()):
            vals += [(i[0], i[1], 38)]
        dims = dict(list(zip(list(dims.keys()), vals)))
        if debug:
            print(dims)

    # Add a time dimension of length n (if add_n_time is an integer)
    if not isinstance(add_n_time_dims, type(None)):
        vals = []
        for i in list(dims.values()):
            vals += [(i[0], i[1], i[2], add_n_time_dims)]
        dims = dict(list(zip(list(dims.keys()), vals)))

    # Dictionary of lon, lat (e.g. for emissions and 2D datasets)
    if just2D:
        vals = []
        for i in list(dims.values()):
            vals += [(i[0], i[1])]
        dims = dict(list(zip(list(dims.keys()), vals)))

    if r_dims:
        if invert == True:
            return {v: k for k, v in list(dims.items())}
        else:
            return dims
    else:
        return dims[res]


def get_latlonalt4res(res=None, centre=True, hPa=False, nest=None,
                      dtype=None, wd=None, filename='ctm.nc',
                      full_vert_grid=False,
                      lat_bounds='latitude_bnds', lon_bounds='longitude_bnds',
                      lon_var='longitude', lat_var='latitude', \
                      #        lon_var=u'lon', lat_var=u'lat',
                      verbose=True, debug=False):
    """
    Get lon, lat, and alt for a given model resolution.

    Parameters
    ----------
    wd (str): Specify the wd to get the results from a run.
    res (str): the resolution if wd not given (e.g. '4x5' )
    debug (bool): legacy debug option, replaced by python logging
    lon_var, lat_var (str): variables names for lon and lat in the NetCDF
    lon_bounds, lat_bounds (str): variables names for lon and lat bounds in the NetCDF
    filename (str): name of NetCDF to use
    dtype (type): type for which data is return as, e.g. np.float64
    nest (str): manual override for retruned variables - vestigle?
    hPa (bool): return altitudes in units of hPa, instead of km
    full_vert_grid (bool): use full vertical grid or reduced (47 vs. 72)

    Returns
    -------
    (list) variables for lon, lat, alt as arrays

    Notes
    -----
     - This function uses an updated version of gchem's variable dictionaries
     - This function replaces most use dictionaries from "gchemgrid"
     - The update to using ctm.nc files has cause an bug linked to the lat
     and lon variabe retrival. Just update to passing a wd with output at the
     correct resolution to fix this.
    """
    logging.info("Calling get_latlonalt4res for res={}".format(res))
    if isinstance(res, type(None)):
        logging.warning("No resolution specified. Assuming 4x5!")
        res = '4x5'
    if isinstance(wd, type(None)):
        # Get AC_tools location, then set example data folder location
        #        this_filename = inspect.getframeinfo(inspect.currentframe()).filename
        #        path = os.path.dirname(os.path.abspath(this_filename))
        AC_tools_dir = os.path.dirname(__file__)
        dwd = AC_tools_dir + '/../data/LM/'
        dir_dict = {
            '4x5': 'LANDMAP_LWI_ctm',
            '2x2.5': 'LANDMAP_LWI_ctm_2x25',
            '1x1': 'work/data/GEOS/HEMCO/EMEP/v2015-03/',\
            # Kludge, use 1x1 for 0.5x0.5 <= remove this
            '0.5x0.5': 'work/data/GEOS/HEMCO/EMEP/v2015-03/',\
            '0.5x0.666': 'LANDMAP_LWI_ctm_05x0666',  \
            '0.25x0.3125': 'LANDMAP_LWI_ctm_025x03125',  \
            '0.25x0.3125_CH': 'LANDMAP_LWI_ctm_025x03125_CH',  \
            '0.25x0.3125_WA': 'LANDMAP_LWI_ctm_025x03125_WA',  \
            # Need to add a 0.5x0.625!
            # Temporary inclusion of local 0.083x0.083 file.
            '0.083x0.083': 'LANDMAP_LWI_ctm_0083x0083', \
            # Temporary inclusion of NASA nature run file
            '0.125x0.125': 'LANDMAP_LWI_ctm_0125x0125',
        }
        try:
            dir = dir_dict[res]
        except KeyError:
            logging.error("{res} not a recognised resolution!".format(res=res))
            raise KeyError
        wd = os.path.join(dwd, dir)
        # Local EMEP files for 1x1 and 0.5x0.5
        if (res == '1x1') or (res == '0.5x0.5'):
            filename = 'EMEP.geos.1x1.nc'
            lat_var = 'lat'
            lon_var = 'lon'
            wd = '/work/data/GEOS/HEMCO/EMEP/v2015-03/'
    # Get the data file name
    data_fname = os.path.join(wd, filename)
    if not os.path.exists(data_fname):
        logging.error("Could not find {fn}".format(fn=data_fname))
        raise IOError("Could not find {fn}".format(fn=data_fname))
    # Error message
    LM_file_msg = "ERROR: are the refernces files in 'AC_tools/data/LM' ?"
    LM_file_msg += "\n (To download just run AC_tools/Scripts/get_data_files.py)"
    if centre:
        try:
            # Extract lat and lon from model output data file
            with Dataset(data_fname, 'r') as d:
                lat = d[lat_var][:]
                lon = d[lon_var][:]
        except:
            try:
                print('WARNING: coord vars not found! -using abrvs.')
                print(('Was using: ', lon_var, lat_var))
                lon_var = 'lon'
                lat_var = 'lat'
                print(('Now using: ', lon_var, lat_var))
                # Extract lat and lon from model output data file
                with Dataset(data_fname, 'r') as d:
                    lat = d[lat_var][:]
                    lon = d[lon_var][:]
            except IOError:
                error = "Could not get {lat}, {lon} from {fn}"\
                    .format(fn=data_fname, lat=lat_var, lon=lon_var)
                logging.error(error)
                if verbose:
                    print(LM_file_msg)
                raise IOError(error)
    # Get edge values
    exception_res = ('1x1', '0.5x0.5', '0.083x0.083', '0.125x0.125')
    if (not centre) and (res not in exception_res):
        # Extract lat and lon from model output data file
        try:
            with Dataset(data_fname, 'r') as d:
                lat = d[lat_bounds][:]
                lon = d[lon_bounds][:]
                # Select lower edge of each bound, and final upper edge
                lat = [i[0] for i in lat]+[lat[-1][1]]
                lon = [i[0] for i in lon]+[lon[-1][1]]
                lat, lon = [np.array(i) for i in (lat, lon)]
        except:
            try:
                print('WARNING: coord vars not found! -using abrvs.')
                print(('Was using: ', lon_var, lat_var))
                lon_var = 'lon'
                lat_var = 'lat'
                print(('Now using: ', lon_var, lat_var))
                # Extract lat and lon from model output data file
                with Dataset(data_fname, 'r') as d:
                    lat = d[lat_var][:]
                    lon = d[lon_var][:]
            except IOError:
                error = "Could not get {lat}, {lon} from {fn}"\
                        .format(fn=data_fname, lat=lat_bounds,
                                lon=lon_bounds)
                if verbose:
                    print(LM_file_msg)
                logging.error(error)
                raise IOError(error)
    # Manually set values for 0.25x0.3125_CH
#    if res=='0.25x0.3125_CH':
#        if centre:
#        lat = np.arange(15., 55.+.25, .25)
#        lon = np.arange(70, 140.+.3125, .3125)
#        else:
#            lat = np.array( [-90]+list(np.arange(-89-0.5, 90+.5, 1))+[90] )
#            lon = np.arange(-180-0.5, 180+.5, 1)
    # Manually set (edge) values for 0.125x0.125:
    if res == '0.125x0.125':
        step_size = 0.125
        if centre:
            # below is correct, but use online values for grid
            pass
#            lat = np.arange(-90, 90+step_size, step_size)
#            lon = np.arange(-180, 180, step_size)
        else:
            lat = np.arange(-90-(step_size/2), 90+(step_size/2), step_size)
            lat = np.append(lat, [89.9375+step_size])
            lon = np.arange(-180-(step_size/2), 180+(step_size/2), step_size)
    # Manually set values for (generic?) 1x1 grid
    if res == '1x1':
        step_size = 1.0
        if centre:
            lat = np.arange(-90, 90+step_size, step_size)
            lon = np.arange(-180, 180, step_size)
        else:
            lat = np.array([-90]+list(np.arange(-89-(step_size/2),
                                                90+(step_size/2), step_size))+[90]
                           )
            lon = np.arange(-180-(step_size/2), 180+(step_size/2), step_size)
    # Manually set values for (generic?) 0.5x0.5 grid
    if res == '0.5x0.5':
        step_size = 0.5
        if centre:
            lat = np.array([-90]+list(np.arange(-89, 90, step_size))+[90])
            lon = np.arange(-180, 180, step_size)
        else:
            lat = np.array([-90]+list(np.arange(-89.75, 90+(step_size/2),
                                                step_size))+[90]
                           )
            lon = np.arange(-180-(step_size/2), 180+(step_size/2), step_size)
    # Manually set values for 0.1x0.1
    # Manually set values for 0.083x0.083
    if res == '0.083x0.083':
        step_size = 0.083333336
        if centre:
            lat = np.arange(-89.95833588, 89.95833588+step_size, step_size)
            lon = np.arange(-179.95832825, 179.95835876, step_size)
            # adjust to center point
            lat = [i+step_size/2 for i in lat[:-1]]
            lon = [i+step_size/2 for i in lon[:-1]]
        else:
            lat = np.arange(-89.95833588, 89.95833588+step_size, step_size)
            lon = np.arange(-179.95832825, 179.95835876, step_size)
    # Get dictionary variable name in Gerrit's GEOS-Chem dimensions list
    # ( now only doing this for alt, as alt values not in model output? )
    if hPa:
        alt = 'c_hPa_geos5'
    else:
        alt = 'c_km_geos5'
    # Use reduced vertical grid? (then add '_r')
    if not full_vert_grid:
        alt += '_r'
    d = gchemgrid(rtn_dict=True)
    alt = d[alt]

    # Also provide high resolution grid if requested from this function all
    if nest == 'high res global':
        lon, lat = np.arange(-180, 180, 0.25), np.arange(-90, 90, 0.25)
        return lon, lat, alt

    if debug:
        print((lon, lat, alt))
    rtn_list = lon, lat, alt
    if not isinstance(dtype, type(None)):
        return [i.astype(dtype) for i in rtn_list]
    else:
        return rtn_list


def hPa_to_Km(input, reverse=False, debug=False):
    """
    convert hPa to km

    Parameters
    ----------
    input (list): list of values (float) to convert
    reverse (bool): Set "reverse" to True to convert Km to hPa
    debug (bool): legacy debug option, replaced by python logging

    Returns
    -------
    (list)
    """
    if reverse:
        return [np.exp(np.float(i) / -7.6)*1013. for i in input]
    else:
        return [-7.6*np.log(float(i) / 1013.) for i in input]


def km2nautical_miles(input):
    """
    Convert km into nautical miles
    """
    return [i*0.539957 for i in input]


def km2degrees(aproximate=True):
    """
    Convert km into degrees
    """
    if aproximate:
        return [i/110. for i in input]
    else:
        print('TODO: setup a full calculation for km from degrees lon')


def m2ft(input):
    """
    Convert m into ft (or km to kft)
    """
    return [i*3.281 for i in input]


def find_nearest_value(array, value):
    """
    Find nearest point.

    Parameters
    ----------
    arrary (np.array): 1D array in which to search for nearest value
    value (float): value to search array for closest point

    Notes
    ----------
     - Uses numpy.argmin, therefore: "In case of multiple occurrences of the
    minimum values, the indices  corresponding to the first occurrence are
    returned."
    """
    # Adapted from (credit:) HappyLeapSecond's Stackoverflow answer.
    # ( http://stackoverflow.com/questions/2566412/  )
    idx = (np.abs(array-value)).argmin()
    return idx


def iGEOSChem_ver(wd, also_return_GC_version=False, verbose=True, debug=False):
    """
    Get iGEOS-Chem verson

    NOTES:
     - These are not GEOSChem versions, but halogen branch versions
    (e.g. iGeosChem 1.1 or 1.2 from dir name ( wd )  )
    """
    # List iGEOSChem versions+ then DataFrame
    versions = [
        '1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '1.6.1', '1.6.2',
        '1.6.3', '1.7', '2.0', '3.0', '4.0', '5.0', '6.0', '6.1', '6.2', '6.3',
        '7.0', '7.0.1', '7.1', '7.1.1', \
        # Also hold a 'NOT FOUND' value to for ease of processing non-halogen code
        'NOT FOUND'
    ]
    df = pd.DataFrame(versions, columns=['Versions'])
    if debug:
        print((wd, versions, df))

    # Which versions are in name?
    def element_in_str(element):
        return (element in wd)
    df['Run Version'] = df['Versions'].apply(element_in_str)

    # Select last value and return as string
    try:
        iGC_ver = df['Versions'][df['Run Version']][-1:].values[0]
    except IndexError:
        err_msq = 'WARNING (i) GEOS-Chem ver. number not found in working dir. path'
        print(('!'*15, err_msq, '!'*15))
        logging.debug(err_msq)
        iGC_ver = 'NOT FOUND'
#        sys.exit()

    if also_return_GC_version:
        # list GEOS-Chem versions (written with dashes and underscores)
        versions = [
            'v11-02', 'v12.0.0',
            'v11-01', 'v11_01', 'v10-01', 'v10_01', 'v9-02', 'v9_02', 'v9-01-03',
            'v9_01_03', 'v9-01-02', 'v9_01_02', 'v9-01-01', 'v9_01_01', 'v8-03-02',
            'v8_03_02', 'v8-03-01', 'v8_03_01', 'v8-02-04', 'v8_02_04', 'v8-02-03',
            'v8_02_03', 'v8-02-02', 'v8_02_02', 'v8-02-01', 'v8_02_01', 'v8-01-04',
            'v8_01_04', 'v8-01-03', 'v8_01_03', 'v8-01-02', 'v8_01_02', 'v8-01-01',
            'v8_01_01', 'v7-04-13', 'v7_04_13', 'v7-04-12', 'v7_04_12'
        ]
        df = pd.DataFrame(versions, columns=['Versions'])
        if debug:
            print((wd, versions, df))
        #
        df['Run Version'] = df['Versions'].apply(element_in_str)
        # selection
        try:
            GC_ver = df['Versions'][df['Run Version']][-1:].values[0]
        except IndexError:
            # map iGEOS-Chem versions to GEOS-Chem versions
            dict_iGC_GC = {
                '1.1': 'v9-2',
                '1.2': 'v9-2',
                '1.3': 'v9-2',
                '1.4': 'v9-2',
                '1.5': 'v9-2',
                '1.6': 'v9-2',
                '1.6.1': 'v9-2',
                '1.6.2': 'v9-2',
                '1.6.3': 'v9-2',
                '1.7': 'v9-2',
                '3.0': 'v10-01',
                '4.0': 'v10-01',
                '5.0': 'v11-01',
                '6.0': 'v11-02',
                '6.1': 'v11-02',
                '6.2': 'v11-02',
                '6.3': 'v11-02',
                '7.0': 'v12.0.0',
                '7.0.1': 'v12.0.0',
                '7.1.1': 'v12.0.0',
                'NOT FOUND': 'NOT FOUND',
            }
            #
            GC_ver = dict_iGC_GC[iGC_ver]
        return iGC_ver, GC_ver

    return iGC_ver


def gchemgrid(input=None, rtn_dict=False, debug=False):
    """
    GeosChem grid lookup table. Can give the latitude, longitude or altitude
    positions of geoschem in units to convert from model units.

    Parameters
    ----------
    rtn_dict (bool): return the lookup dictionary instead
    debug (bool): legacy debug option, replaced by python logging

    Returns
    -------
    (np.array)

    Notes
    -----
    - Detail from original reposotiroy adapted from:
    Updated to dictionary from gchemgrid (credit: Gerrit Kuhlmann ) with
    additional grid adds

    This [function] contains (some) grid coordinates used within GEOS-Chem
    as numpy arrays

    Distribution advice from gchem (credit: Gerrit Kuhlmann):

    Python Script Collection for GEOS-Chem Chemistry Transport Model (gchem)
    Copyright (C) 2012 Gerrit Kuhlmann

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    """
    logging.info(
        'gchemgrid called for:{} (rtn_dict={})'.format(input, rtn_dict))

    if isinstance(input, type(None)) and (rtn_dict == False):
        raise KeyError('gchemgrid requires an input or rtn_dict=True')
    d = {
        # 4x5
        'c_lon_4x5': np.arange(-180, 175+5, 5),
        'e_lon_4x5': np.arange(-182.5, 177.5+5, 5),
        'c_lat_4x5': np.array([-89] + list(range(-86, 90, 4)) + [89]),
        'e_lat_4x5': np.array([-90] + list(range(-88, 92, 4)) + [90]),

        # Altitude - Grid box level edges (eta coordinate):
        'e_eta_geos5_r': np.array([
            1.00179600e+00,   9.86769000e-01,   9.71665000e-01,
            9.56562000e-01,   9.41459000e-01,   9.26356000e-01,
            9.11253000e-01,   8.96152000e-01,   8.81051000e-01,
            8.65949000e-01,   8.50848000e-01,   8.35748000e-01,
            8.20648000e-01,   8.00515000e-01,   7.75350000e-01,
            7.50186000e-01,   7.25026000e-01,   6.99867000e-01,
            6.74708000e-01,   6.36974000e-01,   5.99251000e-01,
            5.61527000e-01,   5.23819000e-01,   4.86118000e-01,
            4.48431000e-01,   4.10759000e-01,   3.73114000e-01,
            3.35486000e-01,   2.85974000e-01,   2.42774000e-01,
            2.06167000e-01,   1.75170000e-01,   1.48896000e-01,
            1.26563000e-01,   1.07578000e-01,   9.14420000e-02,
            7.77260000e-02,   5.58200000e-02,   3.97680000e-02,
            2.80770000e-02,   1.95860000e-02,   9.19100000e-03,
            4.02600000e-03,   1.62500000e-03,   6.01000000e-04,
            1.99000000e-04,   5.50000000e-05,   0.00000000e+00]),

        # Grid box level edges [km]:
        'e_km_geos5_r': np.array([
            6.00000000e-03,   1.35000000e-01,   2.66000000e-01,
            3.99000000e-01,   5.33000000e-01,   6.69000000e-01,
            8.06000000e-01,   9.45000000e-01,   1.08600000e+00,
            1.22900000e+00,   1.37400000e+00,   1.52000000e+00,
            1.66900000e+00,   1.87100000e+00,   2.12800000e+00,
            2.39200000e+00,   2.66300000e+00,   2.94100000e+00,
            3.22800000e+00,   3.67300000e+00,   4.14000000e+00,
            4.63100000e+00,   5.14900000e+00,   5.69800000e+00,
            6.28300000e+00,   6.91000000e+00,   7.58700000e+00,
            8.32400000e+00,   9.41100000e+00,   1.05050000e+01,
            1.15780000e+01,   1.26330000e+01,   1.36740000e+01,
            1.47060000e+01,   1.57310000e+01,   1.67530000e+01,
            1.77730000e+01,   1.98550000e+01,   2.20040000e+01,
            2.42400000e+01,   2.65960000e+01,   3.17160000e+01,
            3.75740000e+01,   4.42860000e+01,   5.17880000e+01,
            5.99260000e+01,   6.83920000e+01,   8.05810000e+01]),

        # Grid box level edges [hPa]:
        'e_hPa_geos5_r': np.array([
            1.01181400e+03,   9.96636000e+02,   9.81382000e+02,
            9.66128000e+02,   9.50874000e+02,   9.35621000e+02,
            9.20367000e+02,   9.05114000e+02,   8.89862000e+02,
            8.74610000e+02,   8.59358000e+02,   8.44107000e+02,
            8.28856000e+02,   8.08522000e+02,   7.83106000e+02,
            7.57690000e+02,   7.32279000e+02,   7.06869000e+02,
            6.81458000e+02,   6.43348000e+02,   6.05247000e+02,
            5.67147000e+02,   5.29062000e+02,   4.90984000e+02,
            4.52921000e+02,   4.14873000e+02,   3.76851000e+02,
            3.38848000e+02,   2.88841000e+02,   2.45210000e+02,
            2.08236000e+02,   1.76930000e+02,   1.50393000e+02,
            1.27837000e+02,   1.08663000e+02,   9.23660000e+01,
            7.85120000e+01,   5.63880000e+01,   4.01750000e+01,
            2.83680000e+01,   1.97920000e+01,   9.29300000e+00,
            4.07700000e+00,   1.65100000e+00,   6.17000000e-01,
            2.11000000e-01,   6.60000000e-02,   1.00000000e-02]),

        # Grid box level centers (eta-coordinates)
        'c_eta_geos5_r': np.array([
            9.94283000e-01,   9.79217000e-01,   9.64113000e-01,
            9.49010000e-01,   9.33908000e-01,   9.18805000e-01,
            9.03703000e-01,   8.88601000e-01,   8.73500000e-01,
            8.58399000e-01,   8.43298000e-01,   8.28198000e-01,
            8.10582000e-01,   7.87933000e-01,   7.62768000e-01,
            7.37606000e-01,   7.12447000e-01,   6.87287000e-01,
            6.55841000e-01,   6.18113000e-01,   5.80389000e-01,
            5.42673000e-01,   5.04968000e-01,   4.67274000e-01,
            4.29595000e-01,   3.91937000e-01,   3.54300000e-01,
            3.10730000e-01,   2.64374000e-01,   2.24471000e-01,
            1.90668000e-01,   1.62033000e-01,   1.37729000e-01,
            1.17070000e-01,   9.95100000e-02,   8.45840000e-02,
            6.67730000e-02,   4.77940000e-02,   3.39230000e-02,
            2.38320000e-02,   1.43890000e-02,   6.60900000e-03,
            2.82500000e-03,   1.11300000e-03,   4.00000000e-04,
            1.27000000e-04,   2.80000000e-05]),

        # Grid box level centers [km]
        'c_km_geos5_r': np.array([
            7.10000000e-02,   2.01000000e-01,   3.32000000e-01,
            4.66000000e-01,   6.01000000e-01,   7.37000000e-01,
            8.75000000e-01,   1.01600000e+00,   1.15700000e+00,
            1.30100000e+00,   1.44700000e+00,   1.59400000e+00,
            1.76900000e+00,   1.99900000e+00,   2.25900000e+00,
            2.52700000e+00,   2.80100000e+00,   3.08400000e+00,
            3.44800000e+00,   3.90400000e+00,   4.38200000e+00,
            4.88600000e+00,   5.41900000e+00,   5.98500000e+00,
            6.59100000e+00,   7.24100000e+00,   7.94700000e+00,
            8.84800000e+00,   9.93800000e+00,   1.10210000e+01,
            1.20860000e+01,   1.31340000e+01,   1.41700000e+01,
            1.51980000e+01,   1.62220000e+01,   1.72430000e+01,
            1.87270000e+01,   2.08360000e+01,   2.30200000e+01,
            2.53070000e+01,   2.86540000e+01,   3.40240000e+01,
            4.01660000e+01,   4.71350000e+01,   5.48340000e+01,
            6.30540000e+01,   7.21800000e+01]),

        # Grid box level centers [hPa]
        'c_hPa_geos5_r': np.array([
            1.00422500e+03,   9.89009000e+02,   9.73755000e+02,
            9.58501000e+02,   9.43247000e+02,   9.27994000e+02,
            9.12741000e+02,   8.97488000e+02,   8.82236000e+02,
            8.66984000e+02,   8.51732000e+02,   8.36481000e+02,
            8.18689000e+02,   7.95814000e+02,   7.70398000e+02,
            7.44984000e+02,   7.19574000e+02,   6.94163000e+02,
            6.62403000e+02,   6.24298000e+02,   5.86197000e+02,
            5.48105000e+02,   5.10023000e+02,   4.71952000e+02,
            4.33897000e+02,   3.95862000e+02,   3.57850000e+02,
            3.13844000e+02,   2.67025000e+02,   2.26723000e+02,
            1.92583000e+02,   1.63661000e+02,   1.39115000e+02,
            1.18250000e+02,   1.00514000e+02,   8.54390000e+01,
            6.74500000e+01,   4.82820000e+01,   3.42720000e+01,
            2.40800000e+01,   1.45420000e+01,   6.68500000e+00,
            2.86400000e+00,   1.13400000e+00,   4.14000000e-01,
            1.39000000e-01,   3.80000000e-02]),

        # GEOS-5 Native Vertical Grid (72 hybrid pressure-sigma levels)
        # http://acmg.seas.harvard.edu/geos/doc/archive/man.v9-01-02/appendix_3.html
        'c_hPa_geos5_bounds': np.array([
            1.00000000e-02,   1.50000000e-02,   2.00000000e-02,
            2.60000000e-02,   3.30000000e-02,   4.00000000e-02,
            4.80000000e-02,   5.70000000e-02,   6.60000000e-02,
            7.80000000e-02,   8.90000000e-02,   1.05000000e-01,
            1.20000000e-01,   1.40000000e-01,   1.59000000e-01,
            1.85000000e-01,   2.11000000e-01,   2.45000000e-01,
            2.79000000e-01,   3.22000000e-01,   3.65000000e-01,
            4.20000000e-01,   4.76000000e-01,   5.46000000e-01,
            6.17000000e-01,   7.06000000e-01,   7.95000000e-01,
            9.07000000e-01,   1.01900000e+00,   1.16000000e+00,
            1.30100000e+00,   1.47600000e+00,   1.65100000e+00,
            1.86800000e+00,   2.08500000e+00,   2.35300000e+00,
            2.62000000e+00,   2.94800000e+00,   3.27600000e+00,
            3.67700000e+00,   4.07700000e+00,   4.56200000e+00,
            5.04700000e+00,   5.63200000e+00,   6.21700000e+00,
            6.91800000e+00,   7.62000000e+00,   8.45600000e+00,
            9.29300000e+00,   1.02850000e+01,   1.12770000e+01,
            1.24600000e+01,   1.36430000e+01,   1.50500000e+01,
            1.64570000e+01,   1.81240000e+01,   1.97920000e+01,
            2.17610000e+01,   2.37300000e+01,   2.60490000e+01,
            2.83680000e+01,   3.10890000e+01,   3.38100000e+01,
            3.69930000e+01,   4.01750000e+01,   4.39100000e+01,
            4.76440000e+01,   5.20160000e+01,   5.63880000e+01,
            6.14960000e+01,   6.66030000e+01,   7.25580000e+01,
            7.85120000e+01,   8.54390000e+01,   9.23660000e+01,
            1.00514000e+02,   1.08663000e+02,   1.18250000e+02,
            1.27837000e+02,   1.39115000e+02,   1.50393000e+02,
            1.63661000e+02,   1.76930000e+02,   1.92587000e+02,
            2.08244000e+02,   2.26745000e+02,   2.45246000e+02,
            2.67087000e+02,   2.88927000e+02,   3.13966000e+02,
            3.39005000e+02,   3.58038000e+02,   3.77070000e+02,
            3.96112000e+02,   4.15155000e+02,   4.34212000e+02,
            4.53269000e+02,   4.72335000e+02,   4.91401000e+02,
            5.10475000e+02,   5.29550000e+02,   5.48628000e+02,
            5.67706000e+02,   5.86793000e+02,   6.05880000e+02,
            6.24967000e+02,   6.44054000e+02,   6.63146000e+02,
            6.82239000e+02,   6.94969000e+02,   7.07699000e+02,
            7.20429000e+02,   7.33160000e+02,   7.45890000e+02,
            7.58621000e+02,   7.71354000e+02,   7.84088000e+02,
            7.96822000e+02,   8.09556000e+02,   8.19743000e+02,
            8.29929000e+02,   8.37570000e+02,   8.45211000e+02,
            8.52852000e+02,   8.60493000e+02,   8.68135000e+02,
            8.75776000e+02,   8.83418000e+02,   8.91059000e+02,
            8.98701000e+02,   9.06342000e+02,   9.13984000e+02,
            9.21626000e+02,   9.29268000e+02,   9.36911000e+02,
            9.44553000e+02,   9.52195000e+02,   9.59837000e+02,
            9.67480000e+02,   9.75122000e+02,   9.82765000e+02,
            9.90408000e+02,   9.98051000e+02,   1.00565000e+03,
            1.01325000e+03]),
        # Kludge - mid point values of bounded  GEOS-5 Native Vertical Grid
        'c_hPa_geos5': np.array([
            2.25000000e-02,   3.65000000e-02,   5.20000000e-02,
            7.05000000e-02,   9.45000000e-02,   1.27500000e-01,
            1.68500000e-01,   2.24000000e-01,   2.96000000e-01,
            3.86500000e-01,   5.04000000e-01,   6.52500000e-01,
            8.39500000e-01,   1.07500000e+00,   1.37150000e+00,
            1.73850000e+00,   2.19350000e+00,   2.75350000e+00,
            3.44000000e+00,   4.27700000e+00,   5.28950000e+00,
            6.50950000e+00,   7.97100000e+00,   9.71150000e+00,
            1.17730000e+01,   1.42345000e+01,   1.71605000e+01,
            2.06260000e+01,   2.47145000e+01,   2.95275000e+01,
            3.51705000e+01,   4.17660000e+01,   4.95110000e+01,
            5.85740000e+01,   6.91565000e+01,   8.14890000e+01,
            9.58295000e+01,   1.12737500e+02,   1.32630500e+02,
            1.56032000e+02,   1.83564500e+02,   2.16072500e+02,
            2.54496500e+02,   2.99847000e+02,   3.51524500e+02,
            3.86586000e+02,   4.24676500e+02,   4.62797500e+02,
            5.00934000e+02,   5.39087500e+02,   5.77245000e+02,
            6.15423500e+02,   6.53597500e+02,   6.91785500e+02,
            7.14064000e+02,   7.39525500e+02,   7.64986500e+02,
            7.90455000e+02,   8.15923000e+02,   8.35022000e+02,
            8.49031500e+02,   8.64313500e+02,   8.79596500e+02,
            8.94879500e+02,   9.10162500e+02,   9.25447000e+02,
            9.40732500e+02,   9.56016000e+02,   9.71301500e+02,
            9.86586500e+02,   1.00187250e+03,   1.01705000e+03]),
        # GEOS-5 Native Vertical Grid (72 hybrid pressure-sigma levels)
        'c_km_geos5': np.array([
            0.0905,   0.2215,   0.3535,   0.4875,   0.623,   0.7605,
            0.899,   1.0395,   1.182,   1.3265,   1.473,   1.6215,
            1.8095,   2.053,   2.3155,   2.5855,   2.862,   3.1465,
            3.552,   4.014,   4.499,   5.0105,   5.5525,   6.1285,
            6.745,   7.4095,   8.1315,   9.1275,  10.22,  11.2995,
            12.3595,  13.404,  14.438,  15.4645,  16.4875,  17.508,
            18.538,  19.582,  20.642,  21.721,  22.8195,  23.944,
            25.098,  26.2835,  27.502,  28.754,  30.0415,  31.3655,
            32.7365,  34.1605,  35.637,  37.1665,  38.7505,  40.388,
            42.078,  43.8205,  45.613,  47.454,  49.3395,  51.271,
            53.2445,  55.2555,  57.299,  59.37,  61.4615,  63.567,
            65.68,  67.8175,  70.0485,  72.496,  75.4755,  79.3635]),
    }
    #
    Temp_arrays = ['c_km_geos5', 'c_hPa_geos5']
    if input in Temp_arrays:
        print(('WARNING array ({}) is temporary! - Please check it!'.format(input)))

    if rtn_dict:
        return d
    else:
        return d[input]


def grids4reses(just_1x1_grids=False):
    """
    Function to store coordinates for various model resolutions (lat and lons)
    """
    # --- LAT, LON dictionary
    use_v_0_0_0_grids = False
    use_v_0_0_1_grids = True
    if use_v_0_0_0_grids:
        d = {
            # BASE
            #  GEOSChem in GEOS5 \citep{Hu2017_ACPD}
            # 0.25 GEOS5
            '0.25x0.25_deg_centre_GENERIC': {
                'lon': np.arange(-180.125, 180.125, 0.25),
                'lat': np.arange(-90.125, 90.125, 0.25)
            },
            # Do a 0.5 degree grid for the LongHurst Provinces
            '0.5x0.5_deg_centre_GENERIC': {
                'lon': np.arange(-180.5, 180.5, 0.5),
                'lat': np.arange(-90.5, 90.5, 0.5)
            },
            # Do a 1x1 degree grid
            '1x1_deg_centre_GENERIC': {
                'lon': np.arange(-180.5, 180.5, 1),
                'lat': np.arange(-90.5, 90.5, 1)
            },
            # Do a 1x1 degree grid (centered on 0.5)
            '1x1_deg_0.5_centre_GENERIC': {
                'lon': np.arange(-180, 181, 1),
                'lat': np.arange(-90, 90, 1),
            },
            # GISS ModelE (Miller et al., 2014)
            '2x2.5_deg_centre_GISS': {
                'lon': np.arange(-178.75, 182.75, 2.5),
                'lat': np.arange(-90, 90, 2)
            },
            # ACCMIP (Lamarque et al., 2013)
            '2x2_deg_centre_ACCMIP': {
                'lon': np.arange(-180, 180, 5),
                'lat': np.arange(-90, 90, 4)
            },
            # GEOSChem (Bey et al., 2001) - 4◦ x5◦
            '2x2.5_deg_centre_GEOSChem': {
                'lon': np.arange(-181.25, 181.25, 2.5),
                'lat': np.arange(-91, 91, 2)
            },
            # UKCA (O’Connor et al., 2014)
            '2x3.75_deg_centre_UKCA': {
                'lon': np.arange(-180, 180, 5),
                'lat': np.arange(-90, 90, 4)
            },
            # GEOSChem (Bey et al., 2001) - 4◦ x5◦
            '4x5_deg_centre_GEOSChem': {
                'lon': np.arange(-182.5, 182, 5),
                'lat': np.arange(-92, 90, 4)
            },
        }
    elif (use_v_0_0_1_grids):
        d = {
            # BASE
            #  GEOSChem in GEOS5 \citep{Hu2017_ACPD}
            # 0.25 GEOS5
            '0.25x0.25_deg_centre_GENERIC': {
                'lon': np.array([-180+(i*0.25) for i in np.arange((360./0.25)-1)]),
                'lat': np.array([-90+(i*0.25) for i in np.arange((180./0.25)+1)])
            },
            # Do a 0.5 degree grid for the LongHurst Provinces
            '0.5x0.5_deg_centre_GENERIC': {
                'lon': np.array([-180+(i*0.5) for i in np.arange((360./0.5)-1)]),
                'lat': np.array([-90+(i*0.5) for i in np.arange((180./0.5)+1)])
            },
            # Do a 1x1 degree grid
            '1x1_deg_centre_GENERIC': {
                'lon': np.array([
                    -179.5, -178.5, -177.5, -176.5, -175.5, -174.5, -173.5, -172.5,
                    -171.5, -170.5, -169.5, -168.5, -167.5, -166.5, -165.5, -164.5,
                    -163.5, -162.5, -161.5, -160.5, -159.5, -158.5, -157.5, -156.5,
                    -155.5, -154.5, -153.5, -152.5, -151.5, -150.5, -149.5, -148.5,
                    -147.5, -146.5, -145.5, -144.5, -143.5, -142.5, -141.5, -140.5,
                    -139.5, -138.5, -137.5, -136.5, -135.5, -134.5, -133.5, -132.5,
                    -131.5, -130.5, -129.5, -128.5, -127.5, -126.5, -125.5, -124.5,
                    -123.5, -122.5, -121.5, -120.5, -119.5, -118.5, -117.5, -116.5,
                    -115.5, -114.5, -113.5, -112.5, -111.5, -110.5, -109.5, -108.5,
                    -107.5, -106.5, -105.5, -104.5, -103.5, -102.5, -101.5, -100.5,
                    -99.5,  -98.5,  -97.5,  -96.5,  -95.5,  -94.5,  -93.5,  -92.5,
                    -91.5,  -90.5,  -89.5,  -88.5,  -87.5,  -86.5,  -85.5,  -84.5,
                    -83.5,  -82.5,  -81.5,  -80.5,  -79.5,  -78.5,  -77.5,  -76.5,
                    -75.5,  -74.5,  -73.5,  -72.5,  -71.5,  -70.5,  -69.5,  -68.5,
                    -67.5,  -66.5,  -65.5,  -64.5,  -63.5,  -62.5,  -61.5,  -60.5,
                    -59.5,  -58.5,  -57.5,  -56.5,  -55.5,  -54.5,  -53.5,  -52.5,
                    -51.5,  -50.5,  -49.5,  -48.5,  -47.5,  -46.5,  -45.5,  -44.5,
                    -43.5,  -42.5,  -41.5,  -40.5,  -39.5,  -38.5,  -37.5,  -36.5,
                    -35.5,  -34.5,  -33.5,  -32.5,  -31.5,  -30.5,  -29.5,  -28.5,
                    -27.5,  -26.5,  -25.5,  -24.5,  -23.5,  -22.5,  -21.5,  -20.5,
                    -19.5,  -18.5,  -17.5,  -16.5,  -15.5,  -14.5,  -13.5,  -12.5,
                    -11.5,  -10.5,   -9.5,   -8.5,   -7.5,   -6.5,   -5.5,   -4.5,
                    -3.5,   -2.5,   -1.5,   -0.5,    0.5,    1.5,    2.5,    3.5,
                    4.5,    5.5,    6.5,    7.5,    8.5,    9.5,   10.5,   11.5,
                    12.5,   13.5,   14.5,   15.5,   16.5,   17.5,   18.5,   19.5,
                    20.5,   21.5,   22.5,   23.5,   24.5,   25.5,   26.5,   27.5,
                    28.5,   29.5,   30.5,   31.5,   32.5,   33.5,   34.5,   35.5,
                    36.5,   37.5,   38.5,   39.5,   40.5,   41.5,   42.5,   43.5,
                    44.5,   45.5,   46.5,   47.5,   48.5,   49.5,   50.5,   51.5,
                    52.5,   53.5,   54.5,   55.5,   56.5,   57.5,   58.5,   59.5,
                    60.5,   61.5,   62.5,   63.5,   64.5,   65.5,   66.5,   67.5,
                    68.5,   69.5,   70.5,   71.5,   72.5,   73.5,   74.5,   75.5,
                    76.5,   77.5,   78.5,   79.5,   80.5,   81.5,   82.5,   83.5,
                    84.5,   85.5,   86.5,   87.5,   88.5,   89.5,   90.5,   91.5,
                    92.5,   93.5,   94.5,   95.5,   96.5,   97.5,   98.5,   99.5,
                    100.5,  101.5,  102.5,  103.5,  104.5,  105.5,  106.5,  107.5,
                    108.5,  109.5,  110.5,  111.5,  112.5,  113.5,  114.5,  115.5,
                    116.5,  117.5,  118.5,  119.5,  120.5,  121.5,  122.5,  123.5,
                    124.5,  125.5,  126.5,  127.5,  128.5,  129.5,  130.5,  131.5,
                    132.5,  133.5,  134.5,  135.5,  136.5,  137.5,  138.5,  139.5,
                    140.5,  141.5,  142.5,  143.5,  144.5,  145.5,  146.5,  147.5,
                    148.5,  149.5,  150.5,  151.5,  152.5,  153.5,  154.5,  155.5,
                    156.5,  157.5,  158.5,  159.5,  160.5,  161.5,  162.5,  163.5,
                    164.5,  165.5,  166.5,  167.5,  168.5,  169.5,  170.5,  171.5,
                    172.5,  173.5,  174.5,  175.5,  176.5,  177.5,  178.5,  179.5]),
                'lat': np.array([
                    -89.5, -88.5, -87.5, -86.5, -85.5, -84.5, -83.5, -82.5, -81.5,
                    -80.5, -79.5, -78.5, -77.5, -76.5, -75.5, -74.5, -73.5, -72.5,
                    -71.5, -70.5, -69.5, -68.5, -67.5, -66.5, -65.5, -64.5, -63.5,
                    -62.5, -61.5, -60.5, -59.5, -58.5, -57.5, -56.5, -55.5, -54.5,
                    -53.5, -52.5, -51.5, -50.5, -49.5, -48.5, -47.5, -46.5, -45.5,
                    -44.5, -43.5, -42.5, -41.5, -40.5, -39.5, -38.5, -37.5, -36.5,
                    -35.5, -34.5, -33.5, -32.5, -31.5, -30.5, -29.5, -28.5, -27.5,
                    -26.5, -25.5, -24.5, -23.5, -22.5, -21.5, -20.5, -19.5, -18.5,
                    -17.5, -16.5, -15.5, -14.5, -13.5, -12.5, -11.5, -10.5,  -9.5,
                    -8.5,  -7.5,  -6.5,  -5.5,  -4.5,  -3.5,  -2.5,  -1.5,  -0.5,
                    0.5,   1.5,   2.5,   3.5,   4.5,   5.5,   6.5,   7.5,   8.5,
                    9.5,  10.5,  11.5,  12.5,  13.5,  14.5,  15.5,  16.5,  17.5,
                    18.5,  19.5,  20.5,  21.5,  22.5,  23.5,  24.5,  25.5,  26.5,
                    27.5,  28.5,  29.5,  30.5,  31.5,  32.5,  33.5,  34.5,  35.5,
                    36.5,  37.5,  38.5,  39.5,  40.5,  41.5,  42.5,  43.5,  44.5,
                    45.5,  46.5,  47.5,  48.5,  49.5,  50.5,  51.5,  52.5,  53.5,
                    54.5,  55.5,  56.5,  57.5,  58.5,  59.5,  60.5,  61.5,  62.5,
                    63.5,  64.5,  65.5,  66.5,  67.5,  68.5,  69.5,  70.5,  71.5,
                    72.5,  73.5,  74.5,  75.5,  76.5,  77.5,  78.5,  79.5,  80.5,
                    81.5,  82.5,  83.5,  84.5,  85.5,  86.5,  87.5,  88.5,  89.5]),
            },
            # Do a 1x1 degree grid (centered on 0.5)
            '1x1_deg_0.5_centre_GENERIC': {
                'lon': np.array([
                    -180, -179, -178, -177, -176, -175, -174, -173, -172, -171, -170,
                    -169, -168, -167, -166, -165, -164, -163, -162, -161, -160, -159,
                    -158, -157, -156, -155, -154, -153, -152, -151, -150, -149, -148,
                    -147, -146, -145, -144, -143, -142, -141, -140, -139, -138, -137,
                    -136, -135, -134, -133, -132, -131, -130, -129, -128, -127, -126,
                    -125, -124, -123, -122, -121, -120, -119, -118, -117, -116, -115,
                    -114, -113, -112, -111, -110, -109, -108, -107, -106, -105, -104,
                    -103, -102, -101, -100,  -99,  -98,  -97,  -96,  -95,  -94,  -93,
                    -92,  -91,  -90,  -89,  -88,  -87,  -86,  -85,  -84,  -83,  -82,
                    -81,  -80,  -79,  -78,  -77,  -76,  -75,  -74,  -73,  -72,  -71,
                    -70,  -69,  -68,  -67,  -66,  -65,  -64,  -63,  -62,  -61,  -60,
                    -59,  -58,  -57,  -56,  -55,  -54,  -53,  -52,  -51,  -50,  -49,
                    -48,  -47,  -46,  -45,  -44,  -43,  -42,  -41,  -40,  -39,  -38,
                    -37,  -36,  -35,  -34,  -33,  -32,  -31,  -30,  -29,  -28,  -27,
                    -26,  -25,  -24,  -23,  -22,  -21,  -20,  -19,  -18,  -17,  -16,
                    -15,  -14,  -13,  -12,  -11,  -10,   -9,   -8,   -7,   -6,   -5,
                    -4,   -3,   -2,   -1,    0,    1,    2,    3,    4,    5,    6,
                    7,    8,    9,   10,   11,   12,   13,   14,   15,   16,   17,
                    18,   19,   20,   21,   22,   23,   24,   25,   26,   27,   28,
                    29,   30,   31,   32,   33,   34,   35,   36,   37,   38,   39,
                    40,   41,   42,   43,   44,   45,   46,   47,   48,   49,   50,
                    51,   52,   53,   54,   55,   56,   57,   58,   59,   60,   61,
                    62,   63,   64,   65,   66,   67,   68,   69,   70,   71,   72,
                    73,   74,   75,   76,   77,   78,   79,   80,   81,   82,   83,
                    84,   85,   86,   87,   88,   89,   90,   91,   92,   93,   94,
                    95,   96,   97,   98,   99,  100,  101,  102,  103,  104,  105,
                    106,  107,  108,  109,  110,  111,  112,  113,  114,  115,  116,
                    117,  118,  119,  120,  121,  122,  123,  124,  125,  126,  127,
                    128,  129,  130,  131,  132,  133,  134,  135,  136,  137,  138,
                    139,  140,  141,  142,  143,  144,  145,  146,  147,  148,  149,
                    150,  151,  152,  153,  154,  155,  156,  157,  158,  159,  160,
                    161,  162,  163,  164,  165,  166,  167,  168,  169,  170,  171,
                    172,  173,  174,  175,  176,  177,  178,  179]),
                'lat': np.array([
                    -90, -89, -88, -87, -86, -85, -84, -83, -82, -81, -80, -79, -78,
                    -77, -76, -75, -74, -73, -72, -71, -70, -69, -68, -67, -66, -65,
                    -64, -63, -62, -61, -60, -59, -58, -57, -56, -55, -54, -53, -52,
                    -51, -50, -49, -48, -47, -46, -45, -44, -43, -42, -41, -40, -39,
                    -38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26,
                    -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13,
                    -12, -11, -10,  -9,  -8,  -7,  -6,  -5,  -4,  -3,  -2,  -1,   0,
                    1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,
                    14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,
                    27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
                    40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,
                    53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,
                    66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,
                    79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90]),
            },
            # GISS ModelE (Miller et al., 2014)
            '2x2.5_deg_centre_GISS': {
                'lon': np.array([
                    -177.5, -175., -172.5, -170., -167.5, -165., -162.5, -160.,
                    -157.5, -155., -152.5, -150., -147.5, -145., -142.5, -140.,
                    -137.5, -135., -132.5, -130., -127.5, -125., -122.5, -120.,
                    -117.5, -115., -112.5, -110., -107.5, -105., -102.5, -100.,
                    -97.5,  -95.,  -92.5,  -90.,  -87.5,  -85.,  -82.5,  -80.,
                    -77.5,  -75.,  -72.5,  -70.,  -67.5,  -65.,  -62.5,  -60.,
                    -57.5,  -55.,  -52.5,  -50.,  -47.5,  -45.,  -42.5,  -40.,
                    -37.5,  -35.,  -32.5,  -30.,  -27.5,  -25.,  -22.5,  -20.,
                    -17.5,  -15.,  -12.5,  -10.,   -7.5,   -5.,   -2.5,    0.,
                    2.5,    5.,    7.5,   10.,   12.5,   15.,   17.5,   20.,
                    22.5,   25.,   27.5,   30.,   32.5,   35.,   37.5,   40.,
                    42.5,   45.,   47.5,   50.,   52.5,   55.,   57.5,   60.,
                    62.5,   65.,   67.5,   70.,   72.5,   75.,   77.5,   80.,
                    82.5,   85.,   87.5,   90.,   92.5,   95.,   97.5,  100.,
                    102.5,  105.,  107.5,  110.,  112.5,  115.,  117.5,  120.,
                    122.5,  125.,  127.5,  130.,  132.5,  135.,  137.5,  140.,
                    142.5,  145.,  147.5,  150.,  152.5,  155.,  157.5,  160.,
                    162.5,  165.,  167.5,  170.,  172.5,  175.,  177.5,  180.]),
                'lat': np.array([
                    -89, -87, -85, -83, -81, -79, -77, -75, -73, -71, -69, -67, -65,
                    -63, -61, -59, -57, -55, -53, -51, -49, -47, -45, -43, -41, -39,
                    -37, -35, -33, -31, -29, -27, -25, -23, -21, -19, -17, -15, -13,
                    -11,  -9,  -7,  -5,  -3,  -1,   1,   3,   5,   7,   9,  11,  13,
                    15,  17,  19,  21,  23,  25,  27,  29,  31,  33,  35,  37,  39,
                    41,  43,  45,  47,  49,  51,  53,  55,  57,  59,  61,  63,  65,
                    67,  69,  71,  73,  75,  77,  79,  81,  83,  85,  87,  89])
            },
            # ACCMIP (Lamarque et al., 2013)
            '2x2_deg_centre_ACCMIP': {
                'lon': np.array([
                    -179, -177, -175, -173, -171, -169, -167, -165, -163, -161, -159,
                    -157, -155, -153, -151, -149, -147, -145, -143, -141, -139, -137,
                    -135, -133, -131, -129, -127, -125, -123, -121, -119, -117, -115,
                    -113, -111, -109, -107, -105, -103, -101,  -99,  -97,  -95,  -93,
                    -91,  -89,  -87,  -85,  -83,  -81,  -79,  -77,  -75,  -73,  -71,
                    -69,  -67,  -65,  -63,  -61,  -59,  -57,  -55,  -53,  -51,  -49,
                    -47,  -45,  -43,  -41,  -39,  -37,  -35,  -33,  -31,  -29,  -27,
                    -25,  -23,  -21,  -19,  -17,  -15,  -13,  -11,   -9,   -7,   -5,
                    -3,   -1,    1,    3,    5,    7,    9,   11,   13,   15,   17,
                    19,   21,   23,   25,   27,   29,   31,   33,   35,   37,   39,
                    41,   43,   45,   47,   49,   51,   53,   55,   57,   59,   61,
                    63,   65,   67,   69,   71,   73,   75,   77,   79,   81,   83,
                    85,   87,   89,   91,   93,   95,   97,   99,  101,  103,  105,
                    107,  109,  111,  113,  115,  117,  119,  121,  123,  125,  127,
                    129,  131,  133,  135,  137,  139,  141,  143,  145,  147,  149,
                    151,  153,  155,  157,  159,  161,  163,  165,  167,  169,  171,
                    173,  175,  177,  179]),
                'lat': np.array([
                    -89, -87, -85, -83, -81, -79, -77, -75, -73, -71, -69, -67, -65,
                    -63, -61, -59, -57, -55, -53, -51, -49, -47, -45, -43, -41, -39,
                    -37, -35, -33, -31, -29, -27, -25, -23, -21, -19, -17, -15, -13,
                    -11,  -9,  -7,  -5,  -3,  -1,   1,   3,   5,   7,   9,  11,  13,
                    15,  17,  19,  21,  23,  25,  27,  29,  31,  33,  35,  37,  39,
                    41,  43,  45,  47,  49,  51,  53,  55,  57,  59,  61,  63,  65,
                    67,  69,  71,  73,  75,  77,  79,  81,  83,  85,  87,  89])
            },
            # GEOSChem (Bey et al., 2001) - 4◦ x5◦
            '2x2.5_deg_centre_GEOSChem': {
                'lon': np.array([
                    -180., -177.5, -175., -172.5, -170., -167.5, -165., -162.5,
                    -160., -157.5, -155., -152.5, -150., -147.5, -145., -142.5,
                    -140., -137.5, -135., -132.5, -130., -127.5, -125., -122.5,
                    -120., -117.5, -115., -112.5, -110., -107.5, -105., -102.5,
                    -100.,  -97.5,  -95.,  -92.5,  -90.,  -87.5,  -85.,  -82.5,
                    -80.,  -77.5,  -75.,  -72.5,  -70.,  -67.5,  -65.,  -62.5,
                    -60.,  -57.5,  -55.,  -52.5,  -50.,  -47.5,  -45.,  -42.5,
                    -40.,  -37.5,  -35.,  -32.5,  -30.,  -27.5,  -25.,  -22.5,
                    -20.,  -17.5,  -15.,  -12.5,  -10.,   -7.5,   -5.,   -2.5,
                    0.,    2.5,    5.,    7.5,   10.,   12.5,   15.,   17.5,
                    20.,   22.5,   25.,   27.5,   30.,   32.5,   35.,   37.5,
                    40.,   42.5,   45.,   47.5,   50.,   52.5,   55.,   57.5,
                    60.,   62.5,   65.,   67.5,   70.,   72.5,   75.,   77.5,
                    80.,   82.5,   85.,   87.5,   90.,   92.5,   95.,   97.5,
                    100.,  102.5,  105.,  107.5,  110.,  112.5,  115.,  117.5,
                    120.,  122.5,  125.,  127.5,  130.,  132.5,  135.,  137.5,
                    140.,  142.5,  145.,  147.5,  150.,  152.5,  155.,  157.5,
                    160.,  162.5,  165.,  167.5,  170.,  172.5,  175.,  177.5]),
                'lat': np.array([
                    -89.5, -88., -86., -84., -82., -80., -78., -76., -74.,
                    -72., -70., -68., -66., -64., -62., -60., -58., -56.,
                    -54., -52., -50., -48., -46., -44., -42., -40., -38.,
                    -36., -34., -32., -30., -28., -26., -24., -22., -20.,
                    -18., -16., -14., -12., -10.,  -8.,  -6.,  -4.,  -2.,
                    0.,   2.,   4.,   6.,   8.,  10.,  12.,  14.,  16.,
                    18.,  20.,  22.,  24.,  26.,  28.,  30.,  32.,  34.,
                    36.,  38.,  40.,  42.,  44.,  46.,  48.,  50.,  52.,
                    54.,  56.,  58.,  60.,  62.,  64.,  66.,  68.,  70.,
                    72.,  74.,  76.,  78.,  80.,  82.,  84.,  86.,  88.,
                    89.5])
            },
            # UKCA (O’Connor et al., 2014)
            '2x3.75_deg_centre_UKCA': {
                'lon': np.array([
                    -178.125, -174.375, -170.625, -166.875, -163.125, -159.375,
                    -155.625, -151.875, -148.125, -144.375, -140.625, -136.875,
                    -133.125, -129.375, -125.625, -121.875, -118.125, -114.375,
                    -110.625, -106.875, -103.125,  -99.375,  -95.625,  -91.875,
                    -88.125,  -84.375,  -80.625,  -76.875,  -73.125,  -69.375,
                    -65.625,  -61.875,  -58.125,  -54.375,  -50.625,  -46.875,
                    -43.125,  -39.375,  -35.625,  -31.875,  -28.125,  -24.375,
                    -20.625,  -16.875,  -13.125,   -9.375,   -5.625,   -1.875,
                    1.875,    5.625,    9.375,   13.125,   16.875,   20.625,
                    24.375,   28.125,   31.875,   35.625,   39.375,   43.125,
                    46.875,   50.625,   54.375,   58.125,   61.875,   65.625,
                    69.375,   73.125,   76.875,   80.625,   84.375,   88.125,
                    91.875,   95.625,   99.375,  103.125,  106.875,  110.625,
                    114.375,  118.125,  121.875,  125.625,  129.375,  133.125,
                    136.875,  140.625,  144.375,  148.125,  151.875,  155.625,
                    159.375,  163.125,  166.875,  170.625,  174.375,  178.125]),
                'lat': np.array([
                    -88.75, -86.25, -83.75, -81.25, -78.75, -76.25, -73.75, -71.25,
                    -68.75, -66.25, -63.75, -61.25, -58.75, -56.25, -53.75, -51.25,
                    -48.75, -46.25, -43.75, -41.25, -38.75, -36.25, -33.75, -31.25,
                    -28.75, -26.25, -23.75, -21.25, -18.75, -16.25, -13.75, -11.25,
                    -8.75,  -6.25,  -3.75,  -1.25,   1.25,   3.75,   6.25,   8.75,
                    11.25,  13.75,  16.25,  18.75,  21.25,  23.75,  26.25,  28.75,
                    31.25,  33.75,  36.25,  38.75,  41.25,  43.75,  46.25,  48.75,
                    51.25,  53.75,  56.25,  58.75,  61.25,  63.75,  66.25,  68.75,
                    71.25,  73.75,  76.25,  78.75,  81.25,  83.75,  86.25,  88.75])
            },
            # GEOSChem (Bey et al., 2001) - 4◦ x5◦
            '4x5_deg_centre_GEOSChem': {
                'lon': np.array([
                    -180, -175, -170, -165, -160, -155, -150, -145, -140, -135, -130,
                    -125, -120, -115, -110, -105, -100,  -95,  -90,  -85,  -80,  -75,
                    -70,  -65,  -60,  -55,  -50,  -45,  -40,  -35,  -30,  -25,  -20,
                    -15,  -10,   -5,    0,    5,   10,   15,   20,   25,   30,   35,
                    40,   45,   50,   55,   60,   65,   70,   75,   80,   85,   90,
                    95,  100,  105,  110,  115,  120,  125,  130,  135,  140,  145,
                    150,  155,  160,  165,  170,  175]),
                'lat': np.array([
                    -89, -86, -82, -78, -74, -70, -66, -62, -58, -54, -50, -46, -42,
                    -38, -34, -30, -26, -22, -18, -14, -10,  -6,  -2,   2,   6,  10,
                    14,  18,  22,  26,  30,  34,  38,  42,  46,  50,  54,  58,  62,
                    66,  70,  74,  78,  82,  86,  89])
            },
        }

    else:
        print('No selection of grid version given')
        sys.exit()

    if just_1x1_grids:
        grids2use = [i for i in d.keys() if '1x1' in i]
        d_orig = d.copy()
        d = {}
        for grid in grids2use:
            d[grid] = d_orig[grid]
        del d_orig
        return d

    else:
        return d


def get_sigfig(x, p=3):
    """
    Return a number with only the significant figures required.

    Parameters
    -------
    x (float): number to convert
    p (int): The number of sig figs you want returned. Default=3

    Returns
    -------
    number with only significant figures.

    Notes
    -----
     - Is this function outdated by get_scientific_number()?
    """

#####################
#    Found at https://github.com/randlet/to-precision/blob/master/to_precision.py
#
#    returns a string representation of x formatted with a precision of p
#    Based on the webkit javascript implementation taken from here:
#    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
##########

    import math
    x = float(x)

    if x == 0.:
        return float("0." + "0"*(p-1))

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e - 1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens - x):
        n = n + 1

    if n >= math.pow(10, p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p - 1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return float("".join(out))

#    output round(x, -int(floor(log10(abs(x)))))
#    return output


def get_scientific_number(number, precision, string=False):
    """
    Gets a number in scientific notation with a given precision.

    Parameters
    -------
    number (float): number you want to change
    precision (Integer): How many significant figures you want
    String (Boolean): Do you want the output returned as a string?
    Output float(default): OR string(if string==True)

    Returns
    -------
    (float) or (str)

    Notes
    -----
     - Returns a rounded number by default, or can be returned as a string
    (Recommended for plotting).
    """
    number = float(number)
    # Special case for 0
    if number == 0.:
        if not string:
            return float("0." + "0"*(precision-1))
        else:
            return ("0." + "0"*(precision-1))

    # If negative prepend with - and use the absolute
    if number < 0:
        sign = "-"
#        precision = precision-1
        number = -number
    else:
        sign = ""

    # Get the exponent
    exponent = int(math.log10(number))
    mantisa = number / math.pow(10, exponent)

    # Fix for leading_zeroes
    if mantisa < 1:
        mantisa = mantisa * 10
        exponent = exponent - 1


#    if exponent < 0:
#        precision = precision+1
    # Get only the sigfigs from the mantisa
    mantisa = round(mantisa, precision)

    # Fix for leading 10
    if mantisa >= 10:
        print("hello from 10")
        mantisa = mantisa/10.0
        exponent = exponent+1

    # Fix for mantisa=1
    if mantisa == 1.0:
        print("mantisa = 1.0")
        mantisa = "1." + "0"*(precision-1)

    # Create the string for the result.
    out = sign + str(mantisa)
    if not exponent == 0:
        out = out + 'E' + str(exponent)

    # Return the result as a float unless asking for a string.
    if not string:
        return float(out)
    else:
        return out
