#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Redundant misc. functions to be eventually removed from AC_tools.
"""

import os
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pandas import DataFrame
# time
import time
import datetime as datetime
# math
from math import radians, sin, cos, asin, sqrt, pi, atan2


def get_arr_edge_indices(arr, res='4x5', extra_points_point_on_edge=None,
                         verbose=True, debug=False):
    """
    Find indices in a lon, lat (2D) grid, where value does not equal a given
    value ( e.g. the edge )
    """
    if verbose:
        print(('get_arr_edge_indices for arr of shape: ', arr.shape))

    # initialise variables
    lon_c, lat_c, NIU = get_latlonalt4res(res=res, centre=True)
    lon_e, lat_e, NIU = get_latlonalt4res(res=res, centre=False)
    lon_diff = lon_e[-5]-lon_e[-6]
    lat_diff = lat_e[-5]-lat_e[-6]
    nn, n, = 0, 0
    last_lat_box = arr[nn, n]
    coords = []
    last_lon_box = arr[nn, n]
    need_lon_outer_edge, need_lat_outer_edge = False, False
    if debug:
        print((lon_e, lat_e))

    # ---- Loop X dimension ( lon )
    for nn, lon_ in enumerate(lon_c):

        # Loop Y dimension ( lat ) and store edges
        for n, lat_ in enumerate(lat_c):

            if debug:
                print((arr[nn, n], last_lat_box, last_lon_box,
                       arr[nn, n] == last_lat_box, arr[nn, n] == last_lon_box))

            if arr[nn, n] != last_lat_box:

                # If 1st lat, selct bottom of box
                point_lon = lon_e[nn]+lon_diff/2
                if need_lat_outer_edge:
                    point_lat = lat_e[n+1]
                else:
                    point_lat = lat_e[n]
                    need_lat_outer_edge = True
                need_lat_outer_edge = False

                # Add mid point to cordinates list
                if isinstance(extra_points_point_on_edge, type(None)):
                    mid_point = [point_lon, point_lat]
                    coords += [mid_point]

                # Add given number of points along edge
                else:
                    coords += [[lon_e[nn]+(lon_diff*i), point_lat] for i in
                               np.linspace(0, 1, extra_points_point_on_edge,
                                           endpoint=True)]

            # temporally save the previous box's value
            last_lat_box = arr[nn, n]

    # ---- Loop Y dimension ( lat )
    for n, lat_ in enumerate(lat_c):

        if debug:
            print((arr[nn, n], last_lat_box, last_lon_box,
                   arr[nn, n] == last_lat_box, arr[nn, n] == last_lon_box))
        # Loop X dimension ( lon ) and store edges
        for nn, lon_ in enumerate(lon_c):

            # If change in value at to list
            if arr[nn, n] != last_lon_box:
                point_lat = lat_e[n]+lat_diff/2

                # Make sure we select the edge lon
                if need_lon_outer_edge:
                    point_lon = lon_e[nn+1]
                else:
                    point_lon = lon_e[nn]
                    need_lon_outer_edge = True
                need_lon_outer_edge = False

                # Add mid point to coordinates list
                if isinstance(extra_points_point_on_edge, type(None)):
                    mid_point = [point_lon, point_lat]
                    coords += [mid_point]

                # Add given number of points along edge
                else:
                    coords += [[point_lon, lat_e[n]+(lat_diff*i)] for i in
                               np.linspace(0, 1, extra_points_point_on_edge,
                                           endpoint=True)]

            # temporally save the previous box's value
            last_lon_box = arr[nn, n]

    return coords


def split_data_by_days(data=None, dates=None, day_list=None,
                       verbose=False, debug=False):
    """
    Takes a list of datetimes and data and returns a list of data and
    the bins ( days )
    """
    if verbose:
        print('split_data_by_days called')

    # Create DataFrame of Data and dates
    df = DataFrame(data, index=dates, columns=['data'])
    # Add list of dates ( just year, month, day ) <= this is mappable, update?
    df['days'] = [datetime.datetime(*i.timetuple()[:3]) for i in dates]
    if debug:
        print(df)

    # Get list of unique days
    if isinstance(day_list, type(None)):
        day_list = sorted(set(df['days'].values))
    # Loop unique days and select data on these days
    data4days = []
    for day in day_list:
        print((day, df[df['days'] == day]))
        data4days += [df['data'][df['days'] == day]]
    # Just return the values ( i.e. not pandas array )
    data4days = [i.values.astype(float) for i in data4days]
    print([type(i) for i in data4days])
#    print data4days[0]
#    sys.exit()

    if debug:
        print(('returning data for {} days, with lengths: '.format(
            len(day_list)), [len(i) for i in data4days]))

    # Return as list of days (datetimes) + list of data for each day
    return data4days, day_list


def obs2grid(glon=None, glat=None, galt=None, nest='high res global',
             sites=None, debug=False):
    """
    values that have a given lat, lon and alt

    Notes
    -------
     - Function flagged for removal
    """
    if isinstance(glon, type(None)):
        glon, glat, galt = get_latlonalt4res(nest=nest, centre=False,
                                             debug=debug)

    # Assume use of known CAST sites... unless others given.
    if isinstance(sites, type(None)):
        loc_dict = get_loc(rtn_dict=True)
        sites = list(loc_dict.keys())

    # Pull out site location indicies
    indices_list = []
    for site in sites:
        lon, lat, alt = loc_dict[site]
        vars = get_xy(lon,  lat, glon, glat)
        indices_list += [vars]
    return indices_list
