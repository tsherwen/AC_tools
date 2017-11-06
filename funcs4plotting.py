#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Generic plotting functions for timeseries/multi-dimensional output.

Use help(<name of function>) to get details on a particular function.

NOTE(S):
 - This module is underdevelopment vestigial/inefficient code is being removed/updated.
 - Where external code is used credit is given.
"""

# ------------------- Section 0 -------------------------------------------
# -------------- Required modules:

# -- Plotting
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import setp
import functools
import matplotlib
#import seaborn as sns

# -- Time
import time
import calendar
import datetime as datetime
from datetime import datetime as datetime_

# -- I/O/Admin...
import gc

# ---  This needs to be updated, imports should be specific and in individual functions
# import tms modules with shared functions
from . funcs_vars import *
from . funcs4generic import *
from . funcs4time import *
from . funcs4pf import *
from . funcs4GEOSC import * # wd2ctms, get_gc_res

# math
from math import log10, floor
import numpy as np
import scipy
# colormaps - Additional maps from Eric Sofen
#from option_c import test_cm as cmc
#from option_d import test_cm as cmd


# ----------------------------- Section 1 ------------------------------------
# -------------- Common Plot Types
#

# ----
# 1.01 - Map plot for given array and resolution (lon, lat
# -----
def map_plot( arr, return_m=False, grid=False, centre=False, cmap=None, \
        no_cb=False, cb=None, rotatecbunits='horizontal',fixcb=None, nticks=10,\
        mask_invalids=False,\
        format='%.2f', adjust_window=0, f_size=20, alpha=1, log=False, \
        set_window=False, res=None, ax=None, case='default', units=None, \
        drawcountries=True,  set_cb_ticks=True, title=None, lvls=None,  \
        interval=1, resolution='c', shrink=0.4, window=False, everyother=1,\
        extend='neither', degrade_resolution=False, discrete_cmap=False, \
        lon_0=None, lon_1=None, lat_0=None, lat_1=None, norm=None,\
        sigfig_rounding_on_cb=2, fixcb_buffered=None, ylabel=True, \
        xlabel=True, wd=None, verbose=True, debug=False, tight_layout=False, \
        axis_titles=False, **Kwargs):
    """
    Plots Global/regional 2D (lon, lat) slices.

    Parameters
    ----------
    adjust_window (int): amount of array entries to remove the edges of array
    alhpa (float): transparency of plotter data
    arr (np.array): input (2D) array
    case (str or int): case for type of plot (vestigle: use log=True of False (default))
    cmap (str): force colormap selection by providing name
    centre (boolean): use centre points of lon/lat grid points for mapping data surface
    drawcountries (boolean): add countries to basemap?
    debug (boolean): legacy debug option, replaced by python logging
    degrade_resolution (boolean): reduce resolution of underlay map detail
    discrete_cmap (boolean): use a discrete instead of conitunous colorbar map
    everyother (int): use "everyother" axis tick (e.g. 3=use every 3rd)
    f_size (float): fontsise
    fixcb (np.array): minimium and maximum to fix colourbar
    fixcb_buffered (array): minimium and maximum to fix colourbar, with buffer space
    format (str): format string for colorbar formating
    grid (boolean): apply a grid over surface plot?
    extend (str): colorbar format settings ( 'both', 'min', 'both' ... )
    interval (int): x/y tick interval in multiples of 15 degrees lat/lon
    lvls (list): manually provide levels for colorbar
    log (boolean): use a log scale for the plot and colorbar
    no_cb (boolean): include a coloubar?
    norm (norm object): normalisation to use for colourbar and array
    nticks (int): number of ticks on colorbar
    mask_invalids (boolean): mask invalid numbers (to allow saving to PDF)
    res (str): GEOS-Chem output configuration resolution ( '4x5' etc... )
    resolution (str): basemasp resolution settings ( 'c' = coarse, 'f' = fine ...  )
    rotatecbunits (str): orientation of colourbar units
    shrink (boolean): colorbar size settings ( fractional shrink )
    set_window (boolean): set the limits of the plotted data (lon_0, lon_1, lat_0, lat_1)
    (for nested boundary conditions )
    sigfig_rounding_on_cb (int): significant figure rounding to use for colourbar
    set_cb_ticks (boolean): mannually set colorbar ticks? (vestigle)
    title (str): plot title (deafult is ==None, therefore no title)
    tight_layout (boolean): use use tight lyaout for figure
    ylabel, xlabel (boolean): label x/y axis?
    units (str): units of given data for plot title
    verbose (boolean): legacy debug option, replaced by python logging
    wd (str): Specify the wd to get the results from a run.
    window (boolean): use window plot settings (fewer axis labels/resolution of map)
    axis_titles (boolean): title X and y axis? (lat and lon)

    Returns
    -------
    optionally returns basemap (return_m==True) and colorbar (no_cb!=True) object

    Notes
    -----
     - Takes a numpy array and the resolution of the output. The plot extent is then set by this output.
    """
    if isinstance(arr, type(None)):
        logging.error("No data given to map_plot!")
        raise AssertionError("No data given to map_plot")
    elif not len(arr.shape)==2:
        logging.error("Input array should be 2D. Got shape {shape}"\
            .format(shape=arr.shape))
    logging.info("map_plot called (array shape={}, res={})".format( \
        arr.shape, res) )

    # Find out what resolution we are using if not specified
    if isinstance(res, type(None)):
        try: # Attempt to extract resolution from wd
            logging.debug("No resolution specified, getting from wd")
            res = get_gc_res(wd)
        except TypeError:  # Assume 4x5 resolution
            # Assume 4x5 resolution
            logging.warning('No resolution specified or found. Assuming 4x5')
            logging.warning('Try specifying the wd or manualy specifying the res')
            res='4x5'

    # Make sure the input data is usable and try to fix it if not.
    (res_lat, res_lon) = get_dims4res(res, just2D=True)
    expected_shape = (res_lon, res_lat)
    print(expected_shape, arr.shape)
    if arr.shape==expected_shape:
        pass
    elif arr.shape!=expected_shape:
        arr = arr.T
        logging.warning("Array was wrong shape and has been transposed!")
        if arr.shape==expected_shape:
            pass
        else:
            err_msg = "Array is the wrong shape. Should be {}. Got {}"\
             .format( str(expected_shape), arr.shape)
            logging.error(err_msg)
            raise AssertionError(err_msg)

    #### Add a invalid warning!
    # Mask for percent arrays containing invalid values ( to allow PDF save )
    if mask_invalids:
        arr = np.ma.masked_invalid( arr )

    # --- Window plot settings
    if window:
        interval = 2   # double interval size
        degrade_resolution=True
    if  res == '0.5x0.666':
        interval,  adjust_window, resolution,shrink  =0.5, 3, 'f', 0.6
    if degrade_resolution:
        resolution = 'l'

    nested_res = ['0.25x0.3125', '0.25x0.3125_CH', '0.25x0.3125_WA']
    if res in nested_res:
        centre=False
        adjust_window = 6

    # Get lons and lats
    lon, lat, NIU = get_latlonalt4res( res, centre=centre, wd=wd )

    if set_window:
        # Convert lats and lons to GC lats and restrict lats, lons, and arr
        if not isinstance( lat_0, type(None) ):
            gclat_0, gclat_1 = [ get_gc_lat(i, res=res) for i in (lat_0, lat_1) ]
            lat = lat[ gclat_0:gclat_1 ]

        if not isinstance( lon_0, type(None) ):
            gclon_0, gclon_1 = [ get_gc_lon(i, res=res) for i in (lon_0, lon_1) ]
            lon = lon[ gclon_0:gclon_1]

    # ----------------  Basemap setup  ----------------
    # Grid/Mesh values
    x, y = np.meshgrid(lon,lat)
    # Set existing axis to current if axis provided
    if not isinstance(ax, type(None)):
        # temporary remove as mpl widget has a bug
        # http://stackoverflow.com/questions/20562582/axes-instance-argument-was-not-found-in-a-figure
#        plt.sca( ax )
        pass

    # ---- Setup map ("m") using Basemap
    m = get_basemap( lat=lat, lon=lon, resolution=resolution, res=res, \
        everyother=everyother, interval=interval, f_size=f_size, ylabel=ylabel, \
        axis_titles=axis_titles,
        xlabel=xlabel, drawcountries=drawcountries )
    # Process data to grid
    x, y = np.meshgrid( *m(lon, lat) )
    # reduce plotted region for nested grids/subregion plots
    if set_window:
        plt.xlim( lon_0, lon_1)
        plt.ylim( lat_0, lat_1 )
    else:
        plt.xlim( lon[0+adjust_window], lon[-1-adjust_window] )
        plt.ylim(lat[0+adjust_window], lat[-1-adjust_window])

    # ----------------  Cases for arrays  ----------------
    # Keep the names for easier reading
####################################################################################################
# Old
#    if ( not isinstance( case, int ) ):
#        case = {
#    'linear':3, 'default':3, 'IO': 1,'limit_to_1_2': 2, 'log': 4,
#        }[case]
####################################################################################################

# New
    if ( isinstance( case, int ) ):
        case = {
                1: 'IO',
                2: 'limit_to_1_2',
                3: 'default',
                4: 'log',
                }[case]

    if case=="IO":
        IO=True
    elif case=="limit_to_1_2":
        limit_to_1_2=True
    elif case=="default":
        default=True
    elif case=="log":
        log=True
    else:
        raise ValueError("Unknown case of {case}".format(case=case))
####################################################################################################

    # -------- colorbar variables...
    # Set cmap range I to limit poly, if not given cmap )
    fixcb_ = fixcb
    # New approach
    if isinstance( fixcb_, type(None) ) or isinstance( cmap, type(None) ):
        fixcb_ = np.array( [ (i.min(), i.max()) for i in [arr ] ][0] )

    if isinstance(cmap, type(None)):
        # Set readable levels for cb, then use these to dictate cmap
        if isinstance(lvls, type(None)):
            lvls = get_human_readable_gradations( vmax=fixcb_[1], vmin=fixcb_[0], \
                nticks=nticks, sigfig_rounding_on_cb=sigfig_rounding_on_cb  )

        # Setup Colormap
        cmap, fixcb_buffered = get_colormap( np.array( fixcb_ ), \
                nticks=nticks, fixcb=fixcb_, buffer_cmap_upper=True )
        # Update colormap with buffer
        cmap = get_colormap( arr=np.array([fixcb_buffered[0], fixcb_buffered[1]]) )

    # Allow function to operate without fixcb_buffered provided
    if isinstance( fixcb_buffered, type(None) ):
        fixcb_buffered = fixcb_
    fixcb_ = fixcb_buffered
    logging.info( 'colorbar variables: ' + str([fixcb_buffered, fixcb, fixcb_, lvls, \
                cmap, lvls]))
    # Use a discrete colour map?
#    if discrete_cmap:
#        if isinstance( fixcb, type(None) ):
#            cmap, norm = mk_discrete_cmap( vmin=arr.min(), vmax=arr.max(), \
#                    nticks=nticks, cmap=cmap )
#        else:
#            cmap, norm = mk_discrete_cmap( vmin=fixcb[0], vmax=fixcb[1], \
#                    nticks=nticks, cmap=cmap

    # --------------  Linear plots -------------------------------
    # standard plot
    linear_cases = ["default",9]
    if case in linear_cases:
        if debug:
            print(fixcb_, arr.shape, [ len(i) for i in (lon, lat) ], norm, cmap)
        poly = m.pcolor( lon, lat, arr, cmap=cmap, norm=norm, alpha=alpha, \
            vmin=fixcb_[0], vmax=fixcb_[1]  )

    # -----------------  Log plots --------------------------------
    if log: # l
        poly = m.pcolor(lon, lat, arr, cmap=cmap, \
            norm=mpl.colors.LogNorm(vmin=fixcb_[0], vmax=fixcb_[1]) )

        if no_cb:
            pass
        else:
            # Get logarithmically spaced integers
            lvls = np.logspace( np.log10(fixcb[0]), np.log10(fixcb[1]), num=nticks)
            # Normalise to Log space
            norm=mpl.colors.LogNorm(vmin=fixcb_[0], vmax=fixcb_[1])

            # Create colourbar instance
            cb = plt.colorbar(poly, ax=m.ax, ticks=lvls, format=format, shrink=shrink, \
                alpha=alpha, norm=norm, extend='min')
        logging.debug(np.ma.min(np.ma.log(arr)), np.ma.max(np.ma.log(arr)), lvls)

    # ----------------  Colorbars  ----------------
    if no_cb:
        pass
    else:
        if isinstance(cb, type(None)):
            # if linear plot without fixcb set, then define here
            ax = plt.gca()

        # Create colourbar instance
        cb = plt.colorbar( poly, ax=ax, shrink=shrink, alpha=alpha, extend=extend )
        # set ylabel tick properties
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(f_size)

        if not isinstance(units, type(None)):
            cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size)

        # Special treatment for log colorbars
        if log:
            round_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))
            tick_locs = [ float('{:.2g}'.format( t )) for t in lvls ]
            # for asectics, round colorbar labels to sig figs given
            for n, lvl in enumerate( lvls ):
                try:
                    lvls[n] = round_to_n( lvl, sigfig_rounding_on_cb)
                except:
                    lvls[n] = lvl
        else:
            tick_locs = np.array( lvls ).copy()

        #make sure the ticks are numbers
#        temp_tick_locs = []
#        for i_tick in tick_locs:
        tick_locs = [float(_tick) for _tick in tick_locs]

        # fix colorbar levels, then provide labels

        cb.set_ticks( np.array(tick_locs) )
        # the format is not correctly being set... - do this manually instead
#        if not isinstance( format, type(None) ):
#            lvls = [ format % (i) for i in lvls ]
        cb.set_ticklabels( lvls )#, format=format )


        #logging.info(tick_locs, lvls, [ type(i) for i in tick_locs, lvls ])
        #logging.info(cb.get_clim(), title, format)

# Set number of ticks
# FIX NEEDED - this currently doesn't doesn't work for log plots
#    if (case != 3) and (not no_cb) and ( case != 4):
#        if set_cb_ticks:
#            tick_locator = ticker.MaxNLocator( nticks=nticks )
#            cb.locator = tick_locator
#            cb.update_ticks()
    # Add grid lines to the plot?
    plt.grid( grid )

    # Add title to plot?
    max_title_len=30

    if tight_layout==True:
        plt.tight_layout()

    if not isinstance( title, type(None) ):
        # Check if the title is too long and if not split it over lines
        if len(title)>max_title_len:
            print("tile takes up multiple lines. Splitting over lines now.")
            import textwrap
            title="\n".join(textwrap.wrap(title,max_title_len))
            print(title)
            # Adjust the top of the plot by 0.05 for every line the title takes
#            plt.subplots_adjust(top=1-0.05*(len(title)%max_title_len))

        plt.title(title, fontsize=f_size*1.5)
    # Setup list of return variables
    return_l = [ plt ]
    if not no_cb:
        return_l += [ cb ]
    if return_m:
        return_l += [ m ]
    return return_l


# ----
# X.XX - UPDATED - Map plot for given array and resolution (lon, lat)
# -----
def plot_map( arr, return_m=False, grid=False, centre=False, cmap=None, no_cb=False, \
        cb=None, rotatecbunits='horizontal',fixcb=None, nticks=10, mask_invalids=False,\
        format='%.2f', adjust_window=0, f_size=20, alpha=1, log=False, \
        set_window=False, res=None, ax=None, case='default', units=None, \
        drawcountries=True,  set_cb_ticks=True, title=None, lvls=None,  \
        interval=15, resolution='c', shrink=0.4, window=False, everyother=1,\
        extend='neither', degrade_resolution=False, discrete_cmap=False, \
        lon_0=None, lon_1=None, lat_0=None, lat_1=None, norm=None,\
        sigfig_rounding_on_cb=2, fixcb_buffered=None, ylabel=True, \
        xlabel=True, wd=None, verbose=True, debug=False, tight_layout=False, \
        **Kwargs):
    """
    Plots Global/regional 2D (lon, lat) slices.

    WARNING - This is an updated version of map_plot (incomplement/develop),
       use map_plot instead!

    Parameters
    ----------
    adjust_window (int): amount of array entries to remove the edges of array
    alhpa (float): transparency of plotter data
    arr (np.array): input (2D) array
    case (str or int): case for type of plot (vestigle: use log=True of False (default))
    cmap (str): force colormap selection by providing name
    centre (boolean): use centre points of lon/lat grid points for mapping data surface
    drawcountries (boolean): add countries to basemap?
    debug (boolean): legacy debug option, replaced by python logging
    degrade_resolution (boolean): reduce resolution of underlay map detail
    discrete_cmap (boolean): use a discrete instead of conitunous colorbar map
    everyother (int): use "everyother" axis tick (e.g. 3=use every 3rd)
    f_size (float): fontsise
    fixcb (np.array): minimium and maximum to fix colourbar
    fixcb_buffered (array): minimium and maximum to fix colourbar, with buffer space
    format (str): format string for colorbar formating
    grid (boolean): apply a grid over surface plot?
    extend (str): colorbar format settings ( 'both', 'min', 'both' ... )
    interval (int): x/y tick interval in degrees lat/lon (default=15)
    lvls (list): manually provide levels for colorbar
    log (boolean): use a log scale for the plot and colorbar
    no_cb (boolean): include a coloubar?
    norm (norm object): normalisation to use for colourbar and array
    nticks (int): number of ticks on colorbar
    mask_invalids (boolean): mask invalid numbers (to allow saving to PDF)
    res (str): GEOS-Chem output configuration resolution ( '4x5' etc... )
    resolution (str): basemasp resolution settings ( 'c' = coarse, 'f' = fine ...  )
    rotatecbunits (str): orientation of colourbar units
    shrink (boolean): colorbar size settings ( fractional shrink )
    set_window (boolean): set the limits of the plotted data (lon_0, lon_1, lat_0, lat_1)
    (for nested boundary conditions )
    sigfig_rounding_on_cb (int): significant figure rounding to use for colourbar
    set_cb_ticks (boolean): mannually set colorbar ticks? (vestigle)
    title (str): plot title (deafult is ==None, therefore no title)
    tight_layout (boolean): use use tight lyaout for figure
    ylabel, xlabel (boolean): label x/y axis?
    units (str): units of given data for plot title
    verbose (boolean): legacy debug option, replaced by python logging
    wd (str): Specify the wd to get the results from a run.
    window (boolean): use window plot settings (fewer axis labels/resolution of map)
    Returns
    -------
    optionally returns basemap (return_m==True) and colorbar (no_cb!=True) object
    Notes
    -----
     - Takes a numpy array and the resolution of the output. The plot extent is then set by this output.
    """
    if isinstance(arr, type(None)):
        logging.error("No data given to map_plot!")
        raise AssertionError("No data given to map_plot")
    elif not len(arr.shape)==2:
        logging.error("Input array should be 2D. Got shape {shape}"\
            .format(shape=arr.shape))
    logging.info("map_plot called")

    # Find out what resolution we are using if not specified
    if isinstance(res, type(None)):
        try: # Attempt to extract resolution from wd
            logging.debug("No resolution specified, getting from wd")
            res = get_gc_res(wd)
        except TypeError:  # Assume 4x5 resolution
            # Assume 4x5 resolution
            logging.warning('No resolution specified or found. Assuming 4x5')
            logging.warning('Try specifying the wd or manualy specifying the res')
            res='4x5'




    # ----------------  Cases for arrays  ----------------
    # Keep the names for easier reading
####################################################################################################
# Old
#    if ( not isinstance( case, int ) ):
#        case = {
#    'linear':3, 'default':3, 'IO': 1,'limit_to_1_2': 2, 'log': 4,
#        }[case]
####################################################################################################

# New
    if ( isinstance( case, int ) ):
        case = {
                1: 'IO',
                2: 'limit_to_1_2',
                3: 'default',
                4: 'log',
                }[case]

    if case=="IO":
        IO=True
    elif case=="limit_to_1_2":
        limit_to_1_2=True
    elif case=="default":
        default=True
    elif case=="log":
        log=True
    else:
        raise ValueError("Unknown case of {case}".format(case=case))
####################################################################################################

    # Make sure the input data is usable and try to fix it if not.
    (res_lat, res_lon) = get_dims4res(res, just2D=True)
    expected_shape = (res_lon, res_lat)

    if arr.shape==expected_shape:
        pass
    elif arr.shape==expected_shape:
        arr = arr.T
        logging.warning("Array was wrong shape and has been transposed!")
    else:
        logging.error("Array is the wrong shape. Should be {}. Got {}"\
         .format( str(expected_shape), arr.shape) )
        raise AssertionError("Incorrect array shape.")

    #### Add a invalid warning!
    # Mask for percent arrays containing invalid values ( to allow PDF save )
    if mask_invalids:
        arr = np.ma.masked_invalid( arr )

    # --- Window plot settings
    if window:
        interval = 30   # double interval size
        degrade_resolution=True
    if  res == '0.5x0.666':
        interval,  adjust_window, resolution,shrink  =0.5, 3, 'f', 0.6
    if degrade_resolution:
        resolution = 'l'

    nested_res = ['0.25x0.3125', '0.25x0.3125_CH', '0.25x0.3125_WA']
    if res in nested_res:
        centre=False
        adjust_window = 6

    # Get lons and lats
    lon, lat, NIU = get_latlonalt4res( res, centre=centre, wd=wd )

    if set_window:
        # Convert lats and lons to GC lats and restrict lats, lons, and arr
        if not isinstance( lat_0, type(None) ):
            gclat_0, gclat_1 = [ get_gc_lat(i, res=res) for i in (lat_0, lat_1) ]
            lat = lat[ gclat_0:gclat_1 ]

        if not isinstance( lon_0, type(None) ):
            gclon_0, gclon_1 = [ get_gc_lon(i, res=res) for i in (lon_0, lon_1) ]
            lon = lon[ gclon_0:gclon_1]

    # ----------------  Basemap setup  ----------------
    # Grid/Mesh values
    x, y = np.meshgrid(lon,lat)
    # Set existing axis to current if axis provided
    if not isinstance(ax, type(None)):
        plt.sca( ax )

    # ---- Setup map ("m") using Basemap
    m = get_basemap( lat=lat, lon=lon, resolution=resolution, res=res, \
        everyother=everyother, interval=interval, f_size=f_size, ylabel=ylabel, \
        xlabel=xlabel, drawcountries=drawcountries )
    # Process data to grid
    x, y = np.meshgrid( *m(lon, lat) )
    # reduce plotted region for nested grids/subregion plots
    if set_window:
        plt.xlim( lon_0, lon_1)
        plt.ylim( lat_0, lat_1 )
    else:
        plt.xlim( lon[0+adjust_window], lon[-1-adjust_window] )
        plt.ylim(lat[0+adjust_window], lat[-1-adjust_window])





####################################################################################################
    # -------- colorbar variables...
    # Set cmap range I to limit poly, if not given cmap )
    fixcb_ = fixcb
    # New approach
    if isinstance( fixcb_, type(None) ) or isinstance( cmap, type(None) ):
        fixcb_ = np.array( [ (i.min(), i.max()) for i in [arr ] ][0] )

    if isinstance(cmap, type(None)):
        # Set readable levels for cb, then use these to dictate cmap
        if isinstance(lvls, type(None)):
            lvls = get_human_readable_gradations( vmax=fixcb_[1], vmin=fixcb_[0], \
                nticks=nticks, sigfig_rounding_on_cb=sigfig_rounding_on_cb  )

        # Setup Colormap
        cmap, fixcb_buffered = get_colormap( np.array( fixcb_ ), \
                nticks=nticks, fixcb=fixcb_, buffer_cmap_upper=True )
        # Update colormap with buffer
        cmap = get_colormap( arr=np.array([fixcb_buffered[0], fixcb_buffered[1]]) )

    # Allow function to operate without fixcb_buffered provided
    if isinstance( fixcb_buffered, type(None) ):
        fixcb_buffered = fixcb_
    fixcb_ = fixcb_buffered
    logging.info( 'colorbar variables: ' + str([fixcb_buffered, fixcb, fixcb_, lvls, \
                cmap, lvls]))
    # Use a discrete colour map?
#    if discrete_cmap:
#        if isinstance( fixcb, type(None) ):
#            cmap, norm = mk_discrete_cmap( vmin=arr.min(), vmax=arr.max(), \
#                    nticks=nticks, cmap=cmap )
#        else:
#            cmap, norm = mk_discrete_cmap( vmin=fixcb[0], vmax=fixcb[1], \
#                    nticks=nticks, cmap=cmap )

    NEW_VERSION=False

    if not NEW_VERSION:
####################################################################################################
            # Old version is here
####################################################################################################
        # --------------  Linear plots -------------------------------
        # standard plot
        linear_cases = ["default",9]
        if case==9 or default:
            if debug:
                print(fixcb_, arr.shape, [ len(i) for i in (lon, lat) ], norm, cmap)
            poly = m.pcolor( lon, lat, arr, cmap=cmap, norm=norm, alpha=alpha, \
                vmin=fixcb_[0], vmax=fixcb_[1]  )

        # -----------------  Log plots --------------------------------
        if log: # l
            poly = m.pcolor(lon, lat, arr, norm=LogNorm(vmin=fixcb_[0], vmax=fixcb_[1]), \
                cmap=cmap)

            if no_cb:
                pass
            else:


                # Get logarithmically spaced integers
                lvls = np.logspace( np.log10(fixcb[0]), np.log10(fixcb[1]), num=nticks)
                # Normalise to Log space
                norm=mpl.colors.LogNorm(vmin=fixcb_[0], vmax=fixcb_[1])

                    # Create colourbar instance
                cb = plt.colorbar(poly, ax=m.ax, ticks=lvls, format=format, shrink=shrink, \
                    alpha=alpha, norm=norm, extend='min')
            logging.debug(np.ma.min(np.ma.log(arr)), np.ma.max(np.ma.log(arr)), lvls)

        # ----------------  Colorbars  ----------------
        if not no_cb:
            if isinstance(cb, type(None)):
                # if linear plot without fixcb set, then define here
                ax = plt.gca()

            # Create colourbar instance
            cb = plt.colorbar( poly, ax=ax, shrink=shrink, alpha=alpha, extend=extend )
            # set ylabel tick properties
            for t in cb.ax.get_yticklabels():
                t.set_fontsize(f_size)

            if not isinstance(units, type(None)):
                cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size)

            # Special treatment for log colorbars
            if log :
                round_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))
                tick_locs = [ float('{:.2g}'.format( t )) for t in lvls ]
                # for asectics, round colorbar labels to sig figs given
                for n, lvl in enumerate( lvls ):
                    try:
                        lvls[n] = round_to_n( lvl, sigfig_rounding_on_cb)
                    except:
                        lvls[n] = lvl
            else:
                tick_locs = np.array( lvls ).copy()

            # Turn tick locations into floats
            tick_locs = [float(_tick) for _tick in tick_locs]

            cb.set_ticks( np.array(tick_locs) )
            # the format is not correctly being set... - do this manually instead
            if not isinstance( format, type(None) ):
                lvls = [ format % (i) for i in lvls ]
            cb.set_ticklabels( lvls )#, format=format )

####################################################################################################

    if NEW_VERSION:

        if case==9 or default:
            vmin=fixcb_[0]
            vmax=fixcb_[1]

        if log:
            norm=LogNorm(vmin=fixcb_[0], vmax=fixcb_[1], cmap=cmap)
            lvls = np.logspace( np.log10(fixcb[0]), np.log10(fixcb[1]), num=nticks)
            # Normalise to Log space
            norm=mpl.colors.LogNorm(vmin=fixcb_[0], vmax=fixcb_[1])
            extend='min'
            ticks=lvls
            ax=m.ax

            rount_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n-1))
            tick_locs = [ float('{:.2g}'.format( t )) for t in lvls ]
            # for asectics, round colorbar labels to sig figs given
            for n, lvl in enumerate( lvls ):
                try:
                    lvls[n] = round_to_n( lvl, sigfig_rounding_on_cb)
                except:
                    lvls[n] = lvl
        else:
            tick_locs = np.array( lvls ).copy()

        # Turn tick locations into floats.
        tick_locs = [float(_tick) for _tick in tick_locs]

        # Create the colormap
        poly = m.pcolor(lon, lat, arr, norm=norm, vmax=vmax, vmin=vmin,
                cmap=cmap, alpha=alpha)

        # Add the colorbar if needed.
        if not no_cb:
            #if linear plot without fixcb set, then define here
            if isinstance(cb, type(None)):
                ax = plt.gca()

            cb = plt.colorbar(poly, ax=ax, ticks=lvls, format=format,
                    shrink=shrink, alpha=alpha, norm=norm, extend=extend)

            # Set ylabel tick properties
            cb.ax.tick_params(labelsize=f_size)

            # Set ylabel units:
            if not isinstance(units, type(None)):
                cb.ax.set_ylabel( units, rotation=rotatecbunits, labelpad=f_size)

            cb.set_ticks (np.array(tick_locs) )
            cb.set_ticklabels( lvls )









        #logging.info(tick_locs, lvls, [ type(i) for i in tick_locs, lvls ])
        #logging.info(cb.get_clim(), title, format)

# Set number of ticks
# FIX NEEDED - this currently doesn't doesn't work for log plots
#    if (not default) and (not no_cb) and ( not log ):
#        if set_cb_ticks:
#            tick_locator = ticker.MaxNLocator( nticks=nticks )
#            cb.locator = tick_locator
#            cb.update_ticks()
    # Add grid lines to the plot?
    plt.grid( grid )

    # Add title to plot?
    max_title_len=30

    if tight_layout==True:
        plt.tight_layout()

    if not isinstance( title, type(None) ):
        # Check if the title is too long and if not split it over lines
        if len(title)>max_title_len:
            print("tile takes up multiple lines. Splitting over lines now.")
            import textwrap
            title="\n".join(textwrap.wrap(title,max_title_len))
            print(title)
            # Adjust the top of the plot by 0.05 for every line the title takes
#            plt.subplots_adjust(top=1-0.05*(len(title)%max_title_len))

        plt.title(title, fontsize=f_size*1.5)
    # Setup list of return variables
    return_l = [ plt ]
    if not no_cb:
        return_l += [ cb ]
    if return_m:
        return_l += [ m ]
    return return_l



# --------
# X.XX - Zonal plot - log or linear
# --------
def zonal_plot( arr, fig, ax=None, title=None, tropics=False, f_size=10, c_off=37, \
        format='%.2f', interval=None, no_cb=False, units=None, shrink=0.4, alpha=1, \
        res='4x5', window=False, cmap=None, log=False, fixcb=None, fixcb_buffered=None, \
        xlimit =None, rotatecbunits='horizontal', extend='neither', ylabel=True, cb=None,\
        lvls=None, sigfig_rounding_on_cb=2, nticks=10, norm=None, set_window=False, \
        lat_0=None, lat_1=None, lat40_2_40=False, xlabel=True, mask_invalids=False, \
        trop_limit=True, verbose=True, debug=False, \
        # redundent?
        lower_limited=False, nlvls=25,
        ):
    """
    Creates a zonal plot from provide 2D array of longditude and latitude.

    Parameters
    ----------
    alhpa (float): transparency of plotter data
    arr (np.array): input (2D) array
    ax (axis instance): axis object to use
    c_off (int): integer cutoff for chemical troposphere
    cmap (str): force colormap selection by providing name
    debug (boolean): legacy debug option, replaced by python logging
    f_size (float): font size
    fixcb (np.array): minimium and maximum to fix colourbar
    fixcb_buffered (array): minimium and maximum to fix colourbar, with buffer space
    fig (figure instance): matplotlib figure instance
    format (str): format string for colorbar formating
    extend (str): colorbar format settings ( 'both', 'min', 'both' ... )
    interval (int): x/y tick interval in multiples of 15 degrees lat/lon
    lvls (list): manually provide levels for colorbar
    log (boolean): use a log scale for the plot and colorbar
    no_cb (boolean): include a coloubar?
    norm (norm object): normalisation to use for colourbar and array
    nticks (int): number of ticks on colorbar
    mask_invalids (boolean): mask invalid numbers (to allow saving to PDF)
    res (str): GEOS-Chem output configuration resolution ( '4x5' etc... )
    rotatecbunits (str): orientation of colourbar units
    shrink (boolean): colorbar size settings ( fractional shrink )
    set_window (boolean): set the limits of the plotted data (lat_0, lat_1)
    (for nested boundary conditions or subregion plots -)
    sigfig_rounding_on_cb (int): significant figure rounding to use for colourbar
    title (str): plot title (deafult is ==None, therefore no title)
    ylabel, xlabel (boolean): label x/y axis?
    units (str): units of given data for plot title
    verbose (boolean): legacy debug option, replaced by python logging
    wd (str): Specify the wd to get the results from a run.
    window (boolean): use window plot settings (fewer axis labels/resolution of map)

    Notes
    -----
     - This function will also apply maskes if set in arguments
     - Input resolution must be provide for non-default (4x5) output
    """
    logging.info( 'zonal plot called for arr of shape: {},{}, {}'.format(
        arr.shape, set_window, xlabel ) )

    # Kludge, for pcent arrays with invalid within them, mask for these.
    if mask_invalids:
        arr = np.ma.masked_invalid( arr )

    # Create axis if not provided <= is this back compatible?
    if isinstance( ax, type(None) ):
        ax = fig.add_subplot(111)

    # Get overall vars
    lon, lat, alt = get_latlonalt4res( res=res )
    alt  = [ i[:len(arr[0,:])] for i in [ alt ]  ][0]

    # === -  limit ars to a given mask, if stipulated
    if tropics:
        mask =  tropics_unmasked(res=res)
        arr = arr * mask[0,:,:c_off+1]
    if lat40_2_40  :
        mask =  tropics_unmasked(res=res)
        arr = arr * mask[0,:,:c_off+1]
#    logging.debug( [lat] + [ np.array(i).shape for i in lat, alt, arr, arr.T ] )

    if set_window:
        arr = arr[ get_gc_lat(lat_0, res=res):get_gc_lat(lat_1, res=res), :]
        lat = lat[ get_gc_lat(lat_0, res=res):get_gc_lat(lat_1, res=res) ]
        logging.debug( 'output post set_window: {}'.format( arr.shape)+
            'lens={}, res={}, mean={}, min={}, max={}'.format(
            [ len(i) for i in (lon,lat,alt)], res,  \
            [(i.mean(), i.min(), i.max()) for i in [arr[: ,:c_off]] ]) )
    min, max = [ (i.min(), i.max()) for i in [ arr[: ,:c_off]  ] ][0]

    # Limit cb to top of (GEOS-Chem chemical) troposphere
    alt = alt[:c_off+1]

    # Plot settings for window plots
    if window:
        if isinstance( interval, type(None) ):
            interval = 3
    else:
        interval = 1
    parallels = np.arange(-90,91,15*interval)

    # Is array reduced to chemistry computed troposphere? - if so limit alt
    if len(arr[0,:] ) != 38:
        arr =arr[:,:38]

    # -------- Colorbar/colomap variables...
    # Set cmap range I to limit poly, if not given cmap )
    fixcb_ = fixcb
    print(fixcb, fixcb_,  [ (i.min(), i.max()) for i in [arr ] ])
    if isinstance( fixcb_, type(None) ) :#and isinstance( cmap, type(None) ):
        fixcb_ = np.array( [ (i.min(), i.max()) for i in [arr ] ][0] )

    if isinstance( cmap, type(None) ):
        # lvls provided, if not calculate these.
        if isinstance( lvls, type(None) ):
            lvls = get_human_readable_gradations( vmax=fixcb_[1], vmin=fixcb_[0], \
                nticks=nticks, sigfig_rounding_on_cb=sigfig_rounding_on_cb  )

        # Setup Colormap
        cmap, fixcb_buffered = get_colormap( np.array( fixcb_ ), nticks=nticks, \
            fixcb=fixcb_, buffer_cmap_upper=True )

        # Update colormap with buffer
        cmap = get_colormap( arr=np.array( [fixcb_buffered[0],fixcb_buffered[1]] ) )
    # Allow function to operate without fixcb_buffered provided
    if isinstance( fixcb_buffered, type(None) ):
        fixcb_buffered = fixcb_
    fixcb_ = fixcb_buffered

    if verbose:
        print('colorbar variables: ', fixcb_buffered, fixcb, fixcb_, lvls, cmap)
    # -----------------  Log plots ---------------------------------------------
    # Create a Log Zonal plot of the given data,  ??? ( updated needed )
    # where fixcb is the min and max of a given array (aka not log )
    # TESTING NEEDED HERE !!! ( Update )
    if log:
        # Normalise to Log space
        if isinstance( norm, type(None) ):
            norm=mpl.colors.LogNorm(vmin=fixcb_[0], vmax=fixcb_[1])
        # Create poly collection
        poly = ax.pcolor( lat, alt, arr.T, norm=norm, cmap=cmap, alpha=alpha)

        if no_cb:
            pass
        else:
            # Get logarithmically spaced integers
            lvls = np.logspace( np.log10(fixcb[0]), np.log10(fixcb[1]), \
                                                num=nticks)
            # Make colorbar
            cb = plt.colorbar(poly, ax=ax, ticks=lvls, format=format, \
                 shrink=shrink, alpha=alpha, norm=norm, extend='min')

        logging.debug('cb min = {}, cb max = {}, lvls = {}'.format( \
            np.ma.min(np.ma.log(arr)), np.ma.max(np.ma.log(arr)), cb_lvls = lvls))

    # --------------  Linear plots -------------------------------
    # standard plot
    else:
        if verbose:
            print(fixcb, fixcb_,  \
                    [ np.array(i).shape for i in (lat, alt, arr, arr.T) ])
        # Create poly collection
        poly = ax.pcolor( lat, alt, arr.T, cmap=cmap, vmin=fixcb_[0],  \
                                vmax=fixcb_[1], norm=norm, alpha=alpha )

    # ----------------  Colorbars  ----------------
    if no_cb:
        pass
    else:
        if isinstance( cb, type(None) ):
            # if linear plot without fixcb set, then define here
            cb = plt.colorbar( poly, ax=ax, shrink=shrink, alpha=alpha,  \
                        format=format, ticks=lvls, norm=norm, \
                        extend=extend )
        # setup colobar labels/ticks
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(f_size)
        if not isinstance( units, type(None) ):
            cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size)

        if log:
            round_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))
            tick_locs = [ float('{:.2g}'.format( t )) for t in lvls ]
            # for asectics, round colorbar labels to sig figs given
            for n, lvl in enumerate( lvls ):
                try:
                    lvls[n] = round_to_n( lvl, sigfig_rounding_on_cb)
                except:
                    lvls[n] = lvl
        else:
            tick_locs = np.array( lvls ).copy()

        # fix colorbar levels, then provide labels
        cb.set_ticks( tick_locs )
        # the format is not correctly being set... - do this manually instead
        if not isinstance( format, type(None) ):
            lvls = [ format % (i) for i in lvls ]
        cb.set_ticklabels( lvls )#, format=format )

        logging.debug( 'variables post log plot setup: ', tick_locs, lvls, \
            [ type(i) for i in (tick_locs, lvls) ], cb.get_clim(), title, format )


    # Setup Y axis
    if (not isinstance( units, type( None) )) and (not no_cb):
        cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size)
    ax.set_ylim( alt[0], alt[-1])
    if ylabel:
        ax.set_ylabel('Altitude (km)', fontsize=f_size*.75)
    else:
        ax.tick_params( axis='y', which='both', labelleft='off')

    if trop_limit:
        ax.set_ylim( 0, 18 )
    if interval != 1:
        ax.set_yticks( ax.get_yticks()[::interval] )
    # Setup X axis
    ax.set_xticks( parallels )
    ax.tick_params( labelsize=f_size*.75 )
    if tropics:
        ax.set_xlim( -26, 26 )
    if lat40_2_40:
        ax.set_xlim( -40, 40 )
    else:
        ax.set_xlim( lat[0], lat[-1] )
    if not isinstance( xlimit, type(None) ):
        ax.set_xlim(xlimit[0], xlimit[1])

    if xlabel:
        ax.set_xlabel('Latitude', fontsize=f_size*.75)
    else:
        ax.tick_params( axis='x', which='both', labelbottom='off')

    # Y axis
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize( fontsize=f_size*.75 )

    # Add title if provided
    if (title != None):
        print([ type( i ) for i in (ax, plt) ])
        print(title)
        ax.set_title(title, fontsize=f_size*1.5)

# --------
# X.XX - Diurnal plot - (UPDATED)
# --------
def BASIC_diurnal_plot( fig=None, ax=None, dates=None, data=None, color='red',\
        title=None, label=None, plt_legend=None, plt_xlabel=True,\
        plt_ylabel=True, show_plt=False, plot_pcent_change_from_max=False, \
        units='ppbv', spec='O3', alt_text=None, loc='lower right', \
        filename2save='diurnal_plot.png', save_plt=False, show_plot=False,\
        stat2plot='50%', alpha = 0.3, time_resolution_str="%H",\
        add_quartiles2plot=True, return_avgs=True, \
        xlabel='Hour of day (UTC)', debug=False ):
    """
    Creates a diurnal plot for given data and dates

    Parameters
    -------
    data (array): numpy array of data
    dates (np.array of datetime.datetime): dates to use for x axis
    fig (fig instance)
    ax (ax instance)
    stat2plot (str): stat (e.g. mean or medium) to use for main line
    xlabel, ylabel (str): label for axis?
    label (str): label for input data
    units (str): unites to label
    time_resolution_str (str): time string to reduce dataframe by
    legend (boolean): add a legend?
    title (str): title for plot
    color (str/hex etc): color for line
    f_size (float): fontsize
    lgnd_f_size (float): fontsize for legend
    loc (str): location for legend
    rotatexlabel (numnber/str): rotation of x axis labels
    pos, posn (int): vestigle(!) location indices for window plots
    ylim (list): min and max y axis limit
    lw (str): linewidth
    return_avgs (boolean): return a list of the average values plotted
    save_plt, show_plot (boolean): show or save plots?
    xlabel (str): label for x axis

    Returns
    -------
    (None)
    """
    debug=True
    if debug:
        prt_str = 'BASIC_diurnal_plot called for {} (data.shape={})'
        print( prt_str.format( spec, data.shape ) )
    import matplotlib
    from matplotlib import dates as d
    import datetime as dt
    #
    if isinstance(fig, type(None)):
        fig = plt.figure()#dpi=Var_rc['dpi'])
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111)#2,2, n_season+1 )
    if debug:
        print( fig, ax)
    # --- Process
    # Form a dataFrame from input numpy arrays. (and remove NaNs... )
    raw_df = pd.DataFrame( {'data':data}, index=dates ).dropna()
    # Add a time coluumn to group by (e.g. "%H:%M" for minutes)
    raw_df['Time'] = raw_df.index.map(lambda x: x.strftime(time_resolution_str))
    df = raw_df.groupby('Time').describe().unstack()
    # Get the labels for time
    if debug:
        print( df.head(), df.index[:5], df.shape)
    time_labels = df['data'][stat2plot].index.values
    # --- Plot
    index = range(len(list(time_labels)))
    if debug:
        print( df['data'][stat2plot])
        print( '!'*20, index, time_labels)
    # Select data for the requested statistic
    avgs = df['data'][stat2plot]
    # Plot the % change from max ?
    if plot_pcent_change_from_max:
        if debug:
            print( avgs, avgs.shape)
        max_ = avgs.max()
        avgs = ( avgs - max_ ) / max_ *100
        if debug:
            print( avgs.shape, max_)
        units='%'
    # - Now plot up
    ax.plot( index, avgs, color=color, linewidth=2.0, label=label )
    # - Add quartiles
    if add_quartiles2plot:
        Q3_data = df['data']['75%']
        # Plot the % change from max ?
        if plot_pcent_change_from_max:
            max_ = Q3_data.max()
            Q3_data = ( Q3_data - max_ ) / max_ *100
        ax.fill_between(index, avgs, Q3_data, alpha=alpha, facecolor=color)
        Q1_data = df['data']['25%']
        # Plot the % change from max ?
        if plot_pcent_change_from_max:
            max_ = Q1_data.max()
            Q1_data = ( Q1_data - max_ ) / max_ *100
        ax.fill_between(index, avgs, Q1_data, alpha=alpha, facecolor=color)
    # Label xticks on axis (if 24 plots, just label every 3rd)
    if index < 6:
        ax.set_xticks(index)
        ax.set_xticklabels(time_labels)
    else:
        ax.set_xticks(index[1::3])
        ax.set_xticklabels(time_labels[1::3])
    xticks = ax.get_xticks()
    if debug:
        print( xticks, ax.get_xticklabels() )
#    ax.set_xticks(np.linspace(3, 21, 7).astype(int))
    # More cosmetic changes...
    if not isinstance(title, type(None)):
        plt.title(title)
    if plt_legend:
        plt.legend( loc=loc )
    if plt_xlabel:
        plt.xlabel(xlabel)
    else:
        ax.tick_params( axis='x', which='both', labelbottom='off')
    if plt_ylabel:
        try:
            spec_ = latex_spec_name(spec)
        except:
            spec_ = spec
        plt.ylabel('{} ({})'.format( spec_, units ) )
    else:
        ax.tick_params( axis='y', which='both', labelright='off')
    # Add alt text if provided.
    if not isinstance(alt_text, type(None)):
        alt_text_y = 0.925
        alt_text_x = 0.975
        plt.text(  alt_text_x, alt_text_y, alt_text, ha='right', va='center', \
                   transform=ax.transAxes )
    # Print or show plot?
    if save_plt: plt.savefig( filename2save  )
    if show_plt: plt.show()
    # Return the average values?
    if return_avgs:
        return np.ma.array(avgs)

# --------
# X.XX - BASIC_seasonal_plot
# --------
def BASIC_seasonal_plot( dates=None, data=None, ax=None,
        title=None, legend=False,  ylabel=None, loc='upper left',  \
        showmeans=False, boxplot=True, return_avgs=False, \
        plt_median=False, plot_Q1_Q3=False,
        ylim=None, xtickrotation=45, alt_text=None, alt_text_x=.925,
        alt_text_y=.925, xlabel=True, rm_yticks=False, log=False, \
        pcent1=25, pcent2=75, color='red' , lw=1, ls='-', label=None,
        debug=False ):
    """ Plot up a basic seasonal plot - adapted from AC_tools """
    if debug:
        print( locals().keys(), ax )
    # get months
    df = pd.DataFrame(data,index=dates )
    # Force use of a standard year
    df['months'] = [i.month for i in dt64_2_dt( df.index.values ) ]
    df.sort_values(by='months', ascending=True, inplace=True)
#    months = list(range(1, 13))
    months = df['months'].values
    df = df.drop('months', axis=1)
    datetime_months = [ datetime.datetime(2009, int(i), 1) for i in set(months)]
    labels = [i.strftime("%b") for i in datetime_months]
    # Get Data by month
    monthly = df.groupby(df.index.month)
    months = list( sorted( monthly.groups.keys() ) )
    monthly = [monthly.get_group(i) for i in months ]
#    [ df[df.index.month==i]  for i in months ]
    # remove nans to allow for percentile calc.
    monthly = [ i.dropna().values  for i in monthly ]
    data_nan = [ i.flatten() for i in monthly ]
    # Plot up median
    medians =  [ np.nanpercentile( i, 50, axis=0 ) for i in data_nan ]
    plt.plot( months, medians, color=color, lw=lw, ls=ls, label=label  )
    # Plot quartiles as shaded area?
    low = [ np.nanpercentile( i, pcent1, axis=0 ) for i in data_nan ]
    high = [ np.nanpercentile( i, pcent2, axis=0 ) for i in data_nan ]
    if isinstance(ax, type(None)):
        ax = plt.gca()
    ax.fill_between( months, low, high, alpha=0.2, color=color   )
    # Beatify plot
    ax.set_xticks( months )
    if xlabel:
        ax.set_xticklabels( labels, rotation=xtickrotation )
    else:
        ax.tick_params( axis='x', which='both', labelbottom='off')
    if debug:
        print('!'*50, alt_text, alt_text_x, alt_text_y)
    if not isinstance( alt_text, type(None) ):
        if debug:
            print('!'*50, alt_text, alt_text_x, alt_text_y)
        plt.text( alt_text_x, alt_text_y, alt_text, ha='right', va='center', \
            transform=ax.transAxes )
    if legend:
        if debug:
            print('>'*500, 'Adding legend', '<'*50, loc)
        plt.legend( loc=loc)#fontsize=f_size*.75, loc=loc )
    if not isinstance( title, type(None) ):
        plt.title( title )
    if not isinstance( ylabel, type(None) ):
        plt.ylabel( ylabel )
#        , fontsize=f_size*0.75 ) # Why is this x0.75?
    else:
        if rm_yticks:
            ax.tick_params( axis='y', which='both', labelleft='off')
    # Log scale?
    if log:
        ax.set_yscale('log')
    else:
        ax.set_yscale('linear')
    if return_avgs:
        return np.ma.array(medians)


# --------
# X.XX - Plot up Sonde data
# --------
def sonde_plot(fig, ax, arr, n=0, title=None, subtitle=None, \
        f_size=10, color=None,  err_bar=False, obs=True, \
        legend=False, units='nmol mol$^{-1}$', stddev=True, \
        c_l=[ 'k', 'red','green', 'blue' , 'purple'], xlimit=None, \
        loc='upper left', c_off=37, label=None, ancillary=True, \
        plt_txt_x=0.5, plt_txt_y=0.94, \
        ylabel=True, xlabel=True, hPa_labels=True,
        debug=False,
        # redundant arguments?
        rasterized=True, tropics=False,
        ):
    """
    Create plot of vertical data for sonde observations/model


    Parameters
    -------
    hPa_labels (boolean): include labels for hPa on axis?
    plt_txt_x, plt_txt_y (float): ordinal locations for alt text
    ancillary (boolean): add ancillary labels etc. to plot
    label (str): label to use for line / legend.
    c_off (int): index to plot model levels too (e.g. 37 for trop. )
    xlimit (float): value to limit x axis to (e.g. 100 ppbv)
    c_l (list):colors to use (index by n )
    uints (str): units to use for axis labels
    obs (boolean): overide plot settings with those for observations.
    err_bar (boolean): apply bar to show quartile range of plots
    stddev (boolean): as above (see err_bar)
    color (str/color instance): color for plotted line
    f_size (float): fontsise
    tropics (boolean): mask for tropics?
    title (str): title for plot
    subtitle (str): subtitle for plot (e.g. location)
    n (int): number of plot within window plot figure
    fig (figure instance): fig. to use
    ax (axis instance): axis to use
    ylabel, xlabel (boolean): include axis labels for plot?


    Returns
    -------

    Notes
    -----
     -  should be re-written (age >3 years) to work with updated AC_tools
    """

    # Get overall vars
    alt, press = [ gchemgrid(i) for i in ('c_km_geos5_r' , 'c_hPa_geos5_r')]

    # Cut off at Tropopause?
    if obs:
        arr = arr[:c_off,:]
    else:
        arr = arr[:c_off]
    alt, press = [ i[:c_off] for i in (alt, press) ]

    # if color  not set, get color
    if isinstance( color, type(None) ):
        color=c_l[n]

    if obs:
#        print [len(i) for i in arr[:,0], alt ]
        ax.plot( arr[:,0] , alt , color=color, label=label )
        # limit cb to top of troposphere
        min, max  =  [ (np.ma.min(i), np.ma.max(i)) for i in [ arr[:,0] ] ][0]
    else:
        ax.plot( arr , alt,  label=label, color=color )
        # limit cb to top of troposphere
        min, max  =  [ (np.ma.min(i), np.ma.max(i)) for i in [ arr ] ][0]

    # Add legend
    if legend:
        plt.legend( loc=loc, fontsize=f_size*.75)

    # Sonde data  = mean, 1st and 3rd Q
    if err_bar:
        if stddev:
            ax.errorbar(arr[:,0] , alt, xerr=arr[:,3], fmt='o', color=color,
                    elinewidth=0.75, markersize=2, alpha=0.75 )
        else:
            ax.errorbar(arr[:,0] , alt, xerr=[arr[:,1], arr[:,2]], fmt='o', \
                color=color, elinewidth=.75, markersize=5, alpha=0.75 )

    # Beautify plot ( e.g. add hPa, units,  etc... )
    if ancillary:
        if not isinstance( title, type(None) ):
            plt.title( title, fontsize = f_size, y=1.0)
            if not isinstance( subtitle, type(None) ):
                plt.text(  plt_txt_x, plt_txt_y, \
                    subtitle, ha='center', va='center', \
                    transform=ax.transAxes, fontsize=f_size*.65)

        if ylabel:
            plt.ylabel('Altitude (km)', fontsize=f_size*.75)
        else:
            plt.tick_params( axis='y', which='both', labelleft='off')
        if xlabel:
            plt.xlabel( units, fontsize=f_size*.75)
        else:
            plt.tick_params( axis='x', which='both', labelbottom='off')

        if xlimit == None:
    #        plt.xlim( min-(min*0.02), max+(max*0.02) )
            plt.xlim( 0.1,  max+(max*0.50) )
        else:
            plt.xlim(0.1, xlimit)
        if hPa_labels:
            ax2 = ax.twinx()
            press = [ myround(i,100) for i in press ][::-1]
            ax2.set_yticks( press[::10] )
            ax2.set_yticklabels( press[::10] )
            majorFormatter = mpl.ticker.FormatStrFormatter('%d')
            ax2.yaxis.set_minor_formatter( majorFormatter )
            ax2.set_ylabel( 'Press. (hPa)', fontsize = f_size*.75)
            ax2.invert_yaxis()
            if ylabel:
                pass
            else:
                ax.tick_params( axis='y', which='both', labelleft='off')
                ax2.tick_params( axis='y', which='both', labelleft='off')



# --------------
# 1.13 - plot up monthly from data provided from DB netCDF
# -------------
def monthly_plot( ax, data, f_size=20, pos=0, posn=1, lw=1,ls='-', color=None, \
        title=None, subtitle=None, legend=False, xrotation=90, \
        window=False, label=None, ylabel=None, xlabel=True, \
        title_loc_y=1.09, plt_txt_x=0.5, plt_txt_y=1.05, \
        plot_Q1_Q3=False, low=None, high=None, loc='upper right' ):
    """
    Plot up seaonal (monthly ) data.

    NOTES:
     - Requires data, and dates in numpy ararry form.
     - Dates must be as datetime.datetime objects.
    """

    # setup color list if not provided
    if isinstance(color, type(None)):
        color = color_list( posn )[ pos ]

    # if this is a window plot, then reduce text size
    if window:
        f_size=int(f_size/2)

    # Plot up provide monthly data
    plt.plot( np.arange(1,len(data)+1), data, color=color, lw=lw, ls=ls, \
        label=label )

    # Also add 5th  and  95th %ile
    if plot_Q1_Q3: # Plot quartiles as shaded area?
            ax.fill_between( np.arange(1,len(data)+1), low, high, alpha=0.2, \
                    color=color   )

    # Beautify
    ax.set_xticklabels( [i.strftime("%b") \
            for i in [datetime.datetime(2009, int(i), 0o1)  \
            for i in np.arange(1,13 ) ] ] )
    plt.xticks(list(range(1,13)),  fontsize=f_size )
    plt.xticks(rotation=xrotation)
    if not xlabel:
        plt.tick_params( axis='x', which='both', bottom='on', top='off',
                                    labelbottom='off')

    plt.xlim(0.5,12.5)
    if ylabel != None:
        plt.ylabel(  ylabel, fontsize=f_size  )
    plt.yticks( fontsize=f_size )
    print(pos, posn-1)
    if not isinstance(title, type(None)):
        t = plt.title(title, fontsize = f_size)
        t.set_y( title_loc_y )
        if not isinstance(subtitle, type(None)):
             plt.text( plt_txt_x, plt_txt_y, subtitle, ha='center', \
                va='center', transform=ax.transAxes, fontsize=f_size*.65)
#    if pos==posn-1:
    if legend:
        plt.legend( loc=loc,  fontsize=int(f_size/1.5) )


# --------------
# 1.14 - plot up monthly timeseries from ...
# -------------
def timeseries_seasonal_plot( ax, dates, data, f_size=20, pos=0, posn=1,  \
        title=None, legend=False, everyother=24,  x_nticks=12, \
        window=False, label=None, ylabel=None, loc='upper left',  \
        lw=1,ls='-', color=None, showmeans=False, boxplot=True, \
        plt_median=False, plot_Q1_Q3=False, pcent1=25, pcent2=75, \
        ylim=None, xtickrotation=45, alt_text=None, alt_text_x=.5,
        alt_text_y=.5, xlabel=None, rm_yticks=False, log=False, \
        debug=False ):
    """
    Plot up timeseries of seasonal data.

    Parameters
    ----------
    ax (axis object): axis to plot onto
    dates (np.array): array of dates (as  datetime.datetime objects)
    data (np.array): 1D array of data
    plot_Q1_Q3 (boolean): plot up quartiles on timeseries plot

    Returns
    -------
    (None)
    """
    # Process data - reduce resolution to daily, and get std
    df = DataFrame(data,index=dates )
    # Force use of a standard year
    months = list(range(1, 13))
    datetime_months = [ datetime.datetime(2009, int(i), 1) for i in months ]
    labels = [i.strftime("%b") for i in datetime_months]
    if debug:
        print(labels)

    # Get Data by month
    monthly =[ df[df.index.month==i]  for i in months ]
    if debug:
        print([i.shape for i in monthly ])

    if boxplot:
        bp = ax.boxplot( monthly, months, showmeans=showmeans )
    else:
        # remove nans to allow for percentile calc.
        data_nan = [ np.ma.array(i).filled( np.nan ) for i in monthly ]
        data_nan = [ i.flatten() for i in data_nan ]

        if plt_median:
            plt.plot( months, \
                [ np.nanpercentile( i, 50, axis=0 ) for i in data_nan ], \
                color=color, lw=lw, ls=ls, label=label  )
            if plot_Q1_Q3: # Plot quartiles as shaded area?
                low = [ np.nanpercentile( i, pcent1, axis=0 ) \
                    for i in data_nan ]
                high = [ np.nanpercentile( i, pcent2, axis=0 ) \
                    for i in data_nan ]
                ax.fill_between( months, low, high, alpha=0.2, \
                    color=color   )
        else:
            plt.plot( months, [i.mean() for i in monthly ], color=color,
            lw=lw, ls=ls, label=label )

    # Beatify plot
    ax.set_xticks( months )
    if xlabel:
        ax.set_xticklabels( labels, rotation=xtickrotation )
    else:
        ax.tick_params( axis='x', which='both', labelbottom='off')
    if not isinstance( ylim, type(None) ):
        ax.set_ylim( ylim )

    if debug:
        print('!'*50, alt_text, alt_text_x, alt_text_y)
    if not isinstance( alt_text, type(None) ):
        if debug:
            print('!'*50, alt_text, alt_text_x, alt_text_y, f_size)
        plt.text(  alt_text_x, alt_text_y, \
            alt_text, ha='center', va='center', \
                    transform=ax.transAxes, fontsize=f_size*.5)

    if legend:
        if debug:
            print('>'*500, 'Adding legend', '<'*50, loc)
        plt.legend( fontsize=f_size*.75, loc=loc )
    if not isinstance( title, type(None) ):
        plt.title( title )
    if not isinstance( ylabel, type(None) ):
        plt.ylabel( ylabel, fontsize=f_size*0.75 ) # Why is this x0.75?
    else:
        if rm_yticks:
            ax.tick_params( axis='y', which='both', labelleft='off')
    # Log scale?
    if log:
        ax.set_yscale('log')
    else:
        ax.set_yscale('linear')

# --------------
# X.XX - plot up daily timeseries from ...
# -------------
def timeseries_daily_plot(fig, ax,  dates, data, pos=1, posn =1,  \
        bin_size=1/24.,widths=0.01, rotatexlabel = 'vertical', \
        white_fill=True, alpha=0.1, linewidth=0.5,xlabel=True, \
        title=None, alt_text=None, f_size=7.5, units='ppbv',
        showmeans=False, plt_median=False, boxplot=True,
        ylabel=True, color='blue', label=None, \
        plot_Q1_Q3=True, pcent1=25, pcent2=75,  debug=False ):
    """
    Plot up daily timeseries of values. Requires data, and dates in numpy
        array form. Dates must be as datetime.datetime objects.

    Notes:
     - Boxplot is the default presentation of data
     - Otherwise a meidan is used
    """

    # get_day_fraction(i)
    dates = np.array( [ get_day_fraction(i) for i in dates ] )

    # bin data
    b_all, bins_used = bin_data( data, dates, bin_size, debug=debug )

    # Plot
    if boxplot:
        bp = ax.boxplot( b_all,  positions=bins_used, widths=widths, \
            showmeans=showmeans, patch_artist=True )
    else:
        # remove nans to allow for percentile calc.
#        print b_all
        data_nan = [ np.ma.array(i).filled( np.nan ) for i in b_all ]
        data_nan = [ i.flatten() for i in data_nan ]

        # Plot average
        if plt_median:
            ln = plt.plot( bins_used, \
                [ np.nanpercentile( i, 50, axis=0 ) for i in data_nan ], \
                color=color, label=None   )
        else:
            ln = plt.plot( bins_used, [i.mean() for i in data_nan], color=color)

        # Plot quartiles as shaded area?
        if plot_Q1_Q3:

            low = [ np.nanpercentile( i, pcent1, axis=0 ) for i in data_nan ]
            high = [ np.nanpercentile( i, pcent2, axis=0 ) for i in data_nan ]
            ax.fill_between( bins_used, low, high, alpha=0.2, color=color   )

    # Beautify
    if not isinstance( title, type(None) ):
        plt.title(title, fontsize=f_size )
    if not isinstance( alt_text, type(None) ):
#        ax.text(x=0.85,y=0.85, s=alt_text, fontsize=f_size*1.5 )
        ax.annotate( alt_text , xy=(0.85, 0.85),  textcoords='axes fraction')

    ax.set_xticklabels( np.arange(0,24,1 )  )
    plt.xticks( np.arange(0,1,1/24. ), fontsize=f_size*.75, \
        rotation=rotatexlabel )
#    plt.xlim(-0.05, 23/24.)
    plt.xlim(-0.05, 1.0)
    if xlabel:
        ax.set_xlabel('Hour of day', labelpad=f_size)
    else:
        ax.tick_params( axis='x', which='both', labelbottom='off')
    # Setup y axis
    if ylabel:
        ax.set_ylabel('{}'.format(units), labelpad=f_size)
    else:
        ax.tick_params( axis='y', which='both', labelleft='off')

    # --- Highlight bins
    bs = np.arange(0, 24, bin_size )#[ bs[0] - bin_size ] + bs
    [ plt.axvline( x=i, color='k', linewidth=linewidth, alpha=alpha,  \
                            linestyle='dashed' ) for i in bs ]




# --------------
# X.XX - plot up timeseries from between two given months
# -------------
def timeseries_month_plot( ax, dates, data, f_size=20, pos=0, posn=1,  \
        title=None, legend=False, everyother=7*24,  x_nticks=12, \
        window=False, label=None, ylabel=None, loc='upper left',  \
        lw=1,ls='-', color=None, start_month=7, end_month=7, \
        boxplot=True, showmeans=False, alt_text=None, r_plt=False, \
        unitrotation=45, color_by_z=False, fig=None,  xlabel=True, \
        second_title='', add_dates2title=True, positive=None, debug=False ):
    """
    Plot up month timeseries of values. Requires data, and dates in numpy
    array form. Dates must be as datetime.datetime objects.

    NOTE(s):
     - This just plot up timeseries, why is the name  "timeseries_*month*_plot"?
     - Update this!
    """

    # Process data - reduce resolution to daily, and get std
    df = DataFrame( data, index=dates )

    # remove dates outside of range (start_month > t  < end_month )
    def get_month(x):
        return x.month
    df[ 'month'] = df.index.map( get_month )
    df = df[ df['month']<=end_month ]
    df = df[ df['month']>=start_month ]
    # remove 'month' column
    df = df.drop('month', 1)

    # label once per week
    days = [i.to_datetime() for i in df.index ]
    labels = [i.strftime("%-d %b") for  i in days ][::everyother]

    # Color in line another provided variables
    if color_by_z:
        if debug:
            print('Coloring line by normalised z values')
        print(df.columns)
        x = df.index
        y, z = [ df[ df.columns[i] ] for i in range(2) ]
        cmap = get_colormap( z.copy(), positive=positive )
        print([ ( i.min(), i.max() ) for i in (x, y, z) ])
        colorline(x, y, z, cmap=cmap, linewidth=lw, ax=ax, \
            norm=plt.Normalize( 0, 360 ), fig=fig )
    else:
        plt.plot( days, df.values, label=label, color=color, ls=ls, lw=lw  )

    # set xticks
    if xlabel:
        plt.xticks( days[::everyother], labels, rotation=unitrotation )
    else:
        plt.tick_params( axis='x', which='both', labelbottom='off')

    # Beatify plot
    if not isinstance( title, type(None) ):
        if add_dates2title:
            title += ' for {}-{} {}'.format( num2month(start_month),\
                num2month(end_month), second_title  )
        plt.title( title)
    if not isinstance( alt_text, type(None) ):
        plt.figtext(x=0.05,y=0.85, s=alt_text, fontsize=f_size*.75 )
    if not isinstance( ylabel, type(None) ):
        plt.ylabel( ylabel, fontsize=f_size*.75 )
    if legend:
        plt.legend( fontsize=f_size*.75, loc=loc )


# --------
# X.XX - North Pole surface plot
# --------
def north_pole_surface_plot( arr, return_m=False, grid=True, centre=False, \
        cmap=None, format='%.2f', m=None, fixcb=False,  nbins=25, \
        res='4x5', ax=None, alpha=1, nticks=10, everyother=1,\
        drawcountries=True, set_cb_ticks=True, title=None, \
#             rotatecbunits='horizontal', extend='neither',
             interval=1, resolution='l', shrink=0.4, window=False, \
        lon_0=0, boundinglat=40, degrade_resolution=False, \
        no_cb=False, cb=None, units=None, f_size=20, \
        debug=False,  **Kwargs):
    """
    Plot up data at north pole as a  2D slice.

    NOTES:
     - Requires data (arr) as numpy array. Arr should be full global size
    (lon, lat) for given resolution
    """

    # ---- Grid/Mesh values for Lat, lon, & alt + cb
    if isinstance( cmap, type(None) ):
        cmap = get_colormap( arr.copy() )

    if debug:
        print('>'*5, [ [ i.min(), i.max(), i.mean(), type(i) ] for i in [arr] ])
    lon, lat, NIU = get_latlonalt4res( res, centre=centre )

    # restrict lat to projection
    lat = lat[ get_gc_lat(boundinglat, res=res): ]

    lon, lat = np.meshgrid(lon, lat)

    # ----------------  Basemap setup  ----------------
    if isinstance( m, type(None) ):
        m = Basemap(projection='npstere', boundinglat=boundinglat,lon_0=lon_0,\
            resolution=resolution, round=True)

        # Beautify basemap plot
        m.drawcoastlines()
        m.drawcountries()
        parallels = np.arange(-90,91,20*interval)
#    meridians = np.arange(-180,181,20*interval)
        meridians = np.arange(-180,180,20*interval)
        m.drawparallels( parallels, labels=[False,False,False,False],\
            fontsize=f_size*.25 )
        m.drawmeridians(meridians, labels=[True,False,True,True],\
            fontsize=f_size*.25 )

    # set x and y
    x,y = m( lon, lat )
    if debug:
        print(1, len(x), len(y))

    if debug:
        print(2, 'len:',  [ len(i) for i in (x,y,lat,lon) ])
        print('>'*5, [ [ i.min(), i.mean(), i.max(), i.shape ] for i in (arr,
            arr[:,get_gc_lat(boundinglat,res=res):])  ])

    # -------- colorbar variables...
    # set cb label sizes
    if not no_cb:
#        if fixcb:
        if not isinstance( fixcb, type(None) ):
            tickmin, tickmax = fixcb[0], fixcb[1]
        else:
            tickmin, tickmax  = arr.min(), arr.max()

    # -----------------------  Linear plots -------------------------------
#    plt.contourf( x, y, arr[:,get_gc_lat(boundinglat,res=res):].T, alpha=alpha)
#    if fixcb:
    polar_arr = arr[:,get_gc_lat(boundinglat, res=res):].T
    if debug:
        print([ (i.min(), i.max(), i.mean()) for i in [polar_arr]  ])
    if not isinstance( fixcb, type(None) ):
        poly = m.pcolor( x, y, polar_arr, cmap=cmap, alpha=alpha,  \
            vmin=fixcb[0], vmax=fixcb[1] )
    else:
        poly = m.pcolor( x, y, polar_arr, cmap=cmap, alpha=alpha )


    # --------------  Log plots ---------------------------------------------


    # ----------------  Colorbars  ----------------
    if not no_cb:
        if isinstance( cb, type(None) ):
            cb = plt.colorbar( poly, ax = m.ax, shrink=shrink, alpha=alpha,  \
                        format=format, ticks= np.linspace(tickmin, tickmax, \
                        nticks ) )
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(f_size)
        if units != None:
            cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size)

    # Set number of ticks
        if set_cb_ticks:
            tick_locator = mpl.ticker.MaxNLocator(nticks=nticks)
            cb.locator = tick_locator
            cb.update_ticks()

    if not isinstance( title, type(None) ):
        plt.title(title, fontsize=f_size*1.5)

    return_list = [ poly ]
    if not no_cb:
        return_list += [ cb ]
    if (return_m):
        return_list += [ m ]
    return return_list


# --------
# X.XX - South Pole surface plot
# --------
def south_pole_surface_plot( arr, return_m=False, grid=True, centre=False,
        cmap=None, format='%.2f', res='4x5', ax=None, alpha=1,  title=None,
        fixcb=False, nbins=25, nticks=10, drawcountries=True, set_cb_ticks=True,
#             rotatecbunits='horizontal', extend='neither',
        interval=1, resolution='l', shrink=0.4, window=False, everyother=1,
        lon_0=0, boundinglat=-40, degrade_resolution=False, no_cb=False, cb=None,
        units=None, f_size=20, debug=False,  **Kwargs):
    """
    Plot up data at south pole 2D slice.

    NOTES:
     - Requires data (arr) as numpy array. Arr should be full global size
    (lon, lat) for given resolution
    """

    # ---- Grid/Mesh values for Lat, lon, & alt + cb
    if isinstance( cmap, type(None) ):
        cmap = get_colormap( arr.copy() )

    if debug:
        print('>'*5, [ [ i.min(), i.max(), i.mean(), type(i) ] for i in [arr] ])
    lon, lat, NIU = get_latlonalt4res( res, centre=centre )

    # restrict lat to projection
    lat = lat[ :get_gc_lat(boundinglat, res=res) +2 ]

    lon, lat = np.meshgrid(lon, lat)

    # ----------------  Basemap setup  ----------------
    m = Basemap(projection='spstere', boundinglat=boundinglat,lon_0=lon_0,
        resolution=resolution, round=True)

    # set x and y
    x,y = m( lon, lat )
    if debug:
        print(1, len(x), len(y))

    # Beautify plot
    m.drawcoastlines()
    m.drawcountries()
    parallels = np.arange(-90,91,20*interval)
#    meridians = np.arange(-180,181,20*interval)
    meridians = np.arange(-180,180,20*interval)
    m.drawparallels( parallels, labels=[False,False,False,False], \
        fontsize=f_size*.25 )
    m.drawmeridians(meridians, labels=[True,False,True,True], \
        fontsize=f_size*.25 )

    if debug:
        print(2, 'len:',  [ len(i) for i in (x,y,lat,lon) ])
        print('>'*5, [ [ i.min(), i.mean(), i.max() ] for i in [arr ] ])

    # -------- colorbar variables...
    # set cb label sizes
#    if fixcb:
    if not isinstance( fixcb, type(None) ):
        tickmin, tickmax = fixcb[0], fixcb[1]
    else:
        tickmin, tickmax  = arr.min(), arr.max()

    # --------------------  Linear plots -------------------------------
#    plt.contourf( x, y, arr[:,get_gc_lat(boundinglat,res=res):].T, alpha=alpha)
#    if fixcb:
    polar_arr = arr[:,:get_gc_lat(boundinglat, res=res)+2].T
    if debug:
        print([ (i.min(), i.max(), i.mean()) for i in [polar_arr]  ])
    if not isinstance( fixcb, type(None) ):
        poly = m.pcolor( x, y, polar_arr, cmap=cmap, alpha=alpha, \
            vmin=fixcb[0], vmax=fixcb[1] )
    else:
        poly = m.pcolor( x, y, polar_arr, cmap=cmap, alpha=alpha )

    # -----------  Log plots ---------------------------------------------


    # ----------------  Colorbars  ----------------
    if not no_cb:
        if isinstance( cb, type(None) ):
            cb = plt.colorbar( poly, ax = m.ax, shrink=shrink, alpha=alpha,  \
                        format=format, ticks= np.linspace(tickmin, tickmax, \
                        nticks ) )
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(f_size)
        if units != None:
            cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size)
    # Set number of ticks
    if  (not no_cb):
        if set_cb_ticks:
            tick_locator = mpl.ticker.MaxNLocator(nticks=nticks)
            cb.locator = tick_locator
            cb.update_ticks()

    if not isinstance( title, type(None) ):
        plt.title(title, fontsize=f_size*1.5)


    if (return_m):
        return plt , cb, m
    elif no_cb:
        return plt
    else:
        return plt , cb


# --------
# X.XX - PDF of monthly surface change plots for given species (takes 5D arr )
# --------
def plot_specs_surface_change_monthly2pdf( arr, res='4x5', dpi=160, \
        no_dstr=True, f_size=20, pcent=True, specs=None, dlist=None, \
        savetitle='', diff=False, extend='neither', column=False, \
        scale=1, units=None, set_window=False, lat_0=None, lat_1=None, \
        mask_invalids=False, debug=False):
    """
    Create multipage PDF with each page containing a 2D (lon,lat) slice
    plot for given species in list of "specs"
    Takes 5D array  ( species ,lon , lat, alt, time)
    """
    logging.info(  'plot_specs_surface_change_monthly2pdf called' )

    # setup pdfs + titles
    if column:
        savetitle = 'Column_by_spec'+savetitle
    else:
        savetitle = 'Surface_by_spec'+savetitle
    pdff = plot2pdfmulti( title=savetitle, open=True, dpi=dpi, no_dstr=no_dstr )
    left=0.05; right=0.9; bottom=0.05; top=0.875; hspace=0.315; wspace=0.1

    # Loop species
    for n, spec in enumerate( specs ):
        if debug:
            print(n, spec)

        # Get units/scale for species + setup fig (allow of
        if isinstance( units, type( None) ):
            if column and (not pcent):
                units, scale = 'DU', 1
            elif pcent:
                units, scale = '%', 1
            else:
                units, scale = tra_unit( spec, scale=True, global_unit=True )

        # setup masking...
        cbarr = arr[n,:,:,0,:].copy() * scale

        if pcent:
            mask_invalids = True

            # mask for changes greater than 500%
            if len( cbarr[ cbarr>500 ]  ) > 0:
                cbarr = np.ma.masked_where( cbarr>500, cbarr )
                extend='max'

            elif len( cbarr[cbarr<-500])>0:
                cbarr = np.ma.masked_where( cbarr<-500, cbarr )
                if  extend == 'max':
                    extend = 'both'
                else:
                    extend = 'min'

            else:
                extend='neither'

        else:
            extend='neither'
            cbarr = cbarr

        # Set the correct title
        ptitle = '{}'.format( latex_spec_name(spec) )
        if column:
            ptitle += ' column'
        else:
            ptitle += ' surface'

        if diff:
            ptitle += ' $\Delta$ concentration'
        else:
            ptitle += ' concentration'
        ptitle += ' ({})'.format( units )

        fig  = plt.figure(figsize=(22, 14), dpi=dpi, facecolor='w', edgecolor='w')

        # set cb ranges for whiole data period
        fixcb  = [( i.min(), i.max() ) for i in [cbarr]][0]

            # Kludge, force max cap at .2
#            if units == 'ratio':
#                fixcb = [ fixcb[0], 0.2 ]
#                extend = 'max'
        if (units == 'pmol mol$^{-1}$ m$^{-3}$') and ( spec == 'AERI' ):
            fixcb = [ fixcb[0], 500 ]
            extend = 'max'

        cmap = get_colormap( fixcb  )

        # Loop thorugh months
        for m, month in enumerate(dlist):
            fig.add_subplot(4,3,m+1)

            # plot up spatial surface change
            map_plot( arr[n,:,:,0,m].T*scale, cmap=cmap, case=9, res=res,\
                no_cb=True, f_size=f_size, fixcb=fixcb, window=True,\
                set_window=set_window,lat_0=lat_0, lat_1=lat_1, \
                mask_invalids=mask_invalids, clevs=clevs, debug=debug)
            plt.title(month.strftime("%b"), fontsize=f_size*2)

        # Add single colorbar
        mk_cb(fig, units=units, left=0.9,  cmap=cmap,vmin=fixcb[0],
            vmax=fixcb[1], f_size=f_size, extend=extend )

        # sort out ascetics -  adjust plots and add title
        fig.subplots_adjust( bottom=bottom, top=top, left=left, \
            right=right,hspace=hspace, wspace=wspace)
        fig.suptitle( ptitle, fontsize=f_size*2, x=.55 , y=.95  )

        # save out figure
        plot2pdfmulti( pdff, savetitle, dpi=dpi, no_dstr=no_dstr )

        # close fig
        plt.clf()
        plt.close()
        del fig

    #  save entire pdf
    plot2pdfmulti( pdff, savetitle, close=True, dpi=dpi, no_dstr=no_dstr )


# --------
# X.XX - PDF of monthly zonal change plots for given species
# --------
def plot_specs_zonal_change_monthly2pdf( Vars, res='4x5', dpi=160, \
        no_dstr=True, f_size=20, pcent=False, specs=None, dlist=None, \
        t_ps=None, savetitle='', diff=False, extend='neither',  \
        set_window=False, lat_0=None, lat_1=None, mask_invalids=False, \
        set_lon=None, units=None, debug=False):
    """
    Create multipage PDF with each page containing a zonalplot for given species in
    list of "specs"
    NOTES:
     - Takes 5D array  ( species ,lon , lat, alt, time)
     - Needs descriptions update.
    """

    savetitle = 'Zonal_by_spec'+savetitle
    pdff = plot2pdfmulti( title=savetitle, open=True, dpi=dpi, no_dstr=no_dstr )
    left=0.05; right=0.9; bottom=0.05; top=0.875; hspace=0.315; wspace=0.2

    # Loop species
    for n, spec in enumerate( specs ):
        if debug:
            print(n, spec, Vars.shape)

        # Get units/scale for species + setup fig
        scale=1
        if isinstance( units, type(None) ):
            units, scale = tra_unit( spec, scale=True, global_unit=True )
        if pcent:
            units, scale = '%', 1
            mask_invalids=True

        # Set the correct title
        ptitle = '{}'.format( latex_spec_name(spec) )
        if diff:
            ptitle += ' zonal $\Delta$ concentration'
        else:
            ptitle += ' zonal concentration'
        ptitle += ' ({})'.format( units )

        fig  = plt.figure(figsize=(22, 14), dpi=dpi, facecolor='w', edgecolor='w')

        cbVars = Vars[n,:,:,:,:].copy()*scale

        # set ranges for whiole data period
        if pcent:
            if len( cbVars[ cbVars>500 ] ) > 0:
                cbVars =  np.ma.masked_where( cbVars>500, cbVars )
                extend='max'
            elif len( cbVars[cbVars<-500])>0:
                cbVars = np.ma.masked_where( cbVars<-500, cbVars )
                if  extend == 'max':
                    extend = 'both'
                else:
                    extend = 'min'
            else:
                extend='neither'

        else:
            extend='neither'

        if set_lon:
            set_lon = get_gc_lon( set_lon, res=res )
            fixcb  = [( i.min(), i.max() ) for i in [ cbVars[set_lon,...] ]][0]
            print('SETTING LON to GC index: ', set_lon)
        else:
            cbVars = cbVars.mean( axis=0 )
        if set_window:
            gclat_0, gclat_1 = [ get_gc_lat(i, res=res ) for i in (lat_0, lat_1) ]
            cbVars = cbVars[...,gclat_0:gclat_1, :]

        fixcb  = [( i.min(), i.max() ) for i in [ cbVars ]][0]

        # Kludge, force max cap at .2
        if units == 'ratio':
            fixcb = [ fixcb[0], 0.2 ]
            extend = 'max'

        cmap = get_colormap( fixcb  )

        # Loop thorugh months
        for m, month in enumerate(dlist):
            axn =[ 3,4,m+1 ]
            ax = fig.add_subplot( *axn )

#                if pcent:
#                    print [ (np.min(i), np.max(i))  \
#                        for i in [ Vars[n,:,:,:,m].mean(axis=0)*scale ] ]
#                    print [ (np.min(i), np.max(i))  \
#                        for i in [ Vars[n,:,:,:,m].median(axis=0)*scale ] ]
            if set_lon:
                arr = Vars[n,set_lon,:,:,m]*scale
            else:
                arr = Vars[n,:,:,:,m].mean(axis=0)*scale

            if set_window:
                arr = arr[...,get_gc_lat(lat_0, res=res): get_gc_lat(lat_1, res=res), :]

            # plot up spatial surface change
            zonal_plot(fig, ax, arr, title=month.strftime("%b"), debug=debug,
                tropics=False, units=units,f_size=f_size, c_off =37, no_cb=True, \
                lat_0=lat_0, lat_1=lat_1, set_window=set_window, fixcb=fixcb, \
                extend=extend, window=True, lower_limited=True, res=res, \
                mask_invalids=mask_invalids, cmap=cmap )

            # only show troposphere
            greyoutstrat( fig, t_ps.mean(axis=0)[:,:,m], axn=axn, res=res )

        # Add single colorbar
        mk_cb(fig, units=units, left=0.915,  cmap=cmap,vmin=fixcb[0],
            vmax=fixcb[1], f_size=f_size, extend=extend )

        # sort out ascetics -  adjust plots and add title
        fig.subplots_adjust( bottom=bottom, top=top, left=left, right=right, \
            hspace=hspace, wspace=wspace)
        fig.suptitle( ptitle, fontsize=f_size*2, x=.55 , y=.95  )

        # save out figure
        plot2pdfmulti( pdff, savetitle, dpi=dpi, no_dstr=no_dstr )

        # close fig
        plt.clf()
        plt.close()
        del fig

    #  save entire pdf
    plot2pdfmulti( pdff, savetitle, close=True, dpi=dpi, no_dstr=no_dstr )


# --------
# X.XX - Change as 2D plot of surface ( for column or surface change )
# --------
def plot_specs_poles_change_monthly2pdf( specs=None, arr=None, res='4x5', \
        dpi=160, no_dstr=True, f_size=20, pcent=False, diff=False, \
        dlist=None, savetitle='', units=None, perspective='north', \
        format=None, extend='neither', boundinglat=50,
        verbose=True, debug=False):
    """
    Takes a 5D np.array ( species, lon, lat, alt, time ) and plots up the
    output by species by month, and saves this as a mulitpage pdf

    Parameters
    -------
    arr (array): 5D np.array ( species, lon, lat, alt, time )
    res (str): the resolution if wd not given (e.g. '4x5' )
    dpi (int): dots per inch of saved output PDF...
    boundinglat (int):
    format (str): formayt of axis labels
    dlist (list): list of dates (datetimes)
    no_dstr (boolean): date string in output filename ("no date string")
    f_size (float): fontsise
    savetitle (str): string to add to filename of PDF
    units (str): units label for colorbar
    pcent (boolean): setup the plot as if the input values were %
    diff (boolean): setup the plot as if the input values were a difference
    boundinglat (int): latitude to show poles until.
    perspective (str): looking at north or south pole?
    extend (str): colorbar format settings ( 'both', 'min', 'both' ... )

    Returns
    -------
    (None)

    Notes
    -----
     - Takes 5D array  ( species ,lon , lat, alt, time)
     - needs update to description/acsetics
    """
    if debug:
        print(arr, no_dstr, f_size, pcent, res, dpi, specs, dlist, savetitle)

    # Setup PDF filename for saving file
    if perspective == 'north':
        savetitle = 'North_'+savetitle
    if perspective == 'south':
        savetitle = 'South_'+savetitle
    # Initialise PDF
    pdff = plot2pdfmulti( title=savetitle, open=True, dpi=dpi, no_dstr=no_dstr)
    # Set ascetics
    left=0.01; right=0.925; bottom=0.025; top=0.95; hspace=0.15; wspace=-0.1

    # Loop species
    for n, spec in enumerate( specs ):
        # Debug print statement?
        if verbose:
            print(n, spec, arr.shape, units, perspective, '<')

        # Get units/scale for species + setup fig
        scale = 1
        if pcent :# and (not units == 'DU') :
            units = '%'
        if isinstance( units, type( None )):
            units, scale = tra_unit( spec, scale=True, global_unit=True )
        parr = arr[n,:,:,0,:]*scale

        if debug:
            print(parr.shape)
            print(n, spec,  units, scale)
            print([ (i.min(), i.max(), i.mean() ) for i in [parr, arr] ])

        # Set the correct title
        ptitle = '{}'.format( latex_spec_name(spec) )

        # Create new figure
        fig = plt.figure(figsize=(22, 14),dpi=dpi,facecolor='w', edgecolor='w')

        # Select north or south polar areas specified to define cb
        if perspective == 'north':
            cbarr=parr[:,get_gc_lat(boundinglat,res=res):,:].copy()
        if perspective == 'south':
            cbarr=parr[:,:get_gc_lat(-boundinglat,res=res),:].copy()

        # Mask above and below 500/-500 % if values in array
        if pcent:
            if len( cbarr[cbarr>500])>0:
                cbarr = np.ma.masked_where( cbarr>500, cbarr )
                extend = 'max'
            elif len( cbarr[cbarr<-500])>0:
                cbarr = np.ma.masked_where( cbarr<-500, cbarr )
                if  extend == 'max':
                    extend = 'both'
                else:
                    extend = 'min'
            else:
                extend = 'neither'
        else:
            extend = 'neither'
        # Setup colormap
        fixcb = np.array([ ( i.min(), i.max() ) for i in [cbarr] ][0])
        if verbose:
            print('fixcb testing ', fixcb, parr.shape)
        # Kludge, force max cap at .2
#        if units == 'ratio':
#            fixcb = [ fixcb[0], 0.2 ]
#                extend = 'max'
        cmap = get_colormap( fixcb )

        # Loop thorugh months
        for m, month in enumerate(dlist):
            ax = fig.add_subplot(3,4,m+1)

            # Plot up spatial surface change
            if perspective == 'north':
                north_pole_surface_plot( parr[:,:,m], no_cb=True, \
                    fixcb=fixcb, diff=diff, pcent=pcent, res=res, \
                    f_size=f_size*2, cmap=cmap, boundinglat=boundinglat)
            if perspective == 'south':
                south_pole_surface_plot( parr[:,:,m], no_cb=True, \
                    fixcb=fixcb, diff=diff, pcent=pcent, res=res, \
                    f_size=f_size*2, cmap=cmap, boundinglat=-boundinglat)

            plt.text( 1.1, -0.01, month.strftime("%b"), fontsize=f_size*2,\
                transform=ax.transAxes, ha='right', va='bottom' )

        # Add single colorbar
        mk_cb(fig, units=units, left=0.915,  cmap=cmap,vmin=fixcb[0],\
            vmax=fixcb[1], f_size=f_size*.75, extend=extend, format=format )

        # Sort out ascetics -  adjust plots and add title
        fig.subplots_adjust( bottom=bottom, top=top, left=left, right=right, \
            hspace=hspace, wspace=wspace)
        fig.suptitle( ptitle, fontsize=f_size*2, x=.475 , y=.975  )

        # Save out figure
        plot2pdfmulti( pdff, savetitle, dpi=dpi,no_dstr=no_dstr )

        # Close fig
        plt.clf()
        plt.close()
        del parr

    # Save entire pdf
    plot2pdfmulti( pdff, savetitle, close=True, dpi=dpi, no_dstr=no_dstr )


# --------
# X.XX - Generic X vs. Y plot
# --------
def X_Y_scatter( x, y, z=None, fig=None, ax=None, vmin=None, vmax=None, \
        left= 0.1, width=0.60, bottom=0.1, height=0.60, widthII=0.2,  \
        lim2std=10,trendline=True, f_size=20, lw=10, title=None, \
        line121=True, X_title=None, Y_title=None ):
    """
    Plot up a X Y scatter plot of x vs. y
    """
#    rect_scatter = [left, bottom, width, height]

    if isinstance( fig, type(None) ):
        fig = plt.figure(1, figsize=(8,8))

    if isinstance( ax, type(None) ):
        ax = plt.axes( )#rect_scatter)

    # Plot up - normalising colors against z  if given
    if isinstance( z, type(None) ):
        plt.scatter( x, y )

    else:
        # Scale colors.
        if isinstance( vmin, type(None) ):
            vmin = float( np.ma.min(z)  )
        if isinstance( vmax, type(None) ):
            vmin = float( np.ma.max(z) )
        s_z = [ cmap( (float(i) - vmin) / ( np.array([vmin,vmax]) ).ptp() ) \
            for i in z ]

        pts = axScatter.scatter(x, y, c=z)


    if lim2std != False:
        stds = [ np.std( i ) for i in (x, y) ]
        means = [ np.mean( i ) for i in (x, y) ]
        print(stds, means)
        mins = [ means[0]-(stds[0]*lim2std), means[1]-(stds[1]*lim2std) ]
        maxs = [ means[0]+(stds[0]*lim2std), means[1]+(stds[1]*lim2std) ]

        # do not let miniums be less than zero.
        ind  = [ n for n, i in enumerate( mins ) if (i < 0 ) ]
        if len(ind)>0:
            for n in ind:
                mins[n] = 0

        plt.xlim( mins[0], maxs[0] )
        plt.ylim( mins[1], maxs[1] )
    else:
        min_, max_  = [  [i.min(), i.max()] for i in  [ np.array([ x,y ]) ] ][0]
        plt.xlim( min_, max_ )
        plt.ylim( min_, max_ )

    if line121:
        xmin, xmax = np.min(x), np.max(x)
        line121 = np.arange( xmin/2, xmax*2 )
        plt.plot( line121, line121, color='red', ls='--', lw=lw )

    # add trendline
    if trendline:
        Trendline( ax, x, y, order =1, intervals= 700, f_size=f_size, lw=lw,
            color='green' )

    if not isinstance( X_title, type( None) ):
        plt.xlabel( X_title, fontsize=f_size )
    if not isinstance( Y_title, type( None) ):
        plt.ylabel( Y_title, fontsize=f_size )
    if not isinstance( title, type( None) ):
        plt.title( title, fontsize=f_size )
    plt.xticks( fontsize=f_size )
    plt.yticks( fontsize=f_size )


# --------
# X.XX - Scatter 3D cube
# --------
def scatter_3D_cube( data, dims=None, res='2x2.5', fig=None, everyother=1, interval=1, \
        f_size=20, cm='RdYlBu', debug=False ):
    """
    Make a 3D scatter cube of given 3D array
    """
    logging.debug('scatter_3D_cube called with input data shape: {}'.format(\
        data.shape) )

    from mpl_toolkits.mplot3d import Axes3D

    if isinstance( dims, type(None) ):
        lon, lat, alt = get_latlonalt4res( res=res )
        (X,Y,Z) = np.mgrid[ \
            lon[0]:lon[-1]:complex( len(lon) ),
            lat[0]:lat[-1]:complex( len(lat) ),
            alt[0]:alt[-1]:complex( len(alt) ) ]

    # Setup Fig + Ax ( is not given )
    cmap = plt.cm.get_cmap(cm)
    if isinstance( fig, type(None) ):
        fig = plt.figure(1)
    fig.clf()
    ax = Axes3D(fig)
    ax.scatter(X,Y,Z, c=data, cmap=cmap)

    # --- Beatify
    # draw meridian lines
    meridians = np.arange(-180,180,20*interval)
    plt.xticks( meridians[::everyother], fontsize = f_size*.75 )

    # draw parrelel lines
    parallels = np.arange(-90,91,20*interval)
    plt.yticks( parallels[::everyother], fontsize = f_size*.75 )

    # limit axis
    ax.set_xlim(  lon[0], lon[-1] )
    ax.set_ylim(  lat[0], lat[-1] )
    ax.set_zlim(  alt[0], alt[-1] )


# --------
# X.XX - Plot up seasonal output from 4D arr (lon, lat, alt, time)
# --------
def get_seasonal_plot( arr, fixcb=None, fig=None, f_size=15, \
        case='linear', format=None, extend='neither', units=None, \
        right=0.9, left=0.05, bottom=0.05, top=0.85, hspace=0.1, \
        wspace=0.1, log=False, title=None, dpi=80, debug=False ):
    """
    Takes any 4D array and plot a 4 subplot window plot by season
    """

    # Split by quater (DJF, MAM, JJA, SON)
    ars, seasons = split_4D_array_into_seasons( arr, annual_plus_seasons=False )

    # create figure
    if isinstance( fig, type(None ) ):
        fig  = plt.figure(figsize=(22, 14), dpi=dpi, facecolor='w', \
                edgecolor='w')

    # fix color mapping
    if isinstance( fixcb, type(None ) ):
        fixcb = [ ( i.min(), i.max() ) for i in [ arr ] ][0]
    cmap = get_colormap( fixcb  )

    # loop seasons
    for n, arr in enumerate( ars ):

        # Plot up  on new axis
        ax = fig.add_subplot( 2, 2, n+1)
        map_plot( arr.T, title=seasons[n], units=None, window=True, \
                        case=case,  f_size=f_size, rotatecbunits='vertical',\
                        no_cb=True, cmap=cmap, fixcb=fixcb, fig=fig )

        # make color bar
        mk_cb(fig, units=units, left=0.925,  cmap=cmap, vmin=fixcb[0],
            vmax=fixcb[1], log=log, f_size=f_size*.5, extend=extend,
            format=format, debug=debug )

        # add title if provided
        if not isinstance( title, type(None) ):
            fig.suptitle( title, fontsize=f_size*2, x=.55 , y=.95  )


    # Adjust figure
    fig.subplots_adjust( bottom=bottom, top=top, left=left,\
        right=right,hspace=hspace, wspace=wspace)


# --------
# X.XX - PDF of annual surface change plots for given species (takes 5D arr )
# --------
def plot_specs_surface_change_annual2pdf( arr, res='4x5', dpi=160, \
        no_dstr=True, f_size=20, pcent=True, specs=None, dlist=None, \
        savetitle='', diff=False, extend='neither', column=False, \
        scale=1, units=None, set_window=False, lat_0=None, lat_1=None, \
        mask_invalids=False, debug=False):
    """
    Create multipage PDF with each page containing a 2D (lon,lat) slice
    plot for given species in list of "specs" Takes 5D array  ( species ,lon , lat, alt,
    time)
    """
    logging.info( 'plot_specs_surface_change_monthly2pdf called' )

    # setup pdfs + titles
    if column:
        savetitle = 'Column_by_spec'+savetitle
    else:
        savetitle = 'Surface_by_spec'+savetitle
    pdff = plot2pdfmulti( title=savetitle, open=True, dpi=dpi, no_dstr=no_dstr )
    left=0.05; right=0.9; bottom=0.05; top=0.875; hspace=0.315; wspace=0.1

    # Loop species
    for n, spec in enumerate( specs ):
#            if debug:
        print(n, spec)

        # Get units/scale for species + setup fig (allow of
        if isinstance( units, type( None) ):
            if column and (not pcent):
                units, scale = 'DU', 1
            elif pcent:
                units, scale = '%', 1
            else:
                units, scale = tra_unit( spec, scale=True, global_unit=True )

        # setup masking...
        cbarr = arr[n,:,:,0,:].mean(axis=-1).copy() * scale

        if pcent:
            mask_invalids = True

            # mask for changes greater than 500%
            if len( cbarr[ cbarr>500 ]  ) > 0:
                cbarr = np.ma.masked_where( cbarr>500, cbarr )
                extend='max'

            elif len( cbarr[cbarr<-500])>0:
                cbarr = np.ma.masked_where( cbarr<-500, cbarr )
                if  extend == 'max':
                    extend = 'both'
                else:
                    extend = 'min'

            else:
                extend='neither'

        else:
            extend='neither'
            cbarr = cbarr

        # Set the correct title
        ptitle = 'Annual Avg. {}'.format( latex_spec_name(spec) )
        if column:
            ptitle += ' column'
        else:
            ptitle += ' surface'

        if diff:
            ptitle += ' $\Delta$ concentration'
        else:
            ptitle += ' concentration'
        ptitle += ' ({})'.format( units )

        fig  = plt.figure(figsize=(22, 14), dpi=dpi, facecolor='w', edgecolor='w')

        # set cb ranges for whiole data period
        fixcb  = [( i.min(), i.max() ) for i in [cbarr]][0]

        # Kludge, force max cap at .2
#            if units == 'ratio':
#                fixcb = [ fixcb[0], 0.2 ]
#                extend = 'max'
        if (units == 'pmol mol$^{-1}$ m$^{-3}$') and ( spec == 'AERI' ):
            fixcb = [ fixcb[0], 500 ]
            extend = 'max'

        cmap = get_colormap( fixcb )

        # Loop thorugh months
#            for m, month in enumerate(dlist):
        fig.add_subplot(111)

        # plot up spatial surface change
        map_plot( arr[n,:,:,0,:].mean(axis=-1).T*scale, cmap=cmap, case=9, \
            res=res, no_cb=True, f_size=f_size, fixcb=fixcb, window=True,\
            set_window=set_window,lat_0=lat_0, lat_1=lat_1, \
            mask_invalids=mask_invalids, debug=debug)

        # Add single colorbar
        mk_cb(fig, units=units, left=0.905,  cmap=cmap,vmin=fixcb[0],
            vmax=fixcb[1], f_size=f_size, extend=extend )

        # sort out ascetics -  adjust plots and add title
        fig.suptitle( ptitle, fontsize=f_size*2, x=.55 , y=.95  )

        # save out figure
        plot2pdfmulti( pdff, savetitle, dpi=dpi, no_dstr=no_dstr )

        # close fig
        plt.clf()
        plt.close()
        del fig

    #  save entire pdf
    plot2pdfmulti( pdff, savetitle, close=True, dpi=dpi, no_dstr=no_dstr )


# --------
# X.XX - PDF of annual zonal change plots for given species
# --------
def plot_specs_zonal_change_annual2pdf( Vars, res='4x5', dpi=160, \
        no_dstr=True, f_size=20, pcent=False, specs=None, dlist=None, \
        t_ps=None, savetitle='', diff=False, extend='neither',  \
        set_window=False, lat_0=None, lat_1=None, mask_invalids=False, \
        set_lon=None, units=None, debug=False):
    """
    Create multipage PDF with each page containing a zonal plot for given species in list
     of "specs"

    NOTES:
     - Takes 5D array  ( species ,lon , lat, alt, time)
    """
    savetitle = 'Annual_Zonal_by_spec'+savetitle
    pdff = plot2pdfmulti( title=savetitle, open=True, dpi=dpi, no_dstr=no_dstr )
    left=0.05; right=0.9; bottom=0.05; top=0.875; hspace=0.315; wspace=0.2

    # Loop species
    for n, spec in enumerate( specs ):
#           if debug:
        print(n, spec, Vars.shape)

        # Get units/scale for species + setup fig
#           scale=1
        if isinstance( units, type(None) ):
            units, scale = tra_unit( spec, scale=True, global_unit=True )
        if pcent:
            units, scale = '%', 1
            mask_invalids=True

        # Set the correct title
        ptitle = 'Annual avg. {}'.format( latex_spec_name(spec) )
        if diff:
            ptitle += ' zonal $\Delta$ concentration'
        else:
            ptitle += ' zonal concentration'
        ptitle += ' ({})'.format( units )

        fig  = plt.figure(figsize=(22, 14), dpi=dpi, facecolor='w', edgecolor='w')

        cbVars = Vars[n,:,:,:,:].mean(axis=-1).copy()*scale

        # set ranges for whiole data period
        if pcent:
            if len( cbVars[ cbVars>500 ] ) > 0:
                cbVars =  np.ma.masked_where( cbVars>500, cbVars )
                extend='max'
            elif len( cbVars[cbVars<-500])>0:
                cbVars = np.ma.masked_where( cbVars<-500, cbVars )
                if  extend == 'max':
                    extend = 'both'
                else:
                    extend = 'min'
            else:
                extend='neither'

        else:
            extend='neither'

        if set_lon:
            set_lon = get_gc_lon( set_lon, res=res )
            fixcb  = [( i.min(), i.max() ) for i in [ cbVars[set_lon,...] ]][0]
            print('SETTING LON to GC index: ', set_lon)
        else:
            cbVars = cbVars.mean( axis=0 )
        if set_window:
            gclat_0, gclat_1 = [ get_gc_lat(i, res=res ) for i in (lat_0, lat_1) ]
            cbVars = cbVars[...,gclat_0:gclat_1, :]

        fixcb  = [( i.min(), i.max() ) for i in [ cbVars ]][0]

        # Kludge, force max cap at .2
        if units == 'ratio':
            fixcb = [ fixcb[0], 0.2 ]
            extend = 'max'

        cmap = get_colormap( fixcb  )

        axn =[ 111 ]
        ax = fig.add_subplot( *axn )

#                if pcent:
#                    print [ (np.min(i), np.max(i))  \
#                        for i in [ Vars[n,:,:,:,m].mean(axis=0)*scale ] ]
#                    print [ (np.min(i), np.max(i))  \
#                        for i in [ Vars[n,:,:,:,m].median(axis=0)*scale ] ]
        if set_lon:
            arr = Vars[n,set_lon,...].mean(axis=-1)*scale
        else:
            arr = Vars[n,...].mean(axis=0).mean(axis=-1)*scale

        if set_window:
            arr = arr[...,get_gc_lat(lat_0, res=res): get_gc_lat(lat_1, res=res), :]

        # plot up spatial surface change
        zonal_plot( arr, fig, ax=ax, title=None, debug=debug, tropics=False, \
            units=units,f_size=f_size, c_off =37, no_cb=True, lat_0=lat_0, lat_1=lat_1, \
            set_window=set_window, fixcb=fixcb, extend=extend, window=True, \
            lower_limited=True, res=res, mask_invalids=mask_invalids, cmap=cmap )

        # only show troposphere
        greyoutstrat( fig, t_ps.mean(axis=0).mean(axis=-1), axn=axn, res=res )

        # Add single colorbar
        mk_cb(fig, units=units, left=0.915,  cmap=cmap,vmin=fixcb[0],
            vmax=fixcb[1], f_size=f_size, extend=extend )

        # sort out ascetics -  adjust plots and add title
#            fig.subplots_adjust( bottom=bottom, top=top, left=left, \
#                right=right,hspace=hspace, wspace=wspace)
        fig.suptitle( ptitle, fontsize=f_size*2, x=.55 , y=.95  )

        # save out figure
        plot2pdfmulti( pdff, savetitle, dpi=dpi, no_dstr=no_dstr )

        # close fig
        plt.clf()
        plt.close()
        del fig

        #  save entire pdf
        plot2pdfmulti( pdff, savetitle, close=True, dpi=dpi, no_dstr=no_dstr )


# --------
# X.XX - Spatial Figure maker ( just provide lon, lat, time,  np array )
# --------
def plot_spatial_figure( arr, fixcb=None, sigfig_rounding_on_cb=2, \
        norm=None, nticks=10, format=None, units=None, extend='neither',
        ax=None, cb_height=0.6, \
        discrete_cmap=False, f_size=15, fig=None, left_cb_pos=0.86, cb_ax=None,\
        bottom=0.005, top=0.95, hspace=0.4, wspace=0.3, left=0.035, right=0.85,\
        dpi=160, res='4x5', show=True, pdf=False, pdftitle=None, title=None, \
        window=False, interval=1, ylabel=True, cb='CMRmap_r', width=0.015,\
        orientation='vertical', rotatecbunits='vertical', title_y=1, \
        title_x=0.5, no_cb=True, return_m=False, log=False, wd=None, \
        resolution='c', lat_min = None, lat_max=None, lon_min=None, \
        lon_max=None, xlabel=True, limit_window=False, axis_titles=False,  \
        bottom_cb_pos=0.2,  figsize=(15, 10), verbose=False, debug=False ):
    """
    Wrapper for map_plot - Creates a "standard" spatial plot with acceptable
    ascethics. Customise with a range of arguements provide during the call
    to function.

    Parameters
    ----------
    adjust_window (int): amount of array entries to remove the edges of array
    alhpa (float): transparency of plotter data
    arr (np.array): input (2D) array
    case (str or int): case for type of plot (vestigle: use log=True of False (default))
    cmap (str): force colormap selection by providing name
    centre (boolean): use centre points of lon/lat grid points for mapping data surface
    drawcountries (boolean): add countries to basemap?
    debug (boolean): legacy debug option, replaced by python logging
    degrade_resolution (boolean): reduce resolution of underlay map detail
    discrete_cmap (boolean): use a discrete instead of conitunous colorbar map
    everyother (int): use "everyother" axis tick (e.g. 3=use every 3rd)
    f_size (float): fontsise
    fixcb (np.array): minimium and maximum to fix colourbar
    fixcb_buffered (array): minimium and maximum to fix colourbar, with buffer space
    format (str): format string for colorbar formating
    grid (boolean): apply a grid over surface plot?
    extend (str): colorbar format settings ( 'both', 'min', 'both' ... )
    interval (int): x/y tick interval in multiples of 15 degrees lat/lon
    lvls (list): manually provide levels for colorbar
    log (boolean): use a log scale for the plot and colorbar
    no_cb (boolean): include a coloubar?
    norm (norm object): normalisation to use for colourbar and array
    nticks (int): number of ticks on colorbar
    mask_invalids (boolean): mask invalid numbers (to allow saving to PDF)
    res (str): GEOS-Chem output configuration resolution ( '4x5' etc... )
    resolution (str): basemasp resolution settings ( 'c' = coarse, 'f' = fine ...  )
    rotatecbunits (str): orientation of colourbar units
    shrink (boolean): colorbar size settings ( fractional shrink )
    set_window (boolean): set the limits of the plotted data (lon_0, lon_1, lat_0, lat_1)
    (for nested boundary conditions )
    sigfig_rounding_on_cb (int): significant figure rounding to use for colourbar
    set_cb_ticks (boolean): mannually set colorbar ticks? (vestigle)
    title (str): plot title (deafult is ==None, therefore no title)
    tight_layout (boolean): use use tight lyaout for figure
    ylabel, xlabel (boolean): label x/y axis?
    units (str): units of given data for plot title
    verbose (boolean): legacy debug option, replaced by python logging
    wd (str): Specify the wd to get the results from a run.
    window (boolean): use window plot settings (fewer axis labels/resolution of map)

    Returns
    -------
    optionally returns basemap (return_m==True) and colorbar (no_cb!=True) object

    NOTES:
        Provide an 3D array of lon, lat, and alt
    """
    # If just lat and lon provided, add a dummy dimension.
    if len(arr.shape) == 2:
        arr = arr[...,None]

    logging.info( 'plot_spatial_figure called, with shape {}, fixcb: {}'.format(\
        arr.shape,  fixcb)+', min: {}, max:{}'.format( arr.min(), arr.max()) )
    logging.debug('@ surface, min: {} and max: {}'.format( arr[...,0].min(), \
        arr[...,0].max()))

    # setup fig if not provided
    if isinstance( fig, type(None) ):
        fig = plt.figure(figsize=figsize, dpi=dpi, facecolor='w',edgecolor='w')
    # setup fig if not provided
    if not isinstance( ax, type(None) ):
        # temporary remove as mpl widget has a bug
        # http://stackoverflow.com/questions/20562582/axes-instance-argument-was-not-found-in-a-figure
#        plt.sca( ax )
        pass

    # Set colourbar limits
    if isinstance( fixcb, type(None) ):
        fixcb = np.array( [ (i.min(), i.max()) for i in [arr[...,0] ] ][0] )

        # Set readable levels for cb, then use these to dictate cmap
        lvls = get_human_readable_gradations( vmax=fixcb[1],  \
            vmin=fixcb[0], nticks=nticks,
            sigfig_rounding_on_cb=sigfig_rounding_on_cb  )

    else:
        if log:
            print('WARNING: Code (to create levels for log colorbar)'+\
                'below needs checking')
            # Get logarithmically spaced integers
            lvls = np.logspace( np.log10(fixcb[0]), np.log10(fixcb[1]), \
                                                 num=nticks)
            # Normalise to Log space
#            norm=mpl.colors.LogNorm(vmin=fixcb_[0], vmax=fixcb_[1])
            if isinstance( norm, type(None) ):
                norm=mpl.colors.LogNorm( vmin=fixcb[0], vmax=fixcb[1] )

        else:
            # Assume numbers provided as fixcb +nticks will allow for creation
            # of human readable levels.
            lvls = np.linspace( fixcb[0], fixcb[1], nticks )

    # Setup Colormap
    cmap, fixcb_buffered = get_colormap( np.array( fixcb ), \
        nticks=nticks, fixcb=fixcb, cb=cb, buffer_cmap_upper=True )

    if discrete_cmap:
        cmap, norm = mk_discrete_cmap( nticks=nticks,\
                vmin=fixcb[0], vmax=fixcb[1], cmap=cmap )
    logging.debug( 'array min, max, mean='.format(  \
        *[ str((i.min(), i.max(), i.mean())) for i in [arr[...,0]]]) )

    # Plot up
    plt_vars = map_plot( arr[...,0].T, format=format, cmap=cmap, ax=ax, \
        fixcb=fixcb, return_m=return_m, log=log, window=window, no_cb=True,
        norm=norm, f_size=f_size*.75,  res=res, wd=wd, resolution=resolution,\
        fixcb_buffered=fixcb_buffered, interval=interval, xlabel=xlabel, \
        ylabel=ylabel, axis_titles=axis_titles, verbose=verbose, debug=debug )

    # if title != None, add to plot
    if not isinstance(title, type(None)):
        try:
            ax.annotate( title , xy=(title_x, title_y), \
                textcoords='axes fraction', fontsize=f_size)
        except:
            logging.info('WARNING! using plt, not axis, for title annotation')
            plt.title( title, fontsize=f_size, y=title_y )
#            plt.text(0.5, title_y, title, fontsize=f_size )

    # limit displayed extent of plot?
#	limit_window=False
    x_y_limits = [ lon_min, lon_max, lat_min, lat_max ]
    x_y_limited = any([ (not isinstance(i, type(None))) for i in x_y_limits ])
    if limit_window or x_y_limited:
        ax = plt.gca()
        #  set axis limits
        ax.set_xlim( lon_min, lon_max )
        ax.set_ylim( lat_min, lat_max )

    # Manually Add colorbar
    logging.debug('colorbar orientation: {}'.format(orientation))
    if no_cb:
        if orientation == 'vertical':
            width = width/2

        cb_ax = mk_cb(fig, units=units, left=left_cb_pos,  cmap=cmap, \
                vmin=fixcb_buffered[0], cb_ax=cb_ax, width=width, \
                height=cb_height,\
                rotatecbunits=rotatecbunits, bottom=bottom_cb_pos, \
                vmax=fixcb_buffered[1], format=format, f_size=f_size*.75, \
                extend=extend, lvls=lvls, log=log, orientation=orientation, \
                sigfig_rounding_on_cb=sigfig_rounding_on_cb, nticks=nticks, \
                norm=norm, discrete_cmap=discrete_cmap, debug=debug )

    # Adjust plot ascetics
    fig.subplots_adjust( bottom=bottom, top=top, left=left, right=right,
                                        hspace=hspace, wspace=wspace)


    # Show/Save as PDF?
    if pdf:
        plot2pdf( title=pdftitle )
    if show:
        plt.show()
    if return_m:
        return [ fig, cmap ] + plt_vars + [ fixcb ] #+= [ cb_ax ]


# --------
# X.XX - Zonal Figure maker ( just provide lon, lat np array )
# --------
def plot_zonal_figure( arr, fixcb=None, sigfig_rounding_on_cb=2, ax=None, \
        norm=None, nticks=10, format=None, units=None, extend='neither', \
        discrete_cmap=False, f_size=15, fig=None, res='4x5', wd=None, t_ps=None, \
        trop_limit=True, axn=None, cb_ax=None, orientation='vertical', \
        rotatecbunits='vertical', width=0.015, height=0.6, \
        bottom=0.1, top=0.925, hspace=0.4, wspace=0.5, left=0.075, right=0.875, \
        cb_bottom=0.125, cb_height=0.825, cb_left=0.885, dpi=160, no_cb=True, \
        region='All', lat_0=None, lat_1=None, pdftitle=None, return_m=False, \
        rtn_plt_vars=False, set_window=False, pdf=False, show=True, log=False, \
        window=False, xlabel=True, ylabel=True, title=None, \
        interval=None, verbose=False, debug=False ):
    """
    Creates zonal zonal plot as formated figure.


    NOTES:
     -
    """
    all_str = 'plot_zonal_figure called ', region, arr.shape, log, units, pdf, \
            show, arr.min(), arr.max()
    logging.info( all_str )
    if verbose:
        print(all_str)

    # If lon, lat, alt array provided then take mean of lon
    if any( [arr.shape[0] ==i for i in (72, 144, 121, 177)] ):
#        arr = arr.mean(axis=0)
        arr = molec_weighted_avg( arr, weight_lon=True, res=res, \
            trop_limit=trop_limit, rm_strat=False, wd=wd)

    # Create figure if not provided
    if isinstance( fig, type(None) ):
        fig = plt.figure(figsize=(15, 10), dpi=dpi,
                    facecolor='w', edgecolor='w')

    # if just plotting over the ocean, remove white space
    if region == 'Oceanic':
        set_window=True
        lat_0, lat_1 = -65, 80
        if isinstance( interval, type(None) ):
            interval = 3

    # Set colourbar limits
    if isinstance( fixcb, type(None) ):
        fixcb = np.array([(i.min(), i.max()) for i in [arr] ][0])

        # Set readable levels for cb, then use these to dictate cmap
        if not log:
            lvls = get_human_readable_gradations( vmax=fixcb[1],  \
                vmin=fixcb[0], nticks=nticks,\
                sigfig_rounding_on_cb=sigfig_rounding_on_cb, \
                verbose=verbose, debug=debug  )
    else:
        # Assume numbers provided as fixcb +nticks will allow for creation
        # of human readable levels.
        lvls = np.linspace( fixcb[0], fixcb[1], nticks )

    # If log plot - overwrite  lvls
    if log:
            # Get logarithmically spaced integers
            lvls = np.logspace( np.log10(fixcb[0]), np.log10(fixcb[1]), \
                                                 num=nticks)
            # Normalise to Log space
#            norm=mpl.colors.LogNorm(vmin=fixcb_[0], vmax=fixcb_[1])
            if isinstance( norm, type(None) ):
                norm=mpl.colors.LogNorm(vmin=fixcb[0], vmax=fixcb[1])

    # Setup Colormap
    cmap, fixcb_buffered = get_colormap( np.array( fixcb ), \
            nticks=nticks, fixcb=fixcb, buffer_cmap_upper=True, \
            verbose=verbose, debug=debug )

    #  Plot
    if isinstance( axn, type( None ) ):
        axn =[ 111 ]
    if isinstance( ax, type( None ) ):
        ax = fig.add_subplot( *axn )
    zonal_plot( arr, fig,  ax=ax, set_window=set_window, log=log,\
            format=format, cmap=cmap, lat_0=lat_0, lat_1=lat_1, \
            fixcb=fixcb, f_size=f_size*.75, res=res, norm=norm, \
            fixcb_buffered=fixcb_buffered, no_cb=True, trop_limit=True,\
            window=window, interval=interval, xlabel=xlabel, ylabel=ylabel,\
            verbose=verbose, debug=debug )

    # Only show troposphere
    if isinstance( t_ps, type(None) ):
        t_ps = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], trop_limit=True )
    greyoutstrat( fig, t_ps.mean(axis=0).mean(axis=-1), axn=axn, res=res )

    if not isinstance( title, type( None ) ):
        plt.title( title, fontsize=f_size*.75 )
#        plt.text(0.5, y_title, title, fontsize=f_size*.75 )

    # Manually Add colorbar
    if no_cb:
        if orientation == 'vertical':
            width = width/2
            height = 0.55

        mk_cb(fig, units=units, left=cb_left,  height=cb_height, \
                bottom=cb_bottom, log=log, orientation=orientation, \
                rotatecbunits=rotatecbunits, width=width, \
                cmap=cmap, vmin=fixcb_buffered[0],\
                vmax=fixcb_buffered[1], format=format, f_size=f_size*.75, \
                extend=extend, lvls=lvls, cb_ax=cb_ax, \
                sigfig_rounding_on_cb=sigfig_rounding_on_cb, nticks=nticks, \
                norm=norm, discrete_cmap=discrete_cmap,
                verbose=verbose, debug=debug )

    # Adjust plot ascetics
    fig.subplots_adjust( bottom=bottom, top=top, left=left, right=right,
                                        hspace=hspace, wspace=wspace)
    # Show/Save as PDF?
    if pdf:
        plot2pdf( title=pdftitle )
    if show:
        plt.show()
    if return_m or rtn_plt_vars:
        return [ fig, cmap, fixcb ]# + plt_vars + [ fixcb ] #+= [ cb_ax ]


# --------
# X.XX - Lat plotter of average + Q1/Q3
# --------
def plot_arr_avg_Q1_Q3( X, Y, ax=None, color='blue', label=None, \
        plt_mean=True, plt_median=False, pcent1=25, pcent2=75, \
        verbose=False, debug=False ):
    """
    Takes a X and Y to plot a mean, Q1 and Q3.

    Notes:
     - Y = 2D array (e.g. lon, lat)
     - X = 1D araray ( e.g. lat )
     - Default fill is Q1 to Q3, but others ranges can be specified
    """
    if isinstance( ax, type(None) ):
        ax= plt.gca()

    # Plot up mean of Y values
    if plt_mean:
        ax = plt.plot( X, Y.mean(axis=0), color=color, label=label )

    # Show quartiles ( 25th 57th  )
    Y_nan = Y.filled( np.nan )
    low = np.nanpercentile( Y_nan, pcent1, axis=0 )
    high = np.nanpercentile( Y_nan, pcent2, axis=0 )
    ax.fill_between( X, low, high, alpha=0.2, color=color   )
    if plt_median:
        ax = plt.plot( X, np.nanpercentile( Y_nan, 50, axis=0 ), \
            color=color, label=label )

    # return axis object
    return ax


# --------------
# 1.35 - Timeseries plotter ( takes datetime + np.array )
# -------------
def timeseries_plot( ax, dates, data, f_size=20, pos=0, posn=1,  \
        title=None, legend=False, everyother=7*24,  x_nticks=12, \
        window=False, label=None, ylabel=None, loc='upper left',  \
        lw=1,ls='-', color=None, start_date=None, end_date=None, \
        boxplot=True, showmeans=False, alt_text=None, r_plt=False, \
        unitrotation=45, color_by_z=False, fig=None,  xlabel=True, \
        positive=None, plt_median=False, add_Q1_Q3=False, pcent1=25, pcent2=75, \
        debug=False ):
    """
    Plot up timeseries of values.

    NOTES:
     - Requires data, and dates in numpy array form.
     - Dates must be as datetime.datetime objects.
    """

    # Process data - reduce resolution to daily, and get std
    df = DataFrame( data, index=dates )

    # Take start and end dates from "dates" if not set in arguments.
    if isinstance( start_date, type(None) ):
        start_date = dates[0]
    if isinstance( end_date, type(None) ):
        end_date = dates[-1]
    df = df[ start_date:end_date]

    # label once per week ( set by "everyother" )
    days = [i.to_datetime() for i in df.index ]
    labels = [i.strftime("%-d %b") for  i in days ][::everyother]

    # Color in line another provided variables
    if color_by_z:
        if debug:
            print('Coloring line by normalised z values')
        print(df.columns)
        x = df.index
        y, z = [ df[ df.columns[i] ] for i in range(2) ]
        cmap = get_colormap( z.copy(), positive=positive )
        print([ ( i.min(), i.max() ) for i in (x, y, z) ])
        colorline(x, y, z, cmap=cmap, linewidth=lw, ax=ax, \
            norm=plt.Normalize( 0, 360 ), fig=fig ) #np.min(z), 1500))
#        colorline(x, y, linewidth=lw, ax=ax)

    else:
        if plt_median: # Plot average
            pass
#            ln = plt.plot( days, np.nanpercentile( df.values, 50, ),
#                color=color, ls=ls, lw=lw, label=None   )
        else: # Plot all
            plt.plot( days, df.values, label=label, color=color, ls=ls, lw=lw  )

    # Plot quartiles as shaded area?
    if plot_Q1_Q3:
        pass
#        low =np.nanpercentile( df.values, pcent1, )
#        high = np.nanpercentile( df.values, pcent2, )
#        ax.fill_between( bins_used, low, high, alpha=0.2, color=color   )

    # Setup X axis
    if xlabel:
        plt.xticks( days[::everyother], labels, rotation=unitrotation, \
            fontsize=f_size )
    else:
        plt.tick_params( axis='x', which='both', labelbottom='off')

    # Beautify plot
    if not isinstance( title, type(None) ):
        plt.title( title + ' for {}-{}'.format( start_date.strftime( \
            '%d/%m/%y' ), end_date.strftime( '%d/%m/%y' ) ),
             fontsize=f_size )
    # Alt text annotate as fig text?
    if not isinstance( alt_text, type(None) ):
        plt.figtext(x=0.05,y=0.85, s=alt_text, fontsize=f_size*.75 )
    # Setup y axis
    plt.yticks( fontsize=f_size )
    if not isinstance( ylabel, type(None) ):
        plt.ylabel( ylabel, fontsize=f_size )
    # Legend?
    if legend:
        plt.legend( fontsize=f_size*.75, loc=loc )


# --------
# 1.36 - Get monthly surface plots for (4D) array
# --------
def plt_4Darray_surface_by_month( arr, res='4x5', dpi=160, \
        no_dstr=True, f_size=10, dlist=None, fixcb=None, format=None, \
        savetitle='', extend='neither',  wd=None, ax=None, fig=None, \
        sigfig_rounding_on_cb=3, nticks=7, discrete_cmap=False, \
        units=None, set_window=False, lat_0=None, lat_1=None, \
        return_m=False, log=False, window=True, interval=3, ylabel=True,\
        norm=None, fig_title=False, pdftitle='',
        pdf=False, show=False, verbose=False, debug=False):
    """
    Create a window plot of surface amp plots from a 4D array
    """
    # Setup local variables + figure
#    left=0.015; right=0.9; bottom=0.05; top=0.95; hspace=0.225; wspace=-0.01
    left=0.015; right=0.87; bottom=0.05; top=0.95; hspace=0.225; wspace=-0.01
    fig  = plt.figure(figsize=(14, 10), dpi=dpi, facecolor='w', edgecolor='w')

    # Get datetime
    if isinstance( dlist, type(None) ):
        dlist = get_gc_datetime( wd=wd )

    # set cb ranges for whole data period
    if isinstance( fixcb, type(None) ):
        fixcb  = [( i.min(), i.max() ) for i in [arr]][0]

    # Create figure if not provided
    if isinstance( fig, type(None) ):
        fig = plt.figure(figsize=(15, 10), dpi=dpi,
                    facecolor='w', edgecolor='w')

    # Set readable levels for cb, then use these to dictate cmap
    lvls = get_human_readable_gradations( vmax=fixcb[1],  \
                    vmin=fixcb[0], nticks=nticks,
                    sigfig_rounding_on_cb=sigfig_rounding_on_cb  )

    # Setup Colormap
    cmap, fixcb_buffered = get_colormap( np.array( fixcb ), \
        nticks=nticks, fixcb=fixcb, buffer_cmap_upper=True )

    if discrete_cmap:
        cmap, norm = mk_discrete_cmap( nticks=nticks,\
                vmin=fixcb[0], vmax=fixcb[1], cmap=cmap )

    if debug:
        print([ (i.min(), i.max(), i.mean()) for i in [ arr.mean(axis=0) ] ])

    # Loop thorugh months
    for m, month in enumerate(dlist):

        # add axis
        axn = [4,3,m+1]
        ax = fig.add_subplot( *axn )

        print(arr[...,m].mean(axis=0).shape, arr.shape)

        # Only show x/y axis on edge plots
        ylabel=False
        xlabel=False
        num_cols=3
        num_rows=4
        if (m in range( len(dlist) )[-num_cols:] ):
            xlabel=True
        if any( [axn[-1]==i for i in range(1, len(dlist)+1)[::num_cols] ]) :
            ylabel=True

        # Plot up
        map_plot( arr[...,0,m].T, format=format, cmap=cmap, ax=ax, \
            fixcb=fixcb, return_m=return_m, log=log, window=window, \
            no_cb=True, norm=norm, f_size=f_size*.75,  res=res, \
            fixcb_buffered=fixcb_buffered, interval=interval,\
            ylabel=ylabel, xlabel=xlabel, verbose=verbose, debug=debug )

        # add month
        plt.title(month.strftime("%b"), fontsize=f_size*1.5)

    # Add single colorbar
#    mk_cb(fig, units=units, left=0.9, cmap=cmap, vmin=fixcb[0], format=format,\
    mk_cb(fig, units=units, left=0.87, cmap=cmap, vmin=fixcb[0], format=format,\
        vmax=fixcb[1], nticks=nticks, f_size=f_size, extend=extend )

    print(nticks, fixcb, lvls)


    # Sort out ascetics -  adjust plots and add title
    if fig_title:
        fig.suptitle(  '{}'.format( latex_spec_name(spec) ), fontsize=f_size*2,
            x=.55 , y=.95  )
        top = 0.9  # allow space for figure title
    fig.subplots_adjust( bottom=bottom, top=top, left=left, \
        right=right, hspace=hspace, wspace=wspace)

    #  save as pdf ?
    if pdf:
        plot2pdf( title=pdftitle )
    if show:
        plt.show()


# --------
# 1.37 - Get monthly surface plots for (4D) array
# --------
def plt_4Darray_zonal_by_month( arr, res='4x5', dpi=160, \
        no_dstr=True, f_size=15, dlist=None, fixcb=None, \
        savetitle='', extend='neither',  wd=None, ax=None, fig=None, \
        sigfig_rounding_on_cb=3, nticks=7, discrete_cmap=False, \
        units=None, set_window=False, lat_0=None, lat_1=None, \
        return_m=False, log=False, window=True, interval=3, ylabel=True,\
        norm=None, fig_title=False, pdftitle='', t_ps=None, xlabel=True, \
        format=None, orientation='vertical', trop_limit=True, region='All',
        pdf=False, show=False, verbose=False, debug=False):
    """
    Create a window plot of surface amp plots from a 4D array
    """

    # Average over lon
    arr = arr.mean(axis=0)

    # Setup local variables + figure
    left=0.075; right=0.875; bottom=0.085; top=0.955; hspace=0.325; wspace=0.1
    fig  = plt.figure(figsize=(7, 7), dpi=dpi, facecolor='w', edgecolor='w')

    # Get datetime
    if isinstance( dlist, type(None) ):
        dlist = get_gc_datetime( wd=wd )

    # set cb ranges for whole data period
    if isinstance( fixcb, type(None) ):
        fixcb  = [( i.min(), i.max() ) for i in [arr]][0]

    # Create figure if not provided
    if isinstance( fig, type(None) ):
        fig = plt.figure(figsize=(15, 10), dpi=dpi,
                    facecolor='w', edgecolor='w')

    # Set readable levels for cb, then use these to dictate cmap
    lvls = get_human_readable_gradations( vmax=fixcb[1],  \
                    vmin=fixcb[0], nticks=nticks,
                    sigfig_rounding_on_cb=sigfig_rounding_on_cb  )

    # Setup Colormap
    cmap, fixcb_buffered = get_colormap( np.array( fixcb ), \
        nticks=nticks, fixcb=fixcb, buffer_cmap_upper=True )

    if discrete_cmap:
        cmap, norm = mk_discrete_cmap( nticks=nticks,\
                vmin=fixcb[0], vmax=fixcb[1], cmap=cmap )

    if debug:
        print([ (i.min(), i.max(), i.mean()) for i in [ arr ] ])

    # Get time in the troposphere diagnostic if not provide as agrument
    if isinstance( t_ps, type(None) ):
        t_ps = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'], trop_limit=True )

    # Loop thorugh months
    for m, month in enumerate(dlist):

        # add axis
        axn =[ 4, 3, m+1 ]
        ax = fig.add_subplot( *axn )

        # set when to use y and x labels
        xlabel=False
        if m in range(12)[-3:]:
            xlabel=True
        ylabel=False
        if m in range(12)[::3]:
            ylabel=True

        # Plot zonally
        zonal_plot( arr[...,m], fig,  ax=ax, set_window=set_window, log=log,\
            format=format, cmap=cmap, lat_0=lat_0, lat_1=lat_1, \
            fixcb=fixcb, f_size=f_size*.75, res=res, norm=norm, \
            fixcb_buffered=fixcb_buffered, no_cb=True, trop_limit=True,\
            window=window, interval=interval, xlabel=xlabel, ylabel=ylabel,\
            verbose=verbose, debug=debug )

        # Only show troposphere
        greyoutstrat( fig, t_ps.mean(axis=0).mean(axis=-1), axn=axn, res=res )

        # add month
        plt.title(month.strftime("%b"), fontsize=f_size*1.5)

    # Add single colorbar
    mk_cb(fig, units=units, left=0.895, cmap=cmap, vmin=fixcb[0], \
        vmax=fixcb[1], nticks=nticks, f_size=f_size*1.25, extend=extend,
         width=0.015*1.5, height=.95, bottom=0.11 )
    if debug:
        print(nticks, fixcb, lvls)

    # sort out ascetics -  adjust plots and add title
    fig.subplots_adjust( bottom=bottom, top=top, left=left, \
        right=right, hspace=hspace, wspace=wspace)
    if fig_title:
        fig.suptitle(  '{}'.format( latex_spec_name(spec) ), fontsize=f_size*2,
            x=.55 , y=.95  )

    #  save as pdf ?
    if pdf:
        plot2pdf( title=pdftitle )
    if show:
        plt.show()


# --------
# 1.38 - Stackplot for variables over X axis
# --------
def X_stackplot( X=None, Y=None, labels=None, baseline='zero', \
        fig=None, ax=None, dpi=160, show=False, f_size=10, legend=False, \
        colors=None, title=None, loc='upper right', ylim=None, xlim=None, \
        lw=8.0, ylabel=None, xlabel=False, log=False, rm_ticks=False, \
        alt_text_x=.15, alt_text_y=0.75, alt_text=None, ncol=1, pcent=False,  \
        stacked=False, verbose=False, debug=False):
    """
    Make a stacked plot (by X axis) for values in Y array.

    Parameters
    -------
    X (list): list of numpy arrays to plot as X
    Y (array): must be a numpy array use use as Y
    labels (list): list of labels for stack data
    baseline (str): if =='zero' then start filling from zero.

    Returns
    -------
    (None)

    Notes
    -----
     - MPL only contains a stackplot function for Y axis, this function is
    based on that code
    (https://github.com/matplotlib/matplotlib/blob/3ae7c2a32ddd9809552315458da1dd70afec1b15/lib/matplotlib/stackplot.py )
    """
    logging.info('X_stackplot called, X[0] & Y[0] shape={}.{}'.format( \
        *[i.shape for i in (X[0], Y) ]) )

    # --- Fig and ax provided? Otherwise create these...
    if isinstance( fig, type(None) ):
        fig = plt.figure( figsize=(8,8), dpi=dpi, facecolor='w', \
            edgecolor='w')
        logging.info('Creating figure' )

    if isinstance( ax, type(None) ):
        ax = fig.add_subplot( 1,1,1  )
        logging.info('Creating ax' )
    logging.debug( '{}'.format( list(zip( labels, [np.sum(i) for i in X] )) )  )

    # --- Stack arrays if not stacked...
    if isinstance( X, np.ndarray ):
        X = np.ma.atleast_2d( X )
        logging.debug('Array passed has been made at least 2D if nessesary')
    else:
        X = np.ma.column_stack( X )
        logging.debug('WARNING - Stacking X data, X shape={}'.format( X.shape))
    # Assume data passed has not been 'stacked', so stack it here.
#    if stacked:
#        stack =  X
#    else:
    stack = np.ma.cumsum( X, axis=1)
    # convert to %?
#    pcent = True
    logging.info('stack shape={}'.format(stack.shape) )
    if pcent:
        max =  np.ma.max( stack, axis=1 ) # get accumulated maximum
        print([ (i.min(), i.max(), i.mean(), i.shape) for i in (stack, max) ])
        stack = np.ma.divide( stack,  max[:,None] ) *100
        print([ (i.min(), i.max(), i.mean(), i.shape) for i in (stack, max) ])
        xlim = [ 0, 100 ]

    if debug:
        print(list(zip( labels, [np.sum(stack[:,n]) for n, i in enumerate(labels) ] )))
        print(list(zip( labels, [np.max(stack[:,n]) for n, i in enumerate(labels) ] )))

    # --- Setup baseline ( can expand to include other options... )
    if baseline == 'zero':
        first_line = np.zeros( stack[:,0].shape)

    # --- Plot by label
    # Get list of colors
    if isinstance( colors, type(None) ):
        colors = color_list( len(stack[0,:]) )
    logging.debug( '{}'.format( list(zip( labels, [ colors[:len(labels)] ] )) ))

    # Color between x = 0 and the first array.
    logging.debug( '{}'.format(stack[:, 0]) )
    logging.debug( 'len colors={}, labels={}, stack={}'.format(
        *[len(i) for i in (colors, labels, stack[0, :]) ]) )
    r =[]
    r += [ ax.fill_betweenx( Y, first_line, stack[:, 0],
                               color=colors[0],
                                label=labels[0]) ]
    # Color between array i-1 and array i
    r +=  [ ax.fill_betweenx(Y, stack[:,i], stack[:,i+1],
                                   color=colors[i+1],
                                   label=labels[i+1] )
                                   for i in range( 0, len(stack[0,:])-1 ) ]

    # Plot transparent lines to get 2D line object to create legend
    [ plt.plot( Y, stack[:,n], alpha=0, color=colors[n], label=i) \
        for n,i in enumerate(labels) ]

    # Log scale?
    if log:
        ax.set_xscale('log')

    # Print maxima
    if debug:
        print(title, [ [ (i.min(), i.max(), i.mean() ) for i in [stack[:,n] ] ]
            for n, label in enumerate(labels) ])

    # --- Beautify plot
    if not isinstance( ylim, type(None) ):
        plt.ylim( ylim )
    if not isinstance( xlim, type(None) ):
        plt.xlim( xlim )
    if not isinstance( title, type(None) ):
        plt.title( title, fontsize=f_size )
    # Add alt text?
    if not isinstance( alt_text, type(None) ):
        ax.annotate( alt_text , xy=(alt_text_x, alt_text_y), \
            textcoords='axes fraction', fontsize=f_size*1.5 )
    if legend:
        # Add legend
        if ncol == 0:
            leg = plt.legend( loc=loc, fontsize=f_size*.75,  )
        else:
            import itertools
            def flip(items, ncol):
                return itertools.chain(*[items[i::ncol] for i in range(ncol)])
            handles, labels = ax.get_legend_handles_labels()
            leg = plt.legend( flip(handles, ncol), flip(labels, ncol), loc=loc,
                ncol=ncol, fontsize=f_size*0.75)

        # ( + update line sizes)
        for legobj in leg.legendHandles:
                legobj.set_linewidth( lw)
                legobj.set_alpha( 1 )

    # Remove tick labels on y axis?
    if ylabel:
        plt.ylabel( ylabel, fontsize=f_size*.75  )
        ax.tick_params( labelsize= f_size*.75 )
    else:
        ax.tick_params( axis='y', which='both', labelleft='off', \
            labelsize= f_size*.75)
    # Remove tick labels on x axis?
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=f_size*.75)
        ax.tick_params(  axis='x', which='both', labelsize= f_size*.75 )
    else:
        if rm_ticks:
            ax.tick_params( axis='x', which='both', labelbottom='off', \
                labelsize= f_size*.75 )
        else:
            ax.tick_params( axis='x', which='both', labelsize= f_size*.75 )

    if show:
        plt.show()

# --------------------------- Section 4 ------------------------------------
# -------------- Plotting Ancillaries
#


# -------------
# X.XX -  color list for rainbow plots
# -------------
def color_list(length, cb='gist_rainbow' ):
    """
    Create a list of colours to generate colors for plots contain multiple datasets
    """
    cm = plt.get_cmap(cb)
    color = [cm(1.*i/length) for i in range(length)]
    return color


# -------------
# X.XX - R squared - credit: Software carpentry
# -------------
def r_squared(x, y):
    """
    Return R^2 where x and y are array-like.

    NOTES:
     - actual_mean = np.mean(actual)
     - ideal_dev = np.sum([(val - actual_mean)**2 for val in ideal])
     - actual_dev = np.sum([(val - actual_mean)**2 for val in actual])
     - return ideal_dev / actual_dev

    """
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return r_value**2


# -------------
# X.XX - setup box plots
# -------------
def set_bp( bp, num, c_list=['k', 'red'], white_fill=True, set_all=True,
        median_color='white', linewidth=2, debug=False ):
    """
    Manually set properties of boxplot ("bp")
    """
    if debug:
        print(num, c_list)
    if set_all:
        setp(bp['boxes'][:], color=c_list[num], linewidth=linewidth*.5)
        setp(bp['caps'][:], color=c_list[num], linewidth=linewidth*.5)
        setp(bp['whiskers'][:], color=c_list[num], linewidth=linewidth*.5)
        setp(bp['fliers'][:], color=c_list[num], linewidth=linewidth*.5)
        setp(bp['medians'][:], color=c_list[num], linewidth=linewidth*.5)
    if white_fill:
        [ box.set( facecolor = 'white') for box in bp['boxes'] ]
    else:
        [ box.set( facecolor = c_list[num] ) for box in bp['boxes'] ]
        setp(bp['medians'][:], color=median_color, linewidth=linewidth)


# -------------
# X.XX - Get all marker types
# -------------
def markers_list( rm_plain_markers=False ):
    """
    Create a list of available markers for use in plots
    """
    from matplotlib.lines import Line2D
    markers = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != ' ':
                markers.append(m)
        except TypeError:
            pass
    # remove the marker similar to a plain line.
    if rm_plain_markers:
        [ markers.pop(i) for i in [3, 5][::-1] ]

    return markers


# -------------
# X.XX linear (polyfit) trendline calculator for X-Y plot (with histograms)
# -------------
def Trendline( ax, X, Y, order=1, intervals=700, f_size=20, color='blue',
        lw=1, debug=False ):
    """
    Adds a trendline and legend instance to existing axis instance.

    Parameters
    -------
    ax (axis instance): (required) axis instance
    X, Y (array): arrays of X,Y values to perform regression on
    order (int): order of line to fit
    intervals (int): number of intervals to use
    f_size (float): font size
    color (str): color of line
    lw (float): line width of plotted trendline
    debug (boolean): legacy debug option, replaced by python logging

    Returns
    -------
    (None)
    """
    # Poly fit data
    params, xp = np.polyfit( X, Y, order  ), \
        np.linspace( min(np.ma.min(X), np.ma.min(Y) ), \
        max(np.ma.max(X), np.ma.max(Y) ), intervals )
    print(params, type(params))
    logging.debug( 'params: {}'.format(str(params)) )
    yp = np.polyval( params, xp )

    # Calculate regression fit
    r_sq = [ r_squared(X, i ) for i in [Y] ]

    # Plot up
    ax.plot( xp, yp, ls='--', lw=lw, color=color,
        label=' (R$^{2}$'+'={0:<,.3f}, y={1:,.3f}x+{2:.3f})'.format(\
        r_sq[0], params[0], params[1] )   )
    # Add legend to plot
    ax.legend()


# -------------
# X.XX - plot_gc_bin_bands - plot up Fast-J  bins  ( 7 longest nm bins)
# -------------
def plot_gc_bin_bands(facecolor='#B0C4DE'):
    """
    Plot highlighted lines of ("tropospheric") GEOS-Chem/Fast-J photolysis bins
    """
    vars = [ GC_var(i) for i in ('FastJ_lower' , 'FastJ_upper')]
    alphas = [ 0.3,0.1 ]*len(vars[0])
    [ plt.axvspan( vars[0][n], vars[1][n], facecolor=facecolor, \
        alpha=alphas[n]) for n in range(len(vars[0])) ]


# --------------
# X.XX - line styles
# -------------
def get_ls(num):
    """
    Get a list of available line styles
    """
    ls =[   ':', '--', '-.', '-', ':', '--', '-.', '-', ':', ':', '--', '-.', '-', ':', '--', '-.', '-', ':'
    ]
    return ls[:num]

# --------
# X.XX - takes time in troposphere diagnostic array (46, 47) overlayes
# --------
def greyoutstrat( fig,  arr, axn=[1,1,1], ax=None, cmap=plt.cm.bone_r, \
        res='4x5', rasterized=True, debug=False):
    """
    Grey out stratosphere in existing zonal plot.

    NOTES:
     - This is used to highlight tropospherically focused work. This function
     requires the array of "time in the troposphere" diagnostic (lon,lat, alt)
    """
    # Get overall vars
    lon, lat, alt = get_latlonalt4res( res=res )
    alt  = [ i[:len(arr[0,:])] for i in [ alt ]  ][0]

    # plot up a grey section for stratosphere
    if isinstance(ax, type(None)):
        ax = fig.add_subplot( *axn )
    arr = np.ma.around(arr)
    arr = np.ma.masked_where( arr >0, arr )
    p = ax.pcolor( lat, alt, arr.T, cmap=cmap)
    if rasterized:
        plt.gcf().set_rasterized(True)

# --------
# X.XX - adjust subplots
# --------
def adjust_subplots( fig, left=None, bottom=None, right=None, top=None, \
        wspace=None, hspace=None ):
    """
    Set subplot adjust in provide figure
    """
    # the left side of the subplots of the figure
    if isinstance( left, type(None)):
        left  = 0.125
    # the right side of the subplots of the figure
    if isinstance( right, type(None)):
        right = 0.9
    # the bottom of the subplots of the figure
    if isinstance( bottom, type(None)):
        bottom = 0.1
    # the top of the subplots of the figure
    if isinstance( top, type(None)):
        top = 0.9
    # the amount of width reserved for blank space between subplots
    if isinstance( wspace, type(None)):
        wspace = 0.2
    # the amount of height reserved for white space between subplots
    if isinstance( hspace, type(None)):
        hspace = 0.5
    # Adjust subplots
    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top, \
        wspace=wspace, hspace=hspace)


# --------
# X.XX - Iodine deposition mask (for sites... e.g Denmark, Germany, Norfolk )
# --------
def mask_not_obs( loc='Denmark', res='4x5', debug=False ):
    """
    provide a mask of all regions apart from the location given
    """

    # Start with all zeros
    arr=np.zeros(get_dims4res(res))

    # Get lats and lons of locations to keep...
    lats, lons  = get_obs_loc( loc )

    # Unmask locations
    lats = [ get_gc_lat(i, res=res) for i in lats]
    lons = [ get_gc_lon(i,res=res) for i in lons ]
    for n, lat in enumerate( lats ):
        arr[lons[n],lat,:] = 1

    return np.ma.not_equal(  arr, 1)

# --------------
# X.XX - Annotate grid
# -------------
def annotate_gc_grid(ax, res='4x5', f_size=6.5, \
        loc_list=[ [-9999,-9999] ], everyother=1, label_gc_grid=True):
    """
    Annotate grid with GEOS-Chem indices
    """

    # Get Vars
    lon, lat, alt = get_latlonalt4res(res=res)
    if  res == '0.5x0.666':
        adjust_window, interval, nticks, nbins, resolution,shrink  = 3, 0.5, \
             int(nticks/3), int(nbins/2), 'l', 0.6
    glat = [get_gc_lat(i, res=res) for i in lat ]
    glon = [get_gc_lon(i, res=res) for i in lon ]

    # set ticks to all
    ax.xaxis.set_ticks(  lon[::everyother] )
    ax.yaxis.set_ticks( lat[::everyother] )
    plt.xticks( fontsize=f_size*2.5  )
    plt.yticks(  fontsize=f_size*2.5 )
    plt.grid(True)

    # label grid boxes by number
    if label_gc_grid:
        for n, la in enumerate( lat ):
            for nn, lo in enumerate( lon ):
                if any( [ [glat[n], glon[nn] ] == i for i in loc_list ] ):
                    color = 'red'
                    fontweight='bold'
                else:
                    color= 'blue'
                    fontweight ='normal'
                ax.text( lo+1, la+1, '{},{}'.format(glon[nn], glat[n]),\
                     fontsize=f_size, color=color )

# --------
# X.XX - Function for centering colorbar
# --------
def shiftedColorMap(cmap, start=0, midpoint=0.5, lower=0, upper=1, \
        start_center=0.5, stop=1.0, maintain_scaling=True, arr=None, \
        name='shiftedcmap', npoints=257, verbose=True, debug=False ):
    """
    ORGINAL DESCRIPTION: Function to offset the "center" of a colormap. Useful
    for data with a negative min and positive max and you want the middle of
    the colormap's dynamic range to be at zero.

    Parameters
    ----------
    cmap (colormap): The matplotlib colormap to be altered
    start (float): Offset from lowest point in the colormap's range.
        Defaults to 0.0 (no lower ofset). Should be between 0.0 and `midpoint`.
    midpoint (float): The new center of the colormap. Defaults to 0.5
    (no shift). Should be between 0.0 and 1.0. In general, this should be
    1 - vmax/(vmax + abs(vmin)). For example if your data range from -15.0 to
    +5.0 and you want the center of the colormap at 0.0, `midpoint` should be
    set to  1 - 5/(5 + 15)) or 0.75
    stop (float): Offset from highets point in the colormap's range. Defaults
    to 1.0 (no upper ofset). Should be between `midpoint` and 1.0.

    Returns
    -------
    (colormap)

    Notes
    -----
     - adapted from stackoverflow to allow for maintaining scaling:
    original: "http://stackoverflow.com/questions/7404116/defining-the-  \
    midpoint-of-a-colormap-in-matplotlib "
    """
    if maintain_scaling:
        # Cut off extra colorbar above point to shift maxoum too...
        vmin, vmax= arr.min(), arr.max()

        # Symetric range
        eq_range = max(  abs(vmin), abs(vmax) ) *2
        # Actual range
        act_range = abs( vmin ) +abs( vmax  )
        cut_off = act_range/eq_range

        if midpoint > 0.5:
            start_center = start_center/cut_off
            lower, upper = start, cut_off
        if midpoint < 0.5:
            start_center = start_center*cut_off
            lower, upper = 1-cut_off, stop
        if midpoint != 0.5:
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list(\
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=lower, \
                b=upper ), cmap(np.linspace(lower, upper, npoints)))
            logging.debug( 'maintain scale for vars: {}'.format( \
                *[str(float(i)) for i in list([ eq_range, act_range, vmin, \
                vmax, abs(vmin), abs(vmax), cut_off, start_center, lower, \
                upper, midpoint, start, stop])] ) )
    cdict = { 'red': [], 'green': [], 'blue': [], 'alpha': [] }

    # Regular index to compute the colors
    reg_index = np.hstack([ np.linspace(start, start_center, 128, \
        endpoint=False), np.linspace(start_center, stop, 129, endpoint=True)  ])

    # Shifted index to match the data
    shift_index = np.hstack([ np.linspace(0.0, midpoint, 128, endpoint=False), \
        np.linspace(midpoint, 1.0, 129, endpoint=True) ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

# --------
# X.XX - Add colorbar to side of plot
# --------
def mk_cb( fig, units=None, left=0.925, bottom=0.2, width=0.015, height=0.6,\
        orientation='vertical', f_size=20, rotatecbunits='vertical', nticks=10, \
        extend='neither', norm=None, log=False, format=None, cmap=None,\
        vmin=0, vmax=10, cb_ax=None, ticklocation='auto', extendfrac=None, \
        sigfig_rounding_on_cb=2, lvls=None, discrete_cmap=False, boundaries=None, \
        verbose=True, debug=False ):
    """
    Create Colorbar based on recieved parameters.

    Parameters
    ----------
    fig (figure instance): figure to add colorbar to
    cb_ax (axis instance): axis for colorbar
    ticklocation (str): how to set tick locations (e.g. 'auto')
    units (str): label for colorbar
    rotatecbunits (str): rotation of colorbar units
    left, bottom, width, height (float): location to place colorbar
    orientation (str): orientation of colorbar
    f_size (float): font size
    nticks (int): number of ticks on colorbar
    extend (str): colorbar format settings ( 'both', 'min', 'both' ... )
    norm (norm object): normalisation to use for colourbar
    log (boolean): use log scale for colorbar
    format (str): format string for colorbar tick formating
    cmap (str): colormap instance
    vmin, vmax (float): miminium and maximum of colorbar
    extendfrac (float): fraction by which to extend "pointers" on colorbar
    sigfig_rounding_on_cb (int): significant figure rounding to use for colourbar
    lvls (list): manually provide levels for colorbar
    boundaries (list):  manually provide levels for discrete colorbar
    discrete_cmap (boolean): use a discrete instead of conitunous colorbar map
    verbose (boolean): legacy debug option, replaced by python logging
    debug (boolean): legacy debug option, replaced by python logging

    Returns
    -------

    Notes
    -----
     - This function allows for avoidance of basemap's issues with spacing
    definitions when conbining with colorbar objects within a plot
    """
    if debug:
        print('mk_cb called with: ', norm, vmin, vmax, log, lvls)

    # Get colormap (by feeding get_colormap array of min and max )
    if isinstance( cmap, type(None) ):
        cmap = get_colormap( arr=np.array( [vmin,vmax]   ) )

    if debug:
        print(left, bottom, width, height)
        print(vmin, vmax, orientation, ticklocation, extend, extendfrac)

    # Make new axis
    if isinstance( cb_ax, type(None) ):
        cb_ax = fig.add_axes([left, bottom, width, height])
        if debug:
            print('>'*5, left, bottom, width, height)

    # Setup normalisation for colorbar
    if isinstance( norm, type(None) ):
        if log:
            # Get logarithmically spaced integers
            lvls = np.logspace( np.log10(vmin), np.log10(vmax) , num=nticks)
            # Normalise to Log space
            norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax)

        else:
            # Normalise to linear space
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    if verbose:
        print(lvls, vmin, vmax, norm)

    if isinstance( lvls, type(None) ):
        # make graduations in colourbar to be human readable
        lvls = get_human_readable_gradations( vmax=vmax,  \
            vmin=vmin, nticks=nticks, \
            sigfig_rounding_on_cb=sigfig_rounding_on_cb  )

    # add an extra decimal place for values under 50,
    # and two for values under 1
    if isinstance( format, type(None) ):
        try:
#            if all([ i<=1 for i in [ len( str( abs( int(i) ) )  ) \
#                    for i in vmax, vmin ] ]):
#                format='%.2f'
#            elif all([ i<=2 for i in [ len( str( abs( int(i) ) )  ) \
#            if all([ i<=2 for i in [ len( str( abs( int(i) ) )  ) \
#                    for i in vmax, vmin ] ]):
#                format='%.1f'
#            if any([ i>=3 for i in [ len( str( abs( int(i) ) )  ) \
#                    for i in vmax, vmin ] ]):
#                format='%.0f'

            # change to scientific notation max+min abs values less than 1
            if ( ( abs( int( vmax)) == 0) and (abs( int( vmax)) == 0) ):
#                format='%.0E'  # hashed prior
                format='%.1E'
                f_size = f_size *.75

        # warning this will need to be restored?! - plotting routines use this
        # This exception is for highly masked arrays
        except:
            format='%.0E'
            # Kludge, set to limits of -500-500 (%) if overly masked.
            vmax, vmin = 500,-500
            extend = 'both'

    if debug:
        print(lvls, norm, extend, format, ticklocation)

    # Make cb with given details
    if discrete_cmap:
        extendfrac=.075  # This was the previous default, still valid? <= update
#        if debug:
#            print 'WARNING: adding extensions to colorbar'
        if log==True:
            print('Will not work as colrbar adjustment not configured'+\
                ' for log scales')
            sys.exit(0)
        else:
            print(lvls)
            if discrete_cmap:
                boundaries = lvls
            else:
                extend_points=(lvls.max()-lvls.min())*extendfrac
                boundaries = [lvls.min()-extend_points]+ list(lvls) + \
                    [extend_points+lvls.max()]

        # the following code isn't active, why is it here? <= update
        # increase the space allowed for the colorbar
        # (the extension  is not summative )
#        if orientation=='horizontal':
#            width += width*extendfrac
#        if orientation=='vertical':
#            height += height*extendfrac

        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, format=format,\
                norm=norm, ticks=lvls, extend=extend, extendfrac=extendfrac,\
                boundaries=boundaries, \
                spacing='proportional',
#                spacing='uniform',
                orientation=orientation, ticklocation=ticklocation)

    # Standard approach below
    else:
        if verbose:
            print(lvls, norm, boundaries, extend, orientation, ticklocation, \
                    cb_ax)
        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, format=format,\
                norm=norm, ticks=lvls, extend=extend, \
                orientation=orientation, ticklocation=ticklocation)


    if log:
        round_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))
        cb.set_ticks( [ float('{:.2g}'.format( t )) for t in lvls ] )
        labels = [ round_to_n( i, sigfig_rounding_on_cb) for i in lvls ]
        cb.set_ticklabels( [ format % i for i in labels] )

    # Set cb label sizes
    if not isinstance( units, type(None) ):
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(f_size)
        if rotatecbunits == 'vertical':
            cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size, \
                fontsize=f_size)
        else:
            cb.set_label( units, fontsize=f_size )
    # set tick sizes regardless whether units (labels) are provided
    cb.ax.tick_params( labelsize=f_size ) #, size=f_size )

    return cb_ax

# --------
# X.XX - Create base map for plotting
# --------
def get_basemap( lat, lon, resolution='l', projection='cyl', res='4x5',\
        everyother=1, f_size=10, interval=1, axis_titles=False, \
        show_grid=True, drawcountries=False, ylabel=True, xlabel=True ):
    """
    Creates a basemap object.

    This should be used for first slide in python animated videos to save computation
    """

    m = Basemap(projection=projection,llcrnrlat=lat[0],urcrnrlat=lat[-1],\
            llcrnrlon=lon[0],urcrnrlon=lon[-1], resolution=resolution )

    if axis_titles:
        plt.ylabel('Latitude', fontsize = f_size*.75)
        plt.xlabel('Longitude',fontsize = f_size*.75)
    if (res == '0.5x0.666') or drawcountries :
        m.drawcountries()
    parallels = np.arange(-90,91,15*interval)
    meridians = np.arange(-180,181,30*interval)
    if (res == '0.25x0.3125') :
        parallels = np.arange(-90,91,15*interval/1.5  )
        meridians = np.arange(-180,181,30*interval/3)

    # use small font size for greater the runs with more axis labele
    f_size = f_size*.75

    # draw meridian lines
    plt.xticks( meridians[::everyother], fontsize = f_size )
    # draw parrelel lines
    plt.yticks( parallels[::everyother], fontsize = f_size )
    m.drawcoastlines()
    plt.grid( show_grid )
    # remove tick labels on y axis
    if not ylabel:
        ax_tmp = ax_tmp = plt.gca()
        ax_tmp.tick_params( axis='y', which='both', labelleft='off')
    if not xlabel:
        ax_tmp = ax_tmp = plt.gca()
        ax_tmp.tick_params( axis='x', which='both', labelbottom='off')

    return m

# --------
# X.XX - Provide an appropriate colormap for given data
# --------
def get_colormap( arr,  center_zero=True, minval=0.15, maxval=0.95, \
        npoints=100, cb='CMRmap_r', maintain_scaling=True, \
        negative=False, positive=False, divergent=False, \
        sigfig_rounding_on_cb=2, buffer_cmap_upper=False, fixcb=None, nticks=10,  \
        verbose=True, debug=False ):
    """
    Create *correct* colormap by checking if it contains just +ve or -ve or
    both then prescribing color map accordingly.

    Parameters
    ----------
    arr (array): array of values to assess colourbar from
    center_zero (boolean): make sure (divergent) colorbar centered around zero
    minval, maxval (float): values to restrict 'gnuplot2' to
    npoints (int): number of points in colormap
    cb (str): name of colorbar string
    maintain_scaling (boolean): maintain scaling for range in color change
    negative (boolean): force colormap to be sequential negative (==True)
    positive (boolean): force colormap to be sequential positive (==True)
    divergent (boolean): force colormap to be divergent (==True)
    sigfig_rounding_on_cb (int): number of sig. figs. to round colourbar ticks
    buffer_cmap_upper (boolean): make sure colorbar has space for maxiumium val.
    fixcb (array): lower and upper values to fix colourmap to.
    nticks (int): number of ticks to use for colorbar
    verbose (boolean): legacy debug option, replaced by python logging
    debug (boolean): legacy debug option, replaced by python logging

    Returns
    -------
    (colormap instance)

    Notes
    -----
     - this function also will can adjust colormaps to fit a given set of ticks
    """
    # Mannual fix maintain scaling to False
#    maintain_scaling=False

    logging.info( 'get_colormap called with fixcb={}'.format(str(fixcb)) )
    # Manually override colourbar?
#    cb='Blues' # Kludge.
#    cb='Reds' # Kludge.

    # Make sure cmap includes range of all readable levels (lvls)
    # i.e head of colormap often rounded for ascetic/readability reasons
    if buffer_cmap_upper:
        lvls, lvls_diff  = get_human_readable_gradations( vmax=fixcb[1],  \
            vmin=fixcb[0], nticks=nticks,  rtn_lvls_diff=True, \
            sigfig_rounding_on_cb=sigfig_rounding_on_cb   )

        # increase maximum value in color by 5% of level diff
        # to allow space for max lvl
        fixcb_ = ( fixcb[0],  lvls[-1]+ ( lvls_diff*0.05 ))
        arr = np.array( fixcb_ )

        # Gotcha - make sure max(lvls)  not increased > 0
        if (max(lvls) <= 0) and ( not (arr.max() < 0 ) ):
            arr = np.array(  [ arr[0], 0 ] )

    # make sure array has a mask
    logging.debug( "arr min and max vals: {}, {} - arr type: {}".format( \
        str(arr.min()), str(arr.max()), type(arr)) )
    if 'ask' not in str(type( arr ) ):
        arr = np.ma.array(arr)
#        s_mask = arr.mask
    logging.debug( 'arr type (post mask check) {}:'.format(type(arr)) )

    # If postive/negative not given, check if +ve/-ve
    if (not positive) or (not negative):
        logging.debug( 'Testing if arr is +ve/-ve, for arr with min' +\
            '{} and max {}'.format( arr.min(), arr.max() ) )

        # --- sequential ?
        # negative?
        arr.mask =False
        arr.mask[arr<=0]=True
        if arr.mask.all():
            negative = True

        # postive?
        arr.mask =False
        arr.mask[arr>=0]=True
        if arr.mask.all():
            positive = True

    # Reverse colourbar if negative
    if negative:
        if cb == 'CMRmap_r':
            cb = cb[:-2]
        else:
            cb = cb+'_r'
    # log colorbar details
    logging.info( 'cmap is: >{}< & data is:'.format( cb ))
    logging.info( '< postive == {}, negative == {}, divergent == {} >'.format(\
        positive, negative, (( not positive) and (not negative)) ))

    # load color map
    cmap = plt.get_cmap( cb )

    # Chop off bottom for 'gnuplot2'
    if (negative and ( 'gnuplot2' in cb)) or (positive and ( 'CMRmap' in cb)):
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list( \
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=0, \
            b=maxval-minval, ), cmap(np.linspace(0, maxval-minval, npoints)))
    if (positive and ( 'gnuplot2' in cb)) or (negative and ( 'CMRmap' in cb)):
#    if positive and ( 'CMRmap' in cb ):
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list( \
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval,\
                b=maxval), cmap(np.linspace(minval, maxval, npoints)))

    # --- divergent
    if ( ( not positive) and (not negative) ) or divergent:
        logging.debug('Data is divergent' )
        arr.mask = False
#        cmap = plt.cm.RdBu_r
        cb = 'RdBu_r'
        cmap = plt.get_cmap( cb )
#        cmap = plt.cm.Spectral
        # Centre color bar around zero
        if center_zero:
            vmin, vmax = arr.min(), arr.max()
            if buffer_cmap_upper:
                pass
            logging.debug('vals. for mid point {}, min {}, max {}'.format(\
                str(1-vmax/(vmax + abs(vmin))), str(vmin), str(vmax)) )
            cmap = shiftedColorMap(cmap, midpoint=1-vmax/(vmax + abs(vmin)), \
                maintain_scaling=maintain_scaling, arr=arr, name='shifted', \
                verbose=verbose, debug=debug )

    # Use updated jet replacements from mpl
    # : https://github.com/bids/colormap
#    cmap  = cmc  # pink, orange, blue alterative...
#    cmap  = cmd  # green - blue alternative

#    arr.mask = s_mask

    logging.debug('colorbar used: {} (center_zero=={})'.format(cb, center_zero))

    if buffer_cmap_upper:
        return cmap, fixcb_
    else:
        return cmap

# --------
# X.XX -Make segments for variable line color plot
# --------
def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


# --------
# X.XX - Make colored line for plot
# --------
def colorline( x, y, z=None, cmap=plt.get_cmap('copper'),  \
        norm=None, linewidth=3, alpha=1.0, ax=None, fig=None, \
        cb_title='Modelled wind direction', convert_x_2datetime=True,
        debug=False ):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    if isinstance( norm, type(None) ):
        norm = plt.Normalize(0.0, 1.0),

    # Special case if a single number:
    # to check for numerical input -- this is a hack
    if not hasattr(z, "__iter__"):
        z = np.array([z])

    z = np.asarray(z)

    if debug:
        print([ type(i[0]) for i in (x, z, y) ])
        print([ i[0] for i in (x, z, y) ])

    # convert to datetime.datetime
    if convert_x_2datetime:
        # --- make sure x axis is float ( not a datetime )
        # convert to dataframe
        df = DataFrame( data=y, index=x)
        x = np.array(df.index.to_pydatetime(), dtype=datetime.datetime)
        if debug:
            print(type(x[0]), x[0])

        # convert in float form
        x = matplotlib.dates.date2num(x)
        if debug:
            print(type(x[0]), x[0])

    # convert x and y into indices
    segments = make_segments(x, y)

    lc = mpl.collections.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)

    if isinstance( ax, type(None) ):
        ax = plt.gca()
    ax.add_collection(lc)

    if not isinstance( fig, type( None) ):
        axcb = fig.colorbar(lc)
        axcb.set_label( cb_title )

    return lc

# --------
# X.XX - Get human readable gradations for plot
# --------
def get_human_readable_gradations( lvls=None, vmax=10, vmin=0, \
        nticks=10, sigfig_rounding_on_cb=2, \
        sigfig_rounding_on_cb_ticks=2, \
        sigfig_rounding_on_cb_lvls=2, rtn_lvls_diff=False, \
        verbose=True, debug=False ):
    """
    Get human readible gradations for ploting ( e.g. colorbars etc ).
    OUTPUT:
    lvls: list of values where the ticks should be
    ticks: list of strings to call the ticks in Sig figs.
    """

    if not np.isfinite(vmin) or not np.isfinite(vmax):
        cb_error ="Colourbar has a NaN in it. vmin={vmin}, vmax={vmax}"\
                .format(vmin=vmin, vmax=vmax)
        logging.error(cb_error)
        raise ValueError(cb_error)

    logging.debug('get_human_readable_gradiations called with the following:')
    logging.debug('vmin = {vmin}, vmax = {vmax}, lvls = {lvls}'\
            .format(vmin=vmin, vmax=vmax, lvls=lvls))

    if isinstance(lvls, type(None)):
        lvls = np.linspace( vmin, vmax, nticks, endpoint=True )

#    verbose=True
    # --- Adjust graduations in colourbar to be human readable
    # in both min and max have absolute values less than 0, then sig figs +1
    # Find the amount of significant figures needed to show a difference
    # between vmin and vmax.
    if vmin==vmax:
        logging.error("There is no difference between vmin and vmax!")
        raise ValueError("There is no difference between the min and max of the data")
    try:
        if ( ( abs( int( vmin)) == 0) and (abs( int( vmax)) == 0) ):
            sigfig_rounding_on_cb += 1
        logging.debug("Significant figures needed for plot is {sf}"\
                .format(sf=sigfig_rounding_on_cb))
    except np.ma.core.MaskError:
        print('Gotcha: numpy.ma.core.MaskError')
        print(lvls, vmin, vmax)


#    # bjn updated sigfig finder
#    # Use logs to find out sig figs needed
#    if vmin==0 or vmax==0:
#        sig_figs_needed = 3
#    else:
#        log_diff = abs( np.log10(abs(vmax)) - np.log10(abs(vmin)) )
#        sig_figs_needed = int(np.ceil(abs(np.log10( log_diff ))))
#
#    sigfig_rounding_on_cb_ticks = sig_figs_needed

    # significant figure ( sig. fig. ) rounding func.
    round_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))
#    round_to_n = lambda x, n: get_sigfig(x,n)

    # --- Get current gradations
#    if debug:
#        print abs(lvls[-4])-abs(lvls[-3]), abs(lvls[-4])-abs(lvls[-3]), lvls,\
#                     sigfig_rounding_on_cb
    try:
        lvls_diff =[ round_to_n( abs( i-l[n+1] ), sigfig_rounding_on_cb_ticks) \
            for n, i in enumerate( l[:-1] ) ]
        lvls_diff = list( set( lvls_diff ) )
        if len(lvls_diff) > 1:
            lvls_diff = max( lvls_diff )
#        lvls_diff = round_to_n( abs(lvls[-3])-abs(lvls[-4]), \
#                                sigfig_rounding_on_cb_ticks)

    # handle if values (2,3) are both negative or abs. of both <0
    except:
        if debug:
            print(abs(lvls[-4])-abs(lvls[-3]), sigfig_rounding_on_cb_ticks)
        try:    # handle if values (2,3) are both negative
            lvls_diff = round_to_n( abs(lvls[-4])-abs(lvls[-3]), \
                                sigfig_rounding_on_cb_ticks)
        except: # If both absolute of vmin and vmax  are <0 ( and +ve )
            if debug:
                print(lvls, lvls[-3], lvls[-4], sigfig_rounding_on_cb_ticks)
            lvls_diff = round_to_n( lvls[-3]-lvls[-4], \
                                sigfig_rounding_on_cb_ticks)


    # ---  Round top of colorbar lvls, then count down from this
    # first get top numer rounded up to nearest 'lvls_diff'
    # caution, this may result in a number outside the cmap,
    # solution: use get_colormap, with buffer_cmap_upper=True
    #  if values are >0,
    if vmax > lvls_diff:
        if debug:
            print(vmax, lvls_diff)
        vmax_rounded = myround( vmax, base=lvls_diff,  integer=False )
        vmax_rounded = round_to_n( vmax_rounded, sigfig_rounding_on_cb)
    else:
        # <= update needed! - add function to round negative numbers
        # ( this method also fails if vmax<lvls_diff )
        vmax_rounded = vmax

#    if debug:
#        print 1, lvls, vmax_rounded, lvls_diff, sigfig_rounding_on_cb_lvls

    lvls = np.array([ vmax_rounded - lvls_diff*i \
            for i in range( nticks ) ][::-1])

    logging.debug("colorbar levels are: {lvls}".format(lvls=lvls))
#    if debug:
#        print lvls, len( lvls )
#        print 2, lvls, vmax_rounded, lvls_diff, sigfig_rounding_on_cb_lvls

    # ensure returned ticks are to a maximum of 2 sig figs
    # ( this only works if all positive ) and are unique
#    try:
#        # Make sure the colorbar labels are not repeated
#        invalid = True
#        while invalid:
#            new_lvls = [ round_to_n( i, sigfig_rounding_on_cb_lvls) \
#                for i in lvls ]
#            if len( set(new_lvls) ) == len(lvls):
#                lvls = new_lvls
#                invalid=False
#            else: # Try with one more sig fig
#                sigfig_rounding_on_cb_lvls += 1
#
#    except:
#        print 'WARNING: unable to round level values to {} sig figs'.format(\
#                   sigfig_rounding_on_cb_lvls  )

    new_lvls = []
    for level in lvls:
        new_lvls.append(get_sigfig(level, sigfig_rounding_on_cb_lvls))

    lvls = new_lvls

#    if any([ (type(i) == str) for i in lvls] ):
#        logging.debug('WARNING str in list of levels, tmp conversion to float'+
#         ' - This is due to error in get_sigfig - seting "0.0" to float' )
#        for n, i in enumerate( lvls ):
#            if i == '0.0':
#                lvls[n] = 0.0

#    if debug:
#        print 3, lvls, vmax_rounded, lvls_diff, sigfig_rounding_on_cb_lvls
#    print [(i, type(i)) for i in vmin, vmax ], lvls
#    print [(i, type(i)) for i in lvls ]


    if rtn_lvls_diff:
        return lvls, lvls_diff
    else:
        return lvls

# --------
# X.XX - mk colourmap discrete
# --------
def mk_discrete_cmap( lvls=None, cmap=None, arr=None,\
        vmin=0, vmax=10, nticks=10, debug=False ):
    """
    Make a discrete colormap from an existing cmap
    NOTE:
     - the data will now need to normalised the range of lvls
    """

    # define bins
    if isinstance( lvls, type(None) ):
        bounds = np.linspace( vmin, vmax, nticks, endpoint=True )
    else:
        bounds = lvls

    # get colormap if not provided
    if isinstance( cmap, type(None) ):
        cmap = get_colormap( np.array([vmin, vmax]) )

    if debug:
        print(lvls, vmin, vmax, nticks)

    # extract colors
#    cmaplist = [ cmap(i) for i in range( len(lvls ) ) ]
    cmaplist = [cmap(i) for i in range(cmap.N)]

    # force the first color entry to be grey
#    cmaplist[0] = (.5,.5,.5,1.0)

    # create the new discrete cmap
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N )

    # make norm... -  define the bins and normalize
    norm = mpl.colors.BoundaryNorm( bounds,  cmap.N)

    return cmap, norm


# --------
# X.XX -
# --------
def show_plot():
    """
    Wrapper for plt.show(). Use to plot to screen.
    """
    plt.show()
    return

# --------
# X.XX -
# --------
def save_plot(title="myplot", location=os.getcwd(),  extensions=['png'], tight=True):
    """
    Save a plot to disk.

    Parameters
    ----------
    title (str): plot title
    location (str): String of directory to save output
    extensions (list):  List of strings of file extensions. (e.g. ['png'])
    tight (boolean): use tight layout - redundent?

    Returns
    -------
    (None)
    """

#    if tight:
#        plt.tight_layout()

    if not os.path.isdir(location):
        os.mkdir(location)
        logging.warning("Plotting location not found")
        logging.warning("Creating dir {folder}".format(folder=location))

    for extension in extensions:
        filename = os.path.join(location, title+"."+extension)
        plt.savefig( filename )
        logging.info("Plot saved to {location}".format(location=filename))
    return




# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# ---------------- Section X -------------------------------------------
# -------------- User specific functions.
# --------------------------------------------------------------------------
#
# NOTE(s):
# (1) These functions should be removed following checks on compatibility
#

# -------------
# X.XX -  print NCAS & York logos in the bottom corners
# -------------
def add_logos_NCAS_york_bottom(fig):
    """
    Add NCAS + York logo for external plots used externally

    NOTE(s):
     - This function is not general enough and needs to be removed from AC_tools
    """
    from PIL import Image

    wd1 = get_dir ('dwd') +'misc/logos/'
    logo_list= [
    'NCAS_national_centre_logo.gif', 'nerclogo1000.gif' , 'york_uni_shield.tif' ,   \
    'york_uni_name.tif'
    ]

    # Setup logo 1
    logo1 = Image.open( wd1 + logo_list[3])  # York name
    # [left, bottom, width, height]
    ax2=fig.add_axes([0.01, 0.015, 0.55, 0.055], frameon=False)
    ax2.imshow(logo1,interpolation="bilinear")
    ax2.axis('off')

    # Setup logo 2
    logo2 = Image.open( wd1 + logo_list[0]) # NCAS logo
    # [left, bottom, width, height]
    ax3 = fig.add_axes([0.7, 0.01, 0.25, 0.1], frameon=False)
    ax3.imshow(logo2,interpolation="bilinear")
    ax3.axis('off')
    return fig


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# ---------------- Section X -------------------------------------------
# -------------- Redundant Functions
# --------------------------------------------------------------------------
#
# NOTE(s):
# (1) These are retained even though they are redundant for back compatibility
# (2) It is not advised to use these.

# --------------
# 1.16 - plot up timeseries from May through Septmeber
# -------------
def timeseries_by_day_plot( ax, dates, data, f_size=20, pos=0, posn=1,  \
        title=None, legend=False, everyother=7,  x_nticks=12, \
        window=False, label=None, ylabel=None, loc='upper right',  \
        lw=1,ls='-', color=None, start_month=5, end_month=9,
        boxplot=True, showmeans=False, debug=False ):
    """
    Timeseries plot.

    ARGUEMENTS:
     - Takes numpy arrays of datetimes and data as arguemnts.

    NOTES:
     - Description needs update.
    """

    # Process data - reduce resolution to daily, and get std
    df = DataFrame(data={'data':data},index=dates )

    # remove dates outside of range (start_month > < end_month )
    def get_month(x):
         return x.month
    df[ 'month'] = df.index.map( get_month )
    df = df[ df['month']<=end_month ]
    df = df[ df['month']>=start_month ]

    # split by day
    def get_day(x):
         return datetime.datetime(*x.timetuple()[:3])
    df[ 'day'] = df.index.map( get_day )
    sel_days = sorted(set(df[ 'day']))
    df = DataFrame(data={'data':df.data, 'day':df.day},index=df.index )
    daily = [ df[df['day']== i].data for i in sel_days ]
    sel_days = [i.to_datetime() for i in sel_days ]

    # label once per week
    labels = [i.strftime("%-d %b") for  i in sel_days ][::everyother]

    if boxplot:
        # output boxplot for range
        bp = ax.boxplot( daily, sel_days, showmeans=showmeans )

        # xticks are numbers, therefore set all to ''
        # except those you want to show dates
        fill_label = ['']*len(sel_days)
        fill_points =list(range(len(sel_days)))[::everyother]
        for n,i in enumerate( fill_points) :
            fill_label[i] = labels[n]

        plt.xticks( list(range(len(sel_days))), fill_label, \
                        rotation='vertical')
    else:
        plt.plot( sel_days, [np.mean(i) for i in daily] )
        plt.xticks( sel_days[::everyother], labels, \
                    rotation='vertical' )

    # Beatify plot
    if not isinstance( title, type(None) ):
        plt.title( title )
    if not isinstance( ylabel, type(None) ):
        plt.ylabel( ylabel )


# -------------
# X.XX - box X-Y plot with histograms
# -------------
def X_Y_hist( x, y, z=None, zlabel=None, fig=None, \
        left= 0.1, width=0.60, bottom=0.1, height=0.60, widthII=0.2, \
        fit = 0.02, binwidth=0.1, cmap=plt.cm.gnuplot2, \
        X_title='X title', Y_title='Y title', f_size=5 , line121=False ):
    """
    Plots a X vs. Y histogram
    NOTES
    ---
    Just use pandas for this or seaborn
     - http://seaborn.pydata.org/generated/seaborn.distplot.html
     - http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.hist.html

    """

    # Use hist code
    nullfmt = mpl.ticker.NullFormatter()         # no labels

    # start with a rectangular Figure
    if isinstance( fig, type(None) ):
        fig = plt.figure(1, figsize=(8,8))

    rect_scatter = [left, bottom, width, height]
    axScatter = plt.axes(rect_scatter)

    if z == None:
        # the scatter plot:
        axScatter.scatter(x, y)

        # definitions for the hist axes & setup.
#        bottom_h = left_h = left+width+fit
#        rect_histx = [left, bottom_h, width, widthII]
#        rect_histy = [left_h, bottom, widthII, height]

    else:
        vmin, vmax = [  ( float(np.ma.min(i)), float(np.ma.max(i)) ) \
            for  i in [ z] ][0]
        s_z = [ cmap(  (float(i) - vmin) / ( np.array([vmin,vmax]) ).ptp() ) \
            for i in z ]
        pts = axScatter.scatter(x, y, c=z)

        # definitions for the hist axes & setup.
#        widthII , height = [ i*.75 for i in widthII, height ]
        widthII  = widthII*.75

    # definitions for the hist axes & setup.
    bottom_h = left_h = left+width+fit
    rect_histx = [left, bottom_h, width, widthII]
    rect_histy = [left_h, bottom, widthII, height]

    # now determine nice limits by hand or prescribe:
    xymax = np.ma.max( [np.ma.max(np.fabs(x)), np.ma.max(np.fabs(y))] )
    lim = ( int(xymax/binwidth) + 1) * binwidth

    # manually?
    lim_min, lim_max  =  0,  max( np.ma.max(x), np.ma.max(y) )
    axScatter.set_xlim( lim_min, lim_max )
    axScatter.set_ylim(  lim_min, lim_max )
    # or by hand?
#    axScatter.set_xlim( (-lim, lim) )
 #   axScatter.set_ylim( (-lim, lim) )

    # add trendline
    Trendline( axScatter, x, y, order =1, intervals= 700 )

    # plot up titles
    plt.xlabel(  X_title )
    plt.ylabel( Y_title )
    if  line121:
        print(lim_min, lim_max)
        plt.plot( np.arange(lim_min, lim_max*1.5 )  , np.arange(lim_min, \
            lim_max*1.5  ), linestyle='--', color=(0.5, 0.5, 0.5))

    # add color bar if using z - [left, bottom, width, height],
    if z != None:
        cbaxes = plt.axes([0.9+fit*.5, bottom, fit*.5, height])
        cb = plt.colorbar(pts, cax = cbaxes)
        if zlabel != None:
            cb.ax.set_ylabel(zlabel, rotation='vertical', labelpad=1 )

    # plot up histograms
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # bin data
    bins = np.arange(-lim, lim + binwidth, binwidth)
    axHistx.hist(x, bins=bins)
    axHisty.hist(y, bins=bins, orientation='horizontal')
    plt.xticks( rotation='vertical' )
#    plt.setp( axHisty.xaxis.get_majorticklabels(), rotation=70 )

    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )

    # make sure scales are linear
    [ i.set_yscale('linear') for i in (axHisty, axScatter) ]
    [ i.set_xscale('linear') for i in (axHistx, axScatter) ]

    plt.rcParams.update({'font.size': f_size})
#    [ i.rcParams.update({'font.size': f_size}) for i in axHisty, axHistx, axScatter ]
#    [ i.xticks( fontsize=f_size ) for i in axHisty, axHistx, axScatter ]
#    [ i.yticks( fontsize=f_size ) for i in axHisty, axHistx, axScatter ]

# --------
# X.XX - Diurnal plot
# --------
def diurnal_plot_df(fig, ax,  dates, data, pos=1, posn =1, color=None, \
        xlabel=True, ylabel=True, label=None, title=None, f_size=10, \
        units='ppbv',lgnd_f_size=None, alpha=0.5, rotatexlabel=45, \
        loc='upper right', time_resolution_str="%H:%M", stat2plot='mean', \
        legend=False, ylim=None, lw=2.0, add_quartiles2plot=True, \
        debug=False ):
    """
    Creates a diurnal plot for given data and dates using pandas Dataframe

    Parameters
    -------
    data (array): numpy array of data
    dates (numpy array of datetime.datetime): dates to use for x axis
    fig (fig instance)
    ax (ax instance)
    stat2plot (str): stat (e.g. mean or medium) to use for main line
    xlabel, ylabel (str): label for axis?
    label (str): label for input data
    units (str): unites to label
    time_resolution_str (str): time string to reduce dataframe by
    legend (boolean): add a legend?
    title (str): title for plot
    color (str/hex etc): color for line
    f_size (float): fontsize
    lgnd_f_size (float): fontsize for legend
    loc (str): location for legend
    rotatexlabel (numnber/str): rotation of x axis labels
    pos, posn (int): vestigle(!) location indices for window plots
    ylim (list): min and max y axis limit
    lw (str): linewidth

    Returns
    -------
    (None)

    Notes
    -----
     - Adapted from David Hagen's example - https://www.davidhagan.me/articles?id=7
    """
    logging.info('diurnal_plot_df called with stat2plot={} '.format(stat2plot))
    # ---  process input data
    # Form a dataFrame from input numpy arrays.
    df = pd.DataFrame( {'data':data}, index=dates )

    # Add a time coluumn to group by (e.g. "%H:%M" for minutes)
    df['Time'] = df.index.map(lambda x: x.strftime(time_resolution_str))

    # Group by time resolution column
    df = df.groupby('Time').describe().unstack()

    # Convert X axis to strings
#    print df
    df.index = pd.to_datetime(df.index.astype(str))
#    df.index = pd.to_datetime(df.index )#.astype(str))
#    df.index = [ datetime.datetime( 2005, 1, 1, i ) for i in range(0, 24) ]

    # --- Aesthetics
    if isinstance(color, type(None)):
        color=color_list(posn)[pos-1]
    # lengend font size
    if isinstance( lgnd_f_size, type(None)):
        lgnd_f_size = f_size

    # --- Plot up average line (median or mean)
    if stat2plot == 'median':
        stat2plot='50%'
    ax.plot(df.index, df['data'][stat2plot], color=color, linewidth=lw, \
        label=label)

    # Add quartiles
    if add_quartiles2plot:
        ax.plot(df.index, df['data']['75%'], color=color, alpha=alpha)
        ax.plot(df.index, df['data']['25%'], color=color, alpha=alpha)

        # And shading for quartiles
        try:
            ax.fill_between(df.index, df['data'][stat2plot], df['data']['75%'],
                alpha=alpha, facecolor=color)
            ax.fill_between(df.index, df['data'][stat2plot], df['data']['25%'],
                alpha=alpha, facecolor=color)
        except:
            logging.info( 'Failed to add percentile shading' )

    # --- Beautify
    # title?
    if not isinstance(title, type(None) ):
        plt.title( title )
    # axis labels?
    if xlabel:
        plt.xlabel('Hour of day', fontsize=f_size )#*.75)
        # Add hourly ticks if labeling xaxis
        plt.xticks( rotation=rotatexlabel, fontsize=f_size )#*.75 )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H') )
    else:
        ax_tmp = ax_tmp = plt.gca()
        ax_tmp.tick_params( axis='x', which='both', labelbottom='off')
    if ylabel:
        plt.ylabel('{}'.format(units), fontsize=f_size)#*.75)
    else:
        ax_tmp = ax_tmp = plt.gca()

    # Add legend?
    if legend:
        plt.legend( fontsize=lgnd_f_size, loc=loc )

    if not isinstance(ylim, type(None)):
        plt.ylim(ylim)


# --------
# X.XX - Diurnal plot
# --------
def diurnal_plot(fig, ax,  dates, data, pos=1, posn =1,  \
        bin_size=2/24.,widths=0.01, rmax=False, \
        ls='-', color=None, fractional=False, diurnal=False, mean=True, \
        xlabel=True, r_avgs=False, marker=None, label=None, \
        markersize=1, title=None, f_size=10, units='ppbv', scale='linear', \
        lw=1,lgnd_f_size=None, alpha=1, debug=False ):
    """
    REDUNDENT!  (use diurnal_plot_df instead)
     Creates a diurnal plot for given data and dates.


    NOTES:
     - Data and dates must be in numpy array form.
     - Dates must also be datetime.datetime objects

     """

    # Convert datetime to fractional day <= set this to map on array at once?
    dates = np.array( [ get_day_fraction(i) for i in dates ] )

    # asectics
    ls =[ls]*5
    if posn> 4:
        ls=get_ls( posn )
    if color == None:
        color=color_list(posn)[pos-1]
    else:
        color=color
    if isinstance( lgnd_f_size, type(None)):
        lgnd_f_size = f_size

    # Bin data
    binned, bins_used = bin_data( data, dates, bin_size, debug=debug )

    # Take average of hourly binned data.
    if mean:
        avgs = np.ma.array([np.ma.mean(i) for i in binned]  )
    else:
        avgs = np.ma.array([np.ma.median(i) for i in binned]  )
    avg = np.ma.mean( avgs )
    max_ = np.ma.max( avgs )
#    print avgs, avg, np.ma.max( avgs )

    if fractional:
#        y = ( avgs - np.ma.max( avgs )  )  / avgs *100
#        y = ( avgs - np.ma.max( avgs )  )  / avg *100
        y = ( avgs - np.ma.max( avgs )  )  / max_*100
        print('test'*100, avg, max_, avgs)

        ymin, ymax =  -0.15, 0.025
        ymin, ymax =  [i*100 for i in (ymin, ymax) ]
        units = '%'
    elif diurnal:
        y =  avgs - np.ma.max( avgs )
        ymin, ymax =  -3.5 , 0.5
    else:
        y = avgs
        ymin, ymax =  None, None#27, 37

    # Plot
    plt.plot( bins_used, y, color=color , label=label, linestyle=ls[pos-1], \
        alpha=alpha, marker=marker, lw=lw, ms=markersize )

    # Beautify
    ax.set_xticklabels( np.arange(0,24,1 )[::2]  )
    plt.xticks( np.arange(0,1,1/24. )[::2], fontsize=f_size )
    plt.xlim(0., 23/24.)

    if ymin != None:
        plt.ylim( ymin, ymax )
    plt.yticks( fontsize=f_size*.75)
    plt.xticks( fontsize=f_size*.75)
    if (title != None):
        plt.title( title )
    if xlabel:
        plt.xlabel('Hour of day', fontsize=f_size*.75)
    plt.ylabel('{}'.format(units), fontsize=f_size*.75)

    # Apply legend to last plot
    if pos == posn:
        plt.legend( fontsize=lgnd_f_size )

    # return max
    if rmax :
        return np.ma.max( avgs )

    if r_avgs:
        return avgs

# --------
# X.XX - Lat plot
# --------
def lat_plot(fig, ax, arr, title=None, f_size=10, units='ppbv', \
        scale='linear', debug=False ):
    """
    Creates a latitude plot on given figure and axis
    """

    NIU, lat, NIU = get_latlonalt4res( res=res )
    del NIU
    plt.plot( lat, arr )
    plt.ylabel(units)
    plt.xlabel('Latitude' )
    if (title != None):
        plt.title(title)
    ax.set_yscale(scale)
    parallels = np.arange(-90,91,15)
    plt.xticks( parallels, fontsize = f_size ) # draw parrelel lines
    plt.rcParams.update({'font.size': f_size})


# --------
# X.XX - Diurnal boxplot
# --------
def diurnal_boxplot(fig, ax,  dates, data, pos=1, posn =1,  bin_size=2/24.,\
        widths=0.01, white_fill=True, alpha=0.1, linewidth=0.5, \
        showmeans=False, title=None, f_size=10, units='ppbv', \
        scale='linear', debug=False ):
    """
    Creates a diurnal plot of boxplots (hourly) for given data and dates.
    NOTES:
     - Data and dates must be in numpy array form.
     - Dates must also be datetime.datetime objects
    """

    # Convert datetime to fractional day <= set this to map on array at once?
    dates = np.array( [ get_day_fraction(i) for i in dates ] )

    # bin data
    binned, bins_used, b_all = avg_n_bin_y_by_x(data, dates, bin_size, \
        binned_data=True, debug=debug)

    # Generate positions
    positions=[ i+(bin_size*.75/posn)+( (bin_size*.75/posn)*float(pos)) \
        for i in bins_used ]

    # Plot
    bp = ax.boxplot( b_all,  positions=positions, widths=widths, \
        showmeans=showmeans, patch_artist=True )
    set_bp( bp, pos, white_fill=white_fill, c_list=color_list(posn) )

    # Beautify
    ax.set_xticklabels( np.arange(0,24,bin_size*24 )[::2]  )
    plt.xticks( np.arange(0,1,bin_size )[::2] )
    plt.xlim(-0.05, 1.05)
    plt.ylabel(units)
    if (title != None):
        plt.title(title)
    plt.xlabel('Hour of day')
    plt.ylabel('{}'.format(units))
    ax.xaxis.grid(True, which='major')
    ax.xaxis.grid(True, which='minor')

    # --- Highlight bins
    bs = np.arange(0, 24, bin_size )#[ bs[0] - bin_size ] + bs
    [ plt.axvline( x=i, color='k', linewidth=linewidth, alpha=alpha, \
         linestyle='dashed' ) for i in bs ]

# --------------
# X.XX - plot up monthly from data provided from DB netCDF
# -------------
def obs_month_plot(data, color=None, title=None, rtn_data=False, plt_day=True, debug=False):
    """
    Plot up seaonal (monthly ) data. Requires data, and dates in numpy
    array form. Dates must be as datetime.datetime objects.

    NOTES:
     - REDUNDENT? (see 1.13)
    """

    # Setup decimal day list
    day_time= [i+float(1/23) for i in range(23) ]
    if debug:
        print(data.shape)
        print(day_time)

    # Plot up all data <= this is inefficient.
    if plt_day:
        for i in range( len( data[0,:] ) ) :
            day_d = data[:,i]  # select day's data and convert to -1
            plt.plot(day_time , day_d , alpha=0.1, color=color)

    # Return data?
    if rtn_data:
        if debug:
            print('rtn_data')
#		day_time,
        data_ = np.ma.mean(data[:,:],axis=1)
        return day_time, data_
    # Plot up daily decimal days
    else :
        if debug:
            print(' not - rtn_data')
        plt.plot( day_time, np.ma.mean(data[:,:],axis=1) , lw=3 , color=color, \
             label=title)
        return plt


# -------------
# X.XX - Input for plotting
# -------------
def get_input_vars(debug=False):
    """
    Get input from command line arguments

    Parameters
    ----------
    debug (boolean): legacy debug option, replaced by python logging

    Returns
    -------
    (list): spec(str), filename(str), category(str), start date(tuple),
    end date  (tuple)
    """
    # If working directory (wd) provided, use command line arguments
    try:
        wd = sys.argv[1]
        print('$>'*5, sys.argv)

    # else prompt for input settings
    except:
        wd = input("ctm.bpch dir: ")
        spec = input("species (default = O3): ")
        if (len(spec) < 1):
            spec = 'O3'
        fn = input("ctm.bpch name (default = ctm.bpch): ")
        if (len(fn) < 1):
            fn = 'ctm.bpch'
        cat_ = input("Category (default = 'IJ-AVG-$'): ")
        if (len(cat_) < 1):
            cat_ = "IJ-AVG-$"
        start = input("Enter start time as 'YYYY,MM,DD' (default = all): ")
        if (len(list(start)) < 1):
            start = None
        else:
            start = datetime.datetime( *tuple( map( int, start.split('-')) ) )
        end = input("Enter end time as 'YYYY,MM,DD' (default = all): ")
        if (len(end) < 1):
            end = None
        else:
            end = datetime.datetime( *tuple( map( int, end.split('-')) ) )
        print(wd, fn, cat_, spec, start, end)

    # Has species file name been provided? ( e.g.  "O3" )
    try:
        spec = sys.argv[2]
    except:
        try:
            spec
        except:
            spec = 'O3'

    # Has ctm.bpch file name been provided? ( e.g.  "ctm.bpch" )
    try:
        fn = sys.argv[3]
    except:
        try:
            fn
        except:
            fn = 'ctm.bpch'

    # Has category been provided? ( e.g.  "IJ-AVG-$" )
    try:
        cat_ = sys.argv[4]
    except:
        try:
            cat_
        except:
            cat_ = "IJ-AVG-$"

    # Has start date been provided? ( in form "YYYY, MM, DD" )
    try:
        print(sys.argv[5].split('-'))
        start = datetime.datetime( *tuple( map( int, sys.argv[5].split('-')) ) )
    except:
        try:
            start
        except:
            start = None

    # Has end date been provided? ( in form "YYYY, MM, DD" )
    try:
        print(sys.argv[6].split('-'))
        end = datetime.datetime( *tuple( map( int, sys.argv[6].split('-')) ) )
    except:
        try:
            end
        except:
            end = None

    if debug:
        print('$>'*5, wd, fn, cat_, spec, start, end)

    # if final character not '/' add this
    if wd[-1]  != '/':
        wd += '/'

    return wd, fn, cat_, spec, start, end


# -------------
# X.XX - weighted average
# -------------
def weighted_average( data, interval , bins_used=False, debug=False):
    """
    Calculate a weighed average of data fro a given interval

    NOTES just use pandas.cut!
    docs: http://pandas-docs.github.io/pandas-docs-travis/generated/pandas.cut.html


    """

    min_v, max_v = int( data.min() ), int( data.max() )
    bins = np.arange(min_v, max_v, interval)
    if debug:
        print([ (i, type(i) ) for i in [ min_v, max_v, bins ] ])
    int_d = np.int_(data)
    binned, bins_used  = [], []
    for bin_ in bins:
        b_mean = np.mean( data[np.where( int_d == bin_ )] )
        if ( b_mean != 0 ):
            binned.append( b_mean )
            bins_used.append( bin_ )
        else:
            print('no data for bin {}'.format(bin_))
            pass
    if debug:
        print([ ( i, len(i) ) for i in [binned , bins] ])
    if (bins_used):
        return binned, bins_used
    else:
        return binned


# --------------
# X.XX - moving average
# -------------
def moving_average(x, n, type='simple'):
    """
    compute an n period moving average.

    NOTE:
     - type is 'simple' | 'exponential'
     - Just use pandas for this?
    (http://pandas.pydata.org/pandas-docs/version/0.17.0/generated/pandas.rolling_mean.html)
    """
    x = np.asarray(x)
    if type=='simple':
        weights = np.ones(n)
    else:
        weights = np.exp(np.linspace(-1., 0., n))

    weights /= weights.sum()

    a =  np.convolve(x, weights, mode='full')[:len(x)]
    a[:n] = a[n]
    return a

# --------------
# X.XX -  Get percentiles
# -------------
def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    NOTES:
     - Credit: {{{ http://code.activestate.com/recipes/511478/ (r1)

    NOTES
    ----
    Just use numpy percentile function or pandas .describe() function?

    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1

# --------
# X.XX - Coupler to select plot type
# --------
def create_plot4case( fig, ax, dates, data, spec, f_size=20, lw=None, ls=None, \
        loc=True, legend=True, units='', boxplot=False, lname='Weyborne', \
        run_name='run', start_month=7, end_month=7, plot_type='daily', \
        case=None, l_dict=None, label=None, alt_text=None, color=None):
    """
    Coupler for timeseries data (dates as datetimes + data) to be
            converted multiple graph forms with just change of case
    """

    # --- select frequency
    # select case
    if isinstance(case, type(None) ):
        case ={
        'diurnal': 1, 'daily': 1, 'monthly':2,'weekly': 3,  \
        'May-Sept':4, 'July' : 5, 'Month': 5,
        'PDF':6
        }[plot_type]

    # --- select labels
    # 'Build' latex run name dictionary is not given.
    if isinstance( l_dict, type(None) ):
        l_dict= GC_var('latex_run_names')

    # create label from dictionary if not provided
    if isinstance( label, type(None) ):
        label= l_dict[ run_name ]

    # --- select line plot details.
    # set linewidth +linestyle default if not given
    if isinstance( lw, type(None) ):
        lw=1
    if isinstance( ls, type(None) ):
        ls='-'

    # --- plot requested timeseries
    if case == 1:
        timeseries_daily_plot( fig, ax, dates, data, f_size=f_size,
            title=latex_spec_name(spec), color=color )

    if case == 2:
        timeseries_seasonal_plot( ax, dates, data,  \
            label=label, ylabel=units, color=color, \
            title=latex_spec_name(spec), \
            ls=ls, lw=lw, loc=loc, legend=legend, f_size=f_size)

    if case == 4:
        timeseries_by_day_plot( ax, dates, data,  \
            label=label, ylabel=units, color=color, \
            title=latex_spec_name(spec), boxplot=boxplot, \
            ls=ls, lw=lw, loc=loc, legend=legend, f_size=f_size)

    if case == 5:
        timeseries_month_plot( ax, dates, data,  \
            label=label, ylabel=units, color=color, alt_text=alt_text, \
            start_month=start_month, end_month=end_month, \
            title=latex_spec_name(spec)+' @ '+lname, boxplot=boxplot, \
            ls=ls, lw=lw, loc=loc, legend=legend, f_size=f_size)

    if case == 6:
        PDF_obs_vs_mod( ax, dates, data, \
            label=label, ylabel=units, color=color, \
            start_month=start_month, end_month=end_month, \
            title=latex_spec_name(spec)+' @ '+lname, alt_text=alt_text,  \
            ls=ls, lw=lw, loc=loc, legend=legend, f_size=f_size)


# --------
# X.XX - Plot up locations (lons and lats) on a map
# --------
def plot_lons_lats_spatial_on_map(lons=None, lats=None, p_size=50, color='red',
        title=None, f_size=15, dpi=320, fig=None, ax=None ):
    """
    Plot a list of lons and lats spatially on a map

    Parameters
    -------
    p_size (int): size of plot location point (lon, lat)
    lons, lats (list): list of locations (in decimal londitude and latitude )
    color (str): color of points on map for locations
    title (str): title for plot
    f_size (float): fontsize
    dpi (int): resolution of figure (dots per sq inch)
    """
    import matplotlib.pyplot as plt
    # --- Setup plot
    if isinstance(fig, type(None)):
        fig = plt.figure(dpi=dpi, facecolor='w', edgecolor='k')
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111)
    # plot up white background  (on a blank basemap plot)
    marker = 'o'
    arr = np.zeros((72, 46))
    plt, m = map_plot(arr.T, return_m=True, cmap=plt.cm.binary,
        f_size=f_size, \
        fixcb=[ 0, 0 ], ax=ax, no_cb=True, resolution='c', \
        ylabel=True, xlabel=True, title=title, axis_titles=True )#
    # Plot up all sites as a scatter plot of points on basmap
    m.scatter( lons, lats, edgecolors=color, c=color, marker=marker, \
        s=p_size, alpha=1)


# --------
# X.XX - Probability distribution plotter
# --------
def PDF_obs_vs_mod( ax, dates, data, f_size=20, pos=0, posn=1,  \
        title=None, legend=False, everyother=7*24,  x_nticks=12, \
        window=False, label=None, ylabel=None, loc='upper left',  \
        lw=1,ls='-', color=None, start_month=1, end_month=12,
        unitrotation=45, alpha=.5, bins=100, alt_text=None, debug=False ):
    """
    Constructs a PDF plot with given dates and data.

    NOTES:
     - If not provided with axis etc, then these are made.
     - data and dates need to be as a np.array, with dates as datetime.datetim objects
    """

    # Process data - reduce resolution to daily, and get std
    df = DataFrame(data={'data':data},index=dates )

    # remove dates outside of range (start_month > < end_month )
    def get_month(x):
         return x.month
    df[ 'month'] = df.index.map( get_month )
    df = df[ df['month']<=end_month ]
    df = df[ df['month']>=start_month ]
    df = DataFrame(data={'data':df['data']}, index=df.index )

    # plot up PDF
    plt.hist( df.values, label=label+' (n={})'.format( len(data) ),
        histtype='stepfilled', color=color, alpha=alpha, bins=bins )

    # Beatify plot
    if not isinstance( title, type(None) ):
        plt.title( title + ' for {}-{}'.format( num2month(start_month),\
            num2month(end_month))  )
    if not isinstance( alt_text, type(None) ):
        plt.figtext(x=0.05,y=0.85, s=alt_text, fontsize=f_size*.75)

    if not isinstance( ylabel, type(None) ):
        plt.ylabel( ylabel )
    if legend:
        plt.legend()


