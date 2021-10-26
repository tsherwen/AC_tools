#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Generic plotting functions for timeseries/multi-dimensional output.

Use help(<name of function>) to get details on a particular function.

NOTE(S):
 - This module is underdevelopment vestigial/inefficient code is being removed/updated.
 - Where external code is used credit is given.
"""
# - Required modules:
import os
import xarray as xr
import inspect
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import setp
import functools
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# Time
import time
import calendar
import datetime as datetime
from datetime import datetime as datetime_
# I/O/Admin...
import gc
# The below imports need to be updated,
# imports should be specific and in individual functions
# import tms modules with shared functions
from . variables import *
from . generic import *
from . AC_time import *
from . planeflight import *
from . GEOSChem_nc import *
from . GEOSChem_bpch import *
# math
from math import log10, floor
import numpy as np
import scipy
# colormaps - Additional maps from Eric Sofen
#from option_c import test_cm as cmc
#from option_d import test_cm as cmd


def quick_map_plot(ds, var2plot=None, extra_str='', projection=ccrs.Robinson,
                   save_plot=True, show_plot=False, savename=None, title=None,
                   LatVar='lat', LonVar='lon', fig=None, ax=None, dpi=320,
                   set_global=True, buffer_degrees=0,
                   verbose=False, debug=False, **kwargs):
    """
    Plot up a quick spatial plot of data using cartopy

    Parameters
    -------
    ds (xr.Dataset): dataset object holding data to plot
    var2plot (str): variable to plot within the dataset
    LatVar, LonVar (str): variables to use for latitude and longitude
    save_plot (bool): save the plot as a .png ?
    show_plot (bool): show the plot on screen
    dpi (int): resolution to use for saved image (dots per square inch)
    savename (str): name to use for png of saved .png
    extra_str (str): extra string to append to save .png
    projection (cartopy.crs obj.):  projection to use
    fig (figure instance): matplotlib figure instance
    ax (axis instance): axis object to use

    Returns
    -------
    (None)

    Notes
    -------
     - pass customisation for matplotlib via **kwargs (to 'plot.imshow')
    """
    # Use the 1st data variable if not variable given
    if isinstance(var2plot, type(None)):
        print("WARNING: No 'var2plot' set, trying 1st data_var")
        var2plot = list(ds.data_vars)[0]
    # Setup figure and axis and plot
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(10, 6))
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111, projection=projection(), aspect='auto')
    # print out the min and max of plotted values
    if verbose:
        Pstr = "In spatial plot of {}, min={} and max={}"
        min_ = float(ds[var2plot].values.min())
        max_ = float(ds[var2plot].values.max())
        print(Pstr.format(var2plot, min_, max_))
    # Call plot via imshow...
    im = ds[var2plot].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                                  transform=ccrs.PlateCarree(),
                                  **kwargs)
    # Beautify the figure/plot
    ax.coastlines()
    if set_global:
        ax.set_global()
    else:
        x0, x1, y0, y1 = ax.get_extent()
        if buffer_degrees != 0:
            lons = ds[LonVar].values
            lats = ds[LatVar].values
            x0 = myround(lons.min()-buffer_degrees, buffer_degrees, )
            x1 = myround(lons.max()+buffer_degrees, buffer_degrees,
                         round_up=True)
            y0 = myround(lats.min()-buffer_degrees, buffer_degrees, )
            y1 = myround(lats.max()+buffer_degrees, buffer_degrees,
                         round_up=True)
        ax.set_extent((x0, x1, y0, y1), projection())
    # Add a generic title if one is not provided
    if isinstance(title, type(None)):
        plt.title('Spatial plot of {}'.format(var2plot))
    else:
        plt.title(title)
    # save the plot?
    if save_plot:
        if isinstance(savename, type(None)):
            savename = 'spatial_plot_{}_{}'.format(var2plot, extra_str)
        savename = rm_spaces_and_chars_from_str(savename)
        plt.savefig(savename+'.png', dpi=dpi)
    if show_plot:
        plt.show()
    return im


def ds2zonal_plot(ds=None, var2plot=None, StateMet=None, AltVar='lev',
                  LatVar='lat', fig=None, ax=None, limit_yaxis2=20,
                  plt_ylabel=True, rm_strat=False,
                  verbose=False, debug=False, **kwargs):
    """
    Make a zonal plot (lat vs. alt) from a xr.dataset object
    """
    # Setup figure and axis if not provided
    if isinstance(fig, type(None)):
        fig = plt.figure()
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(1, 1, 1)
    # Calculate number of molecules
    MolecVar = 'Met_MOLCES'
    StateMet = add_molec_den2ds(StateMet, MolecVar=MolecVar)
    # Remove troposphere
    if rm_strat:
        ds2plot = rm_fractional_troposphere(ds[[var2plot]].copy(),
                                            vars2use=[var2plot],
                                            StateMet=StateMet)
    else:
        ds2plot = ds[[var2plot]].copy()
    # Weight by molecules over lon
    ds2plot = ds2plot * StateMet[MolecVar]
    ds2plot = ds2plot.sum(dim=['lon']) / StateMet[MolecVar].sum(dim=['lon'])
    # Update units of lev (vertical/z axis) to be in km
#    StateMet['Met_PMID']
    LatLonAlt_dict = gchemgrid(rtn_dict=True)
    alt_array = LatLonAlt_dict['c_km_geos5']
    ds2plot = ds2plot.assign_coords({'lev': alt_array[:len(ds.lev.values)]})
    # print out the min and max of plotted values
    if verbose:
        Pstr = "In zonal plot of {}, min={}, max={}"
        min_ = float(ds2plot[var2plot].values.min())
        max_ = float(ds2plot[var2plot].values.max())
        print(Pstr.format(var2plot, min_, max_))
    # Now call plot via xr.dataset
    lat = np.array(ds2plot.lat.values)
    alt = np.array(ds2plot.lev.values)
    im = ax.pcolor(lat, alt, ds2plot[var2plot].values, **kwargs)
#    im = ds2plot[var2plot].plot.imshow(ax=ax, **kwargs)
    # Limit the axis to a value (e.g. 18km to just show tropospheric values)
    if limit_yaxis2:
        plt.ylim(0, limit_yaxis2)

    # TODO
    # plot up second y axis with pressure altitude ?

    # Update axis labels
    ax.set_xlabel('Latitude ($^{\circ}$N)')
    if debug:
        print('plt_ylabel', plt_ylabel)
    if plt_ylabel:
        ax.set_ylabel('Altitude (km)')
    else:
        ax.set_ylabel('')
        ax.tick_params(axis='y', which='both', labelleft='off')
        yticks_labels = ax.get_yticklabels()
        ax.set_yticklabels([None for i in range(len(yticks_labels))])
    return im


def plt_df_X_vs_Y(df=None, x_var='', y_var='', x_label=None, y_label=None,
                  ax=None, fig=None, color=None, plt_all_values=True,
                  alpha=0.5,
                  plot_ODR=True, plot_121=True, save_plot=False, dpi=320,
                  verbose=False):
    """
    Make an generic X vs. Y plot from a dataframe

    Parameters
    -------
    df (pd.DataFrame): a dataframe containing the x_var and y_var data
    x_var, y_var (str): names of the variables to plot for X and Y
    x_label, y_label (str): labels of the variables being plotted
    plot_ODR (bool): plot an orthogonal distance regression line of best fit
    plot_121 (bool): Add a 1:1 line to the X vs. Y plot
    save_plot (bool): save the plot as a .png ?
    dpi (int): resolution to use for saved image (dots per square inch)
    plt_all_values (bool): plot all the values as a scatter plot
    color (str): colour for plotted scatter and regression line
    alpha (float): transparency for plotted (scatter) points

    Returns
    -------
    (None)
    """
    # Setup an figure and axis if not provided
    if isinstance(fig, type(None)):
        fig = plt.figure(dpi=dpi, facecolor='w', edgecolor='k')
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111)
    # if no specific variables are given, then use the df variabel names
    if isinstance(x_label, type(None)):
        x_label = x_var
    if isinstance(y_label, type(None)):
        y_label = y_var
    # Get values to plot
    X = df[x_var].values
    Y = df[y_var].values
    # Get the number of samples (N)
    N = float(df.shape[0])
    # get RMSE
    RMSE = np.sqrt(((Y-X)**2).mean())
    # Plot up all the data underneath as a scatter plot
    if plt_all_values:
        ax.scatter(X, Y, color=color, s=3, facecolor='none', alpha=alpha)
    # Add a 1:1 line
    MinVal = min((min(X), min(Y)))
    MaxVal = max((max(X), max(Y)))
    x_121 = np.arange(MinVal-(MaxVal*0.1), MaxVal+(MaxVal*0.1))
    if plot_121:
        ax.plot(x_121, x_121, alpha=0.5, color='k', ls='--')
    # Add a line for the orthogonal distance regression
    xvalues, Y_ODR = get_linear_ODR(x=X, y=Y, xvalues=x_121,
                                    return_model=False, maxit=10000)
    ODRoutput = get_linear_ODR(x=X, y=Y, xvalues=x_121,
                               return_model=True, maxit=10000)
    if verbose:
        print(x_label, y_label, ODRoutput.beta)
    ax.plot(xvalues, Y_ODR, color=color, label=y_label)
    # Beautify the figure/plot
    ax.set_xlabel(x_label)
    ax.set_xlabel(y_label)
    # Plot N value

    # Save the plotted data as a .png?
    if save_plot:
        png_filename = 'X_vs_Y_{}_vs_{}'.format(x_var, y_var)
        png_filename = rm_spaces_and_chars_from_str(png_filename)
        plt.savefig(png_filename, dpi=dpi)


def plt_df_X_vs_Y_hexbin(x=None, y=None, c=None, xscale='linear',
                         yscale='linear',
                         gridsize=(150, 150),
                         fig=None, ax=None, xlimit=None, ylimit=None, dpi=320,
                         xlabel=None, ylabel=None, clabel=None, cmap=None,
                         vmin=None, vmax=None, save2png=False, show_plot=False,
                         figsize=None, linewidths=0.15):
    """
    Make an generic X vs. Y plot from a dataframe

    Parameters
    -------
    x, y, (np.array): data to plot as X and Y
    xscale, yscale (str): 'log' or 'linear' scale to be used for axis?
    xlabel, ylabel (str): labels for x and y data (e.g. for axis titles)
    xlimit (tuple): limit the extent of the plotted x axis
    ylimit (tuple): limit the extent of the plotted y axis
    vmin, vmax (scalar): set the min and max values for colourbar/plotting of C
    linewidths (float): line widths for the hexbins (decrease if overlap seen)
    clabel (str): Label for the values that the X/Y points are coloured by
    cmap (colormap object): colourmap to use for colouring in c
    save2png (bool): save the resultant plot as a .png
    show_plot (bool): Show the resulting plot of the screen
    dpi (int): resolution of figure (dots per sq inch)
    figsize (tuple): size of the plotted figure in dots per square inch

    Returns
    -------
    (None)
    """
    # Setup an figure and axis if not provided
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=figsize, dpi=dpi,
                         facecolor='w', edgecolor='k')
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111)
    # Set the colourmap if it is not provided
    if not isinstance(cmap, type(None)):
        cmap = get_colormap(x)
    # Plot up the provided data
    mappable = ax.hexbin(x, y, c, gridsize=gridsize, xscale=xscale,
                         yscale=yscale, cmap=cmap,
                         vmin=vmin, vmax=vmax, linewidths=linewidths
                         )
    # Limit X axis to provided values
    if not isinstance(xlimit, type(None)):
        ax.set_xlim(xlimit)
    # Limit Y axis to provided values
    if not isinstance(ylimit, type(None)):
        ax.set_ylim(ylimit)
    # Beautify the figure/plot
    if not isinstance(ylabel, type(None)):
        ax.set_ylabel(ylabel)
    if not isinstance(xlabel, type(None)):
        ax.set_xlabel(xlabel)
    # Add a colourbar?
    if not isinstance(c, type(None)):
        plt.colorbar(mappable, label=clabel)
    # Save?
    if save2png:
        png_filename = 'X_vs_Y_hexbin_{}_vs_{}'.format(xlabel, ylabel)
        if not isinstance(c, type(None)):
            png_filename += '_coloured_by_{}'.format(clabel)
        png_filename = rm_spaces_and_chars_from_str(png_filename)
        plt.savefig(png_filename, dpi=dpi)
    if show_plot:
        plt.show()


def plot_up_diel_by_season(spec='O3', sub_str='UK+EIRE', fig=None,
                           dfs=None, color_dict={'Obs.': 'k', 'Model': 'r'},
                           stat2plot='50%', title=None,
                           dpi=320, plt_legend=True, units=None,
                           show_plot=False, save_plot=False, verbose=False,
                           context='paper', font_scale=0.75,
                           tight_layout=False, use_letters4months=False,
                           debug=False):
    """
    Plot up mulitplot of diel cycle by "season" for given dictionary of DataFrames

    Parameters
    -------
    spec (str): column name in DataFrame for spec to plot up
    dfs (dictionary): dictionary of DataFrames of data+dates (values)
    fig (figure instance): figure to use for to plot subplots onto
    sub_str (str): extra string to use in titles and names of saved files
    stat2plot (str): stat (e.g. mean or medium) to use for main line
    color_dict (dictionary): dictionary of colors for keys in dictionary (dfs)
    plt_legend (bool): Add a single legend to the figure
    units (str): units of the provided data (used in plot labels)
    show_plot (bool): display the plot to screen?
    save_plot (bool): save the plot to disk
    verbose (bool): print verbose statements to screen/error log
    debug (bool): print debugging statements to screen/error log
    dpi (scalar): resolution (in dots per square inch) to use when saving figure

    Returns
    -------
    (None)
    """
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context(context, font_scale=font_scale)
    # Local variables
    seasons = ('DJF', 'MAM', 'JJA',  'SON')
    month2season = np.array([
        None,
        'DJF', 'DJF',
        'MAM', 'MAM', 'MAM',
        'JJA', 'JJA', 'JJA',
        'SON', 'SON', 'SON',
        'DJF'
    ])
    season2text = {
        'DJF': 'Dec-Jan-Feb', 'MAM': 'Mar-Apr-May', 'JJA': 'Jun-Jul-Aug', 'SON': 'Sep-Oct-Nov', None: None,
    }
    if use_letters4months:
        pass
    else:
        seasons = [season2text[i] for i in seasons]
        month2season = np.array([season2text[i] for i in month2season])
    # Split data by season
    for key_ in list(dfs.keys()):
        # Now assign "Seasons"
        dfs[key_]['Season'] = month2season[dfs[key_].index.month]
    # - Loop seasons and Plot data
    # Setup figure and PDF
    if isinstance(fig, type(None)):
        fig = plt.figure(dpi=dpi)
    # setup an inital subplot (to allow axis sharing)
    ax1 = fig.add_subplot(2, 2, 1)
    # Now loop seasons...
    for n_season, season_ in enumerate(seasons):
        if verbose:
            print((n_season, season_, spec, sub_str))
        # Add subplot and set axis labelling
        if (n_season > 0):
            ax = fig.add_subplot(2, 2, n_season+1, sharey=ax1)
        else:
            ax = ax1
        _plt_legend = False
        if (n_season+1) == len(seasons):
            _plt_legend = True
        plt_xlabel = True
        do_not_plt_xlabel_on_subplot = [1, 2]
        if (n_season+1) in do_not_plt_xlabel_on_subplot:
            plt_xlabel = False
        plt_ylabel = True
        do_not_plt_ylabel_on_subplot = [2, 4]
        if (n_season+1) in do_not_plt_ylabel_on_subplot:
            plt_ylabel = False

        # - Loop and plot DataFrames (dfs)
        for n_key, key_ in enumerate(dfs.keys()):
            # Date and dates? ( only select for season )
            tmp_df = dfs[key_]
            tmp_df = tmp_df[tmp_df['Season'] == season_]
            data_ = tmp_df[spec]
            dates_ = pd.to_datetime(tmp_df.index.values)
            # See if color is set in dictionary (color_dict)
            try:
                color = color_dict[key_]
            except KeyError:
                color = None
            # Add a legend to plot?
            legend = False
            if _plt_legend and (n_key == len(list(dfs.keys()))-1):
                legend = True
            if plt_legend == False:
                legend = False
            # Plot up using the basic plotter function
            BASIC_diel_plot(fig=fig, ax=ax, data=data_, units=units,
                            dates=dates_, label=key_, stat2plot=stat2plot,
                            title='{}'.format(season_),
                            plt_xlabel=plt_xlabel, plt_ylabel=plt_ylabel,
                            color=color, spec=spec, plt_legend=legend)
            # Remove tmp data dictionary from memory...
            del tmp_df

    # Show or save if requested (NOTE: default is not to save or show!)
    if isinstance(title, type(None)):
        suptitle = "diel of {} in '{}'"
        # try and use a LaTeX form of the species name
        try:
            specname = latex_spec_name(spec)
        except:
            specname = spec
        fig.suptitle(suptitle.format(specname, sub_str))
    else:
        fig.suptitle(title)
    if tight_layout:
        plt.tight_layout()
    png_filename = 'Seasonal_diel_{}_{}.png'.format(sub_str, spec)
    if save_plot:
        plt.savefig(png_filename, dpi=dpi)
    if show_plot:
        plt.show()


def BASIC_diel_plot(fig=None, ax=None, dates=None, data=None, color='red',
                    title=None, label=None, plt_legend=None, plt_xlabel=True,
                    plt_ylabel=True, show_plt=False,
                    plot_pcent_change_from_max=False,
                    units='ppbv', spec='O3', alt_text=None, loc='best',
                    filename2save='diel_plot.png', save_plt=False,
                    show_plot=False,
                    stat2plot='50%', alpha=0.3, time_resolution_str="%H",
                    add_quartiles2plot=True, return_avgs=True, ls='-', ncol=1,
                    xlabel='Hour of day (UTC)', debug=False, lw=2,
                    force_repeat_of_first_hour_as_last_hour=False,
                    update_xticks=True):
    """
    Creates a diel plot for given data and dates

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
    legend (bool): add a legend?
    title (str): title for plot
    color (str/hex etc): color for line
    f_size (float): fontsize
    lgnd_f_size (float): fontsize for legend
    loc (str): location for legend
    rotatexlabel (numnber/str): rotation of x axis labels
    pos, posn (int): vestigle(!) location indices for window plots
    ylim (list): min and max y axis limit
    lw (str): linewidth
    return_avgs (bool): return a list of the average values plotted
    save_plt, show_plot (bool): show or save plots?
    xlabel (str): label for x axis

    Returns
    -------
    (None)
    """
    prt_str = 'BASIC_diel_plot called for {} (data.shape={})'
    if debug:
        print((prt_str.format(spec, data.shape)))
    import matplotlib
    from matplotlib import dates as d
    import datetime as dt
    #
    if isinstance(fig, type(None)):
        fig = plt.figure()  # dpi=Var_rc['dpi'])
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111)  # 2,2, n_season+1 )
    if debug:
        print((fig, ax))
    # --- Process
    # Form a dataFrame from input numpy arrays. (and remove NaNs... )
    raw_df = pd.DataFrame({'data': data}, index=dates).dropna()
    # Add a time coluumn to group by (e.g. "%H:%M" for minutes)
    raw_df['Time'] = raw_df.index.map(
        lambda x: x.strftime(time_resolution_str))
    # repeat hour zero as hour 24 (to aesthetically loop on plot)
    if force_repeat_of_first_hour_as_last_hour:
        first_hour = raw_df[raw_df['Time'] == '00']
        first_hour.loc[:, 'Time'] = '24'
        raw_df = pd.concat([raw_df, first_hour])
    # Now group
    df = raw_df.groupby('Time').describe().unstack()
    # Get the labels for time
    if debug:
        print((df.head(), df.index[:5], df.shape))
    time_labels = df['data'][stat2plot].index.values
    time_labels = [str(int(i)) for i in time_labels]

    # make sure the values with leading zeros drop these
    index = [float(i) for i in time_labels]
    if debug:
        print((df['data'][stat2plot]))
    if debug:
        print(('!'*20, index, time_labels))
    # Select data for the requested statistic
    avgs = df['data'][stat2plot]
    # Plot the % change from max ?
    if plot_pcent_change_from_max:
        if debug:
            print((avgs, avgs.shape))
        max_ = avgs.max()
        avgs = (avgs - max_) / max_ * 100
        if debug:
            print((avgs.shape, max_))
        units = '%'
    # - Now plot up
    ax.plot(index, avgs, color=color, linewidth=lw, label=label, ls=ls)
    # - Add quartiles
    if add_quartiles2plot:
        Q3_data = df['data']['75%']
        # Plot the % change from max ?
        if plot_pcent_change_from_max:
            max_ = Q3_data.max()
            Q3_data = (Q3_data - max_) / max_ * 100
        ax.fill_between(index, avgs, Q3_data, alpha=alpha, facecolor=color)
        Q1_data = df['data']['25%']
        # Plot the % change from max ?
        if plot_pcent_change_from_max:
            max_ = Q1_data.max()
            Q1_data = (Q1_data - max_) / max_ * 100
        ax.fill_between(index, avgs, Q1_data, alpha=alpha, facecolor=color)
    # Label xticks on axis (if 24 plots, just label every 3rd)
    if update_xticks:
        if len(index) < 6:
            ax.set_xticks(index)
            ax.set_xticklabels(time_labels)
        else:
            ax.set_xticks(index[2::3])
            ax.set_xticklabels(time_labels[2::3])
        xticks = ax.get_xticks()
        if debug:
            print((xticks, ax.get_xticklabels()))
#    ax.set_xticks(np.linspace(3, 21, 7).astype(int))
    # More cosmetic changes...
    if not isinstance(title, type(None)):
        plt.title(title)
    if plt_legend:
        plt.legend(loc=loc, ncol=ncol)
    if plt_xlabel:
        plt.xlabel(xlabel)
    else:
        ax.tick_params(axis='x', which='both', labelbottom='off')
    if plt_ylabel:
        try:
            spec_ = latex_spec_name(spec)
        except:
            spec_ = spec
        plt.ylabel('{} ({})'.format(spec_, units))
    else:
        ax.tick_params(axis='y', which='both', labelleft='off')
    # Add alt text if provided.
    if not isinstance(alt_text, type(None)):
        alt_text_y = 0.925
        alt_text_x = 0.975
        plt.text(alt_text_x, alt_text_y, alt_text, ha='right', va='center',
                 transform=ax.transAxes)
    # Print or show plot?
    if save_plt:
        plt.savefig(filename2save)
    if show_plt:
        plt.show()
    # Return the average values?
    if return_avgs:
        return np.ma.array(avgs)


def BASIC_seasonal_plot(dates=None, data=None, ax=None,
                        title=None, legend=False,  ylabel=None, loc='best',
                        return_avgs=False,
                        plt_median=False, plot_Q1_Q3=False,
                        xtickrotation=45, alt_text=None, alt_text_x=.925,
                        alt_text_y=.925, xlabel=True, rm_yticks=False,
                        log=False,
                        pcent1=25, pcent2=75, color='red', lw=1, ls='-',
                        label=None,
                        ylim=None, debug=False):
    """
    Plot up a basic seasonal plot - adapted from AC_tools
    monthly_plot

    Parameters
    ----------
    dates (nd.array): numpy array of datetime.datetime objects
    data (nd.array): numpy array of data to plot
    loc (str): best location to put legend plot
    legend (bool): Add a legend to the plot
    ylabel (str): label for y axis
    xlabel (bool): label x axis
    return_avgs (bool): return a list of the monthly averages
    plot_Q1_Q3 (bool): plot quartiles on for data?
    title (str): title string for plot
    alt_text (str): subtitle string for plot (insert in plot)
    alt_text_x (float): x axis position for subtitle str
    alt_text_y (float): y axis position for subtitle str
    xtickrotation (float): rotation of x axis ticks
    rm_yticks (bool): remove the y axis ticks
    log (bool): set the y scale to be logarithmic
    ls (str): matplotlibe linestyle to use
    lw (int): linewidth to use
    color (str): colour for line
    label (str): label for line
    pcent1 (int): lower percentile to use (e.g. 25) for shaded region
    pcent2 (int): higher percentile to use (e.g. 75) for shaded region
    ylim (tuple): set limit of the y axis (min, max)
    debug (bool): print debuging statements?

    Returns
    -------
    (None or nd.array)
    """
    if debug:
        print((list(locals().keys()), ax))
    # setup data + dates as a DataFrame
    df = pd.DataFrame(data, index=dates)
    # Force use of a standard year and make sure dates are in order
    df['months'] = [i.month for i in dt64_2_dt(df.index.values)]
    df.sort_values(by='months', ascending=True, inplace=True)
#    months = list(range(1, 13))
    months = df['months'].values
    df = df.drop('months', axis=1)
    datetime_months = [datetime.datetime(2009, int(i), 1) for i in set(months)]
    labels = [i.strftime("%b") for i in datetime_months]
    # Get Data by month
    monthly = df.groupby(df.index.month)
    months = list(sorted(monthly.groups.keys()))
    monthly = [monthly.get_group(i) for i in months]
#    [ df[df.index.month==i]  for i in months ]
    # Remove nans (just from columns!) to allow for percentile calc.
    monthly = [i.dropna(axis=1).values for i in monthly]
    data_nan = [i.flatten() for i in monthly]
    # Plot up median
    medians = [np.nanpercentile(i, 50, axis=0) for i in data_nan]
    plt.plot(months, medians, color=color, lw=lw, ls=ls, label=label)
    # Define an axis variable if one isn't given
    if isinstance(ax, type(None)):
        ax = plt.gca()
    # Plot quartiles as shaded area?
    if plot_Q1_Q3:
        low = [np.nanpercentile(i, pcent1, axis=0) for i in data_nan]
        high = [np.nanpercentile(i, pcent2, axis=0) for i in data_nan]
        ax.fill_between(months, low, high, alpha=0.2, color=color)
    # Beatify plot
    ax.set_xticks(months)
    if xlabel:
        ax.set_xticklabels(labels, rotation=xtickrotation)
    else:
        ax.tick_params(axis='x', which='both', labelbottom='off')
    if debug:
        print(('!'*50, alt_text, alt_text_x, alt_text_y))
    if not isinstance(alt_text, type(None)):
        if debug:
            print(('!'*50, alt_text, alt_text_x, alt_text_y))
        plt.text(alt_text_x, alt_text_y, alt_text, ha='right', va='center',
                 transform=ax.transAxes)
    if legend:
        if debug:
            print(('>'*500, 'Adding legend', '<'*50, loc))
        plt.legend(loc=loc)  # fontsize=f_size*.75, loc=loc )
    if not isinstance(title, type(None)):
        plt.title(title)
    if ylabel:
        plt.ylabel(ylabel)
    else:
        if rm_yticks:
            ax.tick_params(axis='y', which='both', labelleft='off')
    if not isinstance(ylim, type(None)):
        plt.ylim(ylim)
    # Logarithmic scale?
    if log:
        ax.set_yscale('log')
    else:
        ax.set_yscale('linear')
    if return_avgs:
        return np.ma.array(medians)


def binned_boxplots_by_altitude(df=None, fig=None, ax=None, dataset_name=None,
                                num_of_datasets=1, dataset_num=0, label='Obs.',
                                showfliers=False,
                                bins=np.arange(8), binned_var='O3',
                                var2bin_by='ALT',
                                color=None, dpi=320, xlabel=None,
                                ylabel='Altitude (km)',
                                title=None,
                                show_plot=False, widths=0.3, verbose=False,
                                debug=False):
    """
    Plot up dataset (1 or more) as boxplots binned by altitude

    Parameters
    -------
    df (pd.DataFrame): dataframe of alts and variable
    var2bin_by (str): variable to bin by (e.g. altitude)
    num_of_datasets (int): number of dataset for which boxplots will be plotted
    dataset_num (int): number of dataset in order to be order to be plotted
    fig (figure instance): fig. to use
    ax (axis instance): axis to use
    label (str): label for legend
    binned_var (str): variable (in df) to bin by var2bin_by
    bins (np.array): bins to split dataset into
    color (str): color for boxplot
    title (str): title for plot?
    widths (float): width of boxplots
    ylabel, xlabel (bool): include axis labels for plot?

    Returns
    -------
    (None)
    """
    # ----  Local variarbles
    # Setup figure and axis if not provided
    if isinstance(fig, type(None)):
        fig = plt.figure(dpi=dpi)
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111)
    # color provided?
    if isinstance(color, type(None)):
        if 'obs' in dataset_name.lower():
            color = 'k'
        elif 'model' in dataset_name.lower():
            color = 'red'
        else:
            color = 'blue'
    # ---- Process data
    # Make sure there are no NaNs in the arrays that have been passed
    NAN_present = df[[var2bin_by, binned_var]].isnull().values.any()
    assert (not NAN_present), 'bin or binned vars cannots contain NaNs!'
    # Now group by bins and get data as list
    gdf = df.groupby(pd.cut(df[var2bin_by].values, bins))
    data = [i[1][binned_var].values for i in gdf]
    # ----- Now plots
    # - settings for boxplots
    bin_sizes = np.diff(bins)
    mid_points = [(i/2.)+bins[n] for n, i in enumerate(bin_sizes)]
    upper_bin_edge = bins[1:]
    set_of_bin_sizes = list(set(bin_sizes))
    if len(set_of_bin_sizes) == 1:
        if num_of_datasets > 1:
            # what size are the boxplots  (based on number of datasets)
            space4boxplot = set_of_bin_sizes[0]*0.80/num_of_datasets
            # where should the boxplots go?
            positions = mid_points - ((space4boxplot*(num_of_datasets))*0.5)
            positions = positions + \
                (space4boxplot*dataset_num)+(space4boxplot*.5)
            # if there are more than 2 dataset, divide the space up accordingly
            if num_of_datasets > 2:
                widths = widths / num_of_datasets * 1.75
        else:
            positions = mid_points
    else:
        print('Not yet setup to take different bin sizes')
        sys.exit()
    # - Add plots
    bp = ax.boxplot(x=data, positions=positions, patch_artist=True,
                    widths=widths, notch=False, vert=False,
                    showfliers=showfliers)
    # add an (invisible) plot to carry the label
    ax.plot([], [], label=label, color=color)
    # - Beautify
    # set the boxplot format?
    set_bp_style(bp, color=color)  # , linewidth=lw )
    # title and y axis label?
    if not isinstance(ylabel, type(None)):
        ax.set_ylabel(ylabel)
    if not isinstance(xlabel, type(None)):
        ax.set_xlabel(xlabel)
    if not isinstance(title, type(None)):
        plt.title(title)
    # Set y axis range
    if debug:
        print((bins, bins[0], bins[-1]))
    ax.set_ylim(bins[0], bins[-1]+(bin_sizes[-1]*0.5))
    # Set y axis tick labels to be the bins
    ax.set_yticks(upper_bin_edge)
    if all([float(i).is_integer() for i in bins]):
        ax.set_yticklabels([int(i) for i in upper_bin_edge])
    else:
        ax.set_yticklabels([float(i) for i in upper_bin_edge])
    if debug:
        print((ax.get_yticks()))
    if debug:
        print((ax.get_yticks(), bins))
    # - Save or show?
    if show_plot:
        plt.show()


def scatter_3D_cube(data, dims=None, res='2x2.5', fig=None, everyother=1,
                    interval=1,
                    f_size=20, cm='RdYlBu', debug=False):
    """
    Make a 3D scatter cube of given 3D array
    """
    logging.debug('scatter_3D_cube called with input data shape: {}'.format(
        data.shape))

    from mpl_toolkits.mplot3d import Axes3D

    if isinstance(dims, type(None)):
        lon, lat, alt = get_latlonalt4res(res=res)
        (X, Y, Z) = np.mgrid[
            lon[0]:lon[-1]:complex(len(lon)),
            lat[0]:lat[-1]:complex(len(lat)),
            alt[0]:alt[-1]:complex(len(alt))]

    # Setup Fig + Ax ( is not given )
    cmap = plt.cm.get_cmap(cm)
    if isinstance(fig, type(None)):
        fig = plt.figure(1)
    fig.clf()
    ax = Axes3D(fig)
    ax.scatter(X, Y, Z, c=data, cmap=cmap)

    # --- Beatify
    # draw meridian lines
    meridians = np.arange(-180, 180, 20*interval)
    plt.xticks(meridians[::everyother], fontsize=f_size*.75)

    # draw parrelel lines
    parallels = np.arange(-90, 91, 20*interval)
    plt.yticks(parallels[::everyother], fontsize=f_size*.75)

    # limit axis
    ax.set_xlim(lon[0], lon[-1])
    ax.set_ylim(lat[0], lat[-1])
    ax.set_zlim(alt[0], alt[-1])


def plot_zonal_figure(arr, fixcb=None, cb_sigfig=2, ax=None,
                      norm=None, nticks=10, format=None, units=None,
                      extend='neither',
                      discrete_cmap=False, f_size=15, fig=None, res='4x5',
                      wd=None,
                      t_ps=None,
                      trop_limit=True, axn=None, cb_ax=None,
                      orientation='vertical',
                      rotatecbunits='vertical', width=0.015, height=0.6,
                      bottom=0.1, top=0.925, hspace=0.4, wspace=0.5,
                      left=0.075,
                      right=0.875,
                      cb_bottom=0.125, cb_height=0.825, cb_left=0.885, dpi=160,
                      no_cb=True,
                      region='All', lat_0=None, lat_1=None, pdftitle=None,
                      return_m=False,
                      rtn_plt_vars=False, set_window=False, pdf=False,
                      show=True,
                      log=False,
                      window=False, xlabel=True, ylabel=True, title=None,
                      interval=None, verbose=False, debug=False):
    """
    Create a zonal plot as formated figure


    NOTES:
     -
    """
    all_str = 'plot_zonal_figure called ', region, arr.shape, log, units, pdf, \
        show, arr.min(), arr.max()
    logging.info(all_str)
    if verbose:
        print(all_str)

    # If lon, lat, alt array provided then take mean of lon
    if any([arr.shape[0] == i for i in (72, 144, 121, 177)]):
        #        arr = arr.mean(axis=0)
        arr = molec_weighted_avg_BPCH(arr, weight_lon=True, res=res,
                                      trop_limit=trop_limit, rm_strat=False, wd=wd)

    # Create figure if not provided
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(15, 10), dpi=dpi,
                         facecolor='w', edgecolor='w')

    # if just plotting over the ocean, remove white space
    if region == 'Oceanic':
        set_window = True
        lat_0, lat_1 = -65, 80
        if isinstance(interval, type(None)):
            interval = 3

    # Set colourbar limits
    if isinstance(fixcb, type(None)):
        fixcb = np.array([(i.min(), i.max()) for i in [arr]][0])

        # Set readable levels for cb, then use these to dictate cmap
        if not log:
            lvls = get_human_readable_gradations(vmax=fixcb[1],
                                                 vmin=fixcb[0], nticks=nticks,
                                                 cb_sigfig=cb_sigfig,
                                                 verbose=verbose, debug=debug)
    else:
        # Assume numbers provided as fixcb +nticks will allow for creation
        # of human readable levels.
        lvls = np.linspace(fixcb[0], fixcb[1], nticks)

    # If log plot - overwrite  lvls
    if log:
        # Get logarithmically spaced integers
        lvls = np.logspace(np.log10(fixcb[0]), np.log10(fixcb[1]),
                           num=nticks)
        # Normalise to Log space
#            norm=mpl.colors.LogNorm(vmin=fixcb_[0], vmax=fixcb_[1])
        if isinstance(norm, type(None)):
            norm = mpl.colors.LogNorm(vmin=fixcb[0], vmax=fixcb[1])

    # Setup Colormap
    cmap, fixcb_buffered = get_colormap(np.array(fixcb),
                                        nticks=nticks, fixcb=fixcb,
                                        buffer_cmap_upper=True,
                                        verbose=verbose, debug=debug)

    #  Plot
    if isinstance(axn, type(None)):
        axn = [111]
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(*axn)
    zonal_plot(arr, fig,  ax=ax, set_window=set_window, log=log,
               format=format, cmap=cmap, lat_0=lat_0, lat_1=lat_1,
               fixcb=fixcb, f_size=f_size*.75, res=res, norm=norm,
               fixcb_buffered=fixcb_buffered, no_cb=True, trop_limit=True,
               window=window, interval=interval, xlabel=xlabel, ylabel=ylabel,
               verbose=verbose, debug=debug)

    # Only show troposphere
    if isinstance(t_ps, type(None)):
        t_ps = get_GC_output(wd, vars=['TIME_TPS__TIMETROP'], trop_limit=True)
    greyoutstrat(fig, t_ps.mean(axis=0).mean(axis=-1), axn=axn, res=res)

    if not isinstance(title, type(None)):
        plt.title(title, fontsize=f_size*.75)
#        plt.text(0.5, y_title, title, fontsize=f_size*.75 )

    # Manually Add colorbar
    if no_cb:
        if orientation == 'vertical':
            width = width/2
            height = 0.55

        mk_cb(fig, units=units, left=cb_left,  height=cb_height,
              bottom=cb_bottom, log=log, orientation=orientation,
              rotatecbunits=rotatecbunits, width=width,
              cmap=cmap, vmin=fixcb_buffered[0],
              vmax=fixcb_buffered[1], format=format, f_size=f_size*.75,
              extend=extend, lvls=lvls, cb_ax=cb_ax,
              cb_sigfig=cb_sigfig, nticks=nticks,
              norm=norm, discrete_cmap=discrete_cmap,
              verbose=verbose, debug=debug)

    # Adjust plot ascetics
    fig.subplots_adjust(bottom=bottom, top=top, left=left, right=right,
                        hspace=hspace, wspace=wspace)
    # Show/Save as PDF?
    if pdf:
        plot2pdf(title=pdftitle)
    if show:
        plt.show()
    if return_m or rtn_plt_vars:
        return [fig, cmap, fixcb]  # + plt_vars + [ fixcb ] #+= [ cb_ax ]


def plot_arr_avg_Q1_Q3(X, Y, ax=None, color='blue', label=None,
                       plt_mean=True, plt_median=False, pcent1=25, pcent2=75,
                       verbose=False, debug=False):
    """
    Takes a X and Y to plot a mean, Q1 and Q3.

    Notes:
     - Y = 2D array (e.g. lon, lat)
     - X = 1D araray ( e.g. lat )
     - Default fill is Q1 to Q3, but others ranges can be specified
    """
    if isinstance(ax, type(None)):
        ax = plt.gca()

    # Plot up mean of Y values
    if plt_mean:
        ax = plt.plot(X, Y.mean(axis=0), color=color, label=label)

    # Show quartiles ( 25th 57th  )
    Y_nan = Y.filled(np.nan)
    low = np.nanpercentile(Y_nan, pcent1, axis=0)
    high = np.nanpercentile(Y_nan, pcent2, axis=0)
    ax.fill_between(X, low, high, alpha=0.2, color=color)
    if plt_median:
        ax = plt.plot(X, np.nanpercentile(Y_nan, 50, axis=0),
                      color=color, label=label)

    # return axis object
    return ax


def X_stackplot(X=None, Y=None, labels=None, baseline='zero',
                fig=None, ax=None, dpi=160, show=False, f_size=10,
                legend=False,
                colors=None, title=None, loc='upper right', ylim=None,
                xlim=None,
                lw=8.0, ylabel=None, xlabel=False, log=False, rm_ticks=False,
                alt_text_x=.15, alt_text_y=0.75, alt_text=None, ncol=1,
                pcent=False,
                stacked=False, verbose=False, debug=False):
    """
    Make a stacked plot (by X axis) for values in Y array.

    Parameters
    -------
    X (list): list of numpy arrays to plot as X
    Y (array): must be a numpy array use as Y
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
    logging.info('X_stackplot called, X[0] & Y[0] shape={}.{}'.format(
        *[i.shape for i in (X[0], Y)]))
    # - if "fig" and "ax" not provided then create them.
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(8, 8), dpi=dpi, facecolor='w',
                         edgecolor='w')
        logging.info('Creating figure')
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(1, 1, 1)
        logging.info('Creating ax')
    logging.debug('{}'.format(list(zip(labels, [np.sum(i) for i in X]))))
    # - Stack arrays if not stacked...
    if isinstance(X, np.ndarray):
        X = np.ma.atleast_2d(X)
        logging.debug('Array passed has been made at least 2D if nessesary')
    else:
        X = np.ma.column_stack(X)
        logging.debug('WARNING - Stacking X data, X shape={}'.format(X.shape))
    # Assume data passed has not been 'stacked', so stack it here.
#    if stacked:
#        stack =  X
#    else:
    stack = np.ma.cumsum(X, axis=1)
    # convert to %?
#    pcent = True
    logging.info('stack shape={}'.format(stack.shape))
    if pcent:
        max = np.ma.max(stack, axis=1)  # get accumulated maximum
        print([(i.min(), i.max(), i.mean(), i.shape) for i in (stack, max)])
        stack = np.ma.divide(stack,  max[:, None]) * 100
        print([(i.min(), i.max(), i.mean(), i.shape) for i in (stack, max)])
        xlim = [0, 100]
    debug_ptr = (labels, [np.sum(stack[:, n]) for n, i in enumerate(labels)])
    if debug:
        print(list(zip(debug_ptr)))
    # - Setup baseline ( can expand to include other options... )
    if baseline == 'zero':
        first_line = np.zeros(stack[:, 0].shape)
    # - Plot by label
    # Get list of colors
    if isinstance(colors, type(None)):
        colors = color_list(len(stack[0, :]))
    logging.debug('{}'.format(list(zip(labels, [colors[:len(labels)]]))))
    # Color between x = 0 and the first array.
    logging.debug('{}'.format(stack[:, 0]))
    logging.debug('len colors={}, labels={}, stack={}'.format(
        *[len(i) for i in (colors, labels, stack[0, :])]))
    r = []
    r += [ax.fill_betweenx(Y, first_line, stack[:, 0],
                           color=colors[0],
                           label=labels[0])]
    # Color between array i-1 and array i
    r += [ax.fill_betweenx(Y, stack[:, i], stack[:, i+1],
                           color=colors[i+1],
                           label=labels[i+1])
          for i in range(0, len(stack[0, :])-1)]
    # Plot transparent lines to get 2D line object to create legend
    # Not needed. just use the labels from the fill_between calls.
#     [ plt.plot( Y, stack[:,n], alpha=0, color=colors[n], label=i) \
#         for n,i in enumerate(labels) ]
    # Log scale?
    if log:
        ax.set_xscale('log')
    # Print maxima
    debug_ptr = [title]
    debug_ptr += [[(i.min(), i.max(), i.mean()) for i in [stack[:, n]]]
                  for n, label in enumerate(labels)]
    if debug:
        print(debug_ptr)
    # - Beautify plot
    if not isinstance(ylim, type(None)):
        plt.ylim(ylim)
    if not isinstance(xlim, type(None)):
        plt.xlim(xlim)
    if not isinstance(title, type(None)):
        plt.title(title, fontsize=f_size)
    # Add alt text?
    if not isinstance(alt_text, type(None)):
        ax.annotate(alt_text, xy=(alt_text_x, alt_text_y),
                    textcoords='axes fraction', fontsize=f_size*1.5)
    if legend:
        # Add legend
        if ncol == 0:
            leg = plt.legend(loc=loc, fontsize=f_size*.75,)
        else:
            import itertools

            def flip(items, ncol):
                return itertools.chain(*[items[i::ncol] for i in range(ncol)])
            handles, labels = ax.get_legend_handles_labels()
            leg = plt.legend(flip(handles, ncol), flip(labels, ncol), loc=loc,
                             ncol=ncol, fontsize=f_size*0.75)
        # ( + update line sizes)
        for legobj in leg.legendHandles:
            legobj.set_linewidth(lw)
            legobj.set_alpha(1)
    # Remove tick labels on y axis?
    if isinstance(ylabel, type(None)):
        ax.tick_params(axis='y', which='both', labelleft='off',
                       labelsize=f_size*.75)
    else:
        plt.ylabel(ylabel, fontsize=f_size*.75)
        ax.tick_params(labelsize=f_size*.75)
    # Remove tick labels on x axis?
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=f_size*.75)
        ax.tick_params(axis='x', which='both', labelsize=f_size*.75)
    else:
        if rm_ticks:
            ax.tick_params(axis='x', which='both', labelbottom='off',
                           labelsize=f_size*.75)
        else:
            ax.tick_params(axis='x', which='both', labelsize=f_size*.75)
    if show:
        plt.show()


def plot_lons_lats_spatial_on_map_CARTOPY(central_longitude=0,
                                          lats=None, lons=None,
                                          add_background_image=True,
                                          projection=ccrs.PlateCarree,
                                          fig=None, ax=None,
                                          marker='o', s=50, color='red',
                                          show_plot=False,
                                          dpi=320, label=None, alpha=1,
                                          add_gridlines=True,
                                          buffer_degrees=10,
                                          add_detailed_map=False,
                                          map_minor_islands=False,):
    """
    Plot a list of lons and lats spatially on a map (using cartopy)

    Parameters
    -------
    projection (cartopy.crs obj.):  projection to use
    s (int): size of plot location point (lon, lat)
    lons, lats (np.array): list of locations (in decimal londitude and latitude)
    color (str): color of points on map for locations
    dpi (int): resolution of figure (dots per sq inch)
    return_axis (boaol): return the  axis instance
    marker (str): marker style

    Returns
    -------
    (axis instance)

    Notes
    -----
    """
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    # Setup plot
    if isinstance(fig, type(None)):
        fig = plt.figure(dpi=dpi, facecolor='w', edgecolor='k')
    if isinstance(ax, type(None)):
        # setup a cartopy projection for plotting
        ax = fig.add_subplot(111,
                             projection=projection(
                                 central_longitude=central_longitude)
                             )
    # Plot up
    # Add buffer region around plot
#    ax.get_extent()
#     plt.ylim( myround(lats.min()-buffer_degrees, 10, ),
#         myround(lats.max()+buffer_degrees, 10, round_up=True))
#     plt.xlim( myround(lons.min()-buffer_degrees, 10, ),
#         myround(lons.max()+buffer_degrees, 10, round_up=True))
    x0, x1, y0, y1 = ax.get_extent()
    try:
        x0 = myround(lons.min()-buffer_degrees, buffer_degrees, )
        x1 = myround(lons.max()+buffer_degrees, buffer_degrees, round_up=True)
        y0 = myround(lats.min()-buffer_degrees, buffer_degrees, )
        y1 = myround(lats.max()+buffer_degrees, buffer_degrees, round_up=True)
        ax.set_extent((x0, x1, y0, y1), projection())
    except ValueError:
        print('lon and lat buffer not set extent as out of range')

    # Put a background image on for nice sea rendering.
    if add_background_image:
        ax.stock_img()
    if add_detailed_map:
        # Also add minor islands (inc. Cape Verde)?
        if map_minor_islands:
            land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                    edgecolor=None,
                                                    facecolor='none')
            ax.add_feature(land_10m, edgecolor='grey', facecolor='none',
                           zorder=50)
        else:
            ax.add_feature(cfeature.LAN)
            ax.add_feature(cfeature.COASTLINE)
        # Create a feature for States/Admin 1 regions at 1:50m
        # from Natural Earth
        states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            edgecolor='grey',
            facecolor='none')
        ax.add_feature(cfeature.BORDERS, edgecolor='grey')
#        ax.add_feature(states_provinces, edgecolor='gray')

    else:
        ax.coastlines(resolution='110m')
#        ax.drawcountries() # not in cartopy
    # Include gridlines
    if add_gridlines:
        ax.gridlines()

    # Now scatter points on plot
    ax.scatter(lons, lats, color=color, s=s, marker=marker, alpha=alpha,
               transform=projection(), zorder=999, label=label)

    # return  ax (and show plot?)
    if show_plot:
        plt.show()
    return ax


def color_list(length, cb='gist_rainbow'):
    """
    Create a list of colours to generate colors for plots contain multiple datasets
    """
    cm = plt.get_cmap(cb)
    color = [cm(1.*i/length) for i in range(length)]
    return color


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


def set_bp(bp, num, c_list=['k', 'red'], white_fill=True, set_all=True,
           median_color='white', linewidth=2, debug=False):
    """
    Manually set properties of boxplot ("bp")
    """
    if debug:
        print((num, c_list))
    if set_all:
        setp(bp['boxes'][:], color=c_list[num], linewidth=linewidth*.5)
        setp(bp['caps'][:], color=c_list[num], linewidth=linewidth*.5)
        setp(bp['whiskers'][:], color=c_list[num], linewidth=linewidth*.5)
        setp(bp['fliers'][:], color=c_list[num], linewidth=linewidth*.5)
        setp(bp['medians'][:], color=c_list[num], linewidth=linewidth*.5)
    if white_fill:
        [box.set(facecolor='white') for box in bp['boxes']]
    else:
        [box.set(facecolor=c_list[num]) for box in bp['boxes']]
        setp(bp['medians'][:], color=median_color, linewidth=linewidth)


def markers_list(rm_plain_markers=False):
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
        [markers.pop(i) for i in [3, 5][::-1]]

    return markers


def Trendline(ax, X, Y, order=1, intervals=700, f_size=20, color='blue',
              lw=1, debug=False):
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
    debug (bool): legacy debug option, replaced by python logging

    Returns
    -------
    (None)
    """
    # Poly fit data
    params, xp = np.polyfit(X, Y, order), \
        np.linspace(min(np.ma.min(X), np.ma.min(Y)),
                    max(np.ma.max(X), np.ma.max(Y)), intervals)
    print((params, type(params)))
    logging.debug('params: {}'.format(str(params)))
    yp = np.polyval(params, xp)

    # Calculate regression fit
    r_sq = [r_squared(X, i) for i in [Y]]

    # Plot up
    ax.plot(xp, yp, ls='--', lw=lw, color=color,
            label=' (R$^{2}$'+'={0:<,.3f}, y={1:,.3f}x+{2:.3f})'.format(
                r_sq[0], params[0], params[1]))
    # Add legend to plot
    ax.legend()


def plot_gc_bin_bands(facecolor='#B0C4DE'):
    """
    Plot highlighted lines of ("tropospheric") GEOS-Chem/Fast-J photolysis bins
    """
    vars = [GC_var(i) for i in ('FastJ_lower', 'FastJ_upper')]
    alphas = [0.3, 0.1]*len(vars[0])
    [plt.axvspan(vars[0][n], vars[1][n], facecolor=facecolor,
                 alpha=alphas[n]) for n in range(len(vars[0]))]


def get_ls(num):
    """
    Get a list of available line styles
    """
    ls = [
        ':', '--', '-.', '-', ':', '--', '-.', '-', ':', ':', '--',
        '-.', '-', ':', '--', '-.', '-', ':'
    ]
    return ls[:num]


def greyoutstrat(fig,  arr, axn=[1, 1, 1], ax=None, cmap=plt.cm.bone_r,
                 res='4x5', rasterized=True, debug=False):
    """
    Grey out stratosphere in existing zonal plot.

    NOTES
    -------
     - This is used to highlight tropospherically focused work. This function
     requires the array of "time in the troposphere" diagnostic (lon,lat, alt)
    """
    # Get overall vars
    lon, lat, alt = get_latlonalt4res(res=res)
    alt = [i[:len(arr[0, :])] for i in [alt]][0]

    # plot up a grey section for stratosphere
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(*axn)
    arr = np.ma.around(arr)
    arr = np.ma.masked_where(arr > 0, arr)
    p = ax.pcolor(lat, alt, arr.T, cmap=cmap)
    if rasterized:
        plt.gcf().set_rasterized(True)


def adjust_subplots(fig, left=None, bottom=None, right=None, top=None,
                    wspace=None, hspace=None):
    """
    Set subplot adjustment values for provided figure
    """
    # the left side of the subplots of the figure
    if isinstance(left, type(None)):
        left = 0.125
    # the right side of the subplots of the figure
    if isinstance(right, type(None)):
        right = 0.9
    # the bottom of the subplots of the figure
    if isinstance(bottom, type(None)):
        bottom = 0.1
    # the top of the subplots of the figure
    if isinstance(top, type(None)):
        top = 0.9
    # the amount of width reserved for blank space between subplots
    if isinstance(wspace, type(None)):
        wspace = 0.2
    # the amount of height reserved for white space between subplots
    if isinstance(hspace, type(None)):
        hspace = 0.5
    # Adjust subplots
    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                        wspace=wspace, hspace=hspace)


def mask_not_obs(loc='Denmark', res='4x5', debug=False):
    """
    provide a mask of all regions apart from the location given
    """
    # Start with all zeros
    arr = np.zeros(get_dims4res(res))
    # Get lats and lons of locations to keep...
    lats, lons = get_obs_loc(loc)
    # Unmask locations
    lats = [get_gc_lat(i, res=res) for i in lats]
    lons = [get_gc_lon(i, res=res) for i in lons]
    for n, lat in enumerate(lats):
        arr[lons[n], lat, :] = 1
    # Return values
    return np.ma.not_equal(arr, 1)


def annotate_gc_grid(ax, res='4x5', f_size=6.5,
                     loc_list=[[-9999, -9999]], everyother=1,
                     label_gc_grid=True):
    """
    Annotate grid with GEOS-Chem indices
    """
    # Get Vars
    lon, lat, alt = get_latlonalt4res(res=res)
    if res == '0.5x0.666':
        adjust_window, interval, nticks, nbins, resolution, shrink = 3, 0.5, \
            int(nticks/3), int(nbins/2), 'l', 0.6
    glat = [get_gc_lat(i, res=res) for i in lat]
    glon = [get_gc_lon(i, res=res) for i in lon]

    # set ticks to all
    ax.xaxis.set_ticks(lon[::everyother])
    ax.yaxis.set_ticks(lat[::everyother])
    plt.xticks(fontsize=f_size*2.5)
    plt.yticks(fontsize=f_size*2.5)
    plt.grid(True)

    # label grid boxes by number
    if label_gc_grid:
        for n, la in enumerate(lat):
            for nn, lo in enumerate(lon):
                if any([[glat[n], glon[nn]] == i for i in loc_list]):
                    color = 'red'
                    fontweight = 'bold'
                else:
                    color = 'blue'
                    fontweight = 'normal'
                ax.text(lo+1, la+1, '{},{}'.format(glon[nn], glat[n]),
                        fontsize=f_size, color=color)


def shiftedColorMap(cmap, start=0, midpoint=0.5, lower=0, upper=1,
                    start_center=0.5, stop=1.0, maintain_scaling=True,
                    arr=None,
                    name='shiftedcmap', npoints=257, verbose=True,
                    debug=False):
    """
    ORGINAL DESCRIPTION: Function to offset the "centre" of a colormap. Useful
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
        vmin, vmax = arr.min(), arr.max()

        # Symetric range
        eq_range = max(abs(vmin), abs(vmax)) * 2
        # Actual range
        act_range = abs(vmin) + abs(vmax)
        cut_off = act_range/eq_range

        if midpoint > 0.5:
            start_center = start_center/cut_off
            lower, upper = start, cut_off
        if midpoint < 0.5:
            start_center = start_center*cut_off
            lower, upper = 1-cut_off, stop
        if midpoint != 0.5:
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=lower,
                                                    b=upper),
                cmap(np.linspace(lower,
                                 upper, npoints))
            )
            logging.debug('maintain scale for vars: {}'.format(
                *[str(float(i)) for i in list([eq_range, act_range, vmin,
                                               vmax, abs(vmin), abs(
                                                   vmax), cut_off, start_center, lower,
                                               upper, midpoint, start, stop])]))
    cdict = {'red': [], 'green': [], 'blue': [], 'alpha': []}

    # Regular index to compute the colors
    reg_index = np.hstack([np.linspace(start, start_center, 128,
                                       endpoint=False),
                           np.linspace(start_center, stop, 129,
                                       endpoint=True)]
                          )

    # Shifted index to match the data
    shift_index = np.hstack([np.linspace(0.0, midpoint, 128, endpoint=False),
                             np.linspace(midpoint, 1.0, 129, endpoint=True)])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


def mk_cb(fig, units=None, left=0.925, bottom=0.2, width=0.015, height=0.6,
          orientation='vertical', f_size=20, rotatecbunits='vertical',
          nticks=10,
          extend='neither', norm=None, log=False, format=None, cmap=None,
          vmin=0, vmax=10, cb_ax=None, ticklocation='auto', extendfrac=None,
          cb_sigfig=2, lvls=None, discrete_cmap=False,
          boundaries=None, verbose=True, debug=False):
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
    log (bool): use log scale for colorbar
    format (str): format string for colorbar tick formating
    cmap (str): colormap instance
    vmin, vmax (float): miminium and maximum of colorbar
    extendfrac (float): fraction by which to extend "pointers" on colorbar
    cb_sigfig (int): significant figure rounding to use for colourbar
    lvls (list): manually provide levels for colorbar
    boundaries (list):  manually provide levels for discrete colorbar
    discrete_cmap (bool): use a discrete instead of conitunous colorbar map
    verbose (bool): legacy debug option, replaced by python logging
    debug (bool): legacy debug option, replaced by python logging

    Returns
    -------

    Notes
    -----
     - This function allows for avoidance of basemap's issues with spacing
    definitions when conbining with colorbar objects within a plot
     (This function is not specific to matplotlib basemap)
    """
    if debug:
        print(('mk_cb called with: ', norm, vmin, vmax, log, lvls))

    # Get colormap (by feeding get_colormap array of min and max )
    if isinstance(cmap, type(None)):
        cmap = get_colormap(arr=np.array([vmin, vmax]))

    if debug:
        print((left, bottom, width, height))
    if debug:
        print((vmin, vmax, orientation, ticklocation, extend, extendfrac))

    # Make new axis
    if isinstance(cb_ax, type(None)):
        cb_ax = fig.add_axes([left, bottom, width, height])
        if debug:
            print(('>'*5, left, bottom, width, height))

    # Setup normalisation for colorbar
    if isinstance(norm, type(None)):
        if log:
            # Get logarithmically spaced integers
            lvls = np.logspace(np.log10(vmin), np.log10(vmax), num=nticks)
            # Normalise to Log space
            norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)

        else:
            # Normalise to linear space
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    if verbose:
        print((lvls, vmin, vmax, norm))

    if isinstance(lvls, type(None)):
        # make graduations in colourbar to be human readable
        lvls = get_human_readable_gradations(vmax=vmax,
                                             vmin=vmin, nticks=nticks,
                                             cb_sigfig=cb_sigfig)

    # add an extra decimal place for values under 50,
    # and two for values under 1
    if isinstance(format, type(None)):
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
            if ((abs(int(vmax)) == 0) and (abs(int(vmax)) == 0)):
                #                format='%.0E'  # hashed prior
                format = '%.1E'
                f_size = f_size * .75

        # warning this will need to be restored?! - plotting routines use this
        # This exception is for highly masked arrays
        except:
            format = '%.0E'
            # Kludge, set to limits of -500-500 (%) if overly masked.
            vmax, vmin = 500, -500
            extend = 'both'

    if debug:
        print((lvls, norm, extend, format, ticklocation))

    # Make cb with given details
    if discrete_cmap:
        extendfrac = .075  # This was the previous default, still valid? <= update
#        if debug:
#            print 'WARNING: adding extensions to colorbar'
        if log == True:
            print(('Will not work as colrbar adjustment not configured' +
                   ' for log scales'))
            sys.exit(0)
        else:
            print(lvls)
            if discrete_cmap:
                boundaries = lvls
            else:
                extend_points = (lvls.max()-lvls.min())*extendfrac
                boundaries = [lvls.min()-extend_points] + list(lvls) + \
                    [extend_points+lvls.max()]

        # the following code isn't active, why is it here? <= update
        # increase the space allowed for the colorbar
        # (the extension  is not summative )
#        if orientation=='horizontal':
#            width += width*extendfrac
#        if orientation=='vertical':
#            height += height*extendfrac

        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, format=format,
                                       norm=norm, ticks=lvls, extend=extend,
                                       extendfrac=extendfrac,
                                       boundaries=boundaries,
                                       spacing='proportional',
                                       #                spacing='uniform',
                                       orientation=orientation,
                                       ticklocation=ticklocation)

    # Standard approach below
    else:
        if verbose:
            print((lvls, norm, boundaries, extend, orientation, ticklocation,
                   cb_ax))
        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, format=format,
                                       norm=norm, ticks=lvls, extend=extend,
                                       orientation=orientation,
                                       ticklocation=ticklocation)

    if log:
        def round_to_n(x, n): return round(x, -int(floor(log10(x))) + (n - 1))
        cb.set_ticks([float('{:.2g}'.format(t)) for t in lvls])
        labels = [round_to_n(i, cb_sigfig) for i in lvls]
        cb.set_ticklabels([format % i for i in labels])

    # Set cb label sizes
    if not isinstance(units, type(None)):
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(f_size)
        if rotatecbunits == 'vertical':
            cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size,
                             fontsize=f_size)
        else:
            cb.set_label(units, fontsize=f_size)
    # set tick sizes regardless whether units (labels) are provided
    cb.ax.tick_params(labelsize=f_size)  # , size=f_size )

    return cb_ax


def get_colormap(arr,  center_zero=True, minval=0.15, maxval=0.95,
                 npoints=100, cb='CMRmap_r', maintain_scaling=True,
                 negative=False, positive=False, divergent=False,
                 cb_sigfig=2, buffer_cmap_upper=False, fixcb=None, nticks=10,
                 verbose=True, debug=False):
    """
    Create *correct* colormap by checking if it contains just +ve or -ve or
    both then prescribing color map accordingly.

    Parameters
    ----------
    arr (array): array of values to assess colourbar from
    center_zero (bool): make sure (divergent) colorbar centred around zero
    minval, maxval (float): values to restrict 'gnuplot2' to
    npoints (int): number of points in colormap
    cb (str): name of colorbar string
    maintain_scaling (bool): maintain scaling for range in color change
    negative (bool): force colormap to be sequential negative (==True)
    positive (bool): force colormap to be sequential positive (==True)
    divergent (bool): force colormap to be divergent (==True)
    cb_sigfig (int): number of sig. figs. to round colourbar ticks
    buffer_cmap_upper (bool): make sure colorbar has space for maximum val.
    fixcb (array): lower and upper values to fix colourmap to.
    nticks (int): number of ticks to use for colorbar
    verbose (bool): legacy debug option, replaced by python logging
    debug (bool): legacy debug option, replaced by python logging

    Returns
    -------
    (colormap instance)

    Notes
    -----
     - this function also will can adjust colormaps to fit a given set of ticks
    """
    # Mannual fix maintain scaling to False
#    maintain_scaling=False

    logging.info('get_colormap called with fixcb={}'.format(str(fixcb)))
    # Manually override colourbar?
#    cb='Blues' # Kludge.
#    cb='Reds' # Kludge.

    # Make sure cmap includes range of all readable levels (lvls)
    # i.e head of colormap often rounded for ascetic/readability reasons
    if buffer_cmap_upper:
        lvls, lvls_diff = get_human_readable_gradations(vmax=fixcb[1],
                                                        vmin=fixcb[0],
                                                        nticks=nticks,
                                                        rtn_lvls_diff=True,
                                                        cb_sigfig=cb_sigfig
                                                        )

        # increase maximum value in color by 5% of level diff
        # to allow space for max lvl
        fixcb_ = (fixcb[0],  lvls[-1] + (lvls_diff*0.05))
        arr = np.array(fixcb_)

        # Gotcha - make sure max(lvls)  not increased > 0
        if (max(lvls) <= 0) and (not (arr.max() < 0)):
            arr = np.array([arr[0], 0])

    # make sure array has a mask
    logging.debug("arr min and max vals: {}, {} - arr type: {}".format(
        str(arr.min()), str(arr.max()), type(arr)))
    if 'ask' not in str(type(arr)):
        arr = np.ma.array(arr)
#        s_mask = arr.mask
    logging.debug('arr type (post mask check) {}:'.format(type(arr)))

    # If postive/negative not given, check if +ve/-ve
    if (not positive) or (not negative):
        logging.debug('Testing if arr is +ve/-ve, for arr with min' +
                      '{} and max {}'.format(arr.min(), arr.max()))

        # --- sequential ?
        # negative?
        arr.mask = False
        arr.mask[arr <= 0] = True
        if arr.mask.all():
            negative = True

        # postive?
        arr.mask = False
        arr.mask[arr >= 0] = True
        if arr.mask.all():
            positive = True

    # Reverse colourbar if negative
    if negative:
        if cb == 'CMRmap_r':
            cb = cb[:-2]
        else:
            cb = cb+'_r'
    # log colorbar details
    logging.info('cmap is: >{}< & data is:'.format(cb))
    logging.info('< postive == {}, negative == {}, divergent == {} >'.format(
        positive, negative, ((not positive) and (not negative))))

    # load color map
    cmap = plt.get_cmap(cb)

    # Chop off bottom for 'gnuplot2'
    if (negative and ('gnuplot2' in cb)) or (positive and ('CMRmap' in cb)):
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=0,
                                                b=maxval-minval, ),
            cmap(np.linspace(0,
                             maxval-minval,
                             npoints)))

    if (positive and ('gnuplot2' in cb)) or (negative and ('CMRmap' in cb)):
        #    if positive and ( 'CMRmap' in cb ):
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval,
                                                b=maxval),
            cmap(np.linspace(minval,
                             maxval,  npoints))
        )

    # --- divergent
    if ((not positive) and (not negative)) or divergent:
        logging.debug('Data is divergent')
        arr.mask = False
#        cmap = plt.cm.RdBu_r
        cb = 'RdBu_r'
        cmap = plt.get_cmap(cb)
#        cmap = plt.cm.Spectral
        # Centre color bar around zero
        if center_zero:
            vmin, vmax = arr.min(), arr.max()
            if buffer_cmap_upper:
                pass
            logging.debug('vals. for mid point {}, min {}, max {}'.format(
                str(1-vmax/(vmax + abs(vmin))), str(vmin), str(vmax)))
            cmap = shiftedColorMap(cmap, midpoint=1-vmax/(vmax + abs(vmin)),
                                   maintain_scaling=maintain_scaling, arr=arr,
                                   name='shifted',
                                   verbose=verbose, debug=debug)

    # Use updated jet replacements from mpl
    # : https://github.com/bids/colormap
#    cmap  = cmc  # pink, orange, blue alterative...
#    cmap  = cmd  # green - blue alternative

#    arr.mask = s_mask

    logging.debug(
        'colorbar used: {} (center_zero=={})'.format(cb, center_zero))

    if buffer_cmap_upper:
        return cmap, fixcb_
    else:
        return cmap


def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


def colorline(x, y, z=None, cmap=plt.get_cmap('copper'),
              norm=None, linewidth=3, alpha=1.0, ax=None, fig=None,
              cb_title='Modelled wind direction', convert_x_2datetime=True,
              debug=False):
    """
    plot a line (x,y) for plot coloured by variable (z)
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    if isinstance(norm, type(None)):
        norm = plt.Normalize(0.0, 1.0),

    # Special case if a single number:
    # to check for numerical input -- this is a hack
    if not hasattr(z, "__iter__"):
        z = np.array([z])

    z = np.asarray(z)
    # printing for debugging?
    if debug:
        print([type(i[0]) for i in (x, z, y)])
    if debug:
        print([i[0] for i in (x, z, y)])

    # convert to datetime.datetime
    if convert_x_2datetime:
        # --- make sure x axis is float ( not a datetime )
        # convert to dataframe
        df = DataFrame(data=y, index=x)
        x = np.array(df.index.to_pydatetime(), dtype=datetime.datetime)
        if debug:
            print((type(x[0]), x[0]))

        # convert in float form
        x = matplotlib.dates.date2num(x)
        if debug:
            print((type(x[0]), x[0]))

    # convert x and y into indices
    segments = make_segments(x, y)

    lc = mpl.collections.LineCollection(segments, array=z, cmap=cmap,
                                        norm=norm,
                                        linewidth=linewidth, alpha=alpha)

    if isinstance(ax, type(None)):
        ax = plt.gca()
    ax.add_collection(lc)

    if not isinstance(fig, type(None)):
        axcb = fig.colorbar(lc)
        axcb.set_label(cb_title)

    return lc


def get_human_readable_gradations(lvls=None, vmax=10, vmin=0,
                                  nticks=10, cb_sigfig=2,
                                  cb_sigfig_ticks=2,
                                  cb_sigfig_lvls=2, rtn_lvls_diff=False,
                                  verbose=True, debug=False):
    """
    Get human readible gradations for plotting (e.g. colorbars etc)

    OUTPUT:
    lvls: list of values where the ticks should be
    ticks: list of strings to call the ticks in Sig figs.
    """

    if not np.isfinite(vmin) or not np.isfinite(vmax):
        cb_error = "Colourbar has a NaN in it. vmin={vmin}, vmax={vmax}"\
            .format(vmin=vmin, vmax=vmax)
        logging.error(cb_error)
        raise ValueError(cb_error)

    logging.debug('get_human_readable_gradiations called with the following:')
    logging.debug('vmin = {vmin}, vmax = {vmax}, lvls = {lvls}'
                  .format(vmin=vmin, vmax=vmax, lvls=lvls))

    if isinstance(lvls, type(None)):
        lvls = np.linspace(vmin, vmax, nticks, endpoint=True)

#    verbose=True
    # --- Adjust graduations in colourbar to be human readable
    # in both min and max have absolute values less than 0, then sig figs +1
    # Find the amount of significant figures needed to show a difference
    # between vmin and vmax.
    if (vmin == np.NaN) and (vmin == np.NaN):
        err_str = "Both vmin and vmax are NaNs!"
        logging.error(err_str)
        raise ValueError(err_str)
    elif (vmin == vmax):
        err_str = "There is no difference between vmin and vmax!"
        logging.error(err_str)
        raise ValueError(err_str)
    try:
        if ((abs(int(vmin)) == 0) and (abs(int(vmax)) == 0)):
            cb_sigfig += 1
        logging.debug("Significant figures needed for plot is {sf}"
                      .format(sf=cb_sigfig))
    except np.ma.core.MaskError:
        print('Gotcha: numpy.ma.core.MaskError')
        print((lvls, vmin, vmax))


#    # bjn updated sigfig finder
#    # Use logs to find out sig figs needed
#    if vmin==0 or vmax==0:
#        sig_figs_needed = 3
#    else:
#        log_diff = abs( np.log10(abs(vmax)) - np.log10(abs(vmin)) )
#        sig_figs_needed = int(np.ceil(abs(np.log10( log_diff ))))
#
#    cb_sigfig_ticks = sig_figs_needed

    # significant figure ( sig. fig. ) rounding func.

    def round_to_n(x, n): return round(x, -int(floor(log10(x))) + (n - 1))
#    round_to_n = lambda x, n: get_sigfig(x,n)

    # --- Get current gradations
#    if debug:
#        print abs(lvls[-4])-abs(lvls[-3]), abs(lvls[-4])-abs(lvls[-3]), lvls,\
#                     cb_sigfig
    try:
        lvls_diff = [round_to_n(abs(i-l[n+1]), cb_sigfig_ticks)
                     for n, i in enumerate(l[:-1])]
        lvls_diff = list(set(lvls_diff))
        if len(lvls_diff) > 1:
            lvls_diff = max(lvls_diff)
#        lvls_diff = round_to_n( abs(lvls[-3])-abs(lvls[-4]), \
#                                cb_sigfig_ticks)

    # handle if values (2,3) are both negative or abs. of both <0
    except:
        debug_list = (abs(lvls[-4])-abs(lvls[-3]), cb_sigfig_ticks)
        if debug:
            print(debug_list)
        try:    # handle if values (2,3) are both negative
            lvls_diff = round_to_n(abs(lvls[-4])-abs(lvls[-3]),
                                   cb_sigfig_ticks)
        except:  # If both absolute of vmin and vmax  are <0 ( and +ve )
            debug_list = (lvls, lvls[-3], lvls[-4],
                          cb_sigfig_ticks)
            if debug:
                print(debug_list)
            lvls_diff = round_to_n(lvls[-3]-lvls[-4],
                                   cb_sigfig_ticks)

    # ---  Round top of colorbar lvls, then count down from this
    # first get top numer rounded up to nearest 'lvls_diff'
    # caution, this may result in a number outside the cmap,
    # solution: use get_colormap, with buffer_cmap_upper=True
    #  if values are >0,
    if vmax > lvls_diff:
        if debug:
            print((vmax, lvls_diff))
        vmax_rounded = myround(vmax, base=lvls_diff,  integer=False)
        vmax_rounded = round_to_n(vmax_rounded, cb_sigfig)
    else:
        # <= update needed! - add function to round negative numbers
        # ( this method also fails if vmax<lvls_diff )
        vmax_rounded = vmax

#    if debug:
#        print 1, lvls, vmax_rounded, lvls_diff, cb_sigfig_lvls

    lvls = np.array([vmax_rounded - lvls_diff*i
                     for i in range(nticks)][::-1])

    logging.debug("colorbar levels are: {lvls}".format(lvls=lvls))
#    if debug:
#        print lvls, len( lvls )
#        print 2, lvls, vmax_rounded, lvls_diff, cb_sigfig_lvls

    # ensure returned ticks are to a maximum of 2 sig figs
    # ( this only works if all positive ) and are unique
#    try:
#        # Make sure the colorbar labels are not repeated
#        invalid = True
#        while invalid:
#            new_lvls = [ round_to_n( i, cb_sigfig_lvls) \
#                for i in lvls ]
#            if len( set(new_lvls) ) == len(lvls):
#                lvls = new_lvls
#                invalid=False
#            else: # Try with one more sig fig
#                cb_sigfig_lvls += 1
#
#    except:
#        print 'WARNING: unable to round level values to {} sig figs'.format(\
#                   cb_sigfig_lvls  )

    new_lvls = []
    for level in lvls:
        new_lvls.append(get_sigfig(level, cb_sigfig_lvls))

    lvls = new_lvls

#    if any([ (type(i) == str) for i in lvls] ):
#        logging.debug('WARNING str in list of levels, tmp conversion to float'+
#         ' - This is due to error in get_sigfig - seting "0.0" to float' )
#        for n, i in enumerate( lvls ):
#            if i == '0.0':
#                lvls[n] = 0.0

#    if debug:
#        print 3, lvls, vmax_rounded, lvls_diff, cb_sigfig_lvls
#    print [(i, type(i)) for i in vmin, vmax ], lvls
#    print [(i, type(i)) for i in lvls ]

    if rtn_lvls_diff:
        return lvls, lvls_diff
    else:
        return lvls


def mk_discrete_cmap(lvls=None, cmap=None, rtn_norm=False,
                     vmin=None, vmax=None, nticks=10, debug=False):
    """
    Make a discrete colormap from an existing cmap

    Parameters
    ----------
    lvls (np.array): array of the bounding got discrete
    cmap (colormap/str): colour map object to discretise
    nticks (int): number of sections to discretise colourbar into
    vmin/vmax (float): the min and max of the norm object to be created
    rtn_norm (bool): return a norm object as well as colourbar?

    Returns
    -------
    (colormap object, mpl norm object)

    Notes
    -------
     - the data can be optionally to normalised the range of lvls
    """
    # Get colormap if not provided and make sure it is
    if isinstance(cmap, type(None)):
        cmap = get_colormap(np.array([vmin, vmax]))
    if isinstance(cmap, str):
        cmap = matplotlib.cm.get_cmap(cmap)
    if debug:
        print((vmin, vmax, nticks))
    # Extract colors linearly from colormap
    cmaplist = cmap(np.linspace(0, 1, nticks))
    # Create the new discrete colormap object
    cmap_name = '{}_{}'.format(cmap.name, str(nticks))
    cmap = cmap.from_list(cmap_name, cmaplist, nticks)
    # Make a normalisation object... -  define the bins and normalize
    if rtn_norm:
        # Define bins
        if isinstance(lvls, type(None)):
            bounds = np.linspace(vmin, vmax, nticks, endpoint=True)
        else:
            bounds = lvls
        norm = mpl.colors.BoundaryNorm(bounds,  nticks)
        return cmap, norm
    else:
        return cmap


def show_plot():
    """
    Wrapper for plt.show(). Use to plot to screen
    """
    plt.show()
    return


def close_plot():
    """
    Wrapper for plt.close(). Use to closes plots on screen
    """
    plt.close()
    return


def save_plot(title="myplot", location=os.getcwd(),  extensions=['png'],
              tight=False, dpi=320):
    """
    Save a plot to disk

    Parameters
    ----------
    title (str): plot title
    location (str): String of directory to save output
    extensions (list):  List of strings of file extensions. (e.g. ['png'])
    tight (bool): use tight layout - redundant?

    Returns
    -------
    (None)
    """
    if tight:
        plt.tight_layout()
    if not os.path.isdir(location):
        os.mkdir(location)
        logging.warning("Plotting location not found")
        logging.warning("Creating dir {folder}".format(folder=location))

    for extension in extensions:
        filename = os.path.join(location, title+"."+extension)
        plt.savefig(filename, dpi=dpi)
        logging.info("Plot saved to {location}".format(location=filename))
    return


def set_bp_style(bp, color='k', linewidth=1.5, facecolor='none',
                 debug=False):
    """
    Manually set properties of boxplot ("bp")
    """
    from pylab import setp
    setp(bp['boxes'][:], color=color, linewidth=linewidth)
    setp(bp['caps'][:], color=color, linewidth=linewidth)
    setp(bp['whiskers'][:], color=color, linewidth=linewidth)
    setp(bp['fliers'][:], color=color, linewidth=linewidth)
    setp(bp['medians'][:], color=color, linewidth=linewidth)
    [box.set(facecolor=facecolor) for box in bp['boxes']]
    setp(bp['medians'][:], color=color, linewidth=linewidth)


def get_CB_color_cycle():
    """
    Get a list of color blind friednly colors to cycle in plots

    Notes
    -----
     - Credit @ thriveth: https://gist.github.com/thriveth/8560036
    """
    CB_color_cycle = [
        '#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3',
        '#999999', '#e41a1c', '#dede00'
    ]
    return CB_color_cycle


def plot_vertical_fam_loss_by_route(fam='LOx', ref_spec='O3',
                                    wd=None, Mechanism='Halogens',
                                    rm_strat=False,
                                    weight_by_molecs=True, CODE_wd=None,
                                    full_vert_grid=True, dpi=320,
                                    suffix='',
                                    save_plot=True, show_plot=False,
                                    limit_plotted_alititude=True, lw=16,
                                    Ox_fam_dict=None, fontsize=10,
                                    cmap=plt.cm.jet, alt_array=None,
                                    verbose=True, debug=False):
    """
    Plot vertical odd oxygen (Ox) loss via route (chemical family)

    Parameters
    -------
    fam (str): tagged family to track (already compiled in KPP mechanism)
    ref_spec (str): reference species to normalise to
    wd (str): working directory ("wd") of model output
    CODE_wd (str): root of code directory containing the tagged KPP mechanism
    Mechanism (str): name of the KPP mechanism (and folder) of model output
    weight_by_molecs (bool): weight grid boxes by number of molecules
    rm_strat (bool): (fractionally) replace values in statosphere with zeros
    debug, verbose (bool): switches to turn on/set verbosity of output to screen
    alt_array (np.array): array of altitudes for model grid boxes
    full_vert_grid (bool): use the full vertical grid for analysis
    limit_plotted_alititude (bool): limit the plotted vertical extend to troposphere
    suffix (str): suffix in filename for saved plot
    dpi (int): resolution to use for saved image (dots per square inch)
    Ox_fam_dict (dict), dictionary of Ox loss variables/data (from KPP.py)

    Returns
    -------
    (None)

    Notes
    -----
     - AC_tools includes equivlent functions for smvgear mechanisms
    """
    # - Local variables/ Plot extraction / Settings
    # extract variables from data/variable dictionary
    sorted_fam_names = Ox_fam_dict['sorted_fam_names']
    fam_dict = Ox_fam_dict['fam_dict']
    ars = Ox_fam_dict['ars']
    RR_dict_fam_stioch = Ox_fam_dict['RR_dict_fam_stioch']
    RR_dict = Ox_fam_dict['RR_dict']
    tags2_rxn_num = Ox_fam_dict['tags2_rxn_num']
    tags = Ox_fam_dict['tags']
    tags_dict = Ox_fam_dict['tags_dict']
    # Combine to a single array
    arr = np.array(ars)
    if debug:
        print((arr.shape))
    # - Process data for plotting
    fam_tag = [fam_dict[i] for i in tags]
    fam_ars = []
    for fam_ in sorted_fam_names:
        # Get indices for routes of family
        fam_ind = [n for n, i in enumerate(fam_tag) if (i == fam_)]
        if debug:
            print((fam_ind, len(fam_ind)))
        # Select these ...
        fam_ars += [arr[fam_ind, ...]]
    # Recombine and sum by family...
    if debug:
        print(([i.shape for i in fam_ars], len(fam_ars)))
    arr = np.array([i.sum(axis=0) for i in fam_ars])
    if debug:
        print((arr.shape))
    # - Plot up as a stack-plot...
    # Normalise to total and conver to % (*100)
    arr = (arr / arr.sum(axis=0)) * 100
    # Add zeros array to beginning (for stack/area plot )
    arr_ = np.vstack((np.zeros((1, arr.shape[-1])), arr))
    # Setup figure
    fig, ax = plt.subplots(figsize=(9, 6), dpi=dpi,
                           facecolor='w', edgecolor='w')
    # Plot by family
    for n, label in enumerate(sorted_fam_names):
        # Print out some summary stats
        if verbose:
            print(n, label, arr[:n, 0].sum(axis=0),
                  arr[:n+1, 0].sum(axis=0), end=' ')
            print(arr[:n, :].sum(), arr[:n+1, :].sum())
            print([i.shape for i in (alt_array, arr)])
        # Fill between X
        plt.fill_betweenx(alt_array, arr[:n, :].sum(axis=0),
                          arr[:n+1, :].sum(axis=0),
                          color=cmap(1.*n/len(sorted_fam_names)))
        # Plot the line too
        plt.plot(arr[:n, :].sum(axis=0), alt_array, label=label,
                 color=cmap(1.*n/len(sorted_fam_names)), alpha=0,
                 lw=lw,)
    # Beautify the plot
    plt.xlim(0, 100)
    xlabel = '% of total O$_{\\rm x}$ loss'
    plt.xlabel(xlabel, fontsize=fontsize*.75)
    plt.yticks(fontsize=fontsize*.75)
    plt.xticks(fontsize=fontsize*.75)
    plt.ylabel('Altitude (km)', fontsize=fontsize*.75)
    leg = plt.legend(loc='upper center', fontsize=fontsize)
    # Update lengnd line sizes ( + update line sizes)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(lw/2)
        legobj.set_alpha(1)
    plt.ylim(alt_array[0], alt_array[-1])
    # Limit plot y axis to 12km?
    if limit_plotted_alititude:
        plt.ylim(alt_array[0], 12)
    # Show plot or save?
    if save_plot:
        filename = 'Ox_loss_plot_by_vertical_{}_{}'.format(Mechanism, suffix)
        plt.savefig(filename, dpi=dpi)
    if show_plot:
        plt.show()


def plt_box_area_on_global_map(ds=None, var2use='DXYP__DXYP',
                               LatVar='lat', LonVar='lon',
                               x0=-30, x1=-10, y0=0, y1=25,
                               savename=None, cmap=plt.get_cmap('viridis'),
                               projection=ccrs.Robinson(), alpha=0.5,
                               aspect='auto', figsize=(10, 6),
                               dpi=320):
    """
    Plot a global map with a region (x0,x1,y0,y1) outlined with a box

    Parameters
    -------
    x0, x1 (float): longitude min and max to set box extents to
    y0, y1 (float): latitude min and max to set box extents to
    ds (xr.Dataset): dataset to use to for plotting global map
    var2use (str): variable to use from dataset (ds)
    dpi (int): resolution to use for saved image (dots per square inch)
    savename (str): name for file to save image to
    cmap (cmap object): colour map to use for plotting
    aspect (str): aspect ratio for plot
    alpha (float): transparency of plotted box
    figsize (tuple): figure size to use for plot

    Returns
    -------
    (None)
    """
    # Use a generic dataset from AC_tools data files if not provide with one
    if isinstance(ds, type(None)):
        # Get AC_tools location, then set example data folder location
        filename = inspect.getframeinfo(inspect.currentframe()).filename
        path = os.path.dirname(os.path.abspath(filename))
        folder = path+'/../data/LM/LANDMAP_LWI_ctm_0125x0125/'
        # Get coords from LWI 0.125x0.125 data and remove the time dimension
        ds = xr.open_dataset(folder+'ctm.nc')
    # - Select the data
    # Just get an example dataset
    ds = ds[[var2use]]
    # Check input values for lat and lon range to plotting box extent
    assert y0 < y1, 'y0 must be less than y1'
    assert x0 < x1, 'x0 must be less than x1'
    # Set values region
    bool1 = ((ds.lon >= x0) & (ds.lon <= x1)).values
    bool2 = ((ds.lat >= y0) & (ds.lat <= y1)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Set all values to 1
    arr = ds[var2use].values
    arr[:] = 1
    ds[var2use].values = arr
    # Plot the data
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=projection, aspect=aspect,
                         alpha=alpha)
    ds[var2use].plot.imshow(x=LonVar, y=LatVar, ax=ax, cmap=cmap,
                            transform=ccrs.PlateCarree())
    # Beautify the figure/plot
    ax.coastlines()
    # Force global perspective
    ax.set_global()  # this will force a global perspective
    # Remove the colour-bar and force a tighter layout around map
    fig.delaxes(fig.axes[-1])
    plt.tight_layout()
    # Save to png
    if isinstance(savename, type(None)):
        savename = 'spatial_plot_of_region'
    plt.savefig(savename+'.png', dpi=dpi)
