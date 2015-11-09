# ------------------- Section 0 -------------------------------------------
# -------------- Required modules:
#

# -- Plotting                                                                                       
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import FuncFormatter
import matplotlib.font_manager as font_manager
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import matplotlib.ticker
from matplotlib import ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.collections import LineCollection
from pylab import setp
import functools
import matplotlib
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.mplot3d import Axes3D


# -- Time                                                                                           
import time
import calendar
import datetime as datetime
from datetime import datetime as datetime_

# -- I/O/Admin... 
import gc

# tms
from AC_tools.funcs_vars import *
from AC_tools.funcs4generic import *
from AC_tools.funcs4time import *
#from funcs4obs import * #( need: get_CVO_DOAS_obs ... )
from AC_tools.funcs4pf import *
from AC_tools.funcs4GEOSC import * # wd2ctms

# math
from math import log10, floor

# colormaps
#from option_c import test_cm as cmc
#from option_d import test_cm as cmd

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 1 ----- Common Plot Types
# 1.01 - Global/nested region surface plotter ***
# 1.02 - Zonal plotter ***
# 1.05 - Latitudinal binned plot
# 1.06 - Diurnal binned boxplot
# 1.07 - Diurnal binned plot
# 1.11 - Sonde Plot ***
# 1.12 - DB netCDF plotter
# 1.13 - plot up monthly from data provided from DB netCDF
# 1.14 - plot up monthly timeseries from ...
# 1.15 - plot up daily timeseries from ...
# 1.16 - plot up timeseries from May through Septmeber
# 1.17 - plot up timeseries from July 
# 1.18 - North Pole plot
# 1.19 -  South pole plot
# 1.20 - PDF of monthly surface change plots for given species
# 1.21 - PDF of monthly zonal change plots for given species
# 1.22 - PDF of monthly column change plots for given species
# 1.23 - Column change
# 1.24 - Coupler to select plot type 
# 1.25 - Probability distribution plotter
# 1.26 - X vs. Y plot 
# 1.27 - Scatter 3D cube
# 1.28 - Plot up seasonal output from 4D arr (lon, lat, alt, time)
# 1.29 - Hist X vs. Y.
# 1.30 - PDF of annual surface change plots for given species

# 1.31 - PDF of annual zonal change plots for given species
# 1.32 - PDF of annual column change plots for given species

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 4 ----- Plotting Ancillaries 
# 4.01 - Get percentiles
# 4.02 - Moving average 
# 4.03 - Make color list
# 4.04 - Weighted Average
# 4.05 - Get R^2 for data
# 4.06 - Setup box plot ascetics.
# 4.07 - Get Marker types
# 4.08 - Setup Trendline
# 4.09 - Plot up GC bands 
# 4.10 - Get list of linestyles
# 4.21 - Grey out stratosphere
# 4.22 - Adjust subplots
# 4.24 - Setup diurnal plot 
# 4.25 - 
# 4.26 - 
# 4.27 - Print NCAS & York logos in the bottom corners
# 4.28 - Mask all locations apart from given observational site
# 4.29 - 
# 4.30 - Annotate Grid
# 4.31 - Centre colour bar
# 4.35 - Add side colorbar
# 4.36 - build basemap
# 4.37 - Contruct and fill array of with data from ordinals (Lat, Lon) 
# 4.38 - Get colormap ( dif for if all postive , all negative, mix )
# 4.39 - Retrieves color by grouping of sensitivity study
# 4.40 - Make segments for variable line color plot
# 4.41 - Make colored line for plot
# 4.42 - Get human readable gradations for plot
# 4.43 - mk colourmap discrete 
# 4.99 - Get input variables for  plotting


# ----------------------------- Section 1 ------------------------------------
# -------------- Common Plot Types
#

# ----
# 1.01 - Map plot for given array and resolution (lon, lat
# -----
def map_plot( arr, return_m=False, grid=False, gc_grid=False, centre=False,\
            cmap=None, no_cb=False, cb=None, rotatecbunits='horizontal',  \
            fixcb=None, nbins=25, nticks=10, mask_invalids=False,  \
            format='%.2f', adjust_window=0, f_size=20, alpha=1, \
            set_window=False, res='4x5', ax=None, case='default', units=None, \
            drawcountries=True,  set_cb_ticks=True, title=None,   \
            interval=1, resolution='c', shrink=0.4, window=False, everyother=1,\
            extend='neither', degrade_resolution=False, discrete_cmap=False, \
            lon_0=None, lon_1=None, lat_0=None, lat_1=None, norm=None,\
             debug=False, **Kwargs):
    """ Plots Global/regional 2D (lon, lat) slices. Takes a numpy array and the 
        resolution of the output. The plot extent is then set by this output.

        - GEOS-Chem output configuration
            res ( '4x5' etc... )

        - plot settings
            set_window  ( True ... )

        Kwargs:
        - basemasp settings:
            resolution ( 'c' = coarse, 'f' = fine ...  )

        - colorbar settings
            extend ( 'both', 'min', 'both' ... )
            shrink ( size of colorbar )    """

    if debug:
        print [ [ i.min(), i.max(), i.mean(), type(i) ] for i in [arr] ]

    # Kludge, for pcent arrays with invalid within them, mask for these. 
    if mask_invalids:
        arr = np.ma.masked_invalid( arr )

    # --- Window plot settings
    if window:
        interval = 2   # double interval size 
        degrade_resolution=True
#        nticks, nbins, resolution, shrink  =int(nticks/3), int(nbins/2), 'l', 0.2
    if  res == '0.5x0.666':
        interval,  adjust_window, resolution,shrink  =0.5, 3, 'f', 0.6
    if degrade_resolution:
        resolution = 'l'

    if res == '0.25x0.3125':
#        centre=True
        centre=False
        adjust_window = 6
            
    if debug:
        print '>'*5, [ [ i.min(), i.max(), i.mean(), type(i) ] for i in [arr] ]
    lon, lat, NIU = get_latlonalt4res( res, centre=centre )

    if set_window:
        # Convert lats and lons to GC lats and restrict lats, lons, and arr
        if not isinstance( lat_0, type(None) ):
            gclat_0, gclat_1 = [ get_gc_lat(i, res=res) for i in lat_0, lat_1 ]
            lat = lat[ gclat_0:gclat_1 ]

        if not isinstance( lon_0, type(None) ):
            gclon_0, gclon_1 = [ get_gc_lon(i, res=res) for i in lon_0, lon_1 ]
            lon = lon[ gclon_0:gclon_1]

    # ---- Grid/Mesh values for Lat, lon, & alt + cb
    if isinstance( cmap, type(None) ):
        cmap = get_colormap( arr.copy() )
#    if discrete_cmap:
#        if isinstance( fixcb, type(None) ):
#            cmap, norm = mk_discrete_cmap( vmin=arr.min(), vmax=arr.max(), \
#                    nticks=nticks, cmap=cmap )
#        else:
#            cmap, norm = mk_discrete_cmap( vmin=fixcb[0], vmax=fixcb[1], \
#                    nticks=nticks, cmap=cmap )

    # ----------------  Basemap setup  ----------------  
    x, y = np.meshgrid(lon,lat)
    if debug:
        print 1, len(x), len(y)
 
    # ---- Setup map ("m") using Basemap
    m = get_basemap( lat=lat, lon=lon, resolution=resolution, res=res, \
                everyother=everyother, interval=interval, f_size=f_size, \
                drawcountries=drawcountries )

    # Process data to grid
    x, y = np.meshgrid(*m(lon, lat))
    if debug:
        print 2, 'len:',  [ len(i) for i in x,y,lat,lon ]
        print '>'*5, [ [ i.min(), i.mean(), i.max() ] for i in [arr ] ]

    if (set_window):
        plt.xlim( lon_0, lon_1)
        plt.ylim( lat_0, lat_1 )
    else:
        plt.xlim( lon[0+adjust_window], lon[-1-adjust_window] )
        plt.ylim(lat[0+adjust_window], lat[-1-adjust_window])

    # ----------------  Cases for arrays  ----------------  
    if ( not isinstance( case, int ) ):
        case = {
        'linear':3, 'default':3, 'IO': 1,'limit_to_1_2': 2, 'log': 4,  
        }[case]
    if debug:
        print 3, case, [ np.array(i).shape for i in lon, lat, arr ], cmap, alpha
    
    # -------- colorbar variables...
    # set cb label sizes
    if not isinstance( fixcb, type(None) ):
        tickmin, tickmax = fixcb[0], fixcb[1]
    else:
        tickmin, tickmax  = arr.min(), arr.max()
    # --------------  Linear plots -------------------------------
    # standard plot 
    if ( ( case == 3) or ( case== 9 ) ) and ( isinstance( fixcb, type(None) ) ):
        poly = m.pcolor( lon, lat, arr, cmap=cmap, norm=norm, alpha=alpha )

    if ( ( case == 3) or ( case== 9 ) ) and (not isinstance(fixcb,type(None) )):
        poly = m.pcolor( lon, lat, arr, cmap=cmap, norm=norm, 
                        vmin=fixcb[0], vmax=fixcb[1], alpha=alpha )
#        poly = m.pcolormesh( lon, lat, arr, cmap=cmap, clevs=clevs,\
#                        vmin=fixcb[0], vmax=fixcb[1], alpha=alpha )

    # -----------------  Log plots ---------------------------------------------
    if (case == 4 ) and ( isinstance(fixcb,type(None) )):  # log and cmap
        # CAUTION - norm/discrete cmapvariables not passed <= update this
        poly = m.pcolor(lon, lat, arr, norm=LogNorm(vmin=arr.min(), \
            vmax=arr.max()), cmap=cmap)#'PuBu_r')

        if no_cb:
            pass
        else:
            lvls = np.logspace( np.ma.min(np.ma.log(arr)), \
                np.ma.max(np.ma.log(arr)), nbins )

            cb = plt.colorbar(poly, ax = m.ax, ticks=lvls, format=format,\
                shrink=shrink,alpha=alpha, extend=extend)

        if debug:
            print np.ma.min(np.ma.log(arr)), np.ma.max(np.ma.log(arr)), lvls

    if (case == 4 ) and (not isinstance(fixcb,type(None) )):  # log and set min
        # CAUTION - clevs variables not passed <= update this
        poly = m.pcolor(lon, lat, arr, norm=LogNorm(vmin=fixcb[0], \
            vmax=arr.max()), cmap=cmap)#'PuBu_r')

        if no_cb:
            pass
        else:
            # CAUTION - clevs variables not passed <= update this
            # create log with 1000 p
            lvls = np.logspace( np.ma.log(vmin), np.ma.log(vmax), 1E5, 
                endpoint=True )

            # select relevant indices ( linearly spaced integers )
            lvls =  [ lvls[0] ] + list( lvls[ gen_log_space( 1E5,  \
                nticks*15) ] ) + [ lvls[-1] ]
            lvls= np.array( lvls )

            cb = plt.colorbar(poly, ax = m.ax, ticks=lvls, format=format, \
                 shrink=shrink, alpha=alpha, extend='min')

        if debug:
            print np.ma.min(np.ma.log(arr)), np.ma.max(np.ma.log(arr)), lvls

    # ----------------  Colorbars  ----------------  
    if not no_cb:
        # CAUTION - clevs variables not passed <= update this
        if isinstance( cb, type(None) ):
            cb = plt.colorbar( poly, ax = m.ax, shrink=shrink, alpha=alpha,  \
                        format=format, ticks= np.linspace(tickmin, tickmax, \
                        nticks ), extend=extend )        

        for t in cb.ax.get_yticklabels():
            t.set_fontsize(f_size)

        if units != None:
            cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size)  

    # Set number of ticks
    # FIX NEEDED - this currently doesn't doesn't work for log plots
#    if (case != 3) and (not no_cb) and ( case != 4):
#        if set_cb_ticks:
#            tick_locator = ticker.MaxNLocator( nticks=nticks )
#            cb.locator = tick_locator
#            cb.update_ticks()

    plt.grid( grid )

    if not isinstance( title, type(None) ):
        plt.title(title, fontsize=f_size*1.5)

    return_l = [ plt ] 
    if not no_cb:
        return_l += [ cb ]
    if return_m:
        return_l += [ m ]
    return return_l

# --------
# 1.02 - Zonal plot - log or linear
# --------
def zonal_plot( arr, fig, ax=None, title=None, tropics=False, \
    f_size=10, c_off=37, format='%.2f', interval=None, nticks=7,no_cb=False, \
    units=None, shrink=0.4, alpha=1, res='4x5', window=False, cmap=None, \
    log=False, fixcb=None, lower_limited=False, nlvls=25, xlimit =None, \
    rotatecbunits='horizontal', extend='neither', ylabel=True, \
    set_window=False, lat_0=None, lat_1=None, lat40_2_40=False, \
    xlabel=True, mask_invalids=False, debug=False ):
    """ Creates a zonal plot from provide array of CTM output (lon, lat)
        Input resolution must be provide for non-default (4x5) output 
        
        This function will also apply maskes if set in arguments """
    
    # Kludge, for pcent arrays with invalid within them, mask for these. 
    if mask_invalids:
        arr = np.ma.masked_invalid( arr )    
    
    # Create axis if not provided <= is this back comparible?
    if isinstance( ax, type(None) ): 
        ax = fig.add_subplot(111)

    # Get overall vars
    lon, lat, alt = get_latlonalt4res( res=res )
    if isinstance( cmap, type(None) ):
        cmap=get_colormap( arr.copy() )
    alt  = [ i[:len(arr[0,:])] for i in [ alt ]  ][0]

    # === -  limit ars to a given mask, if stipulated
    if tropics:
        mask =  mask_tropics(res=res)
        arr = arr * mask[0,:,:c_off+1] 
    if lat40_2_40  :
        mask =  mask_tropics(res=res)
        arr = arr * mask[0,:,:c_off+1]         

#    print '!'*200, lat, [ np.array(i).shape for i in lat, alt, arr, arr.T ]

    if set_window:
#        arr = arr[ get_gc_lat(lat_0, res=res):get_gc_lat(lat_1, res=res), :]
        lat = lat[ get_gc_lat(lat_0, res=res):get_gc_lat(lat_1, res=res) ]

#    print '!'*200, lat, [ np.array(i).shape for i in lat, alt, arr, arr.T ]
#    sys.exit( 0 )
    
    if debug:
        print arr.shape
        print [ len(i) for i in lon, lat, alt ], res
        print [ (i.mean(), i.min(), i.max()) for i in [ arr[: ,:c_off]  ] ]
    min, max  =  [ (i.min(), i.max()) for i in [ arr[: ,:c_off]  ] ][0] 

    # Limit cb to top of (GEOS-Chem chemical) troposphere
    alt = alt[:c_off+1]

    # Plot settings for window plots
    if window:
        if interval == None:
            interval = 3
    else:
        interval = 1
    parallels = np.arange(-90,91,15*interval)
#    parallels = np.arange(-90,91,5*interval)

    if (lower_limited):
        extend='min'

    # Is array reduced to chemistry computed troposphere? - if so limit alt
    if len(arr[0,:] ) != 38:
        arr =arr[:,:38]

    # Create a Log Zonal plot of the given data, where fixcb is the min and max of a given array (aka not log )
    if (log):
        if (fixcb != None):

            p = ax.pcolor( lat, alt, arr.T, norm=LogNorm(vmin=fixcb[0], 
                        vmax=fixcb[1]), cmap=cmap)        

            l_f = LogFormatter(10, labelOnlyBase=False)
            if fixcb[0] == 0:
                print 'fixcb inputed as: ', fixcb
                print '>'*5, 'USING VALUE OF 0.01 INSTEAD OF 0 FOR LOG', '<'*30
                lvls = np.logspace( 0.01, np.ma.max(np.log(fixcb[1])), nlvls )
                lower_limited= True
#            lvls = np.logspace( np.log(fixcb[0]), np.max(np.log(fixcb[1])), nlvls )
            else:
                lvls = np.logspace( np.ma.log(fixcb[0]), \
                    np.ma.max(np.ma.log(fixcb[1])), nlvls )                        

            if not no_cb:
                cb = plt.colorbar(p, ax=ax, ticks=lvls, extend=extend, \
                    format=format, shrink=shrink , alpha=alpha)                
        else:
            p = ax.pcolor( lat, alt, arr.T, norm=LogNorm(vmin=min, vmax=max), \
                cmap=cmap)
            if not no_cb:
                cb = plt.colorbar(p, ax =ax, extend=extend, format=format, \
                            shrink=shrink, alpha=alpha)
    else:
        if isinstance( fixcb, type( None) ):
            p = ax.pcolor( lat, alt, arr.T, cmap=cmap, vmin=min, vmax=max)
        else:
            print list(lat)
            print fixcb, [ np.array(i).shape for i in lat, alt, arr, arr.T ]
            p = ax.pcolor( lat, alt, arr.T, cmap=cmap, vmin=fixcb[0],  \
                                    vmax=fixcb[1] )

        if not no_cb:
            cb  = fig.colorbar(p, ax=ax, extend=extend, format=format, \
                            alpha=alpha, ticks= np.linspace(np.ma.min(arr), \
                            np.ma.max(arr), nticks ))

    # Setup Y axis    
    if (not isinstance( units, type( None) )) and (not no_cb):
        cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size)                      

    plt.ylim(alt[0], alt[-1])
    if ylabel:
        plt.ylabel('Altitude / km', fontsize=f_size*.75)
    
    # Setup X axis
    plt.xticks( parallels, fontsize=f_size*.75 ) # draw parrelel lines 
    if tropics:
        plt.xlim( -26, 26 )
    if lat40_2_40:
        plt.xlim( -40, 40 )
    else:
        plt.xlim( lat[0], lat[-1] )
    if not isinstance( xlimit, type(None) ):
        plt.xlim(xlimit[0], xlimit[1])
    if xlabel:
        plt.xlabel('Latitude', fontsize=f_size*.75)
    plt.yticks( fontsize=f_size*.75 )

    # Add title if provided    
    if (title != None):
        plt.title(title, fontsize=f_size*1.5)

# --------   
# 1.05 - Lat plot
# --------
def lat_plot(fig, ax, arr, title=None, f_size=10, units='ppbv', \
            scale='linear', debug=False ):
    """ Creates a latitude plot on given figure and axis """

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
# 1.06 - Diurnal boxplot
# --------
def diurnal_boxplot(fig, ax,  dates, data, pos=1, posn =1,  bin_size=2/24.,\
            widths=0.01, white_fill=True, alpha=0.1, linewidth=0.5, \
            showmeans=False, title=None, f_size=10, units='ppbv', \
            scale='linear', debug=False ):
    """ Creates a diurnal plot of boxplots (hourly) for given data and dates. 
            Data and dates must be in numpy array form. 
            Dates must also be datetime.datetime objects """

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


# --------   
# 1.07 - Diurnal plot
# --------
def diurnal_plot(fig, ax,  dates, data, pos=1, posn =1,  \
            bin_size=2/24.,widths=0.01, rmax=False, \
            ls='-', color=None, fractional=False, diurnal=False, mean=True, \
            xlabel = True, r_avgs=False, marker=None, label=None, \
            markersize=1, title=None, f_size=10, units='ppbv', scale='linear', \
            lw=1,lgnd_f_size=None, alpha=1, debug=False ):
    """ Creates a diurnal plot for given data and dates. 
            Data and dates must be in numpy array form. 
            Dates must also be datetime.datetime objects """

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
    avg = np.mean(avgs)
#    print avgs, avg, np.ma.max( avgs )
    
    if fractional:
        y = ( avgs - np.ma.max( avgs )  )  / avgs *100
        ymin, ymax =  -0.15, 0.025 
        ymin, ymax =  [i*100 for i in ymin, ymax ]
        units = '%'
    elif diurnal:        
        y =  avgs - np.ma.max( avgs )
        ymin, ymax =  -3.5 , 0.5
    else:
        y = avgs
        ymin, ymax =  None, None#27, 37
    
    # Plot
    plt.plot( bins_used, y, color=color , label=label, linestyle=ls[pos-1], \
                   alpha=alpha, marker=marker,  lw=lw, ms=markersize )

    # Beautify 
    ax.set_xticklabels( np.arange(0,24,1 )[::2]  )
    plt.xticks( np.arange(0,1,1/24. )[::2], fontsize=f_size )
    plt.xlim(0., 23/24.)

    if ymin != None:    
        plt.ylim( ymin, ymax )
    plt.ylabel(units, fontsize=f_size*.75)
    plt.yticks( fontsize=f_size*.75)
    if (title != None):
        plt.title(title)
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
# 1.11 - Plot up Sonde data
# --------
def sonde_plot(fig, ax, arr, n=0, title=None, subtitle=None, tropics=False, \
            f_size=10, color=None, rasterized=True, err_bar=False, obs=True, \
            legend=False, units='ppbv', hPa_labels=True, stddev=True, \
            c_l=[ 'k', 'red','green', 'blue' , 'purple'], xlimit=None, \
            loc='upper left', c_off = 37, label=None, ancillary=True, \
            debug=False ):
    """ Create plot of vertical data for sonde observations/model """

    # Get overall vars
    alt, press = [ gchemgrid(i) for i in 'c_km_geos5_r' , 'c_hPa_geos5_r']

    # Cut off at Tropopause?
    if obs:
        arr = arr[:c_off,:] 
    else:
        arr = arr[:c_off]         
    alt, press = [ i[:c_off] for i in alt, press ]

    # if color  not set, get color
    if isinstance( color, type(None) ):
        color=c_l[n]
        
    if obs:
#        print [len(i) for i in arr[:,0], alt ]
        ax.plot( arr[:,0] , alt , color=color, label=label)
        # limit cb to top of troposphere
        min, max  =  [ (np.ma.min(i), np.ma.max(i)) for i in [ arr[:,0] ] ][0] 
    else:
        ax.plot( arr , alt,  label=label, color=color )
        # limit cb to top of troposphere
        min, max  =  [ (np.ma.min(i), np.ma.max(i)) for i in [ arr ] ][0] 

    # Add legend
    if legend:
        plt.legend(loc=loc)#, fontsize=f_size*.5)

    # sonde data  = mean, 1st and 3rd Q
    if err_bar:     
        if stddev:
            ax.errorbar(arr[:,0] , alt, xerr=arr[:,3], fmt='o', color=color )
        else:
            ax.errorbar(arr[:,0] , alt, xerr=[arr[:,1], arr[:,2]], fmt='o', \
                color=color )

    # Beautify plot ( e.g. add hPa, units,  etc... )
    if ancillary:
        if (title != None):
            t = plt.title(title, fontsize = f_size*1)
            t.set_y(1.09)
            if not isinstance(subtitle, type(None)):
                plt.text(0.5, 1.05, subtitle, ha='center', va='center', \
                          transform=ax.transAxes, fontsize=f_size*.65)
            plt.subplots_adjust(top=0.86) 
                
        plt.ylabel('Altitude / km', fontsize=f_size*.75)
        plt.xlabel( units, fontsize=f_size*.75)
        if xlimit == None:
    #        plt.xlim( min-(min*0.02), max+(max*0.02) )
            plt.xlim( 0.1,  max+(max*0.50) )
        else:
            plt.xlim(0.1, xlimit)
        if hPa_labels:
            ax2 = ax.twinx()
            press =[ myround(i,100) for i in press ]
            ax2.set_yticks( press[::10])
            majorFormatter = FormatStrFormatter('%d')
            ax2.yaxis.set_minor_formatter(majorFormatter)                        
            ax2.set_ylabel( 'Press. / hPa', fontsize = f_size*.75)

# --------------
# 1.12 - plot up monthly from data provided from DB netCDF
# -------------
def obs_month_plot(data, color=None, title=None, rtn_data=False, \
            plt_day=True, debug=False):
    """ Plot up seaonal (monthly ) data. Requires data, and dates in numpy 
        array form. Dates must be as datetime.datetime objects. 
         - REDUNDENT? (see 1.13) """

    # Setup decimal day list
    day_time= [i+float(1/23) for i in range(23) ]
    if debug:
        print data.shape
        print day_time

    # Plot up all data <= this is inefficient. 
	if plt_day:
	    for i in range( len( data[0,:] ) ) :
	    	day_d = data[:,i]  # select day's data and convert to -1
        	plt.plot(day_time , day_d , alpha=0.1, color=color)

    # Return data?
	if rtn_data:
		if debug:
			print 'rtn_data'
#		day_time, 
		data_ = np.ma.mean(data[:,:],axis=1)
		return day_time, data_
    # Plot up daily decimal days
	else :
		if debug:
			print ' not - rtn_data'
		plt.plot( day_time, np.ma.mean(data[:,:],axis=1) , lw=3 , color=color, \
		     label=title)
		return plt

# --------------
# 1.13 - plot up monthly from data provided from DB netCDF
# -------------
def monthly_plot( ax, data, f_size=20, pos=0, posn=1, lw=1,ls='-', color=None, \
                  title=None, subtitle=None, legend=False, xrotation=90, \
                  window=False, label=None, ylabel=None, loc='upper right' ):
    """ Plot up seaonal (monthly ) data. Requires data, and dates in numpy 
        array form. Dates must be as datetime.datetime objects. """
            
    if color == None:
        color = color_list( posn )[ pos ]

    if window:
        f_size=int(f_size/2)
    # Plot
    plt.plot( np.arange(1,len(data)+1), data, color=color, lw=lw, ls=ls, \
        label=label )

    # Beautify
    ax.set_xticklabels( [i.strftime("%b") for i in [datetime.datetime(2009, \
        int(i), 01) for i in np.arange(1,13 ) ] ] )
    plt.xticks(range(1,13),  fontsize=f_size )
    plt.xlim(0.5,12.5)
    plt.xticks(rotation=xrotation)
    if ylabel != None:
        plt.ylabel(  ylabel, fontsize=f_size  )
    plt.yticks( fontsize=f_size )        
    print pos, posn-1
    if not isinstance(title, type(None)):
        t = plt.title(title, fontsize = f_size)
        t.set_y(1.09)
        if not isinstance(subtitle, type(None)):
             plt.text(0.5, 1.05, subtitle, ha='center', va='center', \
                          transform=ax.transAxes, fontsize=f_size*.65)
#    if pos==posn-1:
    if legend:
        plt.legend( loc=loc,  fontsize=int(f_size/1.5) )


# --------------
# 1.14 - plot up monthly timeseries from ...
# -------------
def timeseries_seasonal_plot( ax, dates, data, f_size=20, pos=0, posn=1,  \
            title=None, legend=False, everyother=24,  x_nticks=12, \
            window=False, label=None, ylabel=None, loc='upper right',  \
            lw=1,ls='-', color=None, showmeans=False, debug=False ):
    """ Plot up timeseries of seasonal data. Requires data, and dates in numpy
        array form. Dates must be as datetime.datetime objects. """

    # Process data - reduce resolution to daily, and get std
    df = DataFrame(data,index=dates )
#    resample = df.resample('M' )#, how='mean')
#    labels = [i.strftime("%b") for i in resample]
    months = range(1, 13) 
    labels = [i.strftime("%b") for i in 
    [ datetime.datetime(2009, int(i), 01) for i in months ]
    ]
    print labels

    # Get Data by month
    monthly =[ df[df.index.month==i]  for i in months ]
    print [i.shape for i in monthly ]
    
    bp = ax.boxplot( monthly, months, showmeans=showmeans ) 
    ax.set_xticklabels( labels )

    # Beatify plot
    if not isinstance( title, type(None) ):
        plt.title( title )
    if not isinstance( ylabel, type(None) ):
        plt.ylabel( ylabel )

# --------------
# 1.15 - plot up daily timeseries from ...
# -------------
def timeseries_daily_plot(fig, ax,  dates, data, pos=1, posn =1,                  
                bin_size=1/24.,widths=0.01, rotatecbunit = 'vertical', 
                white_fill=True, alpha=0.1, linewidth=0.5,xlabel=True, 
                title=None, f_size=7.5, units='ppbv', scale='linear', \
                showmeans=False, debug=False ):
    """ Plot up daily timeseries of values. Requires data, and dates in numpy 
        array form. Dates must be as datetime.datetime objects. """

    # get_day_fraction(i) 
    dates = np.array( [ get_day_fraction(i) for i in dates ] )

    # bin data    
    b_all, bins_used = bin_data( data, dates, bin_size, debug=debug )

    # Plot
    bp = ax.boxplot( b_all,  positions=bins_used, widths=widths, \
        showmeans=showmeans, patch_artist=True )

    # Beautify 
    ax.set_xticklabels( np.arange(0,24,1 )  )
    plt.xticks( np.arange(0,1,1/24. ), fontsize=f_size, rotation=rotatecbunit )
#    plt.xlim(-0.05, 23/24.)
    plt.xlim(-0.05, 1.0)

    plt.ylabel(units, fontsize=f_size)
    plt.yticks( fontsize=f_size)
    if (title != None):
        plt.title(title, fontsize=f_size )
    if xlabel:
        plt.xlabel('Hour of day', fontsize=f_size)
    plt.ylabel('{}'.format(units), fontsize=f_size)


    # --- Highlight bins        
    bs = np.arange(0, 24, bin_size )#[ bs[0] - bin_size ] + bs 
    [ plt.axvline( x=i, color='k', linewidth=linewidth, alpha=alpha,  \
                            linestyle='dashed' ) for i in bs ]


# --------------
# 1.16 - plot up timeseries from May through Septmeber
# -------------
def timeseries_May_Sept_plot( ax, dates, data, f_size=20, pos=0, posn=1,  \
            title=None, legend=False, everyother=7,  x_nticks=12, \
            window=False, label=None, ylabel=None, loc='upper right',  \
            lw=1,ls='-', color=None, start_month=5, end_month=9, 
            boxplot=True, showmeans=False, debug=False ):
    """ Timeseries plot. Takes numpy arrays of datetimes and data as
          arguemnts. """

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
        fill_points =range(len(sel_days))[::everyother] 
        for n,i in enumerate( fill_points) :
            fill_label[i] = labels[n]
    
        plt.xticks( range(len(sel_days)), fill_label, \
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


# --------------
# 1.17 - plot up timeseries from July 
# -------------
def timeseries_month_plot( ax, dates, data, f_size=20, pos=0, posn=1,  \
            title=None, legend=False, everyother=7*24,  x_nticks=12, \
            window=False, label=None, ylabel=None, loc='upper right',  \
            lw=1,ls='-', color=None, start_month=7, end_month=7, \
            boxplot=True, showmeans=False, alt_text=None, r_plt=False, \
            unitrotation=45, color_by_z=False, fig=None,  
            positive=None, debug=False ):
    """ Plot up month timeseries of values. Requires data, and dates in numpy 
        array form. Dates must be as datetime.datetime objects. """

    # Process data - reduce resolution to daily, and get std
    df = DataFrame( data, index=dates )

    # remove dates outside of range (start_month > < end_month )
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

    if color_by_z:
        if debug:
            print 'Coloring line by normalised z values'
        print df.columns
        x = df.index
        y, z = [ df[ df.columns[i] ] for i in range(2) ]
        cmap = get_colormap( z.copy(), positive=positive )
        print [ ( i.min(), i.max() ) for i in x, y, z ]
        colorline(x, y, z, cmap=cmap, linewidth=lw, ax=ax, \
            norm=plt.Normalize( 0, 360 ), fig=fig ) #np.min(z), 1500))
#        colorline(x, y, linewidth=lw, ax=ax)

    else:
        plt.plot( days, df.values, label=label, color=color, ls=ls, lw=lw  )

    # set xticks
    plt.xticks( days[::everyother], labels, rotation=unitrotation )

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

#    if r_plt:
#        return plt

# --------
# 1.18 - North Pole surface plot
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
    """ Plot up data at north pole 2D slice.
        Requires data (arr) as numpy array. Arr should be full global size 
        (lon, lat) for given resolution """

    # ---- Grid/Mesh values for Lat, lon, & alt + cb
    if isinstance( cmap, type(None) ):
        cmap = get_colormap( arr.copy() )

    if debug:
        print '>'*5, [ [ i.min(), i.max(), i.mean(), type(i) ] for i in [arr] ]
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
        print 1, len(x), len(y)


#    debug=True
                    
    if debug:
        print 2, 'len:',  [ len(i) for i in x,y,lat,lon ]
        print '>'*5, [ [ i.min(), i.mean(), i.max(), i.shape ] for i in arr, 
            arr[:,get_gc_lat(boundinglat,res=res):]  ]

    # -------- colorbar variables...
    # set cb label sizes
    if not no_cb:
        if fixcb:
            tickmin, tickmax = fixcb[0], fixcb[1]
        else:
            tickmin, tickmax  = arr.min(), arr.max()

    # -----------------------  Linear plots -------------------------------
#    plt.contourf( x, y, arr[:,get_gc_lat(boundinglat,res=res):].T, alpha=alpha)
    if fixcb:
        poly = m.pcolor( x, y, arr[:,get_gc_lat(boundinglat, \
            res=res):].T, cmap=cmap, alpha=alpha,  vmin=fixcb[0], \
            vmax=fixcb[1] )
    else:
        poly = m.pcolor( x, y, arr[:,get_gc_lat(boundinglat, \
            res=res):].T, cmap=cmap, alpha=alpha )

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
            tick_locator = ticker.MaxNLocator(nticks=nticks)
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
# 1.19- South Pole surface plot
# --------
def south_pole_surface_plot( arr, return_m=False, grid=True, centre=False, 
            cmap=None, format='%.2f', res='4x5', ax=None, alpha=1, 
             fixcb=False, nbins=25, nticks=10, 
             drawcountries=True, set_cb_ticks=True, title=None, 
#             rotatecbunits='horizontal', extend='neither', 
             interval=1, resolution='l', shrink=0.4, window=False, everyother=1, 
            lon_0=0, boundinglat=-40,
             degrade_resolution=False,              
             no_cb=False, cb=None, 
             units=None, f_size=20, 
             debug=False,  **Kwargs):    
    """ Plot up data at south pole 2D slice.
        Requires data (arr) as numpy array. Arr should be full global size 
        (lon, lat) for given resolution """

    # ---- Grid/Mesh values for Lat, lon, & alt + cb
    if isinstance( cmap, type(None) ):
        cmap = get_colormap( arr.copy() )

    if debug:
        print '>'*5, [ [ i.min(), i.max(), i.mean(), type(i) ] for i in [arr] ]
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
        print 1, len(x), len(y)

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
        print 2, 'len:',  [ len(i) for i in x,y,lat,lon ]
        print '>'*5, [ [ i.min(), i.mean(), i.max() ] for i in [arr ] ]

    # -------- colorbar variables...
    # set cb label sizes
    if fixcb:
        tickmin, tickmax = fixcb[0], fixcb[1]
    else:
        tickmin, tickmax  = arr.min(), arr.max()

    # --------------------  Linear plots -------------------------------
#    plt.contourf( x, y, arr[:,get_gc_lat(boundinglat,res=res):].T, alpha=alpha)
    if fixcb:
        poly = m.pcolor( x, y, arr[:,:get_gc_lat(boundinglat, \
            res=res)+2].T, cmap=cmap, alpha=alpha, vmin=fixcb[0],\
             vmax=fixcb[1] )
    else:
        poly = m.pcolor( x, y, arr[:,:get_gc_lat(boundinglat, 
            res=res)+2].T, cmap=cmap, alpha=alpha )

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
            tick_locator = ticker.MaxNLocator(nticks=nticks)
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
# 1.20 - PDF of monthly surface change plots for given species (takes 5D arr )
# --------
def plot_specs_surface_change_monthly2pdf( arr, res='4x5', dpi=160, \
        no_dstr=True, f_size=20, pcent=True, specs=None, dlist=None, \
        savetitle='', diff=False, extend='neither', column=False, \
        scale=1, units=None, set_window=False, lat_0=None, lat_1=None, \
        mask_invalids=False, debug=False):   
        """ Create multipage PDF with each page containing a 2D (lon,lat) slice     
             plot for given species in list of "specs" 
             Takes 5D array  ( species ,lon , lat, alt, time) """
        print 'plot_specs_surface_change_monthly2pdf called'

        # setup pdfs + titles
        if column:
            savetitle = 'Column_by_spec'+savetitle 
        else:
            savetitle = 'Surface_by_spec'+savetitle 
        pdff = plot2pdfmulti( title=savetitle, open=True, \
                              dpi=dpi, no_dstr=no_dstr )    
        left=0.05; right=0.9; bottom=0.05; top=0.875; hspace=0.315; wspace=0.1

        # Loop species
        for n, spec in enumerate( specs ):
            if debug:
                print n, spec

            # Get units/scale for species + setup fig (allow of 
            if isinstance( units, type( None) ):
                if column and (not pcent):
                    units, scale = 'DU', 1        
                elif pcent:
                    units, scale = '%', 1
                else:
                    units, scale = tra_unit( spec, scale=True, \
                        global_unit=True )
            
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

            fig  = plt.figure(figsize=(22, 14), dpi=dpi, 
                facecolor='w', edgecolor='k')

            # set cb ranges for whiole data period
            fixcb  = [( i.min(), i.max() ) for i in [cbarr]][0]

            # Kludge, force max cap at .2
#            if units == 'ratio':
#                fixcb = [ fixcb[0], 0.2 ]
#                extend = 'max'
            if (units == 'pmol / m$^3$') and ( spec == 'AERI' ):
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
                    mask_invalids=mask_invalids, clevs=clevs,\
                    debug=debug)
                plt.title(month.strftime("%b"), fontsize=f_size*2)

            # Add single colorbar
            mk_cb(fig, units=units, left=0.9,  cmap=cmap,vmin=fixcb[0],
                 vmax=fixcb[1], f_size=f_size, extend=extend ) 

            # sort out ascetics -  adjust plots and add title
            fig.subplots_adjust( bottom=bottom, top=top, left=left, \
                right=right,hspace=hspace, wspace=wspace)
            fig.suptitle( ptitle, fontsize=f_size*2, x=.55 , y=.95  )

            # save out figure
            plot2pdfmulti( pdff, savetitle, dpi=dpi,\
                no_dstr=no_dstr )

            # close fig
            plt.clf()
            plt.close()
            del fig

        #  save entire pdf  
        plot2pdfmulti( pdff, savetitle, close=True, dpi=dpi,\
                       no_dstr=no_dstr )
        
# --------
# 1.21 - PDF of monthly zonal change plots for given species
# --------
def plot_specs_zonal_change_monthly2pdf( Vars, res='4x5', dpi=160, \
        no_dstr=True, f_size=20, pcent=False, specs=None, dlist=None, \
        t_ps=None, savetitle='', diff=False, extend='neither',  \
        set_window=False, lat_0=None, lat_1=None, mask_invalids=False, \
        set_lon=None, units=None, debug=False):
        """ Create multipage PDF with each page containing a zonal      
             plot for given species in list of "specs" 
                          Takes 5D array  ( species ,lon , lat, alt, time) """

        savetitle = 'Zonal_by_spec'+savetitle 
        pdff = plot2pdfmulti( title=savetitle, open=True, \
                              dpi=dpi, no_dstr=no_dstr )    
        left=0.05; right=0.9; bottom=0.05; top=0.875; hspace=0.315; wspace=0.2
    
        # Loop species
        for n, spec in enumerate( specs ):
            if debug:
                print n, spec, Vars.shape

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

            fig  = plt.figure(figsize=(22, 14), dpi=dpi, 
                facecolor='w', edgecolor='k')

            cbVars = Vars[n,:,:,:,:].copy()*scale

            # set ranges for whiole data period
            if pcent:
                if len( cbVars[ cbVars>500 ] ) > 0:
                    cbVars =  np.ma.masked_where( cbVars>500, \
                        cbVars )
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
                fixcb  = [( i.min(), i.max() ) \
                    for i in [ cbVars[set_lon,...] ]][0]               
                print 'SETTING LON to GC index: ', set_lon
            else:
                cbVars = cbVars.mean( axis=0 )
            if set_window:
                gclat_0, gclat_1 = [ get_gc_lat(i, res=res ) \
                    for i in lat_0, lat_1 ]
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
                    arr = arr[...,get_gc_lat(lat_0, res=res): \
                        get_gc_lat(lat_1, res=res), :]

                # plot up spatial surface change
                zonal_plot(fig, ax, arr, \
                    title=month.strftime("%b"), debug=debug, tropics=False, \
                    units=units,f_size=f_size, c_off =37, no_cb=True, \
                    lat_0=lat_0, lat_1=lat_1, set_window=set_window,\
                    fixcb=fixcb, extend=extend, window=True, \
                     lower_limited=True, res=res, mask_invalids=mask_invalids, \
                    cmap=cmap )
                
                # only show troposphere
                greyoutstrat( fig, t_ps.mean(axis=0)[:,:,m], axn=axn, \
                    res=res )

            # Add single colorbar
            mk_cb(fig, units=units, left=0.915,  cmap=cmap,vmin=fixcb[0],
                 vmax=fixcb[1], f_size=f_size, extend=extend ) 

            # sort out ascetics -  adjust plots and add title
            fig.subplots_adjust( bottom=bottom, top=top, left=left, \
                right=right,hspace=hspace, wspace=wspace)
            fig.suptitle( ptitle, fontsize=f_size*2, x=.55 , y=.95  )

            # save out figure
            plot2pdfmulti( pdff, savetitle, dpi=dpi,\
                no_dstr=no_dstr )

            # close fig
            plt.clf()
            plt.close()
            del fig

        #  save entire pdf 
        plot2pdfmulti( pdff, savetitle, close=True, dpi=dpi,\
                       no_dstr=no_dstr )

# --------
# 1.22 - 
# --------

# --------
# 1.23 - Change as 2D plot of surface ( for column or surface change )
# --------
def plot_specs_poles_change_monthly2pdf(  specs=None,\
            arr=None, res='4x5', dpi=160, no_dstr=True, f_size=20, pcent=False,\
            diff=False, dlist=None, savetitle='', units=None, \
            perspective='north', format=None,\
            extend='neither', boundinglat=50, debug=False):
        """ Takes a 5D np.array ( species, lon, lat, alt, time ) and plots up
             the output by spcies  by month , and saves this as a mulitpage pdf
                          Takes 5D array  ( species ,lon , lat, alt, time) """

        if debug:
            print arr, no_dstr, f_size, \
                pcent, res, dpi, specs, \
                dlist, savetitle, debug

        if perspective == 'north':
            savetitle = 'North_'+savetitle 
        if perspective == 'south':
            savetitle = 'South_'+savetitle 
        
        pdff = plot2pdfmulti( title=savetitle, open=True, \
                              dpi=dpi, no_dstr=no_dstr )    
        left=0.01; right=0.925; bottom=0.025; top=0.95; hspace=0.15; wspace=-0.1
    
        # Loop species
        for n, spec in enumerate( specs ):

            print n, spec, arr.shape, units, perspective, '<'

            # Get units/scale for species + setup fig
            scale = 1
            if pcent :# and (not units == 'DU') :
                units = '%' 
            if isinstance( units, type( None )):
                units, scale = tra_unit( spec, scale=True, global_unit=True )

            parr = arr[n,:,:,0,:]*scale
            
            if debug:
                print parr.shape
                print n, spec, [ (i.min(), i.max(), i.mean() ) \
                    for i in [parr] ], [ (i.min(), i.max(), i.mean() ) \
                    for i in [arr] ], units, scale

            # Set the correct title        
            ptitle = '{}'.format( latex_spec_name(spec) )
 
            # create new figure            
            fig  = plt.figure(figsize=(22, 14), dpi=dpi, 
                facecolor='w', edgecolor='k')

            # select north or south polar areas specified to define cb
            if perspective == 'north':
                cbarr=parr[:,get_gc_lat(boundinglat,res=res):,:].copy()
            if perspective == 'south':
                cbarr=parr[:,:get_gc_lat(-boundinglat,res=res),:].copy()

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
                    
            fixcb = [ ( i.min(), i.max() ) for i in [cbarr] ][0]
            print 'fixcb testing ', fixcb, parr.shape

            # Kludge, force max cap at .2
#            if units == 'ratio':
#                fixcb = [ fixcb[0], 0.2 ]
#                extend = 'max'

            cmap = get_colormap( fixcb )

            # Loop thorugh months
            for m, month in enumerate(dlist):
                ax = fig.add_subplot(3,4,m+1)

                # plot up spatial surface change
                if perspective == 'north':
                    north_pole_surface_plot( parr[:,:,m], no_cb=True, \
                        fixcb=fixcb, diff=diff, pcent=pcent, res=res,          \
                        f_size=f_size*2,
                        cmap=cmap, boundinglat=boundinglat) 
                if perspective == 'south':
                    south_pole_surface_plot( parr[:,:,m], no_cb=True, \
                        fixcb=fixcb, diff=diff, pcent=pcent, res=res,         \
                        f_size=f_size*2, 
                        cmap=cmap, boundinglat=-boundinglat) 

                plt.text( 1.1, -0.01, month.strftime("%b"), 
                                transform=ax.transAxes,
                                ha='right', va='bottom', 
                                fontsize=f_size*2)

            # Add single colorbar
            mk_cb(fig, units=units, left=0.915,  cmap=cmap,vmin=fixcb[0],\
                 vmax=fixcb[1], f_size=f_size*.75, extend=extend, \
                 format=format ) 

            # sort out ascetics -  adjust plots and add title
            fig.subplots_adjust( bottom=bottom, top=top, left=left, \
                right=right,hspace=hspace, wspace=wspace)
            fig.suptitle( ptitle, fontsize=f_size*2, x=.475 , y=.975  )

            # save out figure
            plot2pdfmulti( pdff, savetitle, dpi=dpi,\
                no_dstr=no_dstr )

            # close fig
            plt.clf()
            plt.close()
            del parr

        #  save entire pdf
        plot2pdfmulti( pdff, savetitle, close=True, dpi=dpi,\
                       no_dstr=no_dstr )

# --------
# 1.24 - Coupler to select plot type 
# --------
def create_plot4case( fig, ax, dates, data, spec, f_size=20, lw=None, ls=None, \
            loc=True, legend=True, units='', boxplot=False, lname='Weyborne', \
            run_name='run', start_month=7, end_month=7, plot_type='daily', \
            case=None, l_dict=None, label=None, alt_text=None, color=None): 
    """ Coupler for timeseries data (dates as datetimes + data) to be 
            converted multiple graph forms with just change of case    """

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
        timeseries_May_Sept_plot( ax, dates, data,  \
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
# 1.25 - Probability distribution plotter
# --------
def PDF_obs_vs_mod( ax, dates, data, f_size=20, pos=0, posn=1,  \
            title=None, legend=False, everyother=7*24,  x_nticks=12, \
            window=False, label=None, ylabel=None, loc='upper right',  \
            lw=1,ls='-', color=None, start_month=1, end_month=12, 
            unitrotation=45, alpha=.5, bins=100, alt_text=None, debug=False ):
    """Constructs a PDF plot with given dates and data.
        If not provided with axis etc, then these are made.
        data and dates need to be as a np.array, 
        with dates as datetime.datetim objects
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


# --------
# 1.26 - Generic X vs. Y plot
# --------
def X_Y_scatter( x, y, z=None, fig=None, ax=None, vmin=None, vmax=None, \
            left= 0.1, width=0.60, bottom=0.1, height=0.60, widthII=0.2,  \
            lim2std=10,trendline=True, f_size=20, lw=10, title=None, \
            line121=True, Trend_line=True, X_title=None, Y_title=None ):
    """ Plot up a X Y scatter plot of x vs. y """
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
        stds = [ np.std( i ) for i in x, y ]
        means = [ np.mean( i ) for i in x, y ]
        print stds, means
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
# 1.27 - Scatter 3D cube
# --------
def scatter_3D_cube( data, dims=None, res='2x2.5', fig=None, 
        everyother=1, interval=1, f_size=20, debug=False ):

    if debug:
        print data.shape

    if isinstance( dims, type(None) ):
        lon, lat, alt = get_latlonalt4res( res=res ) 
        (X,Y,Z) = np.mgrid[ \
            lon[0]:lon[-1]:complex( len(lon) ), 
            lat[0]:lat[-1]:complex( len(lat) ), 
            alt[0]:alt[-1]:complex( len(alt) ) ]

    # Setup Fig + Ax ( is not given )
    if isinstance( fig, type(None) ):
        fig = plt.figure(1)
    fig.clf()
    ax = Axes3D(fig)
    ax.scatter(X,Y,Z, c=data)

    # --- Beatify
    # draw meridian lines
    meridians = np.arange(-180,180,20*interval)
    plt.xticks( meridians[::everyother], fontsize = f_size*.75 ) 

    # draw parrelel lines
    parallels = np.arange(-90,91,20*interval)
    plt.yticks( parallels[::everyother], fontsize = f_size*.75 ) 

    # limit axis
#    plt.xlim( lon[0], lon[-1] )
#    plt.ylim( lat[0], lat[-1] )    
#    plt.zlim( alt[0], alt[-1] )        

# --------
# 1.28 - Plot up seasonal output from 4D arr (lon, lat, alt, time)
# --------
def get_seasonal_plot( arr, fixcb=None, fig=None, f_size=15, \
            case='linear', format=None, extend='neither', units=None, \
           right=0.9, left=0.05, bottom=0.05, top=0.85, hspace=0.1, \
           wspace=0.1, log=False, title=None, dpi=80, debug=False ):
    """ Takes any 4D array and plot a 4 subplot window plot by season"""

    # Split by quater (DJF, MAM, JJA, SON)
    ars, seasons = split_4D_array_into_seasons( arr, annual_plus_seasons=False )
    
    # create figure
    if isinstance( fig, type(None ) ):
        fig  = plt.figure(figsize=(22, 14), dpi=dpi, facecolor='w', \
                edgecolor='k')     

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

# -------------
# 1.29 - box X-Y plot with histograms
# -------------
def X_Y_hist( x, y, z=None, zlabel=None, fig=None, \
            left= 0.1, width=0.60, bottom=0.1, height=0.60, widthII=0.2, \
            fit = 0.02, binwidth=0.1, Trend_line = False, cmap=cm.gnuplot2, \
            X_title='X title'  , Y_title='Y title' , f_size=5 , line121=False ):
    """ Plots a X vs. Y histogram """

    # Use hist code 
    nullfmt   = NullFormatter()         # no labels

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
        print lim_min, lim_max
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
    [ i.set_yscale('linear') for i in axHisty, axScatter ]
    [ i.set_xscale('linear') for i in axHistx, axScatter ]

    plt.rcParams.update({'font.size': f_size})
#    [ i.rcParams.update({'font.size': f_size}) for i in axHisty, axHistx, axScatter ]
#    [ i.xticks( fontsize=f_size ) for i in axHisty, axHistx, axScatter ]
#    [ i.yticks( fontsize=f_size ) for i in axHisty, axHistx, axScatter ]


# --------
# 1.30 - PDF of annual surface change plots for given species (takes 5D arr )
# --------
def plot_specs_surface_change_annual2pdf( arr, res='4x5', dpi=160, \
        no_dstr=True, f_size=20, pcent=True, specs=None, dlist=None, \
        savetitle='', diff=False, extend='neither', column=False, \
        scale=1, units=None, set_window=False, lat_0=None, lat_1=None, \
        mask_invalids=False, debug=False):   
        """ Create multipage PDF with each page containing a 2D (lon,lat) slice     
             plot for given species in list of "specs" 
             Takes 5D array  ( species ,lon , lat, alt, time) """
        print 'plot_specs_surface_change_monthly2pdf called'

        # setup pdfs + titles
        if column:
            savetitle = 'Column_by_spec'+savetitle 
        else:
            savetitle = 'Surface_by_spec'+savetitle 
        pdff = plot2pdfmulti( title=savetitle, open=True, \
                              dpi=dpi, no_dstr=no_dstr )    
        left=0.05; right=0.9; bottom=0.05; top=0.875; hspace=0.315; wspace=0.1

        # Loop species
        for n, spec in enumerate( specs ):
#            if debug:
            print n, spec

            # Get units/scale for species + setup fig (allow of 
            if isinstance( units, type( None) ):
                if column and (not pcent):
                    units, scale = 'DU', 1        
                elif pcent:
                    units, scale = '%', 1
                else:
                    units, scale = tra_unit( spec, scale=True, \
                        global_unit=True )
            
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

            fig  = plt.figure(figsize=(22, 14), dpi=dpi, 
                facecolor='w', edgecolor='k')

            # set cb ranges for whiole data period
            fixcb  = [( i.min(), i.max() ) for i in [cbarr]][0]

            # Kludge, force max cap at .2
#            if units == 'ratio':
#                fixcb = [ fixcb[0], 0.2 ]
#                extend = 'max'
            if (units == 'pmol / m$^3$') and ( spec == 'AERI' ):
                fixcb = [ fixcb[0], 500 ]
                extend = 'max'                

            cmap = get_colormap( fixcb  )
                
            # Loop thorugh months
#            for m, month in enumerate(dlist):
            fig.add_subplot(111)

            # plot up spatial surface change
            map_plot( arr[n,:,:,0,:].mean(axis=-1).T*scale, cmap=cmap, case=9, \
                res=res,\
                no_cb=True, f_size=f_size, fixcb=fixcb, window=True,\
                set_window=set_window,lat_0=lat_0, lat_1=lat_1, \
                mask_invalids=mask_invalids, debug=debug)

            # Add single colorbar
            mk_cb(fig, units=units, left=0.905,  cmap=cmap,vmin=fixcb[0],
                 vmax=fixcb[1], f_size=f_size, extend=extend ) 

            # sort out ascetics -  adjust plots and add title
            fig.suptitle( ptitle, fontsize=f_size*2, x=.55 , y=.95  )

            # save out figure
            plot2pdfmulti( pdff, savetitle, dpi=dpi,\
                no_dstr=no_dstr )

            # close fig
            plt.clf()
            plt.close()
            del fig

        #  save entire pdf  
        plot2pdfmulti( pdff, savetitle, close=True, dpi=dpi,\
                       no_dstr=no_dstr )


# --------
# 1.31 - PDF of annual zonal change plots for given species
# --------
def plot_specs_zonal_change_annual2pdf( Vars, res='4x5', dpi=160, \
        no_dstr=True, f_size=20, pcent=False, specs=None, dlist=None, \
        t_ps=None, savetitle='', diff=False, extend='neither',  \
        set_window=False, lat_0=None, lat_1=None, mask_invalids=False, \
        set_lon=None, units=None, debug=False):
        """ Create multipage PDF with each page containing a zonal      
             plot for given species in list of "specs" 
                          Takes 5D array  ( species ,lon , lat, alt, time) """

        savetitle = 'Annual_Zonal_by_spec'+savetitle 
        pdff = plot2pdfmulti( title=savetitle, open=True, \
                              dpi=dpi, no_dstr=no_dstr )    
        left=0.05; right=0.9; bottom=0.05; top=0.875; hspace=0.315; wspace=0.2
    
        # Loop species
        for n, spec in enumerate( specs ):
#            if debug:
            print n, spec, Vars.shape

            # Get units/scale for species + setup fig        
#            scale=1
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

            fig  = plt.figure(figsize=(22, 14), dpi=dpi, 
                facecolor='w', edgecolor='k')

            cbVars = Vars[n,:,:,:,:].mean(axis=-1).copy()*scale

            # set ranges for whiole data period
            if pcent:
                if len( cbVars[ cbVars>500 ] ) > 0:
                    cbVars =  np.ma.masked_where( cbVars>500, \
                        cbVars )
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
                fixcb  = [( i.min(), i.max() ) \
                    for i in [ cbVars[set_lon,...] ]][0]               
                print 'SETTING LON to GC index: ', set_lon
            else:
                cbVars = cbVars.mean( axis=0 )
            if set_window:
                gclat_0, gclat_1 = [ get_gc_lat(i, res=res ) \
                    for i in lat_0, lat_1 ]
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
                arr = arr[...,get_gc_lat(lat_0, res=res): \
                    get_gc_lat(lat_1, res=res), :]

            # plot up spatial surface change
            zonal_plot( arr, fig, ax=ax, \
                title=None, debug=debug, tropics=False, \
                units=units,f_size=f_size, c_off =37, no_cb=True, \
                lat_0=lat_0, lat_1=lat_1, set_window=set_window,\
                fixcb=fixcb, extend=extend, window=True, \
                    lower_limited=True, res=res, mask_invalids=mask_invalids, \
                    cmap=cmap )
                
            # only show troposphere
            greyoutstrat( fig, t_ps.mean(axis=0).mean(axis=-1), axn=axn, \
                    res=res )

            # Add single colorbar
            mk_cb(fig, units=units, left=0.915,  cmap=cmap,vmin=fixcb[0],
                 vmax=fixcb[1], f_size=f_size, extend=extend ) 

            # sort out ascetics -  adjust plots and add title
#            fig.subplots_adjust( bottom=bottom, top=top, left=left, \
#                right=right,hspace=hspace, wspace=wspace)
            fig.suptitle( ptitle, fontsize=f_size*2, x=.55 , y=.95  )

            # save out figure
            plot2pdfmulti( pdff, savetitle, dpi=dpi,\
                no_dstr=no_dstr )

            # close fig
            plt.clf()
            plt.close()
            del fig

        #  save entire pdf 
        plot2pdfmulti( pdff, savetitle, close=True, dpi=dpi,\
                       no_dstr=no_dstr )

# --------------------------- Section 4 ------------------------------------
# -------------- Plotting Ancillaries 
#

# --------------                                                                                                                                             
# 4.01 -  Get percentiles
# ------------- 
## {{{ http://code.activestate.com/recipes/511478/ (r1)
#import math
#import functools
def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
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

    # median is 50th percentile.
    #median = functools.partial(percentile, percent=0.5)
    ## end of http://code.activestate.com/recipes/511478/ }}}


# --------------                                                                                                                                             
# 4.02 - moving average
# ------------- 
def moving_average(x, n, type='simple'):
    """   compute an n period moving average.
            type is 'simple' | 'exponential' """ 
    x = np.asarray(x)
    if type=='simple':
        weights = np.ones(n)
    else:
        weights = np.exp(np.linspace(-1., 0., n))

    weights /= weights.sum()

    a =  np.convolve(x, weights, mode='full')[:len(x)]
    a[:n] = a[n]
    return a

# -------------
# 4.03 -  color list for rainbow plots
# -------------
def color_list(length, cb='gist_rainbow' ):
    """ Create a list of colours to generate colors for plots contain multiple 
          datasets """
    cm = plt.get_cmap(cb)
    color = [cm(1.*i/length) for i in range(length)]
    return color

# -------------
# 4.04 - weighted average
# -------------
def weighted_average( data, interval , bins_used=False, debug=False):
    """ Calculate a weighed average of data fro a given interval """

    min_v, max_v = int( data.min() ), int( data.max() ) 
    bins = np.arange(min_v, max_v, interval)
    if debug:
        print [ (i, type(i) ) for i in [ min_v, max_v, bins ] ]
    int_d = np.int_(data)
    binned, bins_used  = [], []
    for bin_ in bins:
        b_mean = np.mean( data[np.where( int_d == bin_ )] )
        if ( b_mean != 0 ):            
            binned.append( b_mean )
            bins_used.append( bin_ ) 
        else: 
            print 'no data for bin {}'.format(bin_)
            pass
    if debug:
        print [ ( i, len(i) ) for i in [binned , bins] ]
    if (bins_used):
        return binned, bins_used
    else:
        return binned

# -------------
# 4.05 - R squared - credit: Software carpentry
# -------------
#def r_squared(actual, ideal):
#    actual_mean = np.mean(actual)
#    ideal_dev = np.sum([(val - actual_mean)**2 for val in ideal])
#    actual_dev = np.sum([(val - actual_mean)**2 for val in actual])
#    return ideal_dev / actual_dev
def r_squared(x, y): 
    """ Return R^2 where x and y are array-like."""
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return r_value**2

# -------------
# 4.06 - setup box plots
# -------------
def set_bp( bp, num, c_list=['k', 'red'], white_fill=True, set_all=True,
                    median_color='white', linewidth=5, debug=False ):
    """ Manual set properties of boxplot ("bp") """
    if debug:
        print num, c_list
    if set_all:
        setp(bp['boxes'][:], color=c_list[num])
        setp(bp['caps'][:], color=c_list[num])
        setp(bp['whiskers'][:], color=c_list[num])
        setp(bp['fliers'][:], color=c_list[num])
        setp(bp['medians'][:], color=c_list[num])
        if white_fill:
            [ box.set( facecolor = 'white') for box in bp['boxes'] ]
        else:
            [ box.set( facecolor = c_list[num] ) for box in bp['boxes'] ]
            setp(bp['medians'][:], color=median_color, linewidth=linewidth)                        

# -------------
# 4.07 - Get all marker types
# -------------
def markers_list(rm_plain_markers=False):
    """ Create a list of available markers for use in plots """
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
# 4.08- box X-Y plot with histograms
# -------------
def Trendline( ax, X, Y, order =1, intervals= 700, f_size=20, 
            color='blue', lw=1 ):
        """ Add a trend line to existing plot """
        params, xp = np.polyfit( X, Y, order  ), \
            np.linspace( min(np.ma.min(X), np.ma.min(Y) ), \
            max(np.ma.max(X), np.ma.max(Y) ), intervals )

        yp = np.polyval( params, xp )

        
        print params
#        sys.exit( )

        # Regression 
        r_sq =  [ r_squared(X, i ) for i in [Y]]

        # Plot up
        ax.plot( xp, yp, ls='--', lw=lw, color=color, 
            label = '{0} (R^2 = {1:<,.3f}, m={2:,.3f}, c={3:.3f})'.format(\
             '', r_sq[0], params[0], params[1] )   )

        ax.legend()

# -------------
# 4.09 - plot_gc_bin_bands - plot up Fast-J  bins  ( 7 longest nm bins)
# -------------   
def plot_gc_bin_bands(facecolor='#B0C4DE'):
    """ Plot highlighted lines of GEOS-Chem/Fast-J photolysis bins """
    vars = [ GC_var(i) for i in 'FastJ_lower' , 'FastJ_upper']
    alphas = [ 0.3,0.1 ]*len(vars[0])
    [ plt.axvspan( vars[0][n], vars[1][n], facecolor=facecolor, \
        alpha=alphas[n]) for n in range(len(vars[0])) ]

# --------------
# 4.10 - line styles 
# -------------       
def get_ls(num):
    """ Get a list of available line styles """
    ls =[   ':', '--', '-.', '-', ':', '--', '-.', '-', ':', ':', '--', '-.', '-', ':', '--', '-.', '-', ':'
    ]
    return ls[:num]
                        
# --------   
# 4.21 - takes time in troposphere diagnostic array (46, 47) overlayes 
# --------
def greyoutstrat( fig,  arr, axn=[1,1,1], ax=None, cmap=cm.bone_r, \
            res='4x5', rasterized=True, debug=False):
    """ Grey out stratosphere in existing zonal plot. This is used to highlight 
        tropospherically focused work. This function requires the array of "time 
        in the troposphere" diagnostic (lon,lat, alt) """
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
# 4.22 - adjust subplots
# --------
def adjust_subplots( fig ):
    """ Set subplot adjust in provide figure """

    left  = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9  # the top of the subplots of the figure
    wspace = 0.2 # the amount of width reserved for blank space between subplots
    hspace = 0.5 # the amount of height reserved for white space between subplots

    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, \
        wspace=None, hspace=None)


# ----
# 4.24 - setup diunal
# ----
def setup_diurnal(years, months, f_size=20):
    """ Setup figure for diurnal plot """

    plt.title('Monthly Diurnal Obs (2008-2011) vs. Model ({}{}-{}{}) '.format( \
        years[0], months[0],years[1],months[1]   ) )
    plt.legend(loc='lower left')
    plt.ylabel('Delta O3 (UT 17:00-09:00) /  p.p.b.v')
    plt.xlabel('Hours')
#    plt.ylim(-5, 1)
    plt.ylim(-3, 1)
    plt.xlim(4, 21)
    plt.rcParams.update({'font.size': f_size})

# -------------
# 4.27 -  print NCAS & York logos in the bottom corners
# -------------
def add_logos_NCAS_york_bottom(fig):
    """ Add NCAS + York logo for external plots used externally """
    
    wd1 = get_dir ('dwd') +'misc/logos/'
    logo_list= ['NCAS_national_centre_logo.gif', 'nerclogo1000.gif' , 'york_uni_shield.tif' ,   \
                            'york_uni_name.tif' ]

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

# --------
# 4.28 - Iodine deposition mask (for sites... e.g Denmark, Germany, Norfolk )
# --------
def mask_not_obs( loc='Denmark', res='4x5', debug=False ):
    """ provide a mask of all regions apart from the location given """

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
# 4.30 - Annotate grid 
# -------------
def annotate_gc_grid(ax, res='4x5', f_size=6.5, loc_list=[ [-9999,-9999] ],
                                        everyother=1 , label_gc_grid=True):
    """ Annotate grid with GEOS-Chem indices """

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
# 4.31 - Function for centering colorbar
# --------
def shiftedColorMap(cmap, start=0, midpoint=0.5, lower=0, \
            upper=1,start_center=0.5, stop=1.0, maintain_scaling=True,\
            arr=None, name='shiftedcmap', npoints=257):
    '''
    adapted from stackoverflow to allow for maintaining scaling: 
    original: "http://stackoverflow.com/questions/7404116/defining-the-  \
    midpoint-of-a-colormap-in-matplotlib "
    
    ORGINAL DESCRIPTION:
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    if maintain_scaling:
        # cut off extra colorbar above point to shift maxoum too... 
        vmin, vmax= arr.min(), arr.max()

        # symetric range
        eq_range = max(  abs(vmin), abs(vmax) ) *2 
        # actual range
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

            print 'maintain scale: ', [float(i) for i in list( [ eq_range,\
                act_range, vmin, vmax, abs(vmin), abs(vmax), cut_off, \
                start_center, lower, upper, midpoint, start, stop] )  ]
    cdict = { 'red': [], 'green': [], 'blue': [], 'alpha': [] }

    # regular index to compute the colors
    reg_index = np.hstack([ np.linspace(start, start_center, 128, \
        endpoint=False), np.linspace(start_center, stop, 129, endpoint=True)  ])
    
    # shifted index to match the data
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
# 4.35 - Add colorbar to side of plot
# --------
def mk_cb(fig, units=None, left=0.925, bottom=0.2, width=0.015, height=0.6,\
    orientation='vertical', f_size=20, rotatecbunits='vertical', nticks=10, \
    extend='neither', norm=None, log=False,  format=None, cmap=None,\
    vmin=0, vmax=10, cb_ax=None, ticklocation='auto', extendfrac=.075, \
    sigfig_rounding_on_cb=2, lvls=None, discrete_cmap=False, debug=False):
    """ Create Colorbar. This allows for avoidance of basemap's issues with 
        spacing definitions when conbining with colorbar objects within a plot 
    """
    # Get colormap (by feeding get_colormap array of min and max )
    if isinstance( cmap, type(None) ):
        cmap = get_colormap( arr=np.array( [vmin,vmax]   ) )

    if debug:
        print left, bottom, width, height
        print vmin, vmax, orientation, ticklocation, extend, extendfrac

    # Make new axis
    if isinstance( cb_ax, type(None) ):
        cb_ax = fig.add_axes([left, bottom, width, height])

    print '>'*5, left, bottom, width, height

    # Setup normalisation for colorbar
    if isinstance( norm, type(None) ):
        if log:
            # create log with 1000 p
            lvls = np.logspace( np.ma.log(vmin), np.ma.log(vmax), 1E5, 
                endpoint=True )

            # select relevant indices ( linearly spaced integers )
            print nticks
            lvls =  [ lvls[0] ] + list( lvls[ gen_log_space( 1E5,  \
                nticks*15) ] ) + [ lvls[-1] ]
            lvls= np.array( lvls )

        else:
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    if isinstance( lvls, type(None) ):            
        # make graduations in colourbar to be human readable
        lvls = get_human_readable_gradations( vmax=vmax,  \
            vmin=vmin, nticks=nticks, \
            sigfig_rounding_on_cb=sigfig_rounding_on_cb  )
        print lvls

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
            # Kludge, set to if overly masked. 
            vmax, vmin = 500,-500
            extend = 'both'

    # Make cb with given details 
    if (extend == 'neither') or (log==True):
        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, format=format,\
                norm=norm, ticks=lvls, extend=extend, \
                orientation=orientation, ticklocation=ticklocation)
    else:
        if debug:
            print 'WARNING: adding extensions to colorbar'
        if log==True:
            print 'Will not work as colrbar adjustment not configured'+\
                ' for log scales'
            sys.exit(0)
        else:
            print lvls
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

        print lvls, norm, boundaries

        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, format=format,\
                norm=norm, ticks=lvls, extend=extend, extendfrac=extendfrac,\
                boundaries=boundaries, \
                spacing='proportional',
#                spacing='uniform',
                orientation=orientation, ticklocation=ticklocation)
    
    # set cb label sizes
    if units != None:
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(f_size)
        if rotatecbunits == 'vertical':
            cb.ax.set_ylabel(units, rotation=rotatecbunits, labelpad=f_size, \
                fontsize=f_size)                      
        else:
            cb.set_label( units, fontsize=f_size )
        cb.ax.tick_params(labelsize=f_size) 


# --------
# 4.36 - Create base map for plotting
# --------
def get_basemap( lat, lon, resolution='l', projection='cyl', res='4x5',\
            everyother=1, f_size=10, interval=1, label_y=False, \
            show_grid=True, drawcountries=False ):
    """ Creates a basemap object. 

        This should be used for first slide in python animated 
        videos to save computation """

    m = Basemap(projection=projection,llcrnrlat=lat[0],urcrnrlat=lat[-1],\
                        llcrnrlon=lon[0],\
                        urcrnrlon=lon[-1],\
                        resolution=resolution  )

    if label_y:
        plt.ylabel('Latitude', fontsize = f_size*.75)
        plt.xlabel('Longitude',fontsize = f_size*.75)
    if (res == '0.5x0.666') or drawcountries :
        m.drawcountries()
    parallels = np.arange(-90,91,15*interval)
    meridians = np.arange(-180,181,30*interval)
    if (res == '0.25x0.3125') :
        parallels = np.arange(-90,91,15*interval/3  ) #.25)
        meridians = np.arange(-180,181,30*interval/3)  #.3125)
#        parallels = np.arange(-90,91, .25)
#        meridians = np.arange(-180,181, .3125)

        # use small font size for greater the runs with more axis labele
#        f_size = f_size*.25
#    else:
    f_size = f_size*.75

    # draw meridian lines
    plt.xticks( meridians[::everyother], fontsize = f_size ) 
    # draw parrelel lines
    plt.yticks( parallels[::everyother], fontsize = f_size ) 
    m.drawcoastlines()
    plt.grid( show_grid )

    return m
 
# --------
# 4.37 - Provide an appropriate colormap for given data
# --------
def get_colormap( arr,  center_zero=True, minval=0.15, maxval=0.95, \
            npoints=100, cb='CMRmap_r', maintain_scaling=True, \
            negative=None, positive=None, debug=False ):
    """ Create correct color map for values given array.
        This function checks whether array contains just +ve or -ve or both 
        then prescribe color map  accordingly
    """
    
    # make sure array has a mask
    if debug:
        print arr, [ ( i.min(), i.max() ) for i in [arr] ]
        print '>'*5, ('ask' not in str( type( arr ) )), type( arr )
    if 'ask' not in str(type( arr ) ):
        arr = np.ma.array(arr) 
#        s_mask = arr.mask
    if debug:
        print '>'*5, ('ask' not in str( type( arr ) )), type( arr )
    
    # If postive/negative not given, check if +ve/-ve
    if any( [ not isinstance( i, type(None) ) for i in positive, negative ] ):

        # --- sequential
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

    #  reverse colourbar if negative
    if negative:
        if cb == 'CMRmap_r':
            cb = cb[:-2]
        else:
            cb = cb+'_r'

    print 'cmap is: >{}< & data is:'.format( cb ), 
    print '< postive == {}, negative == {}, divergent == {} >'.format(  \
            positive, negative, (( not positive) and (not negative))   )

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
    if ( not positive) and (not negative):
        if debug:
            print 'diverging data? =={}, for values range: '.format( True ),
        arr.mask = False
        cmap = plt.cm.RdBu_r
        # Centre color bar around zero
        if center_zero:
            vmin, vmax = arr.min(), arr.max()
            print 1-vmax/(vmax + abs(vmin)), vmin, vmax
            cmap = shiftedColorMap(cmap, midpoint=1-vmax/(vmax + abs(vmin)), 
                maintain_scaling=maintain_scaling, arr=arr, name='shifted')

    # Use updated jet replacements from mpl
    # : https://github.com/bids/colormap
#    cmap  = cmc  # pink, orange, blue alterative... 
#    cmap  = cmd  # green - blue alternative

#    arr.mask = s_mask
    if debug:
        print cb, center_zero                  
    return cmap
 
# --------
# 4.38 - Retrieves color by grouping of sensitivity study
# --------
def color4sensitstudy( title=None, rtn_dict=False):
    """ function to define plot colouring used for iodine mechanism 
    paper in ACPD"""
    d =  {
    'I$_{2}$Ox X-sections x2': 'red',
     'I$_{2}$Ox exp. X-sections': 'red',
     'I$_{2}$Ox loss ($\\gamma$) /2': 'deepskyblue',
     'I$_{2}$Ox loss ($\\gamma$) x2': 'deepskyblue',
     'Iodine simulation.': 'purple',
     'Just org. I': 'blue',
     'MBL BrO 2 pptv': 'saddlebrown',
     'Ocean iodide': 'magenta',
     'No I$_{2}$Ox Photolysis': 'red',
     'Sulfate Uptake': 'orange',
     'het. cycle ($\\gamma$) /2': 'orange',
     'het. cycle ($\\gamma$) x2': 'orange',
     'no het. cycle ': 'orange'
     }
    if rtn_dict:
        return d
    else:
         return d[title]
 
# --------
# 4.39 - Retrieves color by grouping of sensitivity study
# --------
def markers4sensitstudy( title=None, rtn_dict=False):
    """ function to define markers used for iodine mechanism
     paper in ACPD"""

    d =  {
    'I$_{2}$Ox X-sections x2': 'd',
     'I$_{2}$Ox exp. X-sections': 'h',
     'I$_{2}$Ox loss ($\\gamma$) /2': '^',
     'I$_{2}$Ox loss ($\\gamma$) x2': 'h',
     'Iodine simulation.': '^',
     'Just org. I': '^',
     'MBL BrO 2 pptv': '^',
     'Ocean iodide': '^',
     'No I$_{2}$Ox Photolysis': '^',
     'Sulfate Uptake': '^',
     'het. cycle ($\\gamma$) /2': 'd',
     'het. cycle ($\\gamma$) x2': 'h',
     'no het. cycle ': '+'
     }
    if rtn_dict:
        return d
    else:
         return d[title]
 
# --------
# 4.40 -Make segments for variable line color plot
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
# 4.41 - Make colored line for plot
# --------
def colorline( x, y, z=None, cmap=plt.get_cmap('copper'),  \
        norm=None, linewidth=3, alpha=1.0, ax=None, fig=None, \
        cb_title='Modelled wind direction', debug=False ):
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
        print [ type(i[0]) for i in x, z, y ]
        print [ i[0] for i in x, z, y ]

    # --- make sure x axis is float ( not a datetime )
    # convert to dataframe
    df = DataFrame( data=y, index=x) 

    # convert to datetime.datetime
    x = np.array(df.index.to_pydatetime(), dtype=datetime.datetime)
    if debug:
        print type(x[0]), x[0]

    # convert in float form
    x = matplotlib.dates.date2num(x) 
    if debug:
        print type(x[0]), x[0]

    # convert x and y into indices    
    segments = make_segments(x, y)

    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)

    if isinstance( ax, type(None) ):                              
        ax = plt.gca()
    ax.add_collection(lc)

    if not isinstance( fig, type( None) ):
        axcb = fig.colorbar(lc)
        axcb.set_label( cb_title )

    return lc

# --------
# 4.42 - Get human readable gradations for plot
# --------
def get_human_readable_gradations( lvls=None, vmax=10, vmin=0, \
            nticks=10, sigfig_rounding_on_cb=2, debug=False ):

    if isinstance( lvls, type(None) ):
        lvls = np.linspace( vmin, vmax, nticks, endpoint=True )

    # --- Adjust graduations in colourbar to be human readable
    if ( ( abs( int( vmax)) == 0) and (abs( int( vmax)) == 0) ):
        sigfig_rounding_on_cb += 1

    # significant figure ( sig. fig. ) rounding func.
    round_to_n = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))
    
    # get current gradations
    if debug:
        print abs(lvls[-2])-abs(lvls[-3])    
    lvls_diff = round_to_n( abs(lvls[-2])-abs(lvls[-3]), \
                                sigfig_rounding_on_cb)

    # round top of colorbar lvls, then count down from this
    vmax_rounded = round_to_n( vmax, sigfig_rounding_on_cb)
    if debug:
        print vmax_rounded,  lvls_diff, nticks
    lvls = np.array([ vmax_rounded - lvls_diff*i \
            for i in range( nticks ) ][::-1])   
        
    return lvls
# --------
# 4.43 - mk colourmap discrete 
# --------
def mk_discrete_cmap( lvls=None, cmap=None, arr=None,\
            vmin=0, vmax=10, nticks=10, debug=False ):
    """ Make a discrete colormap from an existing cmap
    NOTE: the data will now need to normalised the range of lvls """

    debug=True
    
    # define bins
    if isinstance( lvls, type(None) ):
        bounds = np.linspace( vmin, vmax, nticks, endpoint=True )
    else:
        bounds = lvls
    
    # get colormap if not provided
    if isinstance( cmap, type(None) ):
        cmap = get_colormap( np.array([vmin, vmax]) )

    if debug:
        print lvls, vmin, vmax, nticks
    
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

# -------------
# 4.99 - Input for plotting  
# -------------
def get_input_vars(debug=False):
    """ Get input from command line arguments 
        Takes: spec, filename, category, start date , end date  """

    # If working directory (wd) provided, use command line arguments
    try:
        wd = sys.argv[1]
        print '$>'*5, sys.argv

    # else prompt for input settings
    except:
        wd = raw_input("ctm.bpch dir: ")
        spec = raw_input("species (default = O3): ")
        if (len(spec) < 1):
            spec = 'O3'
        fn = raw_input("ctm.bpch name (default = ctm.bpch): ")
        if (len(fn) < 1):
            fn = 'ctm.bpch'
        cat_ = raw_input("Category (default = 'IJ-AVG-$'): ")
        if (len(cat_) < 1):
            cat_ = "IJ-AVG-$"
        start = raw_input("Enter start time as 'YYYY,MM,DD' (default = all): ")
        if (len(list(start)) < 1):
            start = None
        else:
            start = datetime.datetime( *tuple( map( int, start.split('-')) ) )
        end = raw_input("Enter end time as 'YYYY,MM,DD' (default = all): ")
        if (len(end) < 1):
            end = None
        else:
            end = datetime.datetime( *tuple( map( int, end.split('-')) ) )            
        print wd, fn, cat_, spec, start, end

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
        print sys.argv[5].split('-')
        start = datetime.datetime( *tuple( map( int, sys.argv[5].split('-')) ) )
    except:
        try:
            start
        except:
            start = None

    # Has end date been provided? ( in form "YYYY, MM, DD" )
    try:  
        print sys.argv[6].split('-')
        end = datetime.datetime( *tuple( map( int, sys.argv[6].split('-')) ) )
    except:
        try:
            end
        except:
            end = None

    if debug:
        print '$>'*5, wd, fn, cat_, spec, start, end

    # if final character not '/' add this
    if wd[-1]  != '/':
        wd += '/'

    return wd, fn, cat_, spec, start, end
