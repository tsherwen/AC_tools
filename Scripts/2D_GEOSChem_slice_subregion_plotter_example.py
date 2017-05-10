#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Plotter for 2D slices of GEOS-Chem output NetCDFs files.

NOTES
---
 - This is setup for Cly, but many other options (plot/species) are availible 
    by just updating passed variables/plotting function called. 
"""

import AC_tools as AC
import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Basic plotter of NetCDF files using AC_tools
    """
    # --- Local settings hardwired here...  
    fam = 'Cly' # Family to plot
    # print species in family for reference...
    print AC.GC_var(fam)
    
    # --- Get working directory etc from command line (as a dictionary object)
    # (1st argument is fil directory with folder, 2nd is filename)
    Var_rc = AC.get_default_variable_dict()
    # Get details on extracted data (inc. resolution)
    Data_rc = AC.get_shared_data_as_dict( Var_rc=Var_rc )
    
    # --- extract data and units of data for family/species... 
    arr, units = AC.fam_data_extractor( wd=Var_rc['wd'], fam=fam, \
        res=Data_rc['res'], rtn_units=True, annual_mean=False )
    
    # --- Process data (add and extra processing of data here... )
    # take average over time
    print arr.shape
    arr = arr.mean(axis=-1)
    # Select surface values
    print arr.shape
    arr = arr[...,0]
    # convert to pptv
    arr = arr*1E12
    units = 'pptv'
    
    # --- Plot up data... 
    print arr.shape
    #  - Plot a (very) simple plot ... 
#    AC.map_plot( arr.T, res=Data_rc['res'] )
    #  - plot a slightly better plot... 
    # (loads of options here - just type help(AC.plot_spatial_figure) in ipython) 
    # set range for data... 
    fixcb = np.array([ 0., 100. ])
    nticks = 6 # number of ticks on colorbar (make sure the fixcb range divides by this)
    interval = (1/3.) # number of lat/lon labels... (x*15 degrees... )
    # set limits of plot
    lat_min = 5.
    lat_max = 75.
    lon_min = -30.
    lon_max = 60.
    left_cb_pos = 0.85 # set X (fractional) position
    axis_titles = True # add labels for lat and lon
    # title for plot
    title ="Plot of annual average {}".format(fam)
    #  save as pdf (just set to True) or show?
#    figsize = (7,5) # figsize to use? (e.g. square or rectangular plot)

    # call plotter... 
    AC.plot_spatial_figure( arr, res=Data_rc['res'], units=units, fixcb=fixcb,\
        lat_min=lat_min, lat_max=lat_max, lon_min=lon_min, lon_max=lon_max,\
        axis_titles=axis_titles, left_cb_pos=left_cb_pos,\
        nticks=nticks, interval=interval, title=title, show=False )  

    # are the spacings right? - if not just up
    bottom =0.1
    top=0.9
    left=0.1
    right=0.9    
    fig = plt.gcf()
    fig.subplots_adjust( bottom=bottom, top=top, left=left, right=right )

    # show and save as PDF?
    plt.savefig('pete_plot.png')
    AC.show_plot()


if __name__ == "__main__":
    main()