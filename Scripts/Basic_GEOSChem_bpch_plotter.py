#!/usr/bin/python
# modules
from MChem_tools import *
import numpy as np
import sys

# Setup, choose species
species  = 'O3'#'CO2'
RMM_species = 16.*3.
res = '4x5' # ( e.g. '4x5', '2x2.5', '0.5x0.666', '0.25x0.3125'  )
unit, scale = tra_unit( species, scale=True)

# debug/print verbose output?
debug=True

# Only consider GEOS-Chem chemical troposphere
trop_limit=True
calc_burden=False#True

try:    # chcck if a directory was given ad command line
    wd = sys.argv[1]
except: # Otherwise use path below
    wd = '<insert GEOS-Chem run direcotory path here>'

# get data as 4D array ( lon, lat, alt, time ) 
mixing_ratio  = get_GC_output( wd, species=species, category='IJ-AVG-$', \
    trop_limit=trop_limit ) 
print mixing_ratio.shape

# Get data to calculate burden
if calc_burden:

    # Get air mass as numpy array
    air_mass = get_GC_output( wd, vars=['BXHGHT_S__AD'], \
                                  trop_limit=trop_limit )
    # get time in troposphere as fraction from ctm.bpchh diagnostics
    time_in_trop = get_GC_output( wd, vars=['TIME_TPS__TIMETROP'],
                                  trop_limit=trop_limit )
    # print shapes of array to screen
    print [i.shape for i in mixing_ratio, air_mass, time_in_trop ]

    # calc the total & mean mass of a speices - select data you want to calc & print
    # mass (convert to g) / RMM air
    air_moles = ( air_mass *1E3 ) / ( .78*(2.*14.)+.22*(2.*16.) )
    
    # moles * v/v
    species_moles = air_moles * mixing_ratio / scale

    # convert to mass
    species_mass = species_moles * RMM_species              
    print species_mass.shape

    # Get Global Burden, converting to Tg and removing stratosphere
    print np.sum( np.mean( species_mass[:,:,:38,:]*time_in_trop, axis=3 ) ) / 1E12   

#  select data you want to plot - set as mean() of time, but simply remove 'mean(axis=3)' & give number to select desired timestamp
mixing_ratio = mixing_ratio[:,:,:,:].mean( axis=3 ) 

# select surface layer
print mixing_ratio.shape
mixing_ratio = mixing_ratio[:,:,0]
print mixing_ratio.shape

# Transpose to X vs. Y coordinate
mixing_ratio = np.transpose( mixing_ratio )
print mixing_ratio.shape, res

# plot surface
plt, cb = basic_map_plot( mixing_ratio, species, unit, res=res, debug=debug )

# Show plot
plt.show()
