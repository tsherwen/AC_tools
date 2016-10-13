# ------------ Planeflight Plotter - tms -------------------------------------
# --------------  
#!/usr/bin/python
import matplotlib.pyplot as plt 
from MChem_tools import *
import sys

# --------- SET PLOTTING HERE ----------------------------------------
# -------------
# Plot out data?
#debug=True

# Where are the model files? And What Species do you want to plot?
# ---- Get inputs from command line/defaults
try:    # chcck if a directory was given ad command line
    wd    = sys.argv[1]
except: # Otherwise use path below
    wd    = '<insert GEOS-Chem run direcotory path here>'
    
# Which species to plot?  (must be in list form) 
try:    # chcck if a directory was given ad command line
    species_to_plot = sys.argv[2]
except: # Otherwise use path below
    species_to_plot = 'Cl2' #'O3'#'ClNO2'#'Cl'

# What years, months, days to plot?  (must be in list form) 
# set day_to_use by adjusting range
years_to_use, months_to_use= ['2005'], ['12']
days_to_use  = [ "{0:0>2}".format(i)  for i in range(1, 31, 1) ]
print years_to_use, months_to_use, days_to_use 

# Which locations? (where 10 = 1000 Hpa, 08 = 800 Hpa etc ... ) (must be in list form) 
# e.g. locations=['TX1','LA1'] # must be as a list of strings
# must be as a list of strings
#locations=['WE'] 
#locations=['WEY'] 
locations=['CVO'] 
print locations

# Model version
ver='3.0'

# look in the "plane_flight_logs" directory
wd = wd+ '/plane_flight_logs/plane.log.*'
#wd = wd+ '/plane.log.*'
print wd
print sorted(glob.glob(wd))

# Asectics
fontsize = 10

# ----------- START PLOTTING HERE ----------------------------------------
# -------------

# Get species name in TRA_?? (planefligth output) form
species_to_plot=[ what_species_am_i( species_to_plot, ver=ver, invert=True ) ]

# setup figure
fig = plt.figure(figsize=(15,6), dpi=80, facecolor='w', edgecolor='k')

# Loop sites in site list and plot species
for i,site in enumerate(locations):

    # extract data from planeflight (csv) files
    model, names = readfile( sorted(glob.glob(wd) ), site, \
        years_to_use, months_to_use, days_to_use)

    # get species index in list
    k=names.index(species_to_plot[0])

    # plot up extracted data
    plt.plot(  year_to_since_2006(model), model[:,k]*1e9*1E3,color=plt.cm.jet(1.*i/len(locations)),\
        label='{0}'.format(locations[i])  )    

# Beatify plot
plt.xlabel('Time/ CV days', fontsize=fontsize)
plt.ylabel('Conc./ p.p.t.v.', fontsize=fontsize)
plt.legend(loc='upper right',fontsize=fontsize)
plt.grid(b=None, which='major', axis='both', alpha=0.3)
plt.rcParams.update({'font.size': fontsize})

# Show plt
plt.show()
