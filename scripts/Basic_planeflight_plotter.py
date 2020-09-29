#!/usr/bin/python
# ------------ Planeflight Plotter - tms -------------------------------------
# --------------
import glob
import sys
import matplotlib.pyplot as plt
import AC_tools as AC

# --------- SET PLOTTING HERE ----------------------------------------
# -------------
# Plot out data?
# debug=True

# Where are the model files? And What Species do you want to plot?
# ---- Get inputs from command line/defaults
try:    # chcck if a directory was given at command line
    wd = sys.argv[1]
except:  # Otherwise use path below
    wd = '<insert GEOS-Chem run direcotory path here>'

# Which species to plot?  (must be in list form)
try:    # chcck if a directory was given at command line
    species_to_plot = sys.argv[2]
except:  # Otherwise use path below
    species_to_plot = 'Cl2'  # 'O3'#'ClNO2'#'Cl'

# What years, months, days to plot?  (must be in list form)
# set day_to_use by adjusting range
years_to_use, months_to_use = ['2005'], ['05', '06', '07']
days_to_use = ["{0:0>2}".format(i) for i in range(1, 31, 1)]
print((years_to_use, months_to_use, days_to_use))

# Which locations? (where 10 = 1000 Hpa, 08 = 800 Hpa etc ... ) (must be in list form)
# e.g. locations=['TX1','LA1'] # must be as a list of strings
# must be as a list of strings
# locations=['WE']
locations = ['WEY']
# locations=['CVO']
# locations=['BEI']
print(locations)

# Scaling (e.g. pptv or ppbv )
units, scale = 'p.p.t.v.', 1E12
#units, scale = 'p.p.b.v.', 1E9
#units, scale = 's$^{-1}$', 1

# Model version
ver = '3.0'

# look in the "plane_flight_logs" directory
wd = wd + '/plane_flight_logs/plane.log.*'
#wd = wd+ '/plane.log.*'
print(wd)
print((sorted(glob.glob(wd))))

# Asectics
fontsize = 10

# ----------- START PLOTTING HERE ----------------------------------------
# -------------

# Get species name in TRA_?? (planefligth output) form
# if 'TRA' in species_to_plot:
species_to_plot = [AC.what_species_am_i(species_to_plot, ver=ver, invert=True)]
# else:
#    species_to_plot=[species_to_plot]

# setup figure
fig = plt.figure(figsize=(15, 6), dpi=80, facecolor='w', edgecolor='k')

# Loop sites in site list and plot species
for i, site in enumerate(locations):

    # extract data from planeflight (csv) files
    model, names = AC.readfile(sorted(glob.glob(wd)), site,
                               years_to_use, months_to_use, days_to_use)

    # get species index in list
    k = names.index(species_to_plot[0])

    # plot up extracted data
    plt.plot(AC.year_to_since_2006(model), model[:, k]*scale, color=plt.cm.jet(1.*i/len(locations)),
             label='{0}'.format(locations[i]))

# Beatify plot
plt.xlabel('(CVAO days)', fontsize=fontsize)
plt.ylabel('({})'.format(units), fontsize=fontsize)
plt.legend(loc='upper right', fontsize=fontsize)
plt.grid(b=None, which='major', axis='both', alpha=0.3)
plt.rcParams.update({'font.size': fontsize})

# Show plt
plt.show()
