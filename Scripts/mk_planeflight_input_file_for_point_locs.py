#!/usr/bin/python
"""
This programme makes the planeflight*.dat files required to output for specific locations and times in the model. The deafult settting is for hourly output.

NOTES:
 - This programme can be used to produce files to output data for observational sites (e.g. goverment air quality sites)
"""
# --- Packages
import numpy as np
from time import gmtime, strftime
import time
import glob
import AC_tools as AC

# --- Settings
wd = '/work/home/ts551/data/'#IO_obs'
# file of locations?
#pf_loc_dat_file ='ClBrI_ClNO2_PI_O3_PDRA.dat'
pf_loc_dat_file ='Weyborne.dat'
#pf_loc_dat_file ='EU_GRID_0.25x0.3125.dat' 
# Years output is required for?
start_year, end_year = 2014,2016
# debug?
debug = False
# output all reactions?
all_REAs = False
# 
# Which (halogen) code version is being used?
#ver = '1.6' # Iodine simulation in v9-2
#ver = '2.0' # Iodine + Bromine simulation
ver = '3.0' # Cl-Br-I simulation 
# add extra spacing? (needed for large amounts of output, like nested grids)
Extra_spacings =False

# --- Read in site Detail
numbers, lats, lons, pres, locs = AC.readin_gaw_sites( pf_loc_dat_file )
# make sure the format is numpt float 64
lats, lons, pres = [ np.float64(i) for i in lats, lons, pres ]
print lats[0:4]

# --- Set Variables
slist = AC.pf_var('slist', ver=ver )#_REAs_all')

# extra scpaes need for runs with many points
if Extra_spacings:
    pstr = '{:>6}  {:<4} {:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f}'
    endstr = '999999   END  0- 0-   0  0: 0    0.00    0.00    0.00'
else:
    pstr = '{:>5}  {:<3} {:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f}'
    endstr ='99999   END  0- 0-   0  0: 0    0.00    0.00    0.00 '

# remove 2nd reaction following 3 body reactions.
if all_REAs:
    # Get dictionary of rxns for wd/ver
    MUTDwd =  MUTD_runs()[0]
    rdict = AC.rxn_dict_from_smvlog( MUTDwd, ver=ver)

# setup required time range
nvar=len(slist)
yr = range(start_year, end_year )
m  = range(01,13)
da  = range(01,32,1)
h  = range(00,24,1)
mi = range(00,60,1)
minute = 0

# --- loop dates and make Planeflight log files for given points 
for year in yr:
    for b in m: 
        for c in da:

            if debug:
                print year, b, c
            a=open('Planeflight.dat.'+str(year)+'{:0>2}{:0>2}'.format( b, c), 'w')  

            # Print out file headers to pf.dat file
            print >>a, 'Planeflight.dat -- input file for ND40 diagnostic GEOS_FP'
            print >>a, 'Tomas Sherwen'
            print >>a, strftime("%B %d %Y", gmtime())
            print >>a, '-----------------------------------------------'
            print >>a, '{:<4}'.format(nvar),'! Number of variables to be output'
            print >>a, '-----------------------------------------------'

            # Print out species for GEOS-Chem to output to pf.dat file
            for i in range(0,len(slist)):
                print >>a, slist[i]

            # Print out species for GEOS-Chem to output to pf.dat file            
            print >>a, '-------------------------------------------------'
            print >>a, 'Now give the times and locations of the flight'
            print >>a, '-------------------------------------------------'
            print >>a, 'Point   Type DD-MM-YYYY HH:MM     LAT     LON   PRESS'

            # Loop requested dates and times
            counter=0
            for d in h:
                
                for i in range(len(lats)):
                    print >>a, pstr.format( counter, locs[i], c, b,  year, d, \
                                            minute, lats[i], lons[i], pres[i])
                    counter+=1

            # Add footer to pf.dat file
            print >>a, endstr
            a.close()
