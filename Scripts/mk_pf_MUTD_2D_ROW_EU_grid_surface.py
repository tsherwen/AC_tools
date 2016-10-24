# --- Packages
import AC_tools as AC
import numpy as np
from time import gmtime, strftime
import time
import glob

# --- Settings
res='0.25x0.3125'
#res='0.5x0.666'
tpwd = AC.get_dir('tpwd')
#start_year, end_year = 2006,2007
#start_year, end_year = 2012,2013
#start_year, end_year = 2014,2017
start_year, end_year = 2015,2017
debug = False
#ver='1.7'
ver='3.0'
pf_locs_file = 'EU_GRID_{}.dat'.format( res )
#pf_locs_file ='ROW_GRID.dat'

# --- Read in site Detail
location = AC.readin_gaw_sites( pf_locs_file, all=True )
numbers, locs, lats, lons, pres =  [ location[:,i] for i in range(5) ]
lats, lons, pres = [ np.float64(i) for i in lats, lons, pres ]
locs = np.array(locs)
print lats[0:4]

# --- Set Variables
#slist = pf_var( 'slist_v9_2_NREA_red_NOy', ver=ver )#'slist_v9_2_NREA_red' )
#slist = pf_var( 'slist_ClearFlo', ver=ver )
slist = AC.pf_var( 'slist_PEN_WEY_LEI', ver=ver )
# kludge additional of ethanol
slist = slist+ ['TRA_86']

nvar=len(slist)
yr = range(start_year, end_year )
#m  = range(11,13)
m  = range(01,13)
#m  = range(04,05)
#m  = range(06,07)
da  = range(01,32,1)
h  = range(00,24,1)
mi = range(00,60,1)
minute = 0
endstr='999999   END  0- 0-   0  0: 0    0.00    0.00    0.00'
pstr='{:>6}   {:<4}{:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f}'

# --- loop dates and make Planeflight log files for given points
for year in yr:
    for b in m: 
        for c in da:

            if debug:
                print year, b, c
            a=open('Planeflight.dat.'+str(year)+'{:0>2}{:0>2}'.format( b, c), \
                     'w')  

            print >>a, 'Planeflight.dat -- input file for ND40 diagnostic GEOS_FP'
            print >>a, 'Tomas Sherwen'
            print >>a, strftime("%B %d %Y", gmtime())
            print >>a, '-----------------------------------------------'
            print >>a, '{:<4}'.format(nvar)        ,'! Number of variables to be output'
            print >>a, '-----------------------------------------------'

            for i in range(0,len(slist)):
                print >>a, slist[i]
            
            print >>a, '-------------------------------------------------'
            print >>a, 'Now give the times and locations of the flight'
            print >>a, '-------------------------------------------------'
            print >>a, 'Point   Type DD-MM-YYYY HH:MM     LAT     LON   PRESS'
                    
            counter=0
            for d in h:
                
                for i in range(len(lats)):
                    print >>a,pstr.format( counter, locs[i], c, b,  year, d, minute, \
                                           lats[i], lons[i], pres[i])
                    counter+=1
                   
            print >>a, endstr
            a.close()
