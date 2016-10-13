# --- Packages
from AC_tools.funcs4GEOSC import *
from AC_tools.funcs4pf import * 
import numpy as np
from time import gmtime, strftime
import time
import glob

# --- Settings
wd = get_dir( 'dwd' )
#camp = 'planeflight/planeflight_dat_files_Malasapina'#'/CAST_2014_BAE146/'
#camp ='/IO_obs/TORERO/updatedtoreroiodata_surface_data'
camp = 'planeflight/pf_dat_files_CONTRAST_GV_CGV_hv_66_TRA_4_Johan_UPDATE16/'
start_year, end_year = 2013, 2015
debug, tag = True, 'CON' #'TRB'#'MAL'#'ANT'
convert_km_2_hPa, convert_m_2_hPa, convert_2_m = False, False, False 
time_str ='%H:%M' # '%H%M'  # '%H:%M:%S' #    # '%h/%m/%s',
#ver = '1.6' # Iodine simulation in v9-2
#ver = '2.0' # Iodine + Bromine simulation
ver = '3.0' # Cl-Br-I simulation

# --- Set Variables
slist = pf_var('slist', ver=ver)
nvar=len( slist )
yr = range(start_year, end_year )
m  = range(01,13)
da  = range(01,32,1)
# string format for variable print
pstr1 = '{:>5}  {:<3}  {:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f}'
# planeflight.dat footer
pstr2 = '99999   END  0- 0-   0  0: 0    0.00    0.00    0.00'

if debug:
    print yr, m, da

# --- loop files for a given campaign and make Planeflight dat files
for year in yr:
    for b in m: 
        for c in da:

            print year
            if debug:
                print 1, year,b, c , wd + camp + \
                        '/'+'*_{}_{:0>2}_{:0>2}*'.format( str(year),b, c )  
            try:
                fn_ = glob.glob( wd + camp + '/'+'*_{}_{:0>2}_{:0>2}*'.format( \
                                            str(year),b, c )   )[0]
                if debug:
                    print 2, fn_
                
                # --- Extract info from csv.
                times, lats, lons, press = read_in_kml_sites( fn_, debug=debug )
                if debug:
                    print 3, 'esc read_in_kml_sites'
                    print times[:10]
                    print 4, 'times, lats, lons, pres'
                    print 5, [(  i[:3] )  for i in  [times, lats, lons, press] ]
                times = [ time.strptime( i, time_str ) for i in times ]            

                if debug:
                    print 6,times[:10]
                h = [ time.strftime('%H',i) for i in times ]
                if debug:    
                    print 7, h
                minute = [ time.strftime('%M',i) for i in times ]

                times =(times)
                lats, lons, pres = [ np.float64(i) for i in lats, lons, press ]

                if debug:
                    print 8, [len(i) for i in  [times, lats, lons, pres] ]
                # Convert units if requested 
                if convert_2_m: # convert feet to meters
                    pres = [  i*0.3048/1000 for i in pres]  
                if convert_km_2_hPa:
                    pres  = hPa_to_Km(  pres , reverse=True, debug=debug)
                if convert_m_2_hPa:
                    pres = [i/1E3 for i in pres]
                    pres  = hPa_to_Km(  pres, reverse=True, debug=debug)
                if debug:
                    print 9, lats[0:4]                  

                # Create/Open up pf.dat setup     
                a=open( 'Planeflight.dat.'+str(year)+ \
                        '{:0>2}{:0>2}'.format( b, c), 'w')  
                             
                if debug:
                    print a

                # Print out file headers to pf.dat file
                print >>a, 'Planeflight.dat -- input file for ND40 diagnostic GEOS_FP'
                print >>a, 'Tomas Sherwen'
                print >>a, strftime("%d %b %Y", gmtime())
                print >>a, '-----------------------------------------------'
                print >>a, nvar        ,'! Number of variables to be output'
                print >>a, '-----------------------------------------------'

                # Print out species for GEOS-Chem to output to pf.dat file
                for i in range(0,len(slist)):
                    print >>a, slist[i]        

                # Print out locations of GEOS-Chem output to pf.dat file
                print >>a, '-------------------------------------------------'
                print >>a, 'Now give the times and locations of the flight'
                print >>a, '-------------------------------------------------'
                print >>a, 'Point  Type DD-MM-YYYY HH:MM     LAT     LON   PRESS'

                # Loop requested dates and times                    
                for i in range(len(lats)):
                    print >>a, pstr1.format( i,tag, c, b, int(year), int(h[i]),\
                        int(minute[i]), lats[i], lons[i], pres[i])

                # Add footer to pf.dat file
                print >>a, pstr2
            except :
                if debug:
                    print 'no plane flight on: ', b, c, int(year)
