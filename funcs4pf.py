#!/usr/bin/python
""" 
Generic functions for use with output from GEOS-Chem's planeflight module.

NOTE(S):    
 - This module is underdevelopment vestigial/inefficient code is being removed/updated. 
 - Where external code is used credit is given. 
"""

# ----------------------------- Section 0 -----------------------------------
# -------------- Required modules:

# -- I/O functions / Low level          
import sys
import csv
import glob
import pandas as pd
from pandas import HDFStore
from pandas import DataFrame
#from pandas import Series
#from pandas import Panel
import logging

# Math functions
import numpy as np

# -- This needs to be updated, imports should be specific and in individual functions
# import tms modules with shared functions
from funcs4core import *
from funcs4generic import *
from funcs4time import *

# -- Time functions
import datetime as datetime

# -- Variables
from funcs_vars import *

# ------------------------------------------- Section 1 -------------------------------------------
# -------------- Planeflight Setup tools
#
# --------------
# 1.00 - Read in file of GAW sites as lists
# -------------
def readin_gaw_sites(filename, all=False):
    """ 
    Read in list format csv, with details of sites to output planeflight for

    NOTES:
    - REDUNDENT? 
    """

    with open(filename,'rb') as f:
        reader = csv.reader(f, delimiter=',') 
        for row in reader:
            new = row[:]
            try:
                locations.append(new)

            except:
                locations=[new]

        locations=np.array(locations)
        if all:
            return locations
        else:
            numbers = locations[:,0]
        #    IDs = locations[:,1]
            lats = locations[:,2]
            lons = locations[:,3]
            pres = locations[:,4]
            locs = locations[:,5]
    return numbers, lats, lons, pres, locs

# --------------
# 1.01 - Read in file of site lists -  
# -------------
def read_in_kml_sites(filename, limter=10, ind=[0, 3, 1, 2 ], debug=False):
    """ 
    Read in list format csv, with details of sites to output planeflight for.
    NOTES:
     - double up of function 1.00?
    """

    if debug:
        print 'read_in_kml_sites called for ', filename
    reader = csv.reader( open( filename, 'rb' ), delimiter=',' ) 

    if debug:
        print reader
    for row in reader:
        if debug:
            print row
        new = row[:limter]
        try:
            locations.append(new)

        except:
            locations=[new]

    if debug:
        print type( locations)
        print len(locations), locations[:2]
        print len(locations[0]), locations[0]

    if debug:
        print '-'*100
    for ii, i in enumerate( locations ):
        if ( len(i) != 5 ):
            if debug:
                print len(i),  i 
    if debug:
        print '-%'*50

    locations=np.array(locations)

    if debug:
        print type( locations)
        print len(locations), locations[:2]
        print len(locations[0]), locations[0]
        print locations.shape
        
    # extract vars  = 
    times, altitude, lats, lons = [ locations[:,i] for i in ind ]

    if debug:
        print 'complete read of : {}'.format( filename )
        print [ len(i) for i in [times, lats, lons, altitude ]]
        print [ (i[0], i[-1]) for i in [times, lats, lons, altitude ]]
        print 'min max:' , [ (min(i), max(i)) for i in [times, lats, lons, altitude ]]

    return times, lats, lons, altitude

# --------------
# 1.01 - Read files from Dix et al/Volkamer et al
# -------------
def read_TOR_IO_files(filename, debug = False):
    """ 
    Read in csv files from TORERO campaign 
    """
    if debug:
        print 'read_TOR_IO_files called for ', filename

    reader = csv.reader(open(filename,'rb'), delimiter=',') 
    print reader, filename
    data_line = False
    for row in reader:
        new = row[:][0].strip().split()
        print  new[0], (new[0] == 'UTC' ), data_line
        if (data_line):
            try:
                locations.append(new)
            except:
                locations=[new]

        if (new[0] == 'UTC' ):
            data_line = True

    locations=np.array(locations)
    times = locations[:,0]
    altitude = locations[:,1]
    lats = locations[:,2]
    lons = locations[:,3]
    data = locations[:,-2]

    return times, lats, lons, altitude, data


# ------------------------------------------- Section 3 -------------------------------------------
# -------------- Planeflight Extractors
#

# ----
#  3.01 - Get plane flight data for a given location and species
# ----
def wd_pf_2_data( wd, spec, location='TOR', scale=1E12, r_datetime=False,   \
            Kludge_fortan_output=False, ver=None, verbose=True, debug=False):
    """ 
    Read in sites from files in a given working directory

    NOTES:
     - UPDATE NEEDED - the location assignment needs to be more generalised
    """
    if verbose:
        print wd, spec, location, scale

    # Get version if not given as agrument
    if isinstance( ver, type(None) ):
        ver = iGEOSChem_ver( wd )  

    # Extract from planeflight.dat files 
    model, names = readfile_basic( sorted(glob.glob(wd + \
        '/plane_flight_logs/plane.log.2*')), location, \
        Kludge_fortan_output=Kludge_fortan_output,  \
                                                            debug=debug )
    # Set list of surface ( point or ship ) datasets
    surface_data =[ 
    'MAL', 'HAL','GRO', 'CVO', 'ANT', 'M91', 'KEN', 'BTT', 'TRB', \
    'HOU', 'BOL', 'LAC', 'HES', 'SCH', 'TEX', 'CAL', 'PAS', 'WEY', \
    'MON'
    ]
    # And PI surface sites
    surface_data += [
    'HOB', 'ADE','MOV', 'MAU', 'LUA' , 'HOK', 'TOK', 'COI',  'PDM',
    'MVT', 'VIE', 'MON' ]
    # plane data
    plane_data = ['TOR', 'BAE', 'GCV', 'CGV', 'CON' ]

    # Pull values for a variable from data arrays and scale
    try:
        k=names.index( spec )
    except:
        print '>'*30, 'ERROR: failed to find >{}<'.format(spec) + \
            '< in names list ({}), trying planeflight equivilent'. format( \
            ver ), '<'*30, names
		# Translate species
        spec = what_species_am_i( spec, ver=ver, invert=True ) 

        # Get index in list
        k=names.index( spec )
        if debug:
            print k, len(model[:,0] ) , scale

    data = model[:,k]*scale
    
    if debug:
        print location, spec, k, len(data)

    # Provide data and Altitude for vertical measurements 
    # (and ancillaries for if r_datetime == True )
    if  any( [(location == i)  for i in plane_data ]):
        j =names.index( 'PRESS' )
        press = model[:,j]
        if debug:
            print [ ( len(i) , min(i), max(i) )  for i in [data, press]]
        press = [ np.float(i) for i in press ]
        km =  np.array( hPa_to_Km( press )  )
        if r_datetime:
            return [ data, km]  +[ model[:,names.index( i) ] \
                for i in [ 'LAT', 'LON', 'YYYYMMDD', 'HHMM' ]]
        else:
            return data, km 

    # Provide data and extra vars for surface measuremnets
    # " [ 'LAT', 'LON', 'YYYYMMDD', 'HHMM' ]" 
    # data, lat, lon, date, time 
    if any([ (location == i) for i in surface_data] ): 
        vars = [ model[:,names.index( i) ] \
            for i in [ 'LAT', 'LON', 'YYYYMMDD', 'HHMM' ]]
        return [data] + vars

# -------------- 
# 2.02 - Basic planeflight Output reader - mje
# -------------- 
def readfile_basic(files, location, debug=False, \
            Kludge_fortan_output=False, rm_empty=False):
    """ 
    basic readfile function for planeflight output in csv 

    ARUGEMENTS:
     - Kludge_fortan_output (boolean): This is a "kludge" option, for when the numbering 
     is printed out as "***" due to incorrect format statements in planeflight_mod.F
    """
    if debug:
    	print "files: {}, Location: {}".format( files, location )
    for file in files:
        if debug:
            print file
        with open(file,'rb') as f:
            reader=csv.reader(f, delimiter=' ', skipinitialspace = True)
            for row in reader:
                if row[1] == location: 
                    new=row[2:]
                    if debug:
                        print new
                    try:    
                        big.append(new)
                    except:
                        big=[new]
                if row[0] == 'POINT':
                    names = row[2:]

    if Kludge_fortan_output:
        big  = [ row for row in big if (row[4] != '*******' )]
        big  =  [i[:220] for i in big ]

    if debug:
        print np.array(big).shape
    if rm_empty:
        names.pop(-1)
    big = np.float64(big)             
    return big, names

# --------------
# 2.03 - date specific (Year,Month,Day) planeflight output reader - tms
# -------------
def readfile( files, location, years, months, days, plot_all_data=False, 
            debug=False, **kwargs):
    """ 
    Date specific readfile function for planeflight output in csv

    NOTES:
	 - UPDATE NEEDED - The approach be more efficent and readable
    """
    plot_all_data=False
    print 'readfile called, for {} files, with 1st filename: {}'.format(len(files), files[0])
    big, names = [],[]

    # sort for choosen years/months
    for file in files:
        # loop through list on choosen years
        if (not plot_all_data):
            lll = 0
            for year in range(len(years)):
                if (("{0}".format(years[year])) in (file)) :

                    # is it not the last year specificied?
                    if debug:
                        print 'years[year]', years[year], 'years[-1]', years[-1]
                    if (not (years[year] == years[-1])):
                        # just read all years upto point uptyo final year
                        big, names=process_files_to_read(file, location,big, \
                            names, debug=debug)
                        if debug:
                            print 'i got to line 91'

                    # If last year selected, then only plot the given
                    #  months & days
                    if (years[year] == years[-1]):
                        # Plot months exceot last one
                        for month in range(len(months)):                                                                                                  
                            if debug:
                                print 'months[month]', months[month], \
                                    'months[-1]', months[-1], 'months',  \
                                    months, 'type(months)', type(months)
                            if (("{0}{1}".format(years[year], \
                                months[month])) in file) :  
                                if (not (months[month] == months[-1])):
                                    big, names=process_files_to_read(file,\
                                         location,big, names, debug=debug)
                                    if debug:
                                        print 'i got to line 100',  'month',\
                                             month, 'in',len(months),  'year',\
                                              year, 'in' , len(years)
                                if (months[month] == months[-1]):
                                        # For last month, 
                                        # plot days upto last day
                                    for day in range(len(days)):                                                                                          
                                        if (("{0}{1}{2}".format(years[year], \
                                            months[month],days[day])) in file) : 
                                            if debug:
                                                print 'days[day]', days[day], \
                                                    'days[-1]', days[-1]
                                            big, names=process_files_to_read( \
                                                file, location, big, names,\
                                                 debug=debug)
                                            if debug:
                                                print 'i got to line 108'
                                                print 'readfile read big'\
                                                    +' of size: ', len(big)
                                                
        if (plot_all_data):
            big, names=process_files_to_read(files, location,big, names)
            print 'reading all data'

    big=np.float64(big)             
    print 'Read {} pf.log files, contructed arary with dims.: {}'.format( \
        len(files), big.shape )
    if (len(files) == 0) or (big.shape == (0,) ):
        print files, '<'*100, 'No Planeflight.log* files read/extacted,'\
            +' please check dates - dates recieved: ', years, months, days,'>'*100
        sys.exit(0)
    return big, names

# -------------- 
# 2.04 - files sent by a selected read of files
# ----------
def process_files_to_read(files, location, big, names, debug = True):
    """ Function to extract data from GEOS-Chem planeflight csv 
            output to list form  """

    if ( debug ):
        print 'process_files_to_read called'
        print files
    reader=csv.reader(open(files,'rb'), delimiter=' ', skipinitialspace = True)
    for row in reader:
#        print location , (row[0] == 'POINT'), (row[1] == location) ,  len(row) , len(big), len(names)  \ #
#                                        , big.shape#, (row[2:])[0], (row[-1]) , l
        if row[1] == location: 
            new=row[2:]
            try:    
                big.append(new)
            except:
                big=[new]
        if row[0] == 'POINT':
            names = row[2:]
    return big, names


# -------------- 
# 2.06 - Get headers for a given pf file (return var names, and points )
# ----------
def get_pf_headers(file, debug=False):
    """ 
    Extract column headers from a GEOS-Chem planeflight csv file 
    """
    if debug:
        print file
    # open pf file 
    with open(file,'rb') as f:
        reader=csv.reader(f, delimiter=' ', skipinitialspace = True)
        for row in reader:
            if row[0] != 'POINT':
                new=row[1:2][0]
                try:    
                    points.append(new)
                except:
                    points=[new]
            else:
                names = row[2:]
    if debug:
        print names, list( set(points) )
    return names, list( set(points) )

# -------------- 
# 2.07 - pf 2 pandas binary
# ----------
def pf2pandas(wd, files, vars=None, npwd=None, rmvars=None,   \
            debug=False):
    """ 
    Read in GEOS-Chem planeflight output and convert to HDF format

     - Converts date and time columns to datetime format indexes
     - the resultant HDF is in 2D list form 
    ( aka further processing required to 3D /2D output  )
        
    Note: 
     - This function is limited by the csv read speed. for large csv output expect 
     significant processing times or set to automatically run post run
     - Original files are not removed, so this function will double space usage for 
     output unless the original fiels are deleted.
    """

    # Ensure working dorectory string has leading foreward slash
    if wd[-1] != '/':
        wd += '/'

#    pfdate =( re.findall('\d+', file ) )[-1]
    if not isinstance(vars, list ):
        vars, sites = get_pf_headers( files[0], debug=debug )
    if not isinstance(npwd, str ):
        npwd = get_dir('npwd')
    hdf =HDFStore( npwd+ 'pf_{}_{}.h5'.format( wd.split('/')[-3], \
        wd.split('/')[-2], wd.split('/')[-1]  ))
    
    if debug:
        print hdf

    for file in files:
        print file#, pfdate

        # convert planeflight.log to DataFrame
        df = pf_csv2pandas( file, vars )
            
        if file==files[0]:
            hdf.put('d1', df, format='table', data_columns=True)
        else:
            hdf.append('d1', df, format='table', data_columns=True)

        if debug:
            print hdf['d1'].shape, hdf['d1'].index
        del df
    hdf.close()

# -------------- 
# X.0X - converts planeflight.dat file to pandas array
# ----------
def pf_csv2pandas( file=None, vars=None, epoch=False, r_vars=False, \
            debug=False ):
    """ 
    Planeflight.dat CSV reader - used for processor GEOS-Chem PF output
    """

    # Open file
    with open(file, 'rb') as f :

#        if debug:
#            print '>'*30, f
        logging.debug( f )

        # Label 1st column ( + LOC ) if names not in vars
        # ( This effectively means that pandas arrays are the same
        if debug:
            print [ type(i) for i in vars, ['POINT', 'LOC']  ]
        if 'POINT' not in vars:
            names=['POINT', 'LOC'] + vars[:-1]
        else:
            names =vars
        if debug:
            print vars, names 

        # convert to pandas array
        df = pd.read_csv( f, header=None, skiprows=1, \
            delim_whitespace=True, names=names, dtype={'HHMM':str, \
            'YYYYMMDD':str, 'POINT':object}   )

        # convert strings to datetime using pandas mapping
        df = DF_YYYYMMDD_HHMM_2_dt( df, rmvars=None, epoch=epoch )

        if debug:
            print df, df.shape

#    print df.columns   
#    print list( df.columns )

    # return pandas DataFrame
    if r_vars:
        return df, list( df.columns )
    else:
        return df   



# ----
# 2.10 - Extract all pf data for a given site.
# ----
def wd_pf_2_gaw_arr( wd, spec='O3', location='CVO', scale=1E9, debug=False ):
    """ 
    Extract all data rom a GEOS-Chem planeflight csv files in given 
    working directory, returning this in numpy array form 
    """

#    print wd
    files  =  sorted(glob.glob(wd +'/plane_flight_logs/plane.log.2*'))
    model, names = readfile_basic( files, location, debug=debug )
    data = np.float64( model[:,names.index( spec )]*1E9 )
    date = np.int64( model[:,names.index( 'YYYYMMDD' )] )
    time = np.int64( model[:,names.index( 'HHMM' )] )
    print [ ( len(i) , min(i), max(i) )  for i in [data, date, time]]
    return data, date, time

# ----
#  2.11 - Process "raw" csv files from GEOS-Chem planeflight output
# ----
def pro_raw_pf( wd, site='CVO', ext='', run='', frac=False, diurnal=True, \
            res='4x5', debug=False ):
    """ 
    Process csv rom GEOS-Chem planeflight output and save these as numpy memory arrays 

    NOTES:
     -  this is redundent. remove this function.
    """
    np_wd  = get_dir('npwd' )

    # Open & get data
    data, date, time = wd_pf_2_gaw_arr( wd, location=site  )
    print [(type(i),len(i) ) for i in [ time, date, data ] ]

    # get local time and shift time to match, - adjust dates,Cut off adjusted day at end an begining
    time, date, data = adjust_UT_2_lt( time, date, data, site=site)

    # Get O3 max, adjusted diurnal and monthly mean from data and time (int64)
    months, O3_max, O3_min, O3_d = hr_data_2_diurnal(time, date, data, \
        frac, diurnal=diurnal )

    # save out numpy of this data                                                                                                      
    if (diurnal):
#        fp = np.memmap( np_wd+'{}_mean_monthly_diurnal_O3_pf{}_{}.npy'.format(site, ext, run), dtype='float32', mode='w+', shape=(12,24) )
        diurnal_name =  '{}_{}_mean_monthly_diurnal_O3_pf{}_{}_{}.npy'
        fp = np.memmap(np_wd+diurnal_name.format( site, res, ext, run,\
             wd.split('/')[-2] ), dtype='float32', mode='w+', shape=(12,24) )
    else:
        fp = np.memmap( np_wd+'{}_{}_mean_monthly_O3_pf{}_{}.npy'.format(site, \
            res, ext, run), dtype='float32', mode='w+', shape=(12,24) )

    print fp
    fp[:] = months[:]
    np.memmap.flush(fp)

# ----
# 2.12 - Read pf data from 2D NetCDF table file 
# ----
def get_pf_data_from_NetCDF_table( ncfile=None, req_var='TRA_69', spec='IO', \
        loc='CVO', start=None, end=None, ver='1.7', verbose=False, \
        sdate=None, edate=None, debug=False  ):
    """ 
    Extracts data from NetCDF file processed by pf2NetCDF (pandas) converter in PhD_Progs
    """
#    debug=True

    # Convert to plane-flight (pf) variable name ('req_var') if not given
    if isinstance( req_var, type(None) ):
        req_var = what_species_am_i( spec, ver=ver, invert=True ) 

    # --- Open NetCDF within nest, and extract data
    with Dataset( ncfile, 'r' ) as rootgrp:

        # Select only variables for site
        LOC = rootgrp['LOC']
        Epoch = rootgrp['Epoch']
        data =  rootgrp[ req_var ]
        # Convert to numpy array
        LOC, Epoch, data = [ np.array(i) for i in LOC, Epoch, data ]

    # Select just the variable requests
    if debug:
        print LOC
    ind = np.where( LOC == loc  )
    if debug:
        print 'indcies where LOC==loc: ', ind
    Epoch = Epoch[ ind ]
    data = data[ ind ]

    # Covert Epoch to datetime
    if debug:
        print Epoch[0] 
    dates = np.array( [ datetime_.fromtimestamp(i) for i in Epoch ] )
    if debug:
        print dates[0] 

    # Select dates ( <= add this )
    if not isinstance( sdate, type(None) ):
        data = data[ np.where( (dates < edate ) & (dates >= sdate  )  ) ]
        dates = dates[ np.where( (dates < edate ) & ( dates >= sdate  )  ) ]

    return dates, data

# --------------------------------- Section 3  ---------------------------------
# -------------- Planeflight Analysis/Post formating 

# --------------
# 3.01 - Process time/date to datetime equivalent
# -------------
def pf2datetime( model, debug=False ):
    """ 
    Takes planeflight as numpy "model" array and returns list of datetimes 
    """

    # Make sure input is in integer form 
    vars = [ i.astype(np.int64) for i in model[:,0], model[:,1] ] 
    if debug:
        print [ i[:10] for i in vars ], model.shape

    # make pandas dataframe 
    df = DataFrame( data=np.array(vars).T, 
        columns=['YYYYMMDD', 'HHMM' ])

    # convert to datetime and retrun 
    df = DF_YYYYMMDD_HHMM_2_dt( df )

    # convert datetime 64 to datetime <= map this for scalability
    dates =  [ i.to_datetime() for i in df.index ]
    if debug:
        print 2, df.index, type( df.index[0] ), type( dates[0] )
    
    return dates

# --------------                                                                                                                                             
# 3.03 - adjust non UT times to UT
# ------------- 
def pf_UT_2_local_t(time_s, site='CVO', half_hour=None, debug=False):
    """ 
    Adjust pf time in UT to local time
    NOTES:
     - double up of function?
    """
    if debug:
        print 'pf_UT_2_local_t'
    #  UT diff library for sites ....
    UT_diff={'CVO':-1.}                                  

    # look  up diff and adjust to decimal hours
    t_diff= float(  UT_diff[site] )*float(1./24.)     

    # adjust time series to UT
    time_s_adjust = [ i +  t_diff for i in time_s]  

    # cut values to length
    time_s_adjust = [ float(str(i)[:str(i).find('.')+10])   \
            for i  in time_s_adjust  ] 

    # return to as numpy array
    time_s_adjust = np.array(time_s_adjust)  
#    for i,ii in enumerate(time_s):
#        print i , ii,time_s_adjust[i],t_diff
    print t_diff
    return time_s_adjust



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