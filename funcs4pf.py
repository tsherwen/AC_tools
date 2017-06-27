#!/usr/bin/python
# -*- coding: utf-8 -*-
""" 
Generic functions for use with GEOS-Chem's planeflight diagnostic module
"""
# --- compatibility with both python 2 and 3
from __future__ import print_function

# ----------------------------- Section 0 -----------------------------------
# -------------- Required modules:

# -- I/O functions / Low level          
import sys
import csv
import glob
import pandas as pd
import logging

# Math functions
import numpy as np

# -- This needs to be updated, imports should be specific and in individual functions
# import tms modules with shared functions
if __package__ is None:
    from .funcs4core import *
    from .funcs4generic import *
    from .funcs4time import *
    from .funcs_vars import *
else:
    from funcs4core import *
    from funcs4generic import *
    from funcs4time import *
    from funcs_vars import *

# -- Time functions
import datetime as datetime

# ---------------------------------- Section X.X ---------------------------
# -------------- Planeflight Processors
#
# (some example functions using these functions are in AC_tools/Scripts)

# ----
# X.XX - 
# ----
def update_Planeflight_files( wd=None, num_tracers=103, verbose=True ):
    """
    Create new planeflight from old files (with updated # of tracers)

    Parameters
    -------
    wd (str): the working (code) directory to search for files in
    num_tracers (int): the number of tracers (TRA_???) to print in *dat files

    Notes
    -------
     - Used for using existing planeflight output for campaign, but for 
     different output variables
    """
    # --- Local variables
    output_data_str = 'Now give the times and locations of the flight'
    # 
    met_vars = [
    'GMAO_ABSH', 'GMAO_PSFC','GMAO_SURF', 'GMAO_TEMP', 'GMAO_UWND', 'GMAO_VWND'
    ]
    assert isinstance( num_tracers, int), 'num_tracers must be an integer'    
    slist = ['TRA_{:0>3}'.format(i) for i in np.arange(1, num_tracers+1 ) ]
    species  = ['OH', 'HO2']
    slist = slist + species + met_vars
#    slist = pf_var( fill_var_with_zeroes=True, ver=ver )

    # ---  Get files 
    # Get wd 
    if isinstance(wd, type(None)):
        try:
            wd = sys.argv[1]    
        except:
            print('FAIL - Please provide working directory!')
    # vars to use?
    # read in files and extract locations.
    files = glob.glob( wd+'Planeflight.dat*' )
    if verbose:
        print(files)
    
    # --- Loop existing files and extract data
    dfs = []
    for n_file, file in enumerate( files ):
        with open(file, 'rb') as file_:
            # loop variables
            data_from_line_num = 9999    
            data =[]        
            for n, line in enumerate( file_ ):
                # get header
                if ('Point' in line) and ('Type' in line):
                    header = [i.strip().upper() for i in line.split()]
                    data_from_line_num = n
                # Break if passed all data
                elif ('99999' in line) and ('END' in line):   
                    break
                elif (n>data_from_line_num):
                    data += [ [i.strip() for i in line.split()] ]         
                else:
                    pass 
            
            if len(data) > 0: 
                # Add datetime column, then add data to list of dataframes 
                df = pd.DataFrame( np.array(data), columns=header )
                df['datetime'] = df['DD-MM-YYYY'].astype(str)+df['HH:MM'].astype(str)
                df['datetime'] = pd.to_datetime(df['datetime'],\
                format='%d-%m-%Y%H:%M')
                dfs += [ df ] 
            else:
                err_msg = 'WARNING: no data in {}'.format( file )
                logging.info( err_msg )
                print(err_msg)

    # Concatenate extracted data
    if verbose:
        print(dfs[0])
    df = pd.concat(dfs).sort_values('datetime',ascending=True)
    if verbose:
        print('FINAL!!!!' , df)

    # --- Print out new files based on processed DataFrame 
    prt_PlaneFlight_files(df=df, slist=slist, Extra_spacings=False)


# ----
# X.XX - 
# ----
def prt_PlaneFlight_files( df=None, LAT_var='LAT', LON_var='LON', \
        PRESS_var='PRESS', loc_var='TYPE', Username='Tomas Sherwen', \
        Date_var='datetime', slist=None, Extra_spacings=False, \
        verbose=False, debug=False ):
    """
    Takes a dataframe of lats, lons, alts, and times and makes Planeflight.dat.*
    files

    Parameters
    -------
    df (pd.DataFrame): dataframe of data (as floats) indexed by times(datetime)
    wd (str): the working (code) directory to search for files in
    loc_var (str): name for (e.g. plane name), could be more than one. 
    LAT_var, LON_var, PRESS_var (str): name for pressure(HPa),lat and lon in df
    Date_var (str): column name of df containing datetime (UTC) variables
    Username (str): name of the programme's user
    Extra_spacings (boolean): add extra spacing? (needed for large amounts of 
        output, like nested grids)
    slist (list): list of tracers/species to output

    Notes
    -----
     - to get output of a specific frequency for given point locations, just 
    add the times to the axis of the dataframe provided. 
     -  datetime columns is required (as this allows mulitple output loations
     (e.g. sepeerate planes/sites) to be present in input df)
     - This function expects the dataframe to be ordered by datetime
    """
    # --- Packages
    from time import gmtime, strftime
    import time

    # --- Local variables
    # Extra spaces need for runs with many points
    if Extra_spacings:
        pstr = '{:>6}  {:>4} {:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f}'
        endstr = '999999   END  0- 0-   0  0: 0    0.00    0.00    0.00'
    else:
#        pstr = '{:>5}  {:<3} {:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f}'
#        endstr ='99999   END  0- 0-   0  0: 0    0.00    0.00    0.00 '
        pstr = '{:>5}  {:>4} {:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f}'
        endstr ='99999   END 00-00-0000 00:00    0.00    0.00    0.00'
    # Number of variables to output (needed for fortran read of *dat files)
    nvar = len(slist)

    # --- work out how many (UTC) days in the output
    # Get list of unique dates & remove mean from dates
    dates = [datetime.datetime(*i.timetuple()[:3]) for i in df[Date_var] ]
    # Add list of just YYYYMMDD strings to dataframe
    df['YYYYMMDD'] = [i.strftime('%Y%m%d') for i in dates ]
    # Get list of unique days
    dates = np.ma.array( sorted( set(dates) ) )
    
    # --- loop days and create the files 
    for date_ in dates:
    
        # Get data for date
        sub_df = df[ df['YYYYMMDD'].values == date_.strftime('%Y%m%d') ]
        if verbose:
            print('Entries for day ({}): '.format(date_), sub_df.shape)

        # Create/Open up pf.dat setup     
        a=open('Planeflight.dat.'+date_.strftime('%Y%m%d'),'w')  

        # Print out file headers to pf.dat file
        print('Planeflight.dat -- input file for ND40 diagnostic GEOS_FP', file=a)
        print(Username, file=a)
        print(strftime("%B %d %Y", gmtime()), file=a)
        print('-----------------------------------------------', file=a)
        print('{:<4}'.format(nvar),'! Number of variables to be output', file=a)
        print('-----------------------------------------------', file=a)

        # Print out species for GEOS-Chem to output to pf.dat file
        for n in range(0,len(slist)):
            print(slist[n], file=a)

        # Print out species for GEOS-Chem to output to pf.dat file            
        print('-------------------------------------------------', file=a)
        print('Now give the times and locations of the flight', file=a)
        print('-------------------------------------------------', file=a)
        print('Point  Type DD-MM-YYYY HH:MM     LAT     LON   PRESS', file=a)

        # Loop requested times 
        for n, time_ in enumerate( sub_df[Date_var] ):
            # Setup variable list to print
            vars_ = [ n+1, sub_df[loc_var].values[n], time_.day, time_.month ]
            vars_ += [ time_.year, time_.hour, time_.minute ]
            # Extract lat, lon, and pressure for time
            coord_vars = LAT_var, LON_var, PRESS_var
            vars_ += [ float( sub_df[i].values[n] ) for i in coord_vars ]
            # Set formating
            vars_ = pstr.format( *vars_ )
            # Print to file
            print(vars_, file=a)

        # Add footer to pf.dat file
        print(endstr, file=a)
        a.close()


# ---------------------------------- Section X.X ---------------------------
# -------------- Planeflight Extractors
#

# -------------- 
# X.XX - Get headers for a given pf file (return var names, and points )
# ----------
def get_pf_headers(file, debug=False):
    """ 
    Extract column headers from a GEOS-Chem planeflight csv file 

    Parameters
    -------
    file (str): filename to open
    debug (boolean): debug the function?

    Returns
    -------
    (list, list)
    """
    if debug:
        print(file)
    # Open pf file 
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
        print(names, list( set(points) ))
    return names, list( set(points) )


# -------------- 
# X.0X - converts planeflight.dat file to pandas array
# ----------
def pf_csv2pandas(file=None, vars=None, epoch=False, r_vars=False, \
        debug=False):
    """ 
    Planeflight.dat CSV reader - used for processor GEOS-Chem PF output

    Parameters
    -------
    file (str): file name (inc. directory)
    vars (list): vars to extract
    epoch (boolean):
    r_vars (boolean): return list of vars

    Returns
    -------
    (pd.DataFrame)
    """
    # Open file
    with open(file, 'rb') as f :
        logging.debug( f )

        # Label 1st column ( + LOC ) if names not in vars
        # ( This effectively means that pandas arrays are the same )
        if debug:
            print([ type(i) for i in (vars, ['POINT', 'LOC'])  ])
        if 'POINT' not in vars:
            names=['POINT', 'LOC'] + vars[:-1]
        else:
            names =vars
        if debug:
            print(vars, names) 

        # Convert to pandas array
        df = pd.read_csv( f, header=None, skiprows=1, \
            delim_whitespace=True, names=names, dtype={'HHMM':str, \
            'YYYYMMDD':str, 'POINT':object}   )

        # Convert strings to datetime using pandas mapping
        df = DF_YYYYMMDD_HHMM_2_dt( df, rmvars=None, epoch=epoch )
        if debug:
            print(df, df.shape)

    # Return pandas DataFrame
    if r_vars:
        return df, list( df.columns )
    else:
        return df   


# ----
# X.XX - Read pf data from 2D NetCDF table file 
# ----
def get_pf_data_from_NetCDF_table( ncfile=None, req_var='TRA_69', spec='IO', \
        loc='CVO', start=None, end=None, ver='1.7', verbose=False, \
        sdate=None, edate=None, debug=False  ):
    """ 
    Extracts data from NetCDF file processed by pf2NetCDF (pandas) converter in PhD_Progs

    NOTES
    Returns
    -------
    (np.array, np.array)    

    Notes
    -------
     - TODO re-write and update to comments!
    """
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
        LOC, Epoch, data = [ np.array(i) for i in (LOC, Epoch, data) ]

    # Select just the variable requests
    if debug:
        print(LOC)
    ind = np.where( LOC == loc  )
    if debug:
        print('indcies where LOC==loc: ', ind)
    Epoch = Epoch[ ind ]
    data = data[ ind ]

    # Covert Epoch to datetime
    if debug:
        print(Epoch[0]) 
    dates = np.array( [ datetime_.fromtimestamp(i) for i in Epoch ] )
    if debug:
        print(dates[0]) 

    # Select dates ( <= add this )
    if not isinstance( sdate, type(None) ):
        data = data[ np.where( (dates < edate ) & (dates >= sdate  )  ) ]
        dates = dates[ np.where( (dates < edate ) & ( dates >= sdate  )  ) ]

    return dates, data



