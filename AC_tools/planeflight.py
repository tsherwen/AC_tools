#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Generic functions for use with GEOS-Chem's planeflight diagnostic module
"""
# compatibility with both python 2 and 3
from __future__ import print_function

# - Required modules:
# I/O functions / Low level
import sys
import csv
import glob
import pandas as pd
import logging
import numpy as np
# The below list to be updated, imports should be specific and in individual functions
# import tms modules with shared functions
from . core import *
from . generic import *
from . AC_time import *
from . variables import *
# Time
import datetime as datetime


def update_Planeflight_files(wd=None, num_tracers=103, verbose=True):
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
        'GMAO_ABSH', 'GMAO_PSFC', 'GMAO_SURF', 'GMAO_TEMP', 'GMAO_UWND',
        'GMAO_VWND'
    ]
    assert isinstance(num_tracers, int), 'num_tracers must be an integer'
    slist = ['TRA_{:0>3}'.format(i) for i in np.arange(1, num_tracers+1)]
    species = ['OH', 'HO2']
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
    files = glob.glob(wd+'Planeflight.dat*')
    if verbose:
        print(files)

    # --- Loop existing files and extract data
    dfs = []
    for n_file, file in enumerate(files):
        with open(file, 'r') as file_:
            # loop variables
            data_from_line_num = 9999
            data = []
            for n, line in enumerate(file_):
                # get header
                if ('Point' in line) and ('Type' in line):
                    header = [i.strip().upper() for i in line.split()]
                    data_from_line_num = n
                # Break if passed all data
                elif ('99999' in line) and ('END' in line):
                    break
                elif (n > data_from_line_num):
                    data += [[i.strip() for i in line.split()]]
                else:
                    pass

            if len(data) > 0:
                # Add datetime column, then add data to list of dataframes
                df = pd.DataFrame(np.array(data), columns=header)
                df['datetime'] = df['DD-MM-YYYY'].astype(
                    str)+df['HH:MM'].astype(str)
                df['datetime'] = pd.to_datetime(df['datetime'],
                                                format='%d-%m-%Y%H:%M')
                dfs += [df]
            else:
                err_msg = 'WARNING: no data in {}'.format(file)
                logging.info(err_msg)
                print(err_msg)

    # Concatenate extracted data
    if verbose:
        print(dfs[0])
    df = pd.concat(dfs).sort_values('datetime', ascending=True)
    if verbose:
        print('FINAL!!!!', df)

    # --- Print out new files based on processed DataFrame
    prt_PlaneFlight_files(df=df, slist=slist, Extra_spacings=False)


def prt_PlaneFlight_files(df=None, LAT_var='LAT', LON_var='LON',
                          PRESS_var='PRESS', loc_var='TYPE',
                          Username='Tomas Sherwen',
                          Date_var='datetime', slist=None, num_tracers=85,
                          Extra_spacings=False,
                          verbose=False, debug=False):
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
    Extra_spacings (bool): add extra spacing? (needed for large amounts of
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
        endstr = '99999   END 00-00-0000 00:00    0.00    0.00    0.00'
    # Output a general list of species/tracers/met vars if not provided as arguments
    if isinstance(slist, type(None)):
        met_vars = [
            'GMAO_ABSH', 'GMAO_PSFC', 'GMAO_SURF', 'GMAO_TEMP', 'GMAO_UWND', 'GMAO_VWND'
        ]
        assert isinstance(num_tracers, int), 'num_tracers must be an integer'
        slist = ['TRA_{:0>3}'.format(i) for i in np.arange(1, num_tracers+1)]
        species = ['OH', 'HO2']
        slist = slist + species + met_vars
    # Number of variables to output (needed for fortran read of *dat files)
    nvar = len(slist)
    # --- work out how many (UTC) days in the output
    # Get list of unique dates & remove mean from dates
    dates = [datetime.datetime(*i.timetuple()[:3]) for i in df[Date_var]]
    # Add list of just YYYYMMDD strings to dataframe
    df['YYYYMMDD'] = [i.strftime('%Y%m%d') for i in dates]
    # Get list of unique days
    dates = np.ma.array(sorted(set(dates)))
    # --- loop days and create the files
    for date_ in dates:
        # Get data for date
        sub_df = df[df['YYYYMMDD'].values == date_.strftime('%Y%m%d')]
        if verbose:
            print('Entries for day ({}): '.format(date_), sub_df.shape)
        # Create/Open up pf.dat setup
        a = open('Planeflight.dat.'+date_.strftime('%Y%m%d'), 'w')
        # Print out file headers to pf.dat file
        print('Planeflight.dat -- input file for ND40 diagnostic GEOS_FP',
              file=a)
        print(Username, file=a)
        print(strftime("%B %d %Y", gmtime()), file=a)
        print('-----------------------------------------------', file=a)
        print('{:<4}'.format(nvar), '! Number of variables to be output',
              file=a)
        print('-----------------------------------------------', file=a)
        # Print out species for GEOS-Chem to output to pf.dat file
        for n in range(0, len(slist)):
            print(slist[n], file=a)
        # Print out species for GEOS-Chem to output to pf.dat file
        print('-------------------------------------------------', file=a)
        print('Now give the times and locations of the flight', file=a)
        print('-------------------------------------------------', file=a)
        print('Point  Type DD-MM-YYYY HH:MM     LAT     LON   PRESS', file=a)
        # Loop requested times
        for n, time_ in enumerate(sub_df[Date_var]):
            # Setup variable list to print
            vars_ = [n+1, sub_df[loc_var].values[n], time_.day, time_.month]
            vars_ += [time_.year, time_.hour, time_.minute]
            # Extract lat, lon, and pressure for time
            coord_vars = LAT_var, LON_var, PRESS_var
            vars_ += [float(sub_df[i].values[n]) for i in coord_vars]
            # Set formating
            vars_ = pstr.format(*vars_)
            # Print to file
            print(vars_, file=a)
        # Add footer to pf.dat file
        print(endstr, file=a)
        a.close()


def prt_PlaneFlight_files_v12_plus(df=None, LAT_var='LAT', LON_var='LON',
                                   PRESS_var='PRESS', loc_var='TYPE',
                                   OBS_var='OBS',
                                   Date_var='datetime', slist=None,
                                   num_tracers=85, rxn_nums=[],
                                   Extra_spacings=False,
                                   Username='Tomas Sherwen', verbose=False,
                                   debug=False):
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
    Extra_spacings (bool): add extra spacing? (needed for large amounts of
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
        #        pstr = '{:>6}  {:>4} {:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f}'
        #        endstr = '999999   END  0- 0-   0  0: 0    0.00    0.00    0.00'
        print('Extra_spacings not setup for >= v12.0.0')
        sys.exit()
    else:
        #        pstr = '{:>5}  {:<3} {:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f}'
        #        endstr ='99999   END  0- 0-   0  0: 0    0.00    0.00    0.00 '
        pstr = '{:>5}{:>7} {:0>2}-{:0>2}-{:0>4} {:0>2}:{:0>2}  {:>6,.2f} {:>7,.2f} {:>7.2f} {:>10.3f}'
        endstr = '99999   END  00-00-0000 00:00   0.00     0.00    0.00      0.00'
    # Output a general list of species/tracers/met vars if not provided as arguments
    if isinstance(slist, type(None)):
        met_vars = [
            'GMAO_ABSH', 'GMAO_PSFC', 'GMAO_SURF', 'GMAO_TEMP', 'GMAO_UWND',
            'GMAO_VWND', 'GMAO_PRES'
        ]
        assert isinstance(num_tracers, int), 'num_tracers must be an integer'
        slist = ['TRA_{:0>3}'.format(i) for i in np.arange(1, num_tracers+1)]
        species = ['OH', 'HO2']
        slist = slist + species + met_vars
        # Add list of reactions to extract too
        if len(rxn_nums) > 0:
            rxns_are_nums = [(type(i) == int) for i in rxn_nums]
            assert all(rxns_are_nums), 'All rxn numbers must be integers!'
            slist += ['REA_{:0>3}'.format(i) for i in rxn_nums]
    # Number of variables to output (needed for fortran read of *dat files)
    nvar = len(slist)
    # --- Make sure an altitude is defined in df if not provided
    # Updates merged into v12.0.0 mean a OBS altitude is required.
    try:
        df[OBS_var].values
    except KeyError:
        fill_ALT_obs = 99999.00
        df[OBS_var] = fill_ALT_obs
    # --- work out how many (UTC) days in the output
    # Get list of unique dates & remove mean from dates
    dates = [datetime.datetime(*i.timetuple()[:3]) for i in df[Date_var]]
    # Add list of just YYYYMMDD strings to dataframe
    df['YYYYMMDD'] = [i.strftime('%Y%m%d') for i in dates]
    # Get list of unique days
    dates = np.ma.array(sorted(set(dates)))
    # --- loop days and create the files
    for date_ in dates:
        # Get data for date
        sub_df = df[df['YYYYMMDD'].values == date_.strftime('%Y%m%d')]
        if verbose:
            print('Entries for day ({}): '.format(date_), sub_df.shape)
        # Create/Open up pf.dat setup
        a = open('Planeflight.dat.'+date_.strftime('%Y%m%d'), 'w')
        # Print out file headers to pf.dat file
        print('Planeflight.dat -- input file for ND40 diagnostic GEOS_FP',
              file=a)
        print(Username, file=a)
        print(strftime("%B %d %Y", gmtime()), file=a)
        print('-----------------------------------------------', file=a)
        print('{:<4}'.format(nvar), '! Number of variables to be output',
              file=a)
        print('-----------------------------------------------', file=a)
        # Print out species for GEOS-Chem to output to pf.dat file
        for n in range(0, len(slist)):
            print(slist[n], file=a)
        # Print out species for GEOS-Chem to output to pf.dat file
        print('-------------------------------------------------', file=a)
        print('Now give the times and locations of the flight', file=a)
        print('-------------------------------------------------', file=a)
        header = [
            'Point', 'Type', 'DD-MM-YYYY', 'HH:MM', 'LAT', 'LON', 'PRESS',
            'OBS'
        ]
        h_pstr = '{:>5}{:>7} {:>10} {:>5}  {:>6} {:>7} {:>7} {:>10}'

        print(h_pstr.format(*header), file=a)
        # Loop requested times
        for n, time_ in enumerate(sub_df[Date_var]):
            # Setup variable list to print
            vars_ = [n+1, sub_df[loc_var].values[n], time_.day, time_.month]
            vars_ += [time_.year, time_.hour, time_.minute]
            # Extract lat, lon, and pressure for time
            coord_vars = LAT_var, LON_var, PRESS_var, OBS_var
            vars_ += [float(sub_df[i].values[n]) for i in coord_vars]
            # Set formating
            vars_ = pstr.format(*vars_)
            # Print to file
            print(vars_, file=a)
        # Add footer to pf.dat file
        print(endstr, file=a)
        a.close()


def get_pf_headers(file, debug=False):
    """
    Extract column headers from a GEOS-Chem planeflight csv file

    Parameters
    -------
    file (str): filename to open
    debug (bool): debug the function?

    Returns
    -------
    (list, list)
    """
    if debug:
        print(file)
    # Open pf file
    with open(file, 'r') as f:
        Lines = [i for i in f]
        names = Lines[0].strip().split()
        points = [i.strip().split()[0] for i in Lines[1:]]
#         reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
#         for row in f:
#             if row[0] != 'POINT':
#                 new = row[1:2][0]
#                 try:
#                     points.append(new)
#                 except:
#                     points = [new]
#             else:
#                 names = row[2:]
    if debug:
        print(names, list(set(points)))
    return names, list(set(points))


def pf_csv2pandas(file=None, vars=None, epoch=False, r_vars=False,
                  debug=False):
    """
    Planeflight.dat CSV reader - used for processor GEOS-Chem PF output

    Parameters
    -------
    file (str): file name (inc. directory)
    vars (list): vars to extract
    epoch (bool):
    r_vars (bool): return list of vars

    Returns
    -------
    (pd.DataFrame)
    """
    # Open file
    with open(file, 'r') as f:
        logging.debug(f)

        # Label 1st column ( + LOC ) if names not in vars
        # ( This effectively means that pandas arrays are the same )
        if debug:
            print([type(i) for i in (vars, ['POINT', 'LOC'])])
        if 'POINT' not in vars:
            names = ['POINT', 'LOC'] + vars[:-1]
        else:
            names = vars
        if debug:
            print(vars, names)
        # Convert to pandas array
        df = pd.read_csv(f, header=None, skiprows=1,
                         delim_whitespace=True, names=names,
                         dtype={'HHMM': str, 'YYYYMMDD': str, 'POINT': object}
                         )
        # Convert strings to datetime using pandas mapping
        df = DF_YYYYMMDD_HHMM_2_dt(df, rmvars=None, epoch=epoch)
        if debug:
            print(df, df.shape)

    # Return pandas DataFrame
    if r_vars:
        return df, list(df.columns)
    else:
        return df


def get_pf_from_folder(folder='./', dates2use=None, debug=False):
    """
    Get GEOS-Chem planeflight output from folder
    """
    # Which files to use?
    files = list(sorted(glob.glob(folder+'/*plane.*')))
    if debug:
        print(files)
    # Only open dates for certain dates?
    if not isinstance(dates2use, type(None)):
        FileRootsVar = 'FileRoots'
        df = pd.DataFrame(files)
        df = pd.DataFrame({FileRootsVar: files})
        # Which date format to look for in filenames?
        format = '%Y%m%d%H%M'
        # Setup a helper function to extract dates from file strings

        def get_date_from_filename(x, format=format):
            """
            Extract Dates from filenames

            Notes
            -------
             - It is assumed that the date ends the file string before the
             format identifier
            """
            date_str = x.split('.')[-2]
            dt = datetime_.strptime(date_str, format)
            return dt
        dtVar = 'datetime'
        df[dtVar] = df[FileRootsVar].map(get_date_from_filename)
        bool = df[dtVar].isin(dates2use)
        files = list(df.loc[bool, FileRootsVar].values)
    # Get headers
    ALL_vars, sites = get_pf_headers(files[0], debug=debug)
    # Extract dfs
    dfs = [pf_csv2pandas(i, vars=ALL_vars) for i in files]
    # Append the dataframes together
    df = dfs[0].append(dfs[1:])
    return df


def get_pf_data_from_NetCDF_table(ncfile=None, req_var='TRA_69', spec='IO',
                                  loc='CVO', start=None, end=None, ver='1.7',
                                  sdate=None, edate=None,
                                  verbose=False, debug=False):
    """
    Extracts data from NetCDF file processed by pf2NetCDF (pandas) converter in PhD_Progs

    Returns
    -------
    (np.array, np.array)

    Notes
    -------
     - TODO re-write and update to comments!
     (this function should just used xarray)
    """
    # Convert to plane-flight (pf) variable name ('req_var') if not given
    if isinstance(req_var, type(None)):
        req_var = what_species_am_i(spec, ver=ver, invert=True)

    # --- Open NetCDF within nest, and extract data
    with Dataset(ncfile, 'r') as rootgrp:

        # Select only variables for site
        LOC = rootgrp['LOC']
        Epoch = rootgrp['Epoch']
        data = rootgrp[req_var]
        # Convert to numpy array
        LOC, Epoch, data = [np.array(i) for i in (LOC, Epoch, data)]

    # Select just the variable requests
    if debug:
        print(LOC)
    ind = np.where(LOC == loc)
    if debug:
        print('indcies where LOC==loc: ', ind)
    Epoch = Epoch[ind]
    data = data[ind]

    # Covert Epoch to datetime
    if debug:
        print(Epoch[0])
    dates = np.array([datetime_.fromtimestamp(i) for i in Epoch])
    if debug:
        print(dates[0])

    # Select dates ( <= add this )
    if not isinstance(sdate, type(None)):
        data = data[np.where((dates < edate) & (dates >= sdate))]
        dates = dates[np.where((dates < edate) & (dates >= sdate))]

    return dates, data


def mk_planeflight_input4FAAM_flight(folder=None, ds=None,
                                     flight_ID='C216',
                                     folder4csv=None,
                                     PressVar="PS_RVSM",
                                     LonVar='LON_GIN',
                                     LatVar='LAT_GIN', TimeVar='Time',
                                     LocVar='TYPE', LocName='B-146',
                                     DateVar='datetime',
                                     testing_mode=True, csv_suffix='',
                                     num_tracers=203, rxn_nums=[],
                                     Username='Tomas Sherwen',
                                     slist=None,
                                     Extra_spacings=False
                                     ):
    """
    Extract the GEOS-CF model for a given FAAM BAe146 flight

    Parameters
    -------

    Returns
    -------
    (None)
    """
    # Retrieve FAAM BAe146 Core NetCDF files
    if isinstance(ds, type(None)):
        filename = 'core_faam_*_{}_1hz.nc'.format(flight_ID.lower())
        file2use = glob.glob(folder+filename)
        if len(file2use) > 1:
            print('WARNING: more that one file found! (so using latest file)')
            print(file2use)
        ds = xr.open_dataset(file2use[0])
    # Only select the variable of intereest and drop where these are NaNs
    df = ds[[PressVar, LatVar, LonVar, TimeVar]].to_dataframe()
    df = df.dropna()
    # Add a location name (as Type)
    df[LocVar] = LocName
    # remove the index name and add index values to a column
    df.index.name = None
    try:
        df[DateVar]
    except KeyError:
        df['datetime'] = df.index.values
    # If doing a test, then just extract the first 150 points of flight track
    if testing_mode:
        df = df.head(150)
    # Call planeflight maker...
    prt_PlaneFlight_files_v12_plus(df=df, slist=slist,
                                   Extra_spacings=Extra_spacings,
                                   LON_var=LonVar, LAT_var=LatVar,
                                   PRESS_var=PressVar, loc_var=LocVar,
                                   num_tracers=num_tracers, rxn_nums=rxn_nums,
                                   Username=Username,)


def reprocess_split_pf_output_over_2_lines(folder, save_original_file=True):
    """
    Combine planeflight dat file lines where output split over 2 lines
    """
    files2use = list(sorted(glob.glob(folder+'*plane.log*')))
    for file2use in files2use:
        with open(file2use, 'r') as file:
            lines = [i.strip() for i in file]
        file.close()  # Force close
        if save_original_file:
            os.rename(file2use, file2use+'.orig')
        else:
            os.remove(file2use)
        first_part = lines[0::2]
        second_part = lines[1::2]
        a = open(file2use, 'w')
        for n_line, line in enumerate(first_part):
            Str2use = '{}     {}'
            Newline = Str2use.format(first_part[n_line], second_part[n_line])
            print(Newline, file=a)
        a.close()
