#!/usr/bin/python
""" 
Convert planeflight output to NetCDF form en masse. A list of variables can be provided to restrict output 

This function is written to convert any pf output to NetCDF.

NOTES:
 - if GRD output is provided, 3D output can be obtained ( also in NetCDF form )
"""

import glob
import time
import os.path
import sys
import numpy as np
from pandas import DataFrame
import AC_tools as AC

# ---  Master  settings for main call
# Verbose/debug output? (set  debug=True)
DEBUG = True
VERBOSE = False
# Use planeflight files ending with ".out", which have been renumerated
# (Fortran string formatting cuts off output > 5 digits )
renumerated = True  # True#False#True # Use
# make 3D gridded output netCDF?
GRD_input_3D = True  # False#True
# Are there mulitple sites?
Multiple_sites = False  # True  # NOTE: this option is not currently working


def main(wd, vars=None, npwd=None, GRD_input_3D=False, renumerated=False,
         verbose=False, debug=False):
    """
    Driver to process planeflight output from GEOS-Chem

    NOTES:
    ---
     - more details on GEOS-Chem's planeflight diagnostic:
     (http://acmg.seas.harvard.edu/geos/doc/man/chapter_13.html) 
    """
    # Get save directory and set output NC name
    import os
    if not isinstance(npwd, str):
        npwd = get_dir('npwd')
    out_nc = npwd + 'pf_{}_{}.nc'.format(wd.split('/')[-3],
                                         wd.split('/')[-2], wd.split('/')[-1])
    print(("Atempting to append/create file (already exists?:{}): {}".format(
        os.path.isfile(out_nc), out_nc)))

    # Get pf files
    if not os.path.isfile(out_nc):
        files = get_pf_files(wd, renumerated=renumerated)

    # Make NetCDF as table of all pf files.  ( check for file first )
    if not os.path.isfile(out_nc):
        mk_NetCDF_of_pf_files(files, ncfilename=out_nc, debug=debug)

    # If 2D data, make 3D (lon, lat, time) NetCDF file
    if GRD_input_3D:
        make_3D_NetCDF(ncfilename=out_nc, wd=wd, debug=debug)

    # Process multiple sites to "subgrouped" NetCDF file
    # NOTE: this is currently not functioning ... TODO
    if Multiple_sites:
        make_2D_subgroup_NetCDF(ncfilename=out_nc, wd=wd, debug=debug)


def get_pf_files(wd, renumerated=False, debug=False):
    """ 
    Get pf files - edit this give location of planeflight files
    """
    # Ensure working dorectory string has leading foreward slash
    if wd[-1] != '/':
        wd += '/'

    # Get files
    filenames = '/plane_flight_logs/plane*log*'
    if renumerated:
        filenames += '.out'

    files = glob.glob(wd + filenames)
    if debug:
        print((wd, filenames, files))

    return files


def mk_NetCDF_of_pf_files(files, ncfilename=None, debug=False):
    """ 
    Make a table like NetCDF file from to any pf output
    """
    # --- Setup NetCDF file
    ncfile = Dataset(ncfilename, 'w', format='NETCDF4')

    # ---  Loop files, read in and add to NetCDF
    npoint = 1
    for n, file in enumerate(files):

        # If 1st file setup NetCDF
        if n == 0:

            # Get Header infomation from first file
            vars, sites = get_pf_headers(files[0], debug=debug)

            # Extract all points from file
            df, vars = AC.pf_csv2pandas(file=file, vars=vars, epoch=True,
                                        r_vars=True)
            if debug:
                print((df.shape, df.columns))

            # set unlimited data points dimension (POINT)
            POINT = ncfile.createDimension('POINT', None)

            # loop and create variables for each column  (exc. last )
            if debug:
                print(vars)
            [ncfile.createVariable(var, var2type(var), ('POINT'))
             for var in vars]

            # close the file
            ncfile.close()

        else:
            # Extract all points from file
            df, vars = AC.pf_csv2pandas(file=file, vars=vars, epoch=True,
                                        r_vars=True)

        # Open the file in append mode
        ncfile = Dataset(ncfilename, 'a', format='NETCDF4')
        if debug:
            print((df.index))

        # Fill variables for given
        dim_len = len(df.index)
        for var in vars:
            ncfile.variables[var][npoint:npoint+dim_len] = df[var].values

        # Tidy up and count
        npoint += dim_len
        del df
        ncfile.close()


def var2type(var, debug=False):
    """ Insure that strings are i8 type, add additions to list
        for a NetCDF, type must be:
            'f4' (32-bit floating point), 
            'f8' (64-bit floating point), 
            'i4' (32-bit signed integer), 
            'i2' (16-bit signed integer), 
            'i8' (64-bit singed integer), 
            'i1' (8-bit signed integer), 
            'u1' (8-bit unsigned integer), 
            'u2' (16-bit unsigned integer), 
            'u4' (32-bit unsigned integer), 
            'u8' (64-bit unsigned integer), or 
            'S1' (single-character string)
    ... Also:  also a 'S' datatype for variable length strings ( ==numpy object)
    """
    if any([i in var for i in ['TYPE', 'Epoch']]):
        case = 2
    elif any([i in var for i in ['LOC']]):
        case = 3
    else:
        case = 1

    cases = {
        1: 'f8',  # 'f8' (64-bit floating point),
        2: 'i8',  # 'i8' (64-bit singed integer),
        #    3: 'S'   # also a 'S' datatype for variable length strings ( numpy object)
        3: 'S1'
    }

    return cases[case]


def get_3D_vars(vars):
    """ from a list of pf variables, remove known 2D varables """

    known2D = ['Epoch', 'LON', 'LAT', 'YYYYMMDD', 'LOC', 'HHMM', 'POINT']
    return [i for i in vars if (i not in known2D)]


def get_site_description_vars(vars):
    """ from a list of pf variables, get known site specific varables """

    known2D = [
        'Epoch', 'LON', 'LAT', 'YYYYMMDD', 'LOC', 'HHMM', 'POINT', 'PRESS'
    ]
    return [i for i in vars if (i not in known2D)]


def make_3D_NetCDF(ncfilename, wd, debug=False):
    """ Create NetCDF of 3D arrays for all variables in 2D NetCDF file
         Takes a table form NetCDF and build 3D arrays from lat and lon
         in the given file.  """

    # --- Read existing & setup new NetCDF file
    ncfile2D = Dataset(ncfilename, 'r', format='NETCDF4')

    # Setup dimensions
    vars = ncfile2D.variables
    if debug:
        print([i for i in vars])
    lats, lons, Epoch = [ncfile2D[i] for i in ('LAT', 'LON', 'Epoch')]
    if debug:
        print([len(i) for i in (lats, lons, Epoch)])
    lats, lons, Epoch = [np.array(i) for i in (lats, lons, Epoch)]

    lats, lons = [list(sorted(set(i))) for i in (lats, lons)]
    # remove fill value. ( 9.969209968386869e+36 ) <= Improve this approach...
    [i.pop(-1) for i in (lons, lats)]

    # setup 3D NetCDF file
    ncfilename = ncfilename.split('.nc')[0]+'_3D.nc'
    ncfile = Dataset(ncfilename, 'w', format='NETCDF4')
    ncfile.createDimension('lat', len(lats))
    ncfile.createDimension('lon', len(lons))
    ncfile.createDimension('time', None)

    # Define the coordinate variables. They will hold the coordinate
    # information, that is, the latitudes and longitudes.
    time = ncfile.createVariable('time', 'f4', ('time',))
    lat = ncfile.createVariable('lat', 'f4', ('lat',))
    lon = ncfile.createVariable('lon', 'f4', ('lon',))

    # --- Add meta data
    # Assign units attributes to coordinate var data. This attaches a
    # text attribute to each of the coordinate variables, containing the
    # units.
    lat.units = 'degrees_north'
    lat.long_name = 'Latitude'
    lat.standard_name = 'Latitude'
    lat.axis = "Y"

    lon.units = 'degrees_east'
    lon.long_name = 'Longitude'
    lon.standard_name = 'Longitude'
    lon.axis = "X"

    time.units = 'seconds since 1970-01-01 00:00:00'
    time.calendar = "standard"
    time.standard_name = 'Time'
    time.axis = "T"

    # set global varibles
    ncfile.Description = 'planeflight output from '.format(wd)
    ncfile.Contact = 'Tomas Sherwen (ts551@york.ac.uk)'
#    ncfile.History = 'Created {}'.format(  time.ctime(time.time()) )
    ncfile.Grid = 'lat: {}-{}, lon: {}-{}'.format(lats[0], lats[-1],
                                                  lons[0], lons[-1])
#    ncfile.Temp_Res = "Hourly"
#    ncfile.SpatialCoverage='Global'

    # write data to coordinate vars.
    lon[:] = lons
    lat[:] = lats

    # Get unique timesteps
    timesteps = sorted(set(Epoch))
    # masked 1st value ( headers? )
    timesteps.pop(0)

    # set time dimension to timestep values
    time[:] = timesteps

    # select only 3D vars
    vars3D = get_3D_vars(vars)

    # --- Loop 3D species and create variables (with set dimensions)
    for var in vars3D:
        ncfile.createVariable(var, var2type(var), ('time', 'lat', 'lon'), )

    # close NetCDF
    ncfile.close()

    # ---  Loop through timesteps (epoch) and add to NetCDF
    # Loop over timesteps
    for t in timesteps:

        # open NetCDF in append mode
        ncfile = Dataset(ncfilename, 'a', format='NETCDF4')

        # get 1st and last indices for time stamp
        start, end = [(i.min(), i.max())
                      for i in np.where(ncfile2D.variables['Epoch'] == t)][0]

        # Extract Data for timestep & species
        for var in vars3D:

            data_ = ncfile2D.variables[var][start:end]
            lons_ = ncfile2D.variables['LON'][start:end]
            lats_ = ncfile2D.variables['LAT'][start:end]
            if debug:
                print((t, var, [i.shape for i in (data_, lats_, lons_)]))

            # stack data by LAT/LON to 3D array ( using pandas )
            df = DataFrame(data_, index=[lats_, lons_]).unstack()

            # add data to array
            ncfile.variables[var][timesteps.index(t)] = df.values

            # remove from memory
            del df, data_, lats_, lons_

        # Save out final NetCDF file
        ncfile.close()


def make_2D_subgroup_NetCDF(ncfilename, wd, debug=False):
    """ Create NetCDF of 2D arrays for all variables in 2D NetCDF file by
        subgroup ( e.g. multiple sites outputted for at same interval )

         Takes a table form NetCDF and build 2D arrays per unit 'POINT' 
         from sites in the given file.  
    """

    # --- Read existing & setup new NetCDF file
    ncfile2D = Dataset(ncfilename, 'r', format='NETCDF4')

    # Setup dimensions
    vars = ncfile2D.variables
    if debug:
        print([i for i in vars])
    # Get tim in files
    Epoch, LOC = [np.array(ncfile2D[i]) for i in ('Epoch', 'LOC')]
    if debug:
        print([len(i) for i in [Epoch]])

    Epoch, LOC = [list(sorted(set(i))) for i in (Epoch, LOC)]
    # remove 1st empty LOC entry
    LOC = LOC[1:]

    # setup 3D NetCDF file
    ncfilename = ncfilename.split('.nc')[0]+'_2D_by_site.nc'
    ncfile = Dataset(ncfilename, 'w', format='NETCDF4')
    ncfile.createDimension('time', None)

    # Define the coordinate variables. They will hold the coordinate
    # information, that is, the time
    time = ncfile.createVariable('time', 'f4', ('time',))

    # --- Add meta data
    # Assign units attributes to coordinate var data. This attaches a
    # text attribute to each of the coordinate variables, containing the
    # units.
    # Loop each unique variable and add meta data

    time.units = 'seconds since 1970-01-01 00:00:00'
    time.calendar = "standard"
    time.standard_name = 'Time'
    time.axis = "T"

    # set global varibles
#    ncfile.Description = 'planeflight output from '.format( wd )
    ncfile.Contact = 'Tomas Sherwen (ts551@york.ac.uk)'
#    ncfile.History = 'Created {}'.format(  time.ctime(time.time()) )
#    ncfile.Temp_Res = "Hourly"
#    ncfile.SpatialCoverage='Global'

    # Get unique timesteps
    timesteps = sorted(set(Epoch))
    # masked 1st value ( headers? )
    timesteps.pop(0)

    # set time dimension to timestep values
    time[:] = timesteps

    # ---  Add shared variables
    # select only 3D vars
    # ( variables are referred to as 3D as the function was written for grid
    # input. Here the table shape is altered but it remains 2D  )
    vars3D = get_site_description_vars(vars)

    # --- Loop 3D species and create variables (with set dimensions)
    for var in vars3D:
        # Loop sites and add to beginning of variable
        for site in LOC:
            ncfile.createVariable(site+'_'+var, var2type(var), ('time'), )

    # close NetCDF
    ncfile.close()

    # ---  Loop through timesteps (epoch) and add to NetCDF
    # Loop over timesteps in chunks
    chunk_size = 7*24  # Note for 365 days, this does make equally sized chunks
#    chunk_size = 24 # Note for 365 days, this does make equally sized chunks
    ts_chunks = chunks(timesteps, chunk_size)
    # check
    for ts_chunk in ts_chunks:

        # Open NetCDF in append mode
        ncfile = Dataset(ncfilename, 'a', format='NETCDF4')
        if debug:
            print((ts_chunk[0], timesteps[0]))
        # Get 1st and last indices for time stamp chunk
        start = np.where(ncfile2D.variables['Epoch'] == ts_chunk[0])
        start = np.array(start).min()
        end = np.where(ncfile2D.variables['Epoch'] == ts_chunk[-1])
        end = np.array(end).max()
        if debug:
            print((start, end, len(LOC)))

        # Extract Data for timestep & species
        for var in vars3D:
            # Select data for variable
            if debug:
                print(([len(i) for i in (timesteps, ts_chunks, LOC,
                                         ncfile2D.variables['Epoch'])],  LOC))
                print((start, start+len(LOC)+1+chunk_size,
                       start+((end-start)*len(LOC))))
            # the end index is 1 short
            # (as with all of the NetCDF indcies accessed, why?  )
            data_ = ncfile2D.variables[var][start:end+1]
            if debug:
                #                print data_
                print((np.array(data_).shape))
                print((ts_chunk[0], var, data_.shape))  # , data_

            # For each site
            for n, site in enumerate(LOC):
                # Select data for site
                if debug:
                    # , data_
                    print((ts_chunk[0], var, data_.shape,  site, LOC))

                # Add data to array
                if debug:
                    print((timesteps[:200]))
                    print((timesteps.index(ts_chunk[0]),
                           timesteps.index(ts_chunk[-1]), ts_chunk[0], ts_chunk[-1]))
                    print((ts_chunk, len(ts_chunk)))
                    print((np.array(ncfile.variables[site+'_'+var][
                        timesteps.index(ts_chunk[0]):
                        timesteps.index(ts_chunk[-1])+1
                    ]).shape))
                    print((np.array(data_).shape))
                    print((np.array(data_[n::len(LOC)]).shape))

                # Also add an extra index here
                ncfile.variables[site+'_'+var][
                    timesteps.index(ts_chunk[0]):
                    timesteps.index(ts_chunk[-1])+1
                ] = data_[n::len(LOC)]

            # remove from memory
            del data_

    # Save out final NetCDF file
    ncfile.close()


# Call to main with wd.
if __name__ == "__main__":
    # Get wd from command line arguement
    wd = sys.argv[1]
#    wd =' /work/home/ts551/data/all_model_simulations/iodine_runs/iGEOSChem_1.6_G5/run/'

    main(wd=wd, renumerated=renumerated,
         GRD_input_3D=GRD_input_3D,
         verbose=VERBOSE, debug=DEBUG)
