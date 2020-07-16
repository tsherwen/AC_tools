#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Generic functions for use with GEOS-Chem/Data Analysis.

Use help(<name of function>) to get details on a particular function.

NOTE(S):
 - This module is underdevelopment vestigial/inefficient code is being removed/updated.
 - Where external code is used credit is given.
"""
# - Required modules:
# I/O / Low level
import os
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pandas import DataFrame
# time
import time
import datetime as datetime
# math
from math import radians, sin, cos, asin, sqrt, pi, atan2
import xarray as xr

# The below imports need to be updated,
# imports should be specific and in individual functions
# import tms modules with shared functions
from . core import *
from . variables import *


def chunks(l, n):
    """
    Split list in chunks - useful for controlling memory usage
    """
    if n < 1:
        n = 1
    return [l[i:i + n] for i in range(0, len(l), n)]


def file_len(fname):
    """ Get length of file """
    return sum(1 for line in open(fname))


def myround(x, base=5, integer=True, round_up=False):
    """
    Round up values - mappable function for pandas processing
    NOTES:
     - credit: Alok Singhal
    """
#    round_up=True # Kludge - always set this to True

    # round to nearest base number
    rounded = base * round(float(x)/base)

    if round_up:
        # ensure rounded number is to nearest next base
        if rounded < x:
            rounded += base

    if integer:
        return int(rounded)
    else:
        return rounded


def counter_directory_contains_files(model_path, must_contain):
    """
    Count number of files in directory
    """
    model_ouput_file_counter = len(glob.glob1(model_path, must_contain))
    return model_ouput_file_counter


def replace_strs_in_files(wd, input_str, output_str, debug=False):
    """
    replace text in files
    """
    print((wd, input_str, output_str))
    for f in os.listdir(wd):
        if not f.startswith('.'):
            print(f)
            os.rename(wd + f, wd + f.replace(input_str, output_str))
            print((f.replace(input_str, output_str)))


def get_xy(Lon, Lat, lon_edges, lat_edges, debug=False):
    """
    Takes lon,lat values for a point (e.g. station) and arrays of longitude
    and latitudes indicating edges of gridboxes.

    Notes
    -------
     - Returns the index of the gridbox corresponding to this point.
     - Credit: Eric Sofen
     - Could be easily extended to return data values corresponding to points.
    """
    hasobs, lonrange, latrange = np.histogram2d([Lon], [Lat],
                                                [lon_edges, lat_edges]
                                                )
    gridindx, gridindy = np.where(hasobs >= 1)
    if not gridindx:
        if debug:
            pstr = 'Lat, lon outside of x,y range.  Assigning -1 for {} {}'
            print(pstr.format(Lon, Lat))
        return -1, -1
    else:
        # Print Lon, Lat, gridindx, gridindy
        return gridindx[0], gridindy[0]


def plot2pdf(title='new_plot', fig=None, rasterized=True, dpi=320,
             justHH=False, no_dstr=True, save2png=True, wd=None,
             save2eps=False, transparent=True,
             bbox_inches=None, debug=False):
    """
    Save figures (e.g. matplotlib) to pdf file
    """
    # Set save directory ( using default directory dictionary )
    if isinstance(wd, type(None)):
        wd = './'

    # Set pdf name
    if justHH:
        date_str = time.strftime("%y_%m_%d_%H")
    else:
        date_str = time.strftime("%y_%m_%d_%H_%M")
    if no_dstr:
        npdf = wd+title
    else:
        npdf = wd+date_str+'_'+title

    # Setup pdf
    pdf = PdfPages(npdf + '.pdf')

    # Rasterise to save space?
    if rasterized:
        plt.gcf().set_rasterized(True)

    # Save and close
    file_extension = 'PDF'
    pdf.savefig(dpi=dpi, transparent=transparent, bbox_inches=bbox_inches)
    pdf.close()

    # Also export to png and eps?
    if save2png:
        file_extension += '/PDF'
        plt.savefig(npdf+'.png', format='png', dpi=dpi,
                    transparent=transparent, bbox_inches=bbox_inches)
    if save2eps:
        file_extension += '/EPS'
        plt.savefig(npdf+'.eps', format='eps', dpi=dpi,
                    transparent=transparent, bbox_inches=bbox_inches)
    print((file_extension+' saved & Closed as/at: ', npdf))


def plot2pdfmulti(pdf=None, title='new_plot', rasterized=True, wd=None,
                  dpi=320, open=False, close=False, justHH=False,
                  no_dstr=True, tight=False):
    """
    Save figures (e.g. matplotlib) to pdf file with multiple pages
    """
    # Set save directory ( using default directory dictionary )
    if isinstance(wd, type(None)):
        wd = './'
    # Set pdf name
    if justHH and (not no_dstr):
        date_str = time.strftime("%y_%m_%d_%H")
    if not no_dstr:
        date_str = time.strftime("%y_%m_%d_%H_%M")
    if no_dstr:
        npdf = wd+title+'.pdf'
    else:
        npdf = wd+date_str+'_'+title+'.pdf'
    # If 1st call ( open ==True), setup pdf
    if open:
        pdf = PdfPages(npdf)
        print(('pdf opened @: {}'.format(npdf)))
        return pdf
    # Rasterise to save space?
    if rasterized:
        plt.gcf().set_rasterized(True)
    if tight:
        plt.tight_layout()
    # Save and close or keep open to allow additions of plots
    if close:
        pdf.close()
        print(('PDF saved & Closed as/at: ', npdf))
    else:
        pdf.savefig(dpi=dpi)
        print(('pdf is still open @: {}'.format(npdf)))


def sort_sites_by_lat(sites):
    """
    Order given list of GAW sties by latitudes
    """
    # Get info
    vars = [gaw_2_loc(s) for s in sites]  # lat, lon, alt, TZ
    # Sort by lat, index orginal sites list and return
    lats = [i[0] for i in vars]
    slats = sorted(lats)[::-1]
    return [sites[i] for i in [lats.index(ii) for ii in slats]]


def find_nearest(array, value):
    """
    Find nearest number in array to given value.

    Notes
    -------
     - credit: Adapted from HappyLeapSecond's Stackoverflow answer.
   http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    idx = (np.abs(array-value)).argmin()
    return idx


def get_suffix(n):
    """  Add the appropriate suffix (th/st/rd) to any number given  """
    def ordinal(n): return "%d%s" % (
        n, "tsnrhtdd"[(n//10 % 10 != 1)*(n % 10 < 4)*n % 10::4])
    return ordinal(n)


def get_shortest_in(needle, haystack, r_distance=False):
    """
    needle is a single (lat,long) tuple. haystack is a numpy array to find the point in
    that has the shortest distance to needle

    NOTES:
    -------
     - adapted from stackoverflow (Credit: jterrace):
http://stackoverflow.com/questions/6656475/python-speeding-up-geographic-comparison
    """
    # Set Earth's radius
    earth_radius_miles = 3956.0

    # convert to radians
    dlat = np.radians(haystack[:, 0]) - radians(needle[0])
    dlon = np.radians(haystack[:, 1]) - radians(needle[1])

    # get the vector
    a = np.square(np.sin(dlat/2.0)) + cos(radians(needle[0])) * np.cos(
        np.radians(haystack[:, 0])) * np.square(np.sin(dlon/2.0))

    # convert to distance
    great_circle_distance = 2 * np.arcsin(np.minimum(np.sqrt(a),
                                                     np.repeat(1, len(a))))
    d = earth_radius_miles * great_circle_distance

    if r_distance:
        return np.min(d)
    # return the index
    else:
        return list(d).index(np.min(d))


def gen_log_space(limit, n):
    """
    Get logarithmically spaced integers

    Notes
    -------
    - credit: Avaris
 http://stackoverflow.com/questions/12418234/logarithmically-spaced-integers
    """
    result = [1]
    if n > 1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result) < n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # Safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # Problem! same integer. we need to find next_value
            # by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining
            # values will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    return np.array([round(x)-1 for x in result], dtype=np.uint64)


def get_linear_ODR(x=None, y=None, maxit=5000, beta0=(0, 1),
                   xvalues=None, return_model=True, debug=False,
                   verbose=False):
    """
    Wrapper to run orthogonal distance regression (ODR) for arrays of x and y

    Parameters
    ----------
    x (np.array): np.array of x values
    y (np.array): np.array of y values
    maxit (int): integer specifying the maximum number of iterations to perform.
                 For first runs, maxit is the total number of iterations performed
                 and defaults to 50 (in core implimentation, here =1000). For restarts,
                 maxit is the number of additional iterations to perform and defaults
                 to 10.
    beta0 (array_like): starting parameters to use (e.g. y = mx+ c for linear model).
                        default is a 1:1 line between x and y, with an intercept of 0.
    xvalues (array): values of x to predict over with model
    return_model (bool): return the model object
    verbose (bool): print out extra information
    debug (bool): print out debugging information

    Returns
    -------
    (myodr.run or predicted X & Y)

    Notes
    -------
     - The function is setup to either return the ODR model or
     - For original implementation in Fortran, please see the
    https://docs.scipy.org/doc/external/odrpack_guide.pdf
     - Wrapper adapted from example in manual
    https://docs.scipy.org/doc/scipy/reference/odr.html
    """
    import scipy.odr
    # Setup linear model to fit
    def f(B, x):
        '''Linear function y = m*x + b'''
        # B is a vector of the parameters.
        # x is an array of the current x values.
        # x is in the same format as the x passed to Data or RealData.
        #
        # Return an array in the same format as y passed to Data or RealData.
        return B[0]*x + B[1]
    # Create a model
    linear = scipy.odr.Model(f)
    # Create a Data or RealData instance
    mydata = scipy.odr.Data(x, y)
    # Instantiate ODR with your data, model and initial parameter estimate
    myodr = scipy.odr.ODR(mydata, linear, beta0,  maxit=maxit, job=None)
    # Run the fit:
    myoutput = myodr.run()
    # Examine output
    if verbose:
        myoutput.pprint()
    # Return either the entire outputted object, or predicted X and Y values.
    if return_model:
        return myoutput
    else:
        if isinstance(xvalues, type(None)):
            xvalues = np.arange(min(x), max(x), (max(x)-min(x))/100.)
        return xvalues, f(myoutput.beta, xvalues)


def convert_ug_per_m3_2_ppbv(data=None,  spec='O3', rtn_units=False,
                             units='ug m$^{-3}$'):
    """
    Converts units of ugm^-3 to ppbv for a given species assuming SATP
    """
    # --- Get constants
    RMM_air = constants('RMM_air')  # g/mol
    # assume standard air density
    # At sea level and at 15 °C air has a density of approximately 1.225 kg/m3
    # (0.001225 g/cm3, 0.0023769 slug/ft3, 0.0765 lbm/ft3) according to
    # ISA (International Standard Atmosphere).
    AIRDEN = 0.001225  # g/cm3
    # moles per cm3
    #  (1/(g/mol)) = (mol/g) ; (mol/g) * (g/cm3) = mol/cm3
    MOLS = (1/RMM_air) * AIRDEN

    # --- Convert
    # convert ug/m3 to ppbv
    # get g per cm3, ( instead of ug/m3)
    data = data / 1E6 / 1E6
    # get mol/cm3 (mass/ RMM ) = ( ug/m3  /  g/mol )
    data = data/species_mass(spec)
    # convert to ppb
    data = data/MOLS * 1E9
    # update unit string
    units = 'ppbv'
    if rtn_units:
        return data, units
    else:
        return data


def convert_mg_per_m3_2_ppbv(data=None,  spec='O3', rtn_units=False,
                             units='mg m$^{-3}$'):
    """
    Converts units of ugm^-3 to ppbv for a given species assuming SATP
    """
    # --- Get constants
    RMM_air = constants('RMM_air')  # g/mol
    # assume standard air density
    # At sea level and at 15 °C air has a density of approximately 1.225 kg/m3
    # (0.001225 g/cm3, 0.0023769 slug/ft3, 0.0765 lbm/ft3) according to
    # ISA (International Standard Atmosphere).
    AIRDEN = 0.001225  # g/cm3
    # moles per cm3
    #  (1/(g/mol)) = (mol/g) ; (mol/g) * (g/cm3) = mol/cm3
    MOLS = (1/RMM_air) * AIRDEN

    # --- Convert
    # convert mg/m3 to ppbv
    # get g per cm3, ( instead of mg/m3)
    data = data / 1E6 / 1E3
    # get mol/cm3 (mass/ RMM ) =  (mg/m3  /  g/mol )
    data = data/species_mass(spec)
    # convert to ppb
    data = data/MOLS * 1E9
    # update unit string
    units = 'ppbv'
    if rtn_units:
        return data, units
    else:
        return data


def get_2D_solartime_array4_date(date=None, ncfile=None, res='4x5',
                                 lons=None, lats=None, varname='SolarTime',
                                 debug=False):
    """
    Creates 2D (lon,lat) masked (1=Masked) for nighttime for a given list of
    dates

    Parameters
    -------
    date (datetime): date to use (UTC)
    mask_daytime (bool): mask daytime instead of nightime
    ncfile (str): location to netCDF file - not implemented...
    res (str): resolution, if using resolutions listed in get_latlonalt4res
    lons (array): array of longditudes (optional)
    lats (array): array of lattiudes (optional)

    Returns
    -------
    (np.array)

    ncfile (NetCDF file): NetCDF file to extract lat and lon metadata from

    Notes
    -----
     - if ncfile provide programme will work for that grid.
    """
    # Astronomical math
    import ephem
    from ephem import AlwaysUpError, NeverUpError
    # And functions in other AC_tools modules
    from .AC_time import add_days, add_hrs
    logging.info('get_2D_solartime_array4_dates called for {}'.format(date))

    # Profile function...
    if debug:
        start_time = time.time()

    # --- Local variables?
    # reference data for ephem (number of days since noon on 1899 December 31)
    ref_date = datetime.datetime(1899, 12, 31, 12)

    # --- Get LON and LAT variables (if lons/lats not provdided)
    if any([not isinstance(i, type(None)) for i in (lats, lons)]):
        pass
    else:
        if isinstance(ncfile, type(None)):
            # extract from refence files
            lons, lats, NIU = get_latlonalt4res(res=res)
        else:
            # TODO - allow any lat, lon grid to be used by taking input lats and
            # lons from ncfile file/arguments.
            print('Not implemented')
            sys.exit()

    # --- setup function to get solartime based on date, lat and lon
    def _solartime(lon, lat, date=date, rtn_as_epoch=True):
        """
        Get solartime for location (lat, lon) and date

        Returns
        ---
        (float) Epoch time
        """
        # --- get sunrise and sunset for location
        o = ephem.Observer()
        # Set lat (decimal?), lon (decimal?), and date (UTC)
        o.lat = str(lat)
        o.long = str(lon)
        o.date = date
        # Planetary body
        s = ephem.Sun()
        # Compute sun vs observer
        s.compute(o)
        # below code was adapted from stackoverflow (Credit: J.F. Sebastian)
        # http://stackoverflow.com/questions/13314626/local-solar-time-function-from-utc-and-longitude
        # Sidereal time == ra (right ascension) is the highest point (noon)
        hour_angle = o.sidereal_time() - s.ra
        s_time = ephem.hours(
            hour_angle + ephem.hours('12:00')).norm  # norm for 24h
        # ephem.hours is a float number that represents an angle in radians
        # and converts to/from a string as "hh:mm:ss.ff".
        s_time = "%s" % s_time
        if len(s_time) != 11:
            s_time = '0'+s_time
        # return as datetime
        s_time = datetime.datetime.strptime(s_time[:8], '%H:%M:%S')
        if rtn_as_epoch:
            s_time = unix_time(s_time)
        return s_time

    # --- Setup an unstack(ed) pandas dataframe to contain masked values
    # Use list comprehension to setup list of indices for lat and lon
    # Better way of doing this? (e.g. pd.melt?)
    ind_lat_lons_list = [[lon_, lat_] for lat_ in lats for lon_ in lons]
    # Make this into a pd.DataFrame and label columns.
    df = pd.DataFrame(ind_lat_lons_list)
    df.columns = ['lons', 'lats']
    # Apply function to calculate mask value
    df[varname] = df.apply(lambda x: _solartime(x['lons'], x['lats']), axis=1)
    # Re-index by lat and lon
    df = pd.DataFrame(df[varname].values, index=[df['lats'], df['lons']])
    # Unstack and return just as array
    return df.unstack().values


def save_2D_arrays_to_3DNetCDF(ars=None, dates=None, res='4x5', lons=None,
                               lats=None, varname='MASK', Description=None,
                               Contact=None,
                               filename='misc_output', var_type='f8',
                               debug=False):
    """
    makes a NetCDF from a list of dates and list of (lon, lat) arrays

    Parameters
    -------
    ars (list of np.array): list of (lon, lat) arrays
    dates (list of datetime.datetime)
    res (str): reoslution (opitional )
    lons (list): lons for use for NetCDF cooridinate variables
    lats (list): lats for use for NetCDF cooridinate variables
    filename (str): name for output netCDF file
    varname (str): name for variable in NetCDF
    var_type (str): variable type (e.g. 'f8' (64-bit floating point))

    Returns
    -------
    (None)

    Notes
    -----
     - needs updating to take non-standard reoslutions
     (e.g. those in variables)
    """
    logging.info('save_2D_arrays_to_3DNetCDF called')
    from .AC_time import unix_time, dt64_2_dt
    print((locals()))

    # Local Settings
    ncfilename = '{}_{}.nc'.format(filename, res)
    # Get lons and lats...
    if any([isinstance(i, type(None)) for i in (lats, lons)]):
        lons, lats, NIU = get_latlonalt4res(res=res)
        logging.debug('Using offline lons/lats, as either lats/lons==None')
#    else:
#        print 'WARNING: non-standard lats/lons not implemented!!! - TODO. '
#        sys.exit()

    # - Setup new file to save data to
    # write file
    logging.debug('setting up new NetCDF file: {}'.format(ncfilename))
    ncfile = Dataset(ncfilename, 'w', format='NETCDF4')
    ncfile.createDimension('lat', len(lats))
    ncfile.createDimension('lon', len(lons))
    ncfile.createDimension('time', None)

    # Define the coordinate variables. They will hold the coordinate
    # information, that is, the latitudes and longitudes.
    time = ncfile.createVariable('time', 'f4', ('time',))
    lat = ncfile.createVariable('lat', 'f4', ('lat',))
    lon = ncfile.createVariable('lon', 'f4', ('lon',))

    # - Add meta data for coordinates
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

    # - Set global variables
    if not isinstance(Description, type(None)):
        ncfile.Description = Description
    if not isinstance(Contact, type(None)):
        ncfile.Contact = Contact
    ncfile.Grid = 'lat: {}-{}, lon: {}-{}'.format(lats[0], lats[-1],
                                                  lons[0], lons[-1])

    # Write values to coordinate variables (lat, lon)
    lon[:] = lons
    lat[:] = lats

    # - Setup time dimension/variables (as epoch)
    # Convert to Epoch time
    def format(x): return unix_time(x)
    df = pd.DataFrame({'Datetime': dates})
    df['Epoch'] = df['Datetime'].map(format).astype('i8')
    # Store a copy of dates at datetime.datetime, then remove from DataFrame
    dt_dates = dt64_2_dt(df['Datetime'].values.copy())
    del df['Datetime']
    dates = df['Epoch'].values
    # Assign to time variable
    time[:] = dates

    # - Create new NetCDF variable (as f8) with common dimensions
    # (e.g. 'f8' = 64-bit floating point, 'i8'=(64-bit singed integer) )
    ncfile.createVariable(varname, var_type, ('time', 'lat', 'lon'), )
    # Close NetCDF
    ncfile.close()

    # - Now open and add data in append mode
    for n, date in enumerate(dt_dates):
        # Print out debugging text to AC_tools.log
        if debug:
            fmtstr = "%Y%d%m %T"
            logging.debug('saving array date:{}'.format(date.strftime(fmtstr)))
            pcent = (float(n)+1)/float(len(dt_dates))*100
            logging.debug('array #: {} (% complete: {:.1f})'.format(n, pcent))
        # Open NetCDF in append mode
        ncfile = Dataset(ncfilename, 'a', format='NETCDF4')
        # Add data to array
        try:
            ncfile.variables[varname][n] = ars[n]
        except ValueError:
            err_msg = '>total size of new array must be unchanged<'
            err_msg += '(new arr shape {})'.format(str(ars[n].shape))
            print(err_msg)
            logging.info(err_msg)
            sys.exit()
    # Close NetCDF
    ncfile.close()
    if debug:
        logging.debug('saved NetCDF file:{}'.format(ncfilename))


def interpolate_sparse_grid2value(Y_CORDS=None, X_CORDS=None,
                                  X=None, Y=None, XYarray=None, buffer_CORDS=5,
                                  verbose=True, debug=False):
    """
    Get an interpolated value for a location (X,Y) from surrounding array values

    Parameters
    -------
    X_CORDS (np.array): coordinate values for X axis of 2D array (e.g. lon)
    Y_CORDS (np.array): coordinate values for Y axis of 2D array (e.g. lat)
    X (float): value of interest (in same terms at X_CORDS, e.g. lon)
    Y (float): value of interest (in same terms at Y_CORDS, e.g. lat)
    XY (array): array of values with shape (X_CORDS, Y_CORDS)
    buffer_CORDS (int): number of units of X_CORDS and Y_CORDS to interpolate
        around (great this value, greater the cost.)
    verbose (bool): print out extra information
    debug (bool): print out debugging information

    Returns
    -------
    (float)
    """
    # X_CORDS=file_lonc; Y_CORDS=file_latc; X=lon; Y=lat; XYarray=file_data; buffer_CORDS=10
    import scipy.interpolate as interpolate
    # ---  Select a sub grid.
    # WARNING THIS APPRAOCH WILL NOT WORK NEAR ARRAY EDGES!!!
    # (TODO: add functionality for this.)
    # Get indices buffer sub-selection of grid.
#    assert X_CORDS[-1] > X_CORDS[0], 'This -180=>180 and
    if X_CORDS[0] < X_CORDS[-1]:
        s_low_X_ind = find_nearest_value(X_CORDS, X-buffer_CORDS)
        s_high_X_ind = find_nearest_value(X_CORDS, X+buffer_CORDS)
    else:  # Y_CORDS[0] > Y_CORDS[-1]
        s_low_X_ind = find_nearest_value(X_CORDS, X+buffer_CORDS)
        s_high_X_ind = find_nearest_value(X_CORDS, X-buffer_CORDS)
        XYarray[::-1, :]
    # WARNING: This approach assumes grid -90=>90 (but flips slice if not)
    if Y_CORDS[0] < Y_CORDS[-1]:
        s_low_Y_ind = find_nearest_value(Y_CORDS, Y-buffer_CORDS)
        s_high_Y_ind = find_nearest_value(Y_CORDS, Y+buffer_CORDS)
    else:  # Y_CORDS[0] > Y_CORDS[-1]
        s_low_Y_ind = find_nearest_value(Y_CORDS, Y+buffer_CORDS)
        s_high_Y_ind = find_nearest_value(Y_CORDS, Y-buffer_CORDS)
        XYarray[:, ::-1]
    # Print out values in use
    if verbose:
        prt_str = 'Y={}, subrange=({}(ind={}),{}(ind={}))'
        print(prt_str.format(Y, Y_CORDS[s_low_Y_ind], s_low_Y_ind,
                             Y_CORDS[s_high_Y_ind], s_high_Y_ind))
        prt_str = 'X={}, subrange=({}(ind={}),{}(ind={}))'
        print(prt_str.format(X, X_CORDS[s_low_X_ind], s_low_X_ind,
                             X_CORDS[s_high_X_ind], s_high_X_ind))
    # Select sub array and get coordinate axis
    subXY = XYarray[s_low_X_ind:s_high_X_ind, s_low_Y_ind:s_high_Y_ind]
    subX = X_CORDS[s_low_X_ind:s_high_X_ind]
    subY = Y_CORDS[s_low_Y_ind:s_high_Y_ind]
    # Debug (?) by showing 2D grid prior to interpolation
    if debug:
        print(([i.shape for i in (subX, subY, subXY)], XYarray.shape))
        plt.pcolor(subX, subY, subXY.T)
        plt.colorbar()
        plt.show()

    # ---  interpolate over subgrid.
    M = subXY
    rr, cc = np.meshgrid(subX, subY)
    # fill masked values with nans
    M = np.ma.filled(M, fill_value=np.nan)
    # only consider non nan values as values to interpolate with
    vals = ~np.isnan(M)
    if debug:
        print(vals)
    # interpolate
    f = interpolate.Rbf(rr[vals], cc[vals], M[vals], function='linear')
    # extract interpolation...
    interpolated = f(rr, cc)

    # ---  Select value of interest
    # Debug (?) by showing 2D grid post to interpolation
    if debug:
        print((interpolated.shape))
        plt.pcolor(subX, subY, interpolated.T)
        plt.colorbar()
        plt.show()

    # indix value from subgrid?
    Y_ind = find_nearest_value(subY, Y)
    X_ind = find_nearest_value(subX, X)

    return interpolated[X_ind, Y_ind]


def split_NetCDF_by_month(folder=None, filename=None, ext_str='',
                          file_prefix='ts_ctm'):
    """
    Split a NetCDF file by month into new NetCDF files using xarray

    Parameters
    -------
    folder (str): the directory to search for files in
    filename (str): the NetCDF filename (e.g. ctm.nc)
    file_prefix (str): prefix to attach to new saved file
    ext_str (str): extra string for new filenames
    """
    import xarray as xr
    # --- Open data
    ds = xr.open_dataset(folder+filename)
    months = list(sorted(set(ds['time.month'].values)))

    # --- Loop months
    for month_ in months:

        # Custom mask
        def is_month(month):
            return (month == month_)
        # Make sure time is the dimension not module
        time = ds.time
        # Now select for month
        ds_tmp = ds.sel(time=is_month(ds['time.month']))

        # --- Save as NetCDF
        # Name of new NetCDF?
        year_ = list(set(ds_tmp['time.year'].values))[0]
        file2save = '{}_{}_{}_{:0>2}.nc'.format(file_prefix, ext_str, year_,
                                                str(month_))
        logging.info('saving month NetCDF as: {}'.format(file2save))
        # Save the file...
        ds_tmp.to_netcdf(folder+file2save)
        # Delete temporary dataset
        del ds_tmp


def get_2D_df_of_lon_lats_and_time(res='4x5', df_lar_var='lat',
                                   df_lon_var='lon', df_time_var='month',
                                   add_all_months=False,
                                   lons=None, lats=None, verbose=True,
                                   month=9):
    """ stack a 2D table (DataFrame) of lat/lon coords """
    # Get lon and lat resolution (Add other ways to get lat and lon here...
    lons_not_provided = isinstance(lons, type(None))
    lats_not_provided = isinstance(lats, type(None))
    if (lons_not_provided) or (lats_not_provided):
        try:
            assert(type(res) == str), 'Resolution must be a string!'
            lons, lats, NIU = get_latlonalt4res(res=res)
        except:
            print('please provide lons/lats or a res in get_latlonalt4res')
    else:
        pass
    # Make Table like array
    b = np.zeros((len(lons), len(lats)))
    df = pd.DataFrame(b)
    # Use lats and lons as labels for index and columns
    df.index = lons
    df.columns = lats
    # Stack, then reset index to obtain table structure
    df = df.stack()
    df = df.reset_index()
    # Set column headers
    df.columns = df_lon_var, df_lar_var, df_time_var
    # Add time dims...
    if add_all_months:
        dfs = []
        for month in np.arange(1, 13):
            df_tmp = df.copy()
            df_tmp[df_time_var] = month
            dfs.append(df_tmp)
        df = pd.concat(dfs)
    else:  # Just select a single month (September is default )
        print('WARNING: Only September considered')
        df[df_time_var] = month
    return df


def get_vars_from_line_printed_in_txt_file(filename=None, folder=None,
                                           prefix=' HOUR:', var_names=None,
                                           type4var_names=None):
    """
    Get variables from a line printed to non-binary file with a given prefix

    Parameters
    ----------
    filename (Str): name of non binary file (e.g. geos.log)
    folder (str): name of directory where file ("filename") is located
    prefix (str): the string that lines containing data begin with
    var_names (list): optional. names for variables in line (must be # in line)
    type4var_names (dict): dictionary of type to convert variables to

    Returns
    -------
    (pd.DataFrame)

    Notes
    ----------
     - Can be used for debugging or tabulating output. e.g. lines of
     grid index locations and values printed to GEOS-chem's geos.log file
    """
    import pandas as pd
    import sys
    # --- Local vars
    lines_with_var = []
    # --- Open file and loop vars
    with open(folder+filename) as f:

        for line in f:
            #			if prefix in line:
            #            print line
            if line.startswith(prefix):
                var_ = line.split()[1:]
                if not isinstance(var_names, type(None)):
                    try:
                        assert len(var_) == len(var_names)
                    except AssertionError:
                        print((len(var_), var_))
                        sys.exit()
                # Save data until later
                lines_with_var.append(var_)
    # If lines found, return as a DataFrame
    if len(lines_with_var) > 0:
        # Make DataFrame
        df = pd.DataFrame(lines_with_var)
        # Apply provided names to columns?
        if not isinstance(var_names, type(None)):
            df.columns = var_names
        # Covert data type of variable?
        if not isinstance(type4var_names, type(None)):
            for key_ in list(type4var_names.keys()):
                try:
                    df[key_] = df[key_].astype(type4var_names[key_])
                except KeyError:
                    err_str = 'key_ ({}) not in df'.format(key_)
                    logging.info(err_str)

        return df
    else:
        err_str = 'No lines with prefix ({})'.format(prefix)
        logging.info(err_str)


def mk_spatial_dataset_from_longform_df(df=None, LatVar='lat', LonVar='lon',
                                        unstack=True, attrs={},
                                        VarName='New_Variable'):
    """
    Make a xr.dataset from a provided 2D array

    Parameters
    ----------
    df (pd.DataFrame): dataframe containing coordinates and data
    VarName (str): name of the variable that holds the data in the df
    unstack (bool): convert table/long form data into a 2D format
    attrs (dict): attributes to add to the variable (e.g. units)
    LatVar, LonVar (str): variables names for lat and lon in df

    Returns
    -------
    (xr.dataset)
    """
    # Get coordinate values
    lons = df[LonVar].values
    lats = df[LatVar].values
    data = df[VarName].values
    # Convert the dataset from long to 2D form
    if unstack:
        df = pd.DataFrame(data, index=[lats, lons]).unstack()
    # Get coordinate values
    lons4ds = list(df.columns.levels[1])
    lats4ds = list(df.index)
    # Extract the 2D data
    arr = df.values
    # Now make a dataset
    ds = xr.Dataset(data_vars={VarName: (['lat', 'lon', ], arr)},
                    coords={'lat': lats4ds, 'lon': lons4ds, }
                    )
    # Set the attributes for the new variable
    ds[VarName].attrs = attrs
    # Add the coordinate standard information
    ds.lat.attrs = {"axis": 'Y', 'long_name': "latitude",
                    "standard_name": "latitude"}
    ds.lon.attrs = {"axis": 'X', 'long_name': "longitude",
                    "standard_name": "longitude", }
    # Add core attributes
    return ds


def rm_spaces_and_chars_from_str(input_str, remove_slashes=True,
                                 replace_brackets=True, replace_quotes=True,
                                 replace_dots=True, replace_colons=True,
                                 remove_plus=True, swap_pcent=True,
                                 replace_braces=True):
    """
    Remove the spaces and species characters from strings
    """
    input_str = input_str.replace(' ', '_')
    if replace_colons:
        input_str = input_str.replace(':', '_')
    if replace_brackets:
        input_str = input_str.replace('(', '_')
        input_str = input_str.replace(')', '_')
    if replace_braces:
        input_str = input_str.replace('{', '_')
        input_str = input_str.replace('}', '_')
    if replace_quotes:
        input_str = input_str.replace("'", '_')
    if replace_dots:
        input_str = input_str.replace(".", '_')
    if remove_slashes:
        input_str = input_str.replace("\\", '_')
        input_str = input_str.replace("/", '_')
    if remove_plus:
        input_str = input_str.replace("+", '_plus_')
    if swap_pcent:
        input_str = input_str.replace("%", 'pcent')
    return input_str


def mk_folder(folder=None, verbose=False):
    """
    Make a folder at provided location string
    """
    if os.path.isdir(folder):
        if verbose:
            print('WARNING: folder already exists ({})'.format(folder))
    else:
        os.makedirs(folder)
        if verbose:
            print('Created folder for data ({})'.format(folder))


def rm_file(folder=None, filename=None, verbose=False, debug=False):
    """
    Remove a specific file from the working directory or a provided folder
    """
    # Use the current directory if no folder set
    if isinstance(folder, type(None)):
        folder = './'
    # Try and remove the file if present
    try:
        os.remove(folder+filename)
        if verbose:
            print('Removed file: {}'.format(folder+filename))
    except FileNotFoundError:
        pstr = "WARNING: File was not removed as it doesn't exit: {}"
        if verbose:
            print( pstr.format(folder+filename) )


def get_stats_on_files_in_folder_as_dict(folder=None):
    """
    Get statistics on files in a folder
    """
    # Setup a dictionary to hold info on the files
    d = {}
    # list all files in the folder
    files = glob.glob( '{}/*'.format(folder) )
    nfiles = len(files)
    d['#'] = nfiles
    if len(files) == 0:
        print('WARNING: No files found in directory')
        return d
    else:
        # get the sizes of the files
        sizes = [os.path.getsize(i) for i in files]
        # get the extensions of the files
        extensions = [i.split('.')[-1] for i in files]
        extensions = list(set(extensions))
        d['Extensions'] = extensions
        # Smallest sixe
        d['Min. size'] = np.min(sizes)
        # largest size
        d['Max. size'] = np.max(sizes)
        # mean size
        d['Mean size'] = AC.myround(np.mean(sizes), base=0.1, integer=False)
        # mediun size
        d['Medium size'] = np.percentile(sizes, 50)
        return d


def get_avg_2D_conc_of_X_weighted_by_Y(ds, Xvar=None, Yvar='AREA', Yda=None):
    """
    Get average 2D value of X (e.g. O3 v/v) weighted by Y (e.g. surface area)
    """
    # Get Y as data array if not provided
    if isinstance(Yda, type(None)):
        Yda = ds[Yvar]
    # return X as a weighted average of Y
    return float((ds[Xvar]*Yda).sum() / Yda.sum())


def write_lines2txt_file(lines, folder=None, filename=None):
    """
    Write out provided lines to text file
    """
    # Write lines to text file
    with open(folder+filename, "w") as f:
        for line in lines:
            f.write(line)


def read_lines_from_txt_file(folder=None, filename=None):
    """
    Read and return lines from provided text file
    """
    # Write lines to text file
    with open(folder+filename, "r") as f:
        lines = f.readlines()
    return lines


def calc_4D_idx_in_ds(ds=None, df=None, LonVar='lon', LatVar='lat',
                      TimeVar='time', AltVar='hPa', dsAltVar='lev',
                      dsLonVar='lon', dsLatVar='lat', dsTimeVar='time',
                      ds_lat=None, ds_lon=None, ds_hPa=None, ds_time=None,
                      debug=False):
    """
    Calculate the 4D indexes of dataframe locations in a dataset

    Parameters
    ----------
    df (pd.DataFrame): dataframe containing coordinates (lon, lat, alt, time)
    ds (xr.dataset): dataset to calculate indexes within
    LonVar (str): Variable name in DataFrame for Longitude
    dsLonVar (str): Variable name in dataset for Longitude
    LatVar (str): Variable name in DataFrame for latitude
    dsLatVar (str): Variable name in dataset for latitude
    AltVar (str): Variable name in DataFrame for pressure (hPa)
    dsAltVar (str): Variable name in dataset for pressure (hPa)
    TimeVar (str): Variable name in DataFrame for time
    dsTimeVar (str): Variable name in dataset for time
    ds_hPa (np.array): array of pressure coordinate (hPa)
    ds_lat (np.array): array of latitude coordinate in dataset
    ds_lon (np.array):  array of longitude coordinate in dataset
    ds_time (np.array): array of time coordinate in dataset

    Returns
    -------
    (dict)
    """
    # Get arrays of the coordinate variables in the dataset
    if isinstance(ds_lat, type(None)):
        ds_lat = ds[dsLatVar].values
    if isinstance(ds_lon, type(None)):
        ds_lon = ds[dsLonVar].values
    if isinstance(ds_hPa, type(None)):
        ds_hPa = ds[dsAltVar].values
    if isinstance(ds_time, type(None)):
        ds_time = ds[dsTimeVar].values
    # Calculate the index individually by coordinate
    lat_idx = [find_nearest(ds_lat, i) for i in df[LatVar].values]
    lon_idx = [find_nearest(ds_lon, i) for i in df[LonVar].values]
    hPa_idx = [find_nearest(ds_hPa, i) for i in df[AltVar].values]
    time_idx = [find_nearest(ds_time, i) for i in df.index.values]
    # Return a dictionary of the values
    return {LatVar:lat_idx, LonVar:lon_idx, TimeVar:time_idx, AltVar:hPa_idx}


def save_ds2disk_then_reload(ds, savename='TEMP_NetCDF.nc', folder='./',
                             delete_ds=True):
    """
    Save a dataset (subset) to disk locally, then reload into memory
    """
    # Save the dataset to disk
    ds.to_netcdf(folder+savename)
    # Delete the dataset
    if delete_ds:
        del ds
    # Reload the dataset and return it
    return xr.open_dataset(folder+savename)



