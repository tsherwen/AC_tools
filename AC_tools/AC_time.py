#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Time processing functions for use with GEOS-Chem/Data analysis

Use help(<name of function>) to get details on a particular function.

Notes
-----
 - This module is underdevelopment vestigial/inefficient code is being removed/updated.
 - Where external code is used, credit is given.
"""
import logging
import numpy as np
import pandas as pd
import time
import calendar
import datetime as datetime
from datetime import datetime as datetime_
import sys
# Attempt to import ephem if installed
if sys.version_info.major < 3:
    try:
        import ephem
    except ImportError:
        print('ephem package not installed')


def get_day_fraction(date):
    """
    Get day fraction from a datetime object

    Notes
    -----
     - for working with numpy arrays of datetimes, instead of pandas dataframes
    """
    secs = (date.hour * 60.*60.)+(date.minute*60.)+(date.second)
    dsecs = 24.*60.*60.
    return secs/dsecs


def dt64_2_dt(dt64, RtnAsArray=True):
    """
    Convert numpy.datetime64 to datetime.datetime (assuming UTC )

    Parameters
    -----
    dt64 (numpy.datetime64): datetime to convert

    Notes
    -----
     - TODO: Convert this to work as a lamdba function for scalability
    """
    ns = 1e-9  # number of seconds in a nanosecond
    dt = [datetime_.utcfromtimestamp(i.astype(int) * ns) for i in dt64]
    if RtnAsArray:
        return np.array(dt)
    else:
        return dt


def nonISOdate2ISO(ds):
    """
    Convert a non ISO date string to a ISO date string

    Parameters
    -----
    ds(str): date string
    """
    import re
    logging.info('nonISOdate2ISO called')
    regex = re.compile('(\d\d\d\d-\d-\d\d)')
    regexII = re.compile('(.*\s\d:.*)')
    print(ds)
    for d in ds:
        #        print 0, d, len(d[0]), [ re.match(regexII, d[0]) ]
        d = d[0]

        # swap ' ?:00:00' for '00:00:00'
        d = d.replace(' 0:',  ' 00:')
        if re.match(regexII, d):
            d = d[:-7]+'0'+d[-7:]
#        print 1, d, len(d)
        if len(d) != 19:

            # if one digit for day and month
            if len(d) == 17:
                d = d[:5]+'0'+d[5:7]+'0'+d[7:]
#            print 1.1, d, len(d), [ re.match(regex, d) ]

            # if one digit for day
            if (re.match(regex, d)):
                d = d[:5]+'0'+d[5:]
#            print 1.2, d, len(d)

            # if one digit for month
            if len(d) != 19:
                d = d[:8]+'0'+d[8:]
        if len(d) != 19:
            print((1.3, d, len(d[0])))
        d = [d]
        print((2, d, len(d[0])))
    return ds


def nearest(ts, s):
    """
    Find the nearest values (e.g. timestamp)

    Parameters
    -------
    ts (float, int, timestamp): point as object  that nearest to which is being sought
    s (list): list of objects of the same type to be searched

    Returns
    -------
    (timestamp)

    Notes
    -------
     - Credit: Raymond Hettinger
    http://stackoverflow.com/questions/8162379/python-locating-the-closest-timestamp
    """
    # Given a presorted list of timestamps:  s = sorted(index)
    i = bisect_left(s, ts)
    return min(s[max(0, i-1): i+2], key=lambda t: abs(ts - t))


def YYYYMMDD_HHMM_2_datetime(str1=None, str2=None, combined=False,
                             verbose=False, debug=False):
    """
    Mappable converter of strings to datetime.

    Parameters
    -------
    str1 (list): list of strings of times
    str2 (list): list of strings of dates
    combined (bool): if True, then a single list of strings is provided
    debug (bool): print debugging options to screen

    Returns
    -------
    (list)
    """
    # Combined as one string
    if combined:
        dtime = str1
        # Translate from str to datetime
        dtime = [time.strptime(i, '%Y%m%d%H%M') for i in dtime]
        dtime = [datetime_.fromtimestamp(time.mktime(i)) for i in dtime]
    # Combine to one string
    else:
        # Make pandas dataframe
        data = np.array([str1, str2])
        if debug:
            print((data.shape, data[:5, :], [type(i) for i in (str1, str2)]))
        df = pd.DataFrame(data=data.T, columns=['YYYYMMDD', 'HHMM'])
        # Convert to datetime
        dtime = DF_YYYYMMDD_HHMM_2_dt(df=df)
        dtime = dtime.index
    return dtime


def add_months(sourcedate, months):
    """
    Incremental increase of datetime by given months
    """
    month = sourcedate.month - 1 + months
    year = sourcedate.year + month / 12
    month = month % 12 + 1
    day = min(sourcedate.day, calendar.monthrange(year, month)[1])
    return datetime.datetime(year, month, day)


def add_days(sourcedate, days_):
    """
    Incremental increase of  datetime by given days
    """
    sourcedate += datetime.timedelta(days=float(days_))
    return sourcedate


def add_hrs(sourcedate, hrs_, debug=False):
    """
    Incremental increase of datetime by given hours
    """
    if debug:
        print((sourcedate, hrs_))
    sourcedate += datetime.timedelta(hours=float(hrs_))
    return sourcedate


def add_minutes(sourcedate, min_, debug=False):
    """
    Incremental increase of datetime by given minutes
    """
    sourcedate += datetime.timedelta(minutes=float(min_))
    return sourcedate


def add_secs(sourcedate, secs_, debug=False):
    """
    Incremental increase of datetime by given seconds
    """
    sourcedate += datetime.timedelta(seconds=float(secs_))
    return sourcedate


def secs_in_month(months=None, years=None):
    """
    Get number of seconds in a specific month for a specific year (default=2009)
    """
    # Get generica months and year (2009) if not given
    if not isinstance(months, list):
        months = list(range(1, 13))
    if not isinstance(years, list):
        years = [2009] * len(months)
    # Get number of seconds in specific month in year
    # conversion: sec => min => hours => days => months
    ars = []
    for i, m_ in enumerate(months):
        ars += [60*60*24*calendar.monthrange(int(years[i]), int(m_))[1]]
    # Return as a np.array
    return np.array(ars)


def get_dt4run(time_span='year', period=1, startyear=2005, endyear=2005,
               endhour=23, a=None, b=None):
    """
    Make list of datetimes for a given range or between two datetimes

    Parameters
    -------
    a, b (datetime.datetime): dates to create list of dates between (a=first date)
    endhour (int): last hour to  use in list of dates
    startyear, endyear (int): first and last year to output list of dates for
    time_span (str): string of time period (e.g. days)
    period (int): periodicity of returned list of dates (1= 1 hour)

    Returns
    -------
    (list)
    """
    # Set dates
    if isinstance(a, type(None)):
        a = datetime.datetime(startyear, 2, 1, 0, 0)
        if time_span == '3days':
            b = datetime.datetime(endyear, 2, 3, endhour, 0)  # 3 day
        if time_span == 'week':
            b = datetime.datetime(endyear, 2, 7, endhour, 0)  # week
        if time_span == 'month':
            b = datetime.datetime(endyear, 3, 1, endhour, 0)  # one month
        if time_span == '6months':
            b = datetime.datetime(endyear, 8, 1, endhour, 0)  # 6 month(s)
        if time_span == 'year':
            endyear = 2006  # Kludge as Data ran from Feb to Feb
            b = datetime.datetime(endyear, 1, 31, endhour, 0)  # full year

    # Make list of dates to view (hourly intervals between a and b)
    dates = dt_hrs_a2b(a, b)
    return dates


def dt_hrs_a2b(a, b, period=1, debug=False):
    """
    Returns list of hour spaced datetimes between two given datetimes

    Parameters
    -------
    a, b (datetime.datetime): dates to create list of dates between (a=first date)
    period (int): periodicity of returned list of dates (1= 1 hour)

    Returns
    -------
    (list)
    """
    dates = [a]
    if debug:
        print((dates, a, b, period))
    while dates[-1] < b:
        dates += [add_hrs(dates[-1], period)]
    if debug:
        print((dates[0], dates[-1]))
    return dates


# def normalise2dailymax(dates, data, debug=False):
#     """
#     Normalise data to daily maximiun.
#
#     ARGUMENTS:
#      - list of dates as datetime.datetime objects.
#      - list of of
#     """
#     logging.info('normalise2dailymax called')
#     if debug:
#         logging.debug([(type(i), i.shape) for i in (data, dates)])
#
#     # Get list of unique dates & remove mean from dates
#     dates = np.ma.array([datetime.datetime(*i.timetuple()[:3]) for i in dates])
#     idates = np.ma.array((sorted(set(dates))))
#
#     if debug:
#         logging.debug([(np.min(i), np.max(i), np.mean(i)) for i in [data]])
#     for s in idates:
#         #        print s, len(data[np.ma.where( dates == s) ]),  np.ma.max(data[np.ma.where( dates == s )] )
#         data[np.ma.where(dates == s)] = data[np.ma.where(
#             dates == s)] - np.ma.max(data[np.ma.where(dates == s)])
#     if debug:
#         logging.debug([(np.min(i), np.max(i), np.mean(i)) for i in [data]])
#     return data


def time2datetime(dates):
    """
    Convert time object to datetime object
    """
    assert type(dates) == list, 'Please provide a list of times to unc'
    return [datetime_.fromtimestamp(time.mktime(i)) for i in dates]


def num2month(input=None, reverse=False, rtn_dict=False):
    """
    Convert number (1-12) to abbreviated name of month

    Parameters
    -------
    reverse (bool): invert dictionary if reverse==True.
    rtn_dict (bool): return the entire dictionary instead of a value for a key

    Notes
    -------
     - input is either a 3 character month string or an integer 1=>12
    """
    d = {
        1: 'Jan',
        2: 'Feb',
        3: 'Mar',
        4: 'Apr',
        5: 'May',
        6: 'Jun',
        7: 'Jul',
        8: 'Aug',
        9: 'Sep',
        10: 'Oct',
        11: 'Nov',
        12: 'Dec'
    }

    if reverse:
        d = {v: k for k, v in list(d.items())}

    if rtn_dict:
        return d
    else:
        return d[input]


def DF_YYYYMMDD_HHMM_2_dt(df, date_header='YYYYMMDD', time_header='HHMM',
                          rmvars=None, epoch=False):
    """
    Convert times to datetime from time strings of HHMM and YYYYMMDD

    Parameters
    -------
    df (pd.DataFrame): dataframe containing columns of datetimes in string format
    time_header, date_header (str): column titles for time and date (?_header)
    rmvars (list): list of variables to remove from dataframe
    epoch (bool): return the values in terms of epoch (unix) time

    Returns
    -------
    (pd.DataFrame)
    """
    # Function to map integer to 4 char str
    def format(x): return '{:0>4}'.format(int(x))
    # Use mapped function for speed.
    df[time_header] = df[time_header].map(format)
    # Combine to make datetime.
    # ( go via integer for dates, to ensure no floating zeros appear )
    df['Datetime'] = df[date_header].astype(int).astype(str) + \
        df[time_header].astype(str)
    logging.debug('1st 10 dates: '.format(logging.debug(df['Datetime'][:10])))
    df['Datetime'] = pd.to_datetime(df['Datetime'], format='%Y%m%d%H%M')
    # Remove variables if list provided as "rmvars"
    if isinstance(rmvars, list):
        [df.drop(i, 1) for i in rmvars]
    # Convert to Epoch if requested
    if epoch:
        def format(x): return unix_time(x)
        df['Epoch'] = df['Datetime'].map(format).astype('i8')
        del df['Datetime']
    else:
        df.index = df['Datetime']
    return df


def get_TZ4loc(lat=50, lon=0):
    """
    Get the UTC offset (timezone/TZ) in hours for a given location

    Parameters
    -------
    lon (float): longitude in decimal degrees North
    lat (float): latitude in decimal degrees east

    Notes
    -------
     - Original file with timezone boundaries can be found here  http://ftp.as.harvard.edu/gcgrid/geos-chem/data/ExtData/HEMCO/TIMEZONES/v2015-02/
     - This may not include all locations (e.g. Cape Verde) and at somepoint should be updated to use the latest best data (linked below)
    https://github.com/evansiroky/timezone-boundary-builder

    Notes
    -------
    (float)
    """
    import os
    import xarray as xr
    import inspect
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    path = os.path.dirname(os.path.abspath(filename))
    folder = path+'/../data/'
    filename = 'timezones_voronoi_1x1.nc'
    folder = '/users/ts551/scratch/data/TIMEZONES/v2015-02/'
    ds = xr.open_dataset(folder+filename)
    UTC_offset = ds.sel(lat=lat, lon=lon, method='nearest').squeeze()
    UTC_offset = UTC_offset['UTC_OFFSET'].values.astype('timedelta64[h]')
    return UTC_offset.astype(float)


def unix_time(dt):
    """
    Convert datetime to Unix time.

    Parameters
    -------
    dt (datetime.datetime): Single datetime object

    Notes
    -------
     - epoch is counted from a reference time of:
    datetime.datetime(1970, 1, 1, 0, 0)
    """
    epoch = datetime.datetime.utcfromtimestamp(0)
    delta = dt - epoch
    return delta.days*86400+delta.seconds+delta.microseconds/1e6


def dt_days_a2b(a, b, period=1, debug=False):
    """
    Calculate days between two dattime.datetime format dates

    Parameters
    -------
    a, b (datetime.datetime): dates to create list of dates between (a=first date)
    period (int): periodicity of returned list of dates (1= 1 hour)

    Returns
    -------
    (list)
    """
    dates = [a]
    if debug:
        print((dates, a, b, period))
    while dates[-1] < b:
        dates += [add_days(dates[-1], period)]
    if debug:
        print((dates[0], dates[-1]))
    return dates


def get_nighttime_values(dates=None, data=None, select_nighttime=True,
                         select_daytime=False,
                         daybreak=datetime.datetime(1970, 1, 1, 6),
                         dayend=datetime.datetime(1970, 1, 1, 18)):
    """
    Calculate nighttime values using dates array and pandas
    """
    # use dataframe to map daytime boolean
    df = pd.DataFrame(np.array(dates))
    print(df)
    df.columns = ['Datetime']
    # function to generate boolean for daytime

    def is_daytime(input, daybreak=daybreak, dayend=dayend):
        """
        Takes datetime.datetime and retruns True (bool) if daytime
        """
        daytime = False
        # after daybreak
        if (input.hour >= daybreak.hour):
            daytime = True
        # ... and after nightfall
        if (input.hour > dayend.hour):
            daytime = False
        return daytime
    df['ind'] = df.index.values
    df['daytime'] = df['Datetime'].map(is_daytime)
    # Just select nighttime or daytime
    if select_nighttime:
        df = df[df['daytime'] == False]
    if select_daytime:  # select daytime
        df = df[df['daytime'] == True]
    # Select just indexed values
    data = np.array(data)[df['ind'].values, ...]
    dates = np.array(dates)[df['ind'].values]

    return data, dates


def get_daily_maximum(dates=None, data=None):
    """
    Calculate daily maximum values using dates array and pandas
    """
    # Use dataframe to hold dates and name column datetime
    df = pd.DataFrame(np.array(dates))
    df.columns = ['Datetime']
    # Add column of index numbers to allow for later indexing...
    df['ind'] = df.index.values
    # Add column for days

    def convert_datetime2days(input):
        return datetime.datetime(*input.timetuple()[:3])
    df['days'] = df['Datetime'].map(convert_datetime2days)

    # - loop days
    daily_max_data = []
    # Make sure data is a numpy array
    data = np.array(data)
    for day in sorted(set(df['days'])):
        print((day, df['days'][:5]))
        # Select data for day
        a_day_ind = df[df['days'] == day]
        # Select data for day
        a_day_data = data[a_day_ind['ind'].values, ...]
        print([i.shape for i in (a_day_data, a_day_ind, data)])
        # Get daily maximum
        daily_max_data += [a_day_data.max(axis=0)]
    # Get average daily maximum
    avg_data = np.array(daily_max_data).mean(axis=0)
    return avg_data


def get_8hr_rolling_mean(df, window=8):
    """
    Get 8 hour rolling mean of pandas dataframe/series.

    Parameters
    -------
    df (pd.DataFrame):
    window (int): the window (hrs) over which to calculate mean (default=8 hrs)

    Returns
    -------
    (pd.DataFrame)
    """
    # loop columns if Dataframe
    dfs = []
    try:
        for col in df.columns:
         # apply mean
            dfs += [df[col].rolling(window=window, center=False).mean()]
    # Just process values if Series
    except AttributeError:
        df = df.rolling(window=window, center=False).mean()
    # Combine dataframes
    if len(dfs) > 1:
        # concatenate
        df = pd.concat(dfs, axis=1)
    return df


def solartime(observer, sun=None):
    """
    Get Solartime  for location of 'observer' relative to 'sun'

    Parameters
    -------
    observer (ephem observer object): Location of the observer
    sun (ephem sun object): Which dun to use? (default: our sun)

    Returns
    -------
    (float)

    Notes
    -------
     - Credit: J.F. Sebastian
    http://stackoverflow.com/questions/13314626/local-solar-time-function-from-utc-and-longitude
    """
    import ephem
    if isinstance(sun, type(None)):
        sun = ephem.Sun()
    # Astronomical math - compute the angle between the sun and observe
    sun.compute(observer)
    # sidereal time == ra (right ascension) is the highest point (noon)
    hour_angle = observer.sidereal_time() - sun.ra
    return ephem.hours(hour_angle + ephem.hours('12:00')).norm  # norm for 24h
