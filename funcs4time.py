#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Time processing functions for use with GEOS-Chem/Data analysis.

Use help(<name of function>) to get details on a particular function.

NOTE(S):
 - This module is underdevelopment vestigial/inefficient code is being removed/updated.
 - Where external code is used credit is given.
"""
# ----------------------------- Section X.X -----------------------------------
# -------------- Required modules:

import logging

# -- Math/Analysis
import numpy as np
#from time import mktime
#from pandas import DataFrame
import pandas as pd

# -- Time
import time
import calendar
import datetime as datetime
from datetime import datetime as datetime_
# Attempt to import ephem if installed
try:
    import ephem
except:
    print( 'ephem package not installed')

# ----------------------- Section X.X -------------------------------------------
# -------------- Time Processing
#

# --------------
# X.XX - Datetime to fractional day
# --------------
def get_day_fraction(date):
    """
    Get day fraction.

    NOTES:
     - for working with numpy arrays of datetimes, instead of pandas dataframes
    """
    secs = (date.hour *60.*60.)+(date.minute*60.)+(date.second)
    dsecs = 24.*60.*60.
    return  secs/dsecs


# --------------
# X.XX - numpy.datetime64 to datetime.datetime (assuming UTC )
# -------------
def dt64_2_dt( dt64 ):
    """
    Convert numpy.datetime64 to datetime.datetime (assuming UTC )

    NOTES:
     - ACTION NEEDED: Convert this to work as a lamdba function for
        scalability
    """
    ns = 1e-9 # number of seconds in a nanosecond
    return  [ datetime_.utcfromtimestamp(i.astype(int) * ns) for i in dt64 ]


# --------------
# X.XX - processor for non-ISO dates
# -------------
def nonISOdate2ISO( ds ):
    """
    Convert a non ISO date string to a ISO date string
    NOTES:
     - ds: date string
    """
    import re

    logging.info( 'nonISOdate2ISO called' )
    regex = re.compile( '(\d\d\d\d-\d-\d\d)')
    regexII = re.compile( '(.*\s\d:.*)')
    print(ds)
    for d in ds:
#        print 0, d, len(d[0]), [ re.match(regexII, d[0]) ]
        d = d[0]

        # swap ' ?:00:00' for '00:00:00'
        d= d.replace(' 0:',  ' 00:')
        if re.match(regexII, d):
            d= d[:-7]+'0'+d[-7:]
#        print 1, d, len(d)
        if len(d) != 19:

            # if one digit for day and month
            if len(d) == 17:
                d= d[:5]+'0'+d[5:7]+'0'+d[7:]
#            print 1.1, d, len(d), [ re.match(regex, d) ]

            # if one digit for day
            if (re.match(regex, d)) :
                d= d[:5]+'0'+d[5:]
#            print 1.2, d, len(d)

            # if one digit for month
            if len(d) != 19:
                d= d[:8]+'0'+d[8:]
        if len(d) != 19:
            print(1.3, d, len(d[0]))
        d = [d]
        print(2, d, len(d[0]))
    return ds


# --------------
# X.XX - Find nearest timestamp - Credit: Raymond Hettinger
# -------------
def nearest(ts, s):
    """
    Find nearest timestamp.

    ARGUEMENTS:
     - ts: point as object (float, int, timestamp) that nearset to which is being sought
     - s: list of objects of the same type to be searched
    NOTES:
     - Credit: Raymond Hettinger  -
    http://stackoverflow.com/questions/8162379/python-locating-the-closest-timestamp
    """
    # Given a presorted list of timestamps:  s = sorted(index)
    i = bisect_left(s, ts)

    return min(s[max(0, i-1): i+2], key=lambda t: abs(ts - t))


# --------------
# X.XX -  convert two lists/numpy arrays for date andtime to a datetime numpy
# -------------
def YYYYMMDD_HHMM_2_datetime( str1=None, str2=None, conbined=False,  \
            verbose=False, debug=False ):
    """
    Mappable converter of strings to datetime.

    ARGUEMENTS:
     - list of strings of dates (str1) and times (str2)
     - cobined (boolean): if True, then a single list of strings is provided
     - shared functionality with "DF_YYYYMMDD_HHMM_2_dt" ?
    """
    # combined as one string
    if conbined :
        dtime = str1

        # translate from str to datetime
        dtime = [ time.strptime( i, '%Y%m%d%H%M' ) for i in dtime ]
        dtime = [ datetime_.fromtimestamp(time.mktime(i)) for i in dtime ]

    # combine to one string
    else:

        # make pandas dataframe
        data = np.array( [str1,str2] )
        if debug:
            print(data.shape, data[:5,:], [ type(i) for i in (str1,str2) ])
        df = pd.DataFrame( data=data.T, columns=['YYYYMMDD', 'HHMM'] )

        # convert to datetime
        dtime = DF_YYYYMMDD_HHMM_2_dt( df=df )
        dtime = dtime.index

    return dtime


# -------------
# X.XX - Incremental increase datetime by given months - credit: Dave Webb
# -------------
def add_months(sourcedate,months):
    """
    Incremental increase of datetime by given months
    """
    month = sourcedate.month - 1 + months
    year = sourcedate.year + month / 12
    month = month % 12 + 1
    day = min(sourcedate.day,calendar.monthrange(year,month)[1])
    return datetime.datetime(year,month,day)


# TODO? Combine these functions?
# ( by adding specification of timedelta/case approach )
# -------------
# X.XX - Incremental increase datetime by given days
# -------------
def add_days(sourcedate,days_):
    """
    Incremental increase of  datetime by given days
    """
    sourcedate += datetime.timedelta(days=days_)
    return sourcedate

# -------------
# X.XX - Incremental increase datetime by given hours
# -------------
def add_hrs(sourcedate,hrs_, debug=False):
    """
    Incremental increase of datetime by given hours
    """
    if debug:
        print(sourcedate, hrs_)
    sourcedate += datetime.timedelta(hours=hrs_)
    return  sourcedate


# -------------
# X.XX - Incremental increase datetime by given minutes
# -------------
def add_minutes( sourcedate, min_, debug=False ):
    """
    Incremental increase of datetime by given minutes
    """
    sourcedate += datetime.timedelta(minutes=min_)
    return sourcedate


# -------------
# X.XX - Incremental increase datetime by given seconds
# -------------
def add_secs( sourcedate, secs_, debug=False ):
    """
    Incremental increase of datetime by given seconds
    """
    sourcedate += datetime.timedelta(seconds=secs_)
    return sourcedate


# --------------
# X.XX - day adjust -  seconds to months
# -------------
def d_adjust( months=None, years=None ):
    """
    Get adjustment values to convert an array of per second to per month
    """
    # Get months and years if not given
    if not isinstance(months, list):
        months = list(range(1,13))
    if not isinstance(years, list):
        years  = [2009]* len(months)

    # ajust to months ( => min => hours => days => months )
    return  np.array( [ 60*60*24*calendar.monthrange(int( years[i] ), \
        int( m_ ))[1] for i, m_ in enumerate(months) ]  )

# --------------
# X.XX - Get datetimes for run period
# -------------
def get_dt4run(time_span='year', period=1, startyear=2005, endyear=2005, \
        endhour=23, a=None, b=None ):
    """
    Get list of datetimes for a given range or between two provided
    datetimes  ( "a" and "b" )

    ARGUMENTS:
     - time_span : string of time period (e.g. days)
     - period: periodicty (1= 1 hour)
     - endhour: last hour of datetime list
     - first (startyear) and last year (endyear) requied
    """
    # Set dates
    if isinstance(a, type(None) ):
        a = datetime.datetime(startyear,2,1, 0, 0)
        if time_span == '3days':
            b = datetime.datetime(endyear,2,3, endhour, 0)  # 3 day
        if time_span == 'week':
            b = datetime.datetime(endyear,2,7, endhour, 0)  # week
        if time_span == 'month':
            b = datetime.datetime(endyear,3,1, endhour, 0)  # one month
        if time_span == '6months':
            b = datetime.datetime(endyear,8,1, endhour, 0)  # 6 month(s)
        if time_span == 'year':
            endyear=2006 # Kludge as Data ran from Feb to Feb
            b = datetime.datetime(endyear,1,31, endhour, 0)  # full year

    # --- Make list of dates to view (hourly intervals between a and b )
    dates = dt_hrs_a2b(a, b)
    return dates

# --------------
# X.XX - Datetime hours between datetime a and datetime b
# -------------
def dt_hrs_a2b( a, b, period=1, debug=False ) :
    """
    Returns list of hour spaced datetimes between two given datetimes

    ARGUMENTS:
     - two dates, one before (a) the other (b)
     - periodicty (1= 1 hour)
    """
    dates = [a]
    if debug:
        print(dates, a, b, period)
    while dates[-1] < b:
        dates += [ add_hrs(dates[-1], period) ]
    if debug:
        print(dates[0], dates[-1])
    return dates



# --------------
#  X.XX - Normalise data to daily maximiun.
# --------------
def normalise2dailymax(dates, data, debug=False ):
    """
    Normalise data to daily maximiun.

    ARGUMENTS:
     - list of dates as datetime.datetime objects.
     - list of of
    """
    logging.info( 'normalise2dailymax called' )
    if debug:
        logging.debug( [( type(i), i.shape ) for i in (data, dates) ] )

    # Get list of unique dates & remove mean from dates
    dates = np.ma.array([ datetime.datetime(*i.timetuple()[:3]) for i in dates ] )
    idates =np.ma.array((sorted(set( dates ) ) ))

    if debug:
        logging.debug( [(np.min(i), np.max(i), np.mean(i) ) for i in [data] ] )
    for s in idates:
#        print s, len(data[np.ma.where( dates == s) ]),  np.ma.max(data[np.ma.where( dates == s )] )
        data[np.ma.where( dates == s) ]  = data[np.ma.where( dates == s) ] - np.ma.max(data[np.ma.where( dates == s )] )
    if debug:
        logging.debug(  [(np.min(i), np.max(i), np.mean(i) ) for i in [data] ] )
    return data


# ----
#  X.XX - Translate from time to datetime
# ----
def time2datetime( dates ):
    """
    Convert time object to datetime object
    """
    return [ datetime_.fromtimestamp(time.mktime(i)) for i in dates ]


# ----
#  X.XX - return abbreviated month for a given month number or vice versa
# ----
def num2month(input=None, reverse=False, rtn_dict=False):
    """
    Convert number (1-12) to month in year

    ARGUMENTS:
     - reverse (boolean): invert dictionary if reverse==True.
     - input is either a 3 character month string or an integer 1=>12
    """
    d={
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

# --------------
#  X.XX - convert times to datetime from HHMM and YYYYMMDD
# ----------
def DF_YYYYMMDD_HHMM_2_dt(df, date_header='YYYYMMDD',
        time_header='HHMM', rmvars=None, epoch=False,
        verbose=False, debug=False):
    """
    Convert times to datetime from time strings of HHMM and YYYYMMDD

    ARGUMENTS:
     - column titles for time (time_header) and date (date_header)
     - rmvars: list of variables to remove from dataframe

    NOTES:
     - Use pandas DataFrame to allow for converting date and time strings
    by mapped functions for speed.
    """

    # --- Process time and dates
    # Map integer to 4 char str
    format = lambda x: '{:0>4}'.format( int( x ) )

    # Use mapped function for speed.
    df[time_header] = df[time_header].map( format )

     # Combine to make datetime.
     # ( go via integer for dates, to ensure no floating zeros appear )
    df['Datetime'] = df[date_header].astype(int).astype(str) + \
                                    df[time_header].astype(str)
    logging.debug('1st 10 dates: '.format(logging.debug( df['Datetime'][:10])))
    df['Datetime']  = pd.to_datetime( df['Datetime'], format='%Y%m%d%H%M' )

    # remove stated variables.
#    if not isinstance(rmvars, list ):
#        rmvars =['POINT','LAT', 'LON', 'PRESS', 'HHMM','YYYYMMDD'  ]
    if isinstance(rmvars, list ):
        [ df.drop(i, 1) for i in  rmvars ]

    # Convert to Epoch if requested
    if epoch:
        format = lambda x: unix_time(x)
        df['Epoch'] = df['Datetime'].map( format ).astype('i8')
        del df['Datetime']

    else:
        df.index=df['Datetime']

    return df

# --------------
#  X.XX - Convert datetime.datetine to  Unix time
# ----------
def unix_time(dt):
    """
    Convert datetime to Unix time.

    ARGUMENTS:
     - Single datetime object

    NOTES:
     - epoch = datetime.datetime(1970, 1, 1, 0, 0)
    """
    epoch = datetime.datetime.utcfromtimestamp(0)
    delta = dt - epoch
#    return delta.total_seconds()
    return delta.days*86400+delta.seconds+delta.microseconds/1e6

# --------------
# X.XX - Datetime days between datetime a and datetime b
# -------------
def dt_days_a2b( a, b, period=1, debug=False ) :
    """
    Calculate days between two dattime.datetime format dates

    ARGUMENTS:
     - two dates, one before (a) the other (b)
     - periodicty (1= 1day)
    """
    dates = [a]
    if debug:
        print(dates, a, b, period)
    while dates[-1] < b:
        dates += [ add_days(dates[-1], period) ]
    if debug:
        print(dates[0], dates[-1])
    return dates

# --------------
# X.XX - Select just nighttime values from an array.
# -------------
def get_nighttime_values( dates=None, data=None, select_nighttime=True,\
        select_daytime=False,
        daybreak=datetime.datetime(1970,0o1,0o1,0o6),
        dayend=datetime.datetime(1970,0o1,0o1,18) ):
    """
    Calculate nighttime values using dates array and pandas
    """
    # use dataframe to map daytime boolean
    df = pd.DataFrame( np.array(dates) )
    print(df)
    df.columns = ['Datetime']
    # function to generate boolean for daytime
    def is_daytime( input, daybreak=daybreak, dayend=dayend ):
        """
        Takes datetime.datetime and retruns True (boolean) if daytime
        """
        daytime=False
        # after daybreak
        if (input.hour >= daybreak.hour):
            daytime = True
        # ... and after nightfall
        if (input.hour > dayend.hour):
            daytime = False
        return daytime
    df['ind'] = df.index.values
    df['daytime'] = df['Datetime'].map( is_daytime )
    # just select nighttime or daytime
    if select_nighttime:
        df = df[ df['daytime'] == False ]
    if select_daytime: # select daytime
        df = df[ df['daytime'] == True ]

    # Select just indexed values
    data = np.array(data)[ df['ind'].values, ...  ]
    dates = np.array(dates)[ df['ind'].values ]

    return data, dates

# --------------
# X.XX - Get daily maximum
# -------------
def get_daily_maximum( dates=None, data=None ):
    """
    Calculate daily maximum values using dates array and pandas
    """
    # Use dataframe to hold dates and name column datetime
    df = pd.DataFrame( np.array(dates) )
    df.columns = ['Datetime']

    # Add column of index numbers to allow for later indexing...
    df['ind'] = df.index.values

    # Add column for days
    def convert_datetime2days(input):
        return datetime.datetime(*input.timetuple()[:3])
    df['days'] = df['Datetime'].map(convert_datetime2days)

    # --- loop days
    daily_max_data = []
    # make sure data is a numpy array
    data = np.array(data)
    for day in sorted(set( df['days']) ):

        print(day, df['days'][:5])
        # select data for day
        a_day_ind = df[ df['days'] == day ]
        # select data for day
        a_day_data = data[a_day_ind['ind'].values, ... ]
        print([i.shape for i in (a_day_data, a_day_ind, data) ])
        # get daily maximum
        daily_max_data += [ a_day_data.max(axis=0) ]

    # Get average daily maximum
    avg_data = np.array( daily_max_data ).mean( axis=0 )

    return avg_data

# --------------
# X.XX - Get 8hour rolling average
# -------------
def get_8hr_rolling_mean( df ):
    """
    Get 8 hour rolling mean of pandas dataframe/series.

    Parameters
    -------
    df (DataFrame)

    Returns
    -------
    (DataFrame)
    """

    # loop columns if Dataframe
    dfs = []
    try:
        for col in df.columns:
         # apply mean
            dfs += [ df[col].rolling(window=8,center=False).mean() ]
    # our just process values if Series
    except AttributeError:
        df = df.rolling(window=8,center=False).mean()

    #  combine dataframes
    if len(dfs) > 1:
        # concatenate
        df = pd.concat(dfs, axis=1)

    return df


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# ---------------- Section X -------------------------------------------
# -------------- User specific functions.
# --------------------------------------------------------------------------
#
# NOTE(s):
# (1) These functions should be removed following checks on compatibility
#


# --------------
# X.XX - returns data as a mean monthly value
# -------------
def data2monthly( data, dates ):
    """
    Resample list of numpy data by month (taken from given list of dates )

    NOTES:
        - Why is this a seperate function?
    """
    df = pd.DataFrame(data,index=dates )
    df['Month'] =  [ i.month for i in dates ]
    grouped = df.groupby('Month')
    totals = grouped.mean()
    return  totals.values, totals.index

# --------------
# X.XX - adjust non UT times to UT
# -------------
def gaw_lc_2_UT(time_s, site, half_hour=None, debug=False):
    """
    Adjust list of local times to UTC for a given GAW site

    NOTEs:
        - This function is redundent.
    """
    if debug:
        print('gaw_lc_2_UT called')
    UT_diff={'CVO':-1.}                    #  UT diff library for sites ....
    t_diff= float(  UT_diff[site] )*float(1./24.)     # look  up diff and adjust to decimal hours
    time_s_adjust = [ i + t_diff for i in time_s]  # adjust time series to UT
    time_s_adjust = [ float(str(i)[:str(i).find('.')+10])  for i  in time_s_adjust  ] # cut values to length
    time_s_adjust = np.array(time_s_adjust)  # return to as numpy array
    debug = True
    if debug:
        for i,ii in enumerate(time_s):
            print(i, ii, time_s_adjust[i])
    return time_s_adjust


# ----
# X.XX - Adjust to lt from UT
# ----
def adjust_UT_2_lt( time, date, data, site='CVO', dUTC=None, debug=False ):
    """
    Adjust list of UTC times to local time a given GAW site

    NOTEs:
     - This function is redundent. It is reverse of function "gaw_lc_2_UT"
    """
    from .funcs4generic import chunks

    if (dUTC ==  None ):
        dUTC   = gaw_2_loc(site)[-1]
    if debug:
        print('adjust_UT_2_lt called and dUTC = {}'.format(dUTC))

    # np.roll works in reverse to UTC change (reversed date also)
    time = np.array([np.roll( i, -1*dUTC ) for i in chunks(time,24) ]).flatten()
    date = np.roll(date, -1*dUTC )
    dUTC = -1*dUTC
    if (dUTC >= 0 ):
        print([ ( len(i),len(i)/24. ) for i in [time, data, date] ])
        time, data, date = [i[dUTC:-24+dUTC] for i in [time, data, date] ]
        print([ ( len(i),len(i)/24. ) for i in [time, data, date] ])
    else:
        print([ ( len(i),len(i)/24. ) for i in [time, data, date] ])
        time, data, date = [i[24+dUTC:dUTC] for i in [time, data, date] ]
        print([ ( len(i),len(i)/24. ) for i in [time, data, date] ])
    return time,date, data


# --------------
# X.XX - returns data as a mean monthly value
# -------------
def data2daily( data, dates ):
    """
    resample list of numpy data by day (e.g. taken from given list of dates )

    NOTES:
     - Redundent? Why is this a seperate function?
    """
    df = pd.DataFrame(data,index=dates )
    totals = df.resample('D', how='mean')
    return  totals.values, totals.index

# ----
# X.XX - Get interval dates assuming month gap in output.
# -----
def get_int_btwn(start, end, months=False, years=False ):
    """
    Get interval dates assuming month gap in output.

    NOTES:
     - Is this function redundent?
    """
    # list of dates
    m_,i = [], 0
    while(add_months(start,i) != end ):  # -1 so that final (end) value included
        m_.append( add_months(start,i) )
        i += 1

    if months:     # if months reqested...
        return [ i.strftime('%m') for i in m_]
    elif years:      # if years reqested...
        return [ i.strftime('%Y') for i in m_]
    else:
        return m_
#    else:
#        print 'State whether years or months are required as boolean arg (e.g. months=True)'

# ----
# X.XX -
# -----
def solartime(observer, sun=None):
    """
    """
    if isinstance(sun, type(None)):
        ephem.Sun()
    # Astronomical math
    import ephem
    # Credit: J.F. Sebastian
    # http://stackoverflow.com/questions/13314626/local-solar-time-function-from-utc-and-longitude
    sun.compute(observer)
    # sidereal time == ra (right ascension) is the highest point (noon)
    hour_angle = observer.sidereal_time() - sun.ra
    return ephem.hours(hour_angle + ephem.hours('12:00')).norm  # norm for 24h



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


# --------------
# X.XX - Process time/date to CV days equivilent - mje
# -------------
def year_to_since_2006(model):
    """
    Converts planeflight output date and time to unit Cape Verde (CVAO) years.

    ARGUEMTNS:
     -  extracted model data from funcs4pf.readfile_basic

    NOTES:
     - Credit MJE
    """
    year=(model[:,0]//10000)
    month=((model[:,0]-year*10000)//100)
    day=(model[:,0]-year*10000-month*100)
    hour=model[:,1]//100
    min=(model[:,1]-hour*100)
    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]
    since2006=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
    return since2006