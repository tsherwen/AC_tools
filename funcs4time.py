"""
    Time processing functions for use with GEOS-Chem/Data analysis
"""

# --------------- 
# ---- Section 1 ----- Time processing
# 1.01 - Datetime to fractional day 
# 1.02 - numpy.datetime64 to datetime.datetime 
# 1.03 - Non ISO date to ISO (sonde strings )
# 1.04 - Find nearest timestamp in list (datetime)
# 1.05 - Process HHMM and YYYYMM strings to datetime
# 1.06 - Increase datetime by given months
# 1.07 - Increase datetime by given days
# 1.08 - Increase datetime by given hours
# 1.09 - Increase datetime by given minutes
# 1.10 - Increase datetime by given seconds
# 1.11 - convert dates to mod years (Kludge for obs. comparisons)
# 1.12 - Get days in months for a given year ( to adjust monhtly averages )
# 1.13 - non UT to UTC time corrector
# 1.14 - non UTC to UTC ... (redundent?) 
# 1.15 - Sort list of ctms into chronological order
# 1.16 - Avg. data by month
# 1.17 - Takes monthly outputs and sort to chronological order
# 1.18 - Get datetimes for period of run... 
# 1.19 - adjust data to daily
# 1.20 - Datetime hours between datetime a and datetime b
# 1.21 - Get interval dates assuming month gap in output
# 1.22 - Normalise data to daily mean. 
# 1.23 - Normalise data to daily min. 
# 1.24 - Normalise data to daily max. 
# 1.25 - hourly data to diurnal ... 
# 1.26 - Translate from time to datetime
# 1.27 - translate abrev. month name to num ( or vice versa ) 
# 1.28 - convert times to datetime from HHMM and YYYYMMDD
# 1.29 - Convert datetime.datetine to  Unix time
# 1.30 - Process time/date to CV days equivilent 

# - Math/Analysis                                                                                   
import numpy as np
from time import mktime
import scipy.stats as stats
from pandas import HDFStore
from pandas import DataFrame
from pandas import Series
from pandas import Panel
import math

# -- Time                                                                                           
import time
import calendar
import datetime as datetime
from datetime import datetime as datetime_

# --- Extras
from AC_tools.funcs4core import *
from AC_tools.funcs_vars import *
from AC_tools.funcs4generic import *

# ----------------------- Section 1 -------------------------------------------
# -------------- Time Processing
#

# --------------
# 1.01 - Datetime to fractional day 
# --------------
def get_day_fraction(date):
    secs = (date.hour *60.*60.)+(date.minute*60.)+(date.second)
    dsecs = 24.*60.*60.
    return  secs/dsecs
    
# --------------
# 1.02 - numpy.datetime64 to datetime.datetime (assuming UTC )
# ------------- 
def dt64_2_dt( dt64 ):
    """  ACTION NEEDED: Convert this to work as a lamdba function for 
            scalability"""
    ns = 1e-9 # number of seconds in a nanosecond
    return  [ datetime_.utcfromtimestamp(i.astype(int) * ns) for i in dt64 ]

# --------------
# 1.03 - numpy.datetime64 to datetime.datetime (assuming UTC )
# ------------- 
def nonISOdate2ISO( ds ):
    print 'nonISOdate2ISO'
    regex = re.compile( '(\d\d\d\d-\d-\d\d)')
    regexII = re.compile( '(.*\s\d:.*)')    
    print ds
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
            print 1.3, d, len(d[0])                        
        d = [d] 
        print 2, d, len(d[0])
    return ds

# --------------
# 1.04 - Find nearest timestamp
# -------------
def nearest(ts, s):

    """ Credit: Raymond Hettinger  -      http://stackoverflow.com/questions/8162379/python-locating-the-closest-timestamp              
    """

    # Given a presorted list of timestamps:  s = sorted(index)
    i = bisect_left(s, ts)

    return min(s[max(0, i-1): i+2], key=lambda t: abs(ts - t))

# --------------
# 1.05 -  convert two lists/numpy arrays for date andtime to a datetime numpy 
# -------------
def YYYYMMDD_HHMM_2_datetime( str1=None, str2=None,  conbined=False,  \
            verbose=False, debug=False ):    

    # combined as one string
    if conbined : 
        dtime = str1

        # translate from str to datetime
        dtime = [ time.strptime( i, '%Y%m%d%H%M' ) for i in dtime ]
        dtime = [ datetime_.fromtimestamp(mktime(i)) for i in dtime ]

    # combine to one string
    else: 

        # make pandas dataframe 
        data = np.array( [str1,str2] )
        if debug:
            print data.shape, data[:5,:], [ type(i) for i in str1,str2 ]
        df = DataFrame( data=data.T, columns=['YYYYMMDD', 'HHMM'] )

        # convert to datetime 
        dtime = DF_YYYYMMDD_HHMM_2_dt( df=df )
        dtime = dtime.index

    return dtime

# -------------
# 1.06 - Incremental increase datetime by given months - credit: Dave Webb
# -------------
def add_months(sourcedate,months):
    month = sourcedate.month - 1 + months
    year = sourcedate.year + month / 12
    month = month % 12 + 1
    day = min(sourcedate.day,calendar.monthrange(year,month)[1])
    return datetime.datetime(year,month,day)

# <= Combine 1.06- 1.10 , by adding specification of timedelta/case approach
# -------------
# 1.07 - Incremental increase datetime by given days 
# -------------
def add_days(sourcedate,days_):
    sourcedate += datetime.timedelta(days=days_)
    return sourcedate

# -------------
# 1.08 - Incremental increase datetime by given hours
# -------------
def add_hrs(sourcedate,hrs_, debug=False):
    if debug:
        print sourcedate, hrs_
    sourcedate += datetime.timedelta(hours=hrs_)    
    return  sourcedate

# -------------
# 1.09 - Incremental increase datetime by given minutes
# -------------
def add_minutes(sourcedate,min_, debug=False):
    sourcedate += datetime.timedelta(minutes=min_)
    return sourcedate

# -------------
# 1.10 - Incremental increase datetime by given seconds
# -------------
def add_secs(sourcedate,secs_, debug=False):
    sourcedate += datetime.timedelta(seconds=secs_)
    return sourcedate


# -------------
# 1.11 - convert_days_2_mod_year - Kludge for comparing obs and model data from different years
# ------------- 
# <=remove this (1.11). 
#def convert_days_2_mod_year( days, run1213=True ):
    
    # Set year. - standard years run x/07 to y/07
#    if run1213:
#        y = [2012, 2012, 2012, 2012, 2012, 2012, 2013, 2013, 2013, 2013, 2013, 2013]
#    else:
#        y = [2005]*12        
#    m = [7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6]
#    modd = dict(zip(m, y))

    # get input year, month, day
#    id = [ [i.year, i.month, i.day] for i  in days ]
#    print id

    # change years to mod year if month == input month
#    for v in id:
#        v[0] = modd[ v[1]  ]

#    print id
#    days = [ datetime.datetime( *v) for v in id]
#    return days    

# --------------
# 1.12 - day adjust -  seconds to months
# -------------
def d_adjust( months=None, years=None ):

    # Get months and years if not given
    if not isinstance(months, list):
        months = range(1,13)
    if not isinstance(years, list):
        years  = [2009]* len(months)

    # ajust to months ( => min => hours => days => months )
    return  np.array( [ 60*60*24*calendar.monthrange(int( years[i] ), \
        int( m_ ))[1] for i, m_ in enumerate(months) ]  )  

# --------------                                                                                                                                             
# 1.13 - adjust non UT times to UT
# ------------- 
def gaw_lc_2_UT(time_s, site, half_hour=None, debug=False):
    if debug:
        print 'gaw_lc_2_UT called'
    UT_diff={'CVO':-1.}                    #  UT diff library for sites ....
    t_diff= float(  UT_diff[site] )*float(1./24.)     # look  up diff and adjust to decimal hours
    time_s_adjust = [ i + t_diff for i in time_s]  # adjust time series to UT
    time_s_adjust = [ float(str(i)[:str(i).find('.')+10])  for i  in time_s_adjust  ] # cut values to length
    time_s_adjust = np.array(time_s_adjust)  # return to as numpy array
    debug = True
    if debug:
        for i,ii in enumerate(time_s):
            print i, ii, time_s_adjust[i]
    return time_s_adjust

# ----  
# 1.14 - Adjust to lt from UT 
# ----  
def adjust_UT_2_lt( time, date, data, site='CVO', dUTC=None, debug=False ):
    if (dUTC ==  None ):
        dUTC   = gaw_2_loc(site)[-1]
    if debug:
        print 'adjust_UT_2_lt called and dUTC = {}'.format(dUTC)

    # np.roll works in reverse to UTC change (reversed date also)  
    time = np.array([np.roll( i, -1*dUTC ) for i in chunks(time,24) ]).flatten()
    date = np.roll(date, -1*dUTC )
    dUTC = -1*dUTC
    if (dUTC >= 0 ):
        print [ ( len(i),len(i)/24. ) for i in [time, data, date] ]
        time, data, date = [i[dUTC:-24+dUTC] for i in [time, data, date] ] 
        print [ ( len(i),len(i)/24. ) for i in [time, data, date] ]
    else:
        print [ ( len(i),len(i)/24. ) for i in [time, data, date] ]
        time, data, date = [i[24+dUTC:dUTC] for i in [time, data, date] ] 
        print [ ( len(i),len(i)/24. ) for i in [time, data, date] ]
    return time,date, data


# --------------
# 1.16 - returns data as a mean monthly value
# -------------
def data2monthly( data, dates ):
    df = DataFrame(data,index=dates )
    df['Month'] =  [ i.month for i in dates ]
    grouped = df.groupby('Month')
    totals = grouped.mean()
    return  totals.values, totals.index
    
# --------------
# 1.18 - Get datetimes for run period
# -------------
def get_dt4run(time_span='year', period=1, startyear=2005,endyear=2005, 
                endhour=23, a=None, b=None  ): 

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
# 1.19 - returns data as a mean monthly value
# -------------
def data2daily( data, dates ):
    df = DataFrame(data,index=dates )
    totals = df.resample('D', how='mean')
    return  totals.values, totals.index

# --------------
# 1.20 - Datetime hours between datetime a and datetime b
# -------------
def dt_hrs_a2b( a, b, period=1, debug=False ) :
    dates = [a]
    if debug:
        print dates, a, b, period
    while dates[-1] < b:
        dates += [ add_hrs(dates[-1], period) ]
    if debug:
        print dates[0], dates[-1]
    return dates
    
# ----
# 1.21 - Get interval dates assuming month gap in output.
# -----
def get_int_btwn(start, end, months=False, years=False ):
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

# --------------
# 1.22 - Normalise data to daily mean.  - mv'd to funcs4time
# --------------
def normalise2dailymean(dates, data, debug=False ):
    if debug:
        print [( type(i), i.shape ) for i in data, dates ]#,  [(i.shape, np.ma.max(i), np.ma.min(i)) for i in [ data ] ]

    # Get list of unique dates & remove mean from dates
    dates = np.ma.array([datetime.datetime(*i.timetuple()[:3]) for i in dates ])
    idates =np.ma.array((sorted(set( dates ) ) ))
    for s in idates:
        data[np.ma.where( dates == s) ]  = data[np.ma.where( dates == s) ] - np.ma.mean(data[np.ma.where( dates == s )] )
    return data
    
# --------------
# 1.23 - Normalise data to daily minimum.   - mv'd to funcs4time
# --------------
def normalise2dailymin(dates, data, debug=False ):
    if debug:
        print [( type(i), i.shape ) for i in data, dates ]#,  [(i.shape, np.ma.max(i), np.ma.min(i)) for i in [ data ] ]

    # Get list of unique dates & remove mean from dates
    dates = np.ma.array(
    [ datetime.datetime(*i.timetuple()[:3]) for i in dates ] 
    )                                    
    idates =np.ma.array((sorted(set( dates ) ) )) 
    for s in idates:
        data[np.ma.where( dates == s) ]  = data[np.ma.where( dates == s) ] - np.ma.min(data[np.ma.where( dates == s )] )
    return data

# --------------
#  1.24 - Normalise data to daily maximiun.   - mv'd to funcs4time
# --------------
def normalise2dailymax(dates, data, debug=False ):
    if debug:
        print [( type(i), i.shape ) for i in data, dates ]#,  [(i.shape, np.ma.max(i), np.ma.min(i)) for i in [ data ] ]

    # Get list of unique dates & remove mean from dates
    dates = np.ma.array([ datetime.datetime(*i.timetuple()[:3]) for i in dates ] )
    idates =np.ma.array((sorted(set( dates ) ) ))

    if debug:
        print [(np.min(i), np.max(i), np.mean(i) ) for i in [data] ]
    for s in idates:
#        print s, len(data[np.ma.where( dates == s) ]),  np.ma.max(data[np.ma.where( dates == s )] )
        data[np.ma.where( dates == s) ]  = data[np.ma.where( dates == s) ] - np.ma.max(data[np.ma.where( dates == s )] )
    if debug:
        print [(np.min(i), np.max(i), np.mean(i) ) for i in [data] ]
    return data

# ----  
#  1.25 - Adjust hr data to diurnal
# ----  
def hr_data_2_diurnal(time, date, data, frac=False, diurnal=True, debug=False ):
    if debug:
        print 'hr_data_2_diurnal called'

    # get O3 max (09:00 local time) (9+UTadjust::24) & min (17:00 local time) (9+UTadjust::24
#    O3_max = data[np.where( time==np.int64(900)  )]
#    O3_min = data[np.where( time==np.int64(1700) )
    O3_max = data[np.where( time==np.int64(800)  )]
    O3_min = data[np.where( time==np.int64(1600) )]
    O3_d   = O3_max - O3_min
    print [ i[5000:5010] for i in  [ O3_d, O3_max, O3_min ]]

    # Adjust to diurnal
    print [ ( type(i), i.shape, np.min(i), np.max(i) ) for i in [data] ]
    if (diurnal):    
        if (frac):
            data = np.array( [ ( (data[np.where(date==day)]-O3_max[ii] )/data[np.where(date==day)] )*100
                               for ii, day in enumerate( sorted( set(date) )) ] ).flatten()
        else:
            data = np.array( [ data[ np.where(date==day) ] -O3_max[ii] for ii, day in enumerate( sorted( set(date) )) ] ).flatten()
    print [ ( type(i), i.shape, np.min(i), np.max(i) ) for i in [data]]

    # split by month for all years and provide daily data and diurnal value
    months = []

    # Remove data with a delta greater than 30 ppbv, as this is not chemical destruction but meteorology 
    cap    = 30
    for m_ in range(1,13,1):
        print m_, '{:0>2}'.format(m_ ) 
        m_all = data[ [ int(ii) for ii, i in enumerate(date) if ( str(i)[4:6] == '{:0>2}'.format(m_) )  ] ]
        print [ ( len(i),len(i)/24. ) for i in [m_all ]]
        z = np.array( chunks(m_all, 24) )
        print [ ( type(i), i.shape, np.min(i), np.max(i) ) for i in [z]]
        z = np.ma.masked_where( (z <= -cap), z )
        print [ ( type(i), i.shape, np.min(i), np.max(i) ) for i in [z]]
        z = np.ma.masked_where( (z >= cap), z )
        print [ ( type(i), i.shape, np.min(i), np.max(i) ) for i in [z]]
        z = np.mean( z, axis=0 )
        print [ ( type(i), i.shape, np.min(i), np.max(i) ) for i in [z]]
        months.append( z ) 
    months = np.array( months )
    return months, O3_max,O3_min,O3_d 


# ----  
#  1.26 - Translate from time to datetime
# ----  
def time2datetime( dates ):
    return [ datetime_.fromtimestamp(mktime(i)) for i in dates ]
    
# ----  
#  1.27 - return abbreviated month for a given month number or vice versa
# ----  
def num2month(input=None, reverse=False, rtn_dict=False):

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
        d= {
        'Jan' : 1,
        'Feb' : 2,
        'Mar' : 3,
        'Apr' : 4,
        'May' : 5,
        'Jun' : 6,
        'Jul' : 7,
        'Aug' : 8,
        'Sep' : 9, 
        'Oct' : 10,
        'Nov' : 11,
        'Dec' : 12 
        }

    if rtn_dict:
        return d 
    else:
        return d[input]
    
# -------------- 
#  1.28 - convert times to datetime from HHMM and YYYYMMDD
# ----------
def DF_YYYYMMDD_HHMM_2_dt(df, date_header='YYYYMMDD', 
        time_header='HHMM', rmvars=None, epoch=False,  
        verbose=False, debug=False):
    """ Use pandas DataFrame to allow for converting date and time strings
        by mapped functions for speed. """

    # --- Process time and dates
    # Map integer to 4 char str
    format = lambda x: '{:0>4}'.format( int( x ) )

    # Use mapped function for speed. 
    df[time_header] = df[time_header].map( format )

     # Combine to make datetime.
     # ( go via integer for dates, to ensure no floating zeros appear )
    df['Datetime'] = df[date_header].astype(int).astype(str) + \
                                    df[time_header].astype(str) 
    if debug:
        print df['Datetime'][:10]
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
#  1.29 - Convert datetime.datetine to  Unix time
# ----------
def unix_time(dt):
    """ time from epoch - datetime.datetime(1970, 1, 1, 0, 0) """
    epoch = datetime.datetime.utcfromtimestamp(0)
    delta = dt - epoch
#    return delta.total_seconds()
    return delta.days*86400+delta.seconds+delta.microseconds/1e6
    
    
# --------------
# 1.30 - Process time/date to CV days equivilent - mje
# -------------
""" translate year to "since2006" function """
def year_to_since_2006(model):
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
    
# --------------
# 1.31 - Datetime hours between datetime a and datetime b
# -------------
def dt_days_a2b( a, b, period=1, debug=False ) :
    dates = [a]
    if debug:
        print dates, a, b, period
    while dates[-1] < b:
        dates += [ add_days(dates[-1], period) ]
    if debug:
        print dates[0], dates[-1]
    return dates