"""
    Creates K-USGS Indices from H and D time-series data.
"""

import numpy as np
import obspy.core
import obspy.core.utcdatetime as UTC
from Algorithm import Algorithm
from geomagio import TimeseriesFactoryException

MINUTESPERHOUR = 60
MINUTESPERDAY = 1440  # 24h * 60min
ONEDAY = np.timedelta64(1, 'D')
ONEMINUTE = np.timedelta64(1, 'm')
ONEHOUR = np.timedelta64(1, 'h')

class KUSGSAlgorithm(Algorithm):
    """
        Algorithm for creating K-USGS indices.

        Parameters
        ----------
    """

    def __init__(self):
        Algorithm.__init__(self, inchannels=['H'], outchannels=['H','H','H','H'])

    def check_stream(self, timeseries, channels):
        """
            Checks a stream to make certain all the required channels exist.

            Parameters
            ----------
            timeseries: obspy.core.Stream
                Stream to be checked.
            channels: array_like
                Channels that are expected in stream.
        """
        for channel in channels:
            if len(timeseries.select(channel=channel)) == 0:
                print 'Channel %s not found in input' % channel
                return False
        return True

    def process(self, timeseries):
        """
            Process the data to calculate K-USGS indices.

            Parameters
            ----------
            out_stream: obspy.core.Stream
                new stream object containing the converted coordinates.
        """
        out_stream = timeseries

        clean_MHVs(timeseries)

        return out_stream


def clean_MHVs(timeseries):
    """
        SR-Curve uses 24 Mean Hourly Values (MHVs) for 1 entire calendar
        day, plus the last 2 MHVs from the previous day and the first 2
        MHVs from the following day. Thus the data is cleaned in daily
        chuncks using 28 MHVs. MHVs are excluded if they contain minute
        values that have a large range, or they fall in the tails of the
        monthly MHV distribution. MHVs that are excluded or don't exist
        are replaced with a daily or monthly mean.
    """
    # type = <class 'obspy.core.trace.Trace'>
    trace = timeseries.select(channel='H')[0]
    statistics(trace)

    totalMinutes = trace.stats.npts
    # This algorithm operates on entire calendar days of 1-Minute values.
    if (totalMinutes % MINUTESPERDAY) != 0:
        raise TimeseriesFactoryException(
                'Entire calendar days of minute data required for K.')

    starttime = trace.stats.starttime
    endtime = trace.stats.endtime

    # type = <type 'list'>
    days = get_days(starttime, endtime)

    months = get_months(days)

    monthlyStats = []
    monthly_slices(months, monthlyStats, trace)

    dailyStats = []
    daily_slices(days, dailyStats, trace)

    hours = []
    for day in days:
        hours.append(get_hours(day))

    hourlyStats = []
    hourly_slices(hours, hourlyStats, trace)

    clean_hours(hourlyStats)

    print_months(monthlyStats)
    # print_days(dailyStats)
    # print_hours(hourlyStats)

    print_all(trace.stats)

def clean_hours(hourlyStats):
    print "#####  Cleaning MHVs  #####"

    nanMhvCount = 0
    month = 0        # Numerical month value so we can see when it changes.
    monthCount = -1  # Reference the same month until we're done with it.
    monthCounts = [] # List of months ('month') with counts of NaN's.
    for hour in hourlyStats:
        if np.isnan(hour.stats.statistics['average']):
            nanMhvCount += 1

            if hour.stats.starttime.month != month:
                month = hour.stats.starttime.month
                monthCount += 1
                monthCounts.append([month, 1])
            else:
                monthCounts[monthCount][1] += 1

        else:
            nanMinuteCount = 0
            for minute in hour:
                # Keep track of minutes that are NaN for eliminating hours.
                if np.isnan(minute):
                    nanMinuteCount += 1

            # If half of the minute values are bad, the entire hour is bad.
            if nanMinuteCount >= 30:
                hour.stats.statistics['average'] = np.nan
                nanMhvCount += 1

                if hour.stats.starttime.month != month:
                    month = hour.stats.starttime.month
                    monthCount += 1
                    monthCounts.append([month, 1])
                else:
                    monthCounts[monthCount][1] += 1

    print str(nanMhvCount) + " NaN MHVs found."
    for mhv in monthCounts:
        if mhv[1] >= 30:
            print "A MHV needs to be replaced in month " + str(mhv[0])

    print "#####  MHVs Cleaned  #####\n"

def hourly_slices(hours, hourlyStats, trace):
    for day in hours:
        for hour in day:
            end = np.datetime64(hour) + ONEHOUR - ONEMINUTE
            # TODO Look into using the raw time value instead of a string
            end = UTC.UTCDateTime(str(end))

            thisHour = trace.slice(hour, end)
            if thisHour.stats.npts != MINUTESPERHOUR:
                raise TimeseriesFactoryException(
                        '1 Hour should have 60 minutes.')

            statistics(thisHour)
            hourlyStats.append(thisHour)

def daily_slices(days, dailyStats, trace):
    for day in days:
        end = np.datetime64(day) + ONEDAY - ONEMINUTE
        # TODO Look into using the raw time value instead of a string
        end = UTC.UTCDateTime(str(end))

        thisDay = trace.slice(day, end)
        if thisDay.stats.npts != MINUTESPERDAY:
            raise TimeseriesFactoryException(
                    'Entire calendar days of minute data required for K.')

        statistics(thisDay)
        dailyStats.append(thisDay)

def monthly_slices(months, monthlyStats, trace):
    for month in months:
        # Numpy doesn't know how to add a month...so work-around.
        date = np.datetime64(month, timezone='UTC')
        monthNum = int(np.datetime_as_string(date)[5:7])

        year = np.datetime_as_string(date, timezone='UTC')[:5]

        if monthNum < 10:
            monthNum = "0" + str(monthNum+1)
        elif monthNum < 12:
            monthNum = str(monthNum + 1)
        else:
            monthNum = str("01")
            year = str(int(year[:4])+1) + "-"

        end = year + monthNum + "-01T00:00:00.000000Z"
        end = np.datetime64(end, timezone='UTC')

        endtime = np.datetime_as_string(end - ONEMINUTE, timezone='UTC')
        # end work-around

        starttime = np.datetime_as_string(date, timezone='UTC')

        starttime = UTC.UTCDateTime(str(starttime))
        endtime = UTC.UTCDateTime(str(endtime))
        thisMonth = trace.slice(starttime, endtime)

        statistics(thisMonth)
        monthlyStats.append(thisMonth)

def statistics(trace):
    """
        Calculate average, standard deviation, minimum and maximum on given
        trace, add them to a 'statistics' object and attach them to the trace.
    """
    H = trace.data

    mean = np.nanmean(H)
    # Skip some calculations if this entire trace is NaN's.
    if not (np.isnan(mean)):
        # Ignoring any NaN's for these calculations.
        trace.stats.statistics = {
            'average': mean,
            'maximum': np.nanmax(H),
            'minimum': np.nanmin(H),
            'standarddeviation': np.nanstd(H)
        }
    else:
        trace.stats.statistics = {
            'average': np.nan,
            'maximum': np.nan,
            'minimum': np.nan,
            'standarddeviation': np.nan
        }

def get_months(days):
    """
        Get months between (inclusive) starttime and endtime.

        Returns
        -------
        array_like
            List of times, one per month, for all days between and including
            ``starttime`` and ``endtime``.
    """
    months = []

    month = 0
    for day in days:
        if day.month != month:
            date = np.datetime64(day)
            date = UTC.UTCDateTime(str(date))
            months.append(date)
            month = day.month

    return months

def get_hours(day):
    """
        Get all of the hours in the given day.

        Returns
        -------
        array_like
            List of times, one per month, for all days between and including
            ``starttime`` and ``endtime``.
    """
    hours = []
    date = np.datetime64(day)

    for i in range(0, 24):
        hour = date + i * ONEHOUR
        hour = UTC.UTCDateTime(str(hour))

        hours.append(hour)

    return hours

def get_days(starttime, endtime):
    """
        Get days between (inclusive) starttime and endtime.

        Parameters
        ----------
        starttime : obspy.core.UTCDateTime
            The start time
        endtime : obspy.core.UTCDateTime
            The end time

        Returns
        -------
        array_like
            List of times, one per day, for all days between and including
            ``starttime`` and ``endtime``.

        Raises
        ------
        TimeseriesFactoryException
            If starttime is after endtime
    """
    if starttime > endtime:
        raise TimeseriesFactoryException('starttime must be before endtime')

    days = []
    day = starttime
    lastday = (endtime.year, endtime.month, endtime.day)

    while True:
        days.append(day)
        if lastday == (day.year, day.month, day.day):
            break
        # move to next day
        day = obspy.core.UTCDateTime(day.timestamp + 86400)

    return days

def print_all(stats):
    statistics = stats.statistics
    print "Start Time            : " + str(stats.starttime)
    print "End Time              : " + str(stats.endtime)
    print "Total # of Minutes    : " + str(stats.npts) + " (" + str(stats.npts/MINUTESPERDAY) + " * " + str(MINUTESPERDAY) + ")"
    print "Average of all Minutes: " + str(statistics['average']) + "nT"
    print "Std Dev of all Minutes: " + str(statistics['standarddeviation']) + "nT"
    print "Range of all Minutes  : " + str(statistics['maximum'] - statistics['minimum']) + "nT"

def print_months(monthlyStats):
    ### Example output ###
    #  Month          : 2013-12-31T00:00:00.000000Z
    #  Monthly Average: 20894.2173562
    #  Monthly Std Dev: 9.39171243572
    #  Monthly Range  : 44.319
    for month in monthlyStats:
        stats = month.stats
        statistics = stats.statistics
        print "  Month          : " + str(stats.starttime)
        print "  Monthly Average: " + str(statistics['average'])
        print "  Monthly Std Dev: " + str(statistics['standarddeviation'])
        print "  Monthly Range  : " + str(statistics['maximum'] - statistics['minimum']) + "\n"

def print_days(dailyStats):
    ### Example output ###
    #    Day          : 2013-12-31T00:00:00.000000Z
    #    Daily Average: 20894.2173562
    #    Daily Std Dev: 9.39171243572
    #    Daily Range  : 44.319
    for day in dailyStats:
        stats = day.stats
        statistics = stats.statistics
        print "    Day          : " + str(stats.starttime)
        print "    Daily Average: " + str(statistics['average'])
        print "    Daily Std Dev: " + str(statistics['standarddeviation'])
        print "    Daily Range  : " + str(statistics['maximum'] - statistics['minimum']) + "\n"

def print_hours(hourlyStats):
    ### Example output ###
    #      Hour          : 2014-01-02T17:00:00.000000Z
    #      Hourly Average: 20855.7571167
    #      Hourly Std Dev: 10.1907743067
    #      Hourly Range  : 36.883
    for hour in hourlyStats:
        stats = hour.stats
        statistics = stats.statistics
        print "      Hour          : " + str(stats.starttime)
        print "      Hourly Average: " + str(statistics['average'])
        print "      Hourly Std Dev: " + str(statistics['standarddeviation'])
        print "      Hourly Range  : " + str(statistics['maximum'] - statistics['minimum']) + "\n"
