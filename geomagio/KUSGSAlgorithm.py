"""
    Creates K-USGS Indices from H and D time-series data.
"""

import numpy as np
# import obspy.core
from obspy.core import Trace, Stats, Stream, UTCDateTime
import matplotlib
import matplotlib.pyplot as plot
from matplotlib.dates import DateFormatter, WeekdayLocator, DayLocator, MONDAY
from Algorithm import Algorithm
from geomagio import TimeseriesFactoryException
import copy
# import pandas as pd

ONEMINUTE = 60
ONEHOUR = 60 * ONEMINUTE
ONEDAY = 24 * ONEHOUR

MINUTESPERHOUR = 60
MINUTESPERDAY = 24 * MINUTESPERHOUR

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
                New stream object containing the converted coordinates.
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
    trace.stats.statistics = statistics(trace.data)

    totalMinutes = trace.stats.npts
    # This algorithm operates on entire calendar days of 1-Minute values.
    if (totalMinutes % MINUTESPERDAY) != 0:
        raise TimeseriesFactoryException(
                'Entire calendar days of minute data required for K.')

    starttime = trace.stats.starttime
    endtime = trace.stats.endtime

    # type = <type 'list'>
    # hours = get_hours(starttime, endtime)
    # days = get_days(starttime, endtime)
    # months = get_months(days)

    hours = []
    rawHours = []
    days = []
    rangeLimit = 1.0
    months = get_traces(trace, 'months')
    print len(months), "months"
    print Stream(months)
    for month in months:
        print month.stats.starttime
        stdD = month.stats.statistics['standarddeviation']
        maxRange = stdD * rangeLimit

        monthDays = get_traces(month, 'days')
        print len(monthDays), "days in month"
        print Stream(monthDays)
        for day in monthDays:
            print day.stats.starttime
            dayHours = get_traces(day, 'hours')
            print len(dayHours), "hours in day"
            for dayHour in dayHours:
                print dayHour.stats.starttime
                hour = clean_range(dayHour, maxRange)
                # clean again, later
                hours.append(hour)
                rawHours.append(dayHour)
            days.append(day)

    plot_hours(rawHours, hours, 'Before cleaning', 'After cleaning')
    # monthBefore = copy.deepcopy(months)

    # days = get_traces(trace, 'days')
    # dayBefore = copy.deepcopy(days)

    # hours = get_traces(trace, 'hours')
    # hourBefore = copy.deepcopy(hours)

    # cleanedHours = clean_range(hours, months)

    # hourBefore = copy.deepcopy(hours)
    # clean_distribution(hourBefore, hours, months)

    # plot_days(dayBefore, days, 'Input data', 'After')
    # print_days(days, 'wide')

    # timeseries.plot() # This doesn't show anything
    # trace.plot()      # This also shows nothing...
    # print_months(months)
    # plot_months(monthBefore, months)

    # print_all(trace.stats)
    # plot_all(monthBefore, months, dayBefore, days, hourBefore, hours)

def clean_distribution(hourBefore, hours, months, exclude=1.0):
    """
        Elminiate any MHVs within a month window that fall in the tails of the
        monthly distribution, which is defined as a number of standard
        deviations away from the monthly mean.

        Parameters
        ----------
        hourBefore : List <obspy.core.trac.Trace>
            List of hourly statistics before eliminating points
        hours : List <obspy.core.trac.Trace>
            List of hourly traces with statistics
        months : List <obspy.core.trac.Trace>
            List of monthly traces statistics
        exclude : Float
            Number of standard deviations (from the monthly statistics) to
            use as the acceptable range of MHVs within the month. Default is
            3 Standard Deviations, which should eliminate ~0.3% of the data,
            assuming a normal distribution, which I'm not sure is a valid
            assumptions, but it's what the papers say about the algorithm.
    """
    print "cleaning distribution"

    fig = plot.figure("Montly Statistics")

    dist_plot(fig, hours, months, 2, 0)

    prevMonth = -1
    for hour in hours:
        hourMonth = hour.stats.starttime.month
        for month in months:
            if (month == hourMonth):
                stddev = month.stats.statistics['standarddeviation']
                avg = month.stats.statistics['average']
                maximum = avg + stddev * exclude
                minimum = avg + stddev * exclude
                if (hour.stats.statistics['average'] > maximum) or (hour.stats.statistics['average'] < minimum):
                    print "Elminiate point"

    dist_plot(fig, hours, months, 2, 1)

    # plot.subplots_adjust(hspace=0.23, wspace=0.01)
    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.show()

    print "done cleaning distribution"
    # Uncomment to see hour data printed and/or plotted for evalutating.
    plot_hours(hourBefore, hours, 'Raw input data, entire data range',
        'After removing MHVs at tails of distribution, entire data range')
    # print_hours(hours, 'wide')

def clean_range(hour, maxRange):
    """Replace any hours within each month that have a range of minutes that
        is too extreme with NaN. Too extreme is defined by standard deviations
        of all minutes minutes within a month. Also eliminates any hours that
        have less than 30 minutes of valid (non-NaN) data.

    Parameters
    ----------
        hours : List <obspy.core.trac.Trace>
            List of hourly traces statistics
        maxRange : Float
            Number of standard deviations (from the monthly statistics) to
            use as the acceptable range of minutes within each hour in the
            month. Default is 1 Standard Deviation.
    """
    # print "#####  Cleaning MHVs  #####"
    # TODO: simplify this method by not passing in hourBefore, leave the
    # hours alone return the array.
    # TODO: Arrays should be hours instead of hourlyStat, months instead of
    #      months, days instead of dailyStats
    # TODO: "Here is an array of traces and an acceptable range, clear average
    #      for any traces that fall outside of this range"
    # nanMhvCount = 0
    trace = []

    clearAvg = False

    # stdD = month.stats.statistics['standarddeviation']
    # allowedRange = stdD * rangeLimit

    maximum = hour.stats.statistics['maximum']
    minimum = hour.stats.statistics['minimum']
    hourRange = maximum - minimum

    # if np.isnan(hour.stats.statistics['average']):
    #     nanMhvCount += 1

    if hourRange > maxRange:
        clearAvg = True
        # hour.stats.statistics['average'] = np.nan
        # nanMhvCount += 1

    else:
        nanMinuteCount = len(hour.data[hour.data == np.nan])
        # for minute in hour:
        #     # Keep track of minutes that are NaN for eliminating hours.
        #     if np.isnan(minute):
        #         nanMinuteCount += 1

        # If half of the minute values are bad, the entire hour is bad.
        if nanMinuteCount >= 30:
            clearAvg = True
            # hour.stats.statistics['average'] = np.nan
            # nanMhvCount += 1

    if clearAvg:
        stats = Stats(hour.stats)
        stats.statistics['average'] = np.nan
        return Trace(hour.data, stats)

    return hour
    # print str(nanMhvCount) + " Total NaN MHVs found."

    # print "#####  MHVs Cleaned  #####\n"
    # Uncomment to see hour data printed and/or plotted for evalutating.
    # plot_hours(hourBefore, hours, 'Raw input data',
    #     'After removing large minute ranges')
    # print_hours(hours, 'wide')

def clean_hours_other(hours, months, rangeLimit=2.0):
    """
        rangeLimit as number of standard devations from the monthly mean.
            Default is 2 standard deviations on top and bottom of the monthly
            average. Decimal values are allowed.
    """
    print "#####  Cleaning MHVs  #####"

    nanMhvCount = 0
    month = 0        # Numerical month value so we can see when it changes.
    monthCount = -1  # Reference the same month until we're done with it.
    monthCounts = [] # List of months ('month') with counts of NaN's.
    for hour in hours:
        maximum = hour.stats.statistics['maximum']
        minimum = hour.stats.statistics['minimum']
        hourRange = maximum - minimum

        hourMonth = hour.stats.starttime.month
        for monthStat in months:
            if monthStat.stats.starttime.month == hourMonth:
                stdD = monthStat.stats.statistics['standarddeviation']
                allowedRange = 2 * stdD * rangeLimit

        if np.isnan(hour.stats.statistics['average']):
            nanMhvCount += 1

            if hourMonth != month:
                month = hourMonth
                monthCount += 1
                monthCounts.append([month, 1])
            else:
                monthCounts[monthCount][1] += 1

        elif hourRange > allowedRange:
            hour.stats.statistics['average'] = np.nan
            nanMhvCount += 1

            if hourMonth != month:
                month = hourMonth
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

    print str(nanMhvCount) + " Total NaN MHVs found."
    for mhv in monthCounts:
        print "  " + str(mhv[1]) + " NaN MHVs in month: " + str(mhv[0])
        if mhv[1] >= 30:
            print "A MHV needs to be replaced in month " + str(mhv[0])

    print "#####  MHVs Cleaned  #####\n"

def dist_plot(fig, hours, months, sets=1, offset=0):
    means = []
    times = []

    prevMonth = -1
    monthCount = 0

    for hour in hours:
        means.append(hour.stats.statistics['average'])
        times.append(hour.stats.starttime)

        hourMonth = hour.stats.starttime.month
        if (prevMonth != hourMonth):
            if len(times) > 1:
                dist_subplot(fig, months, monthCount, prevTitle,
                    times, means, sets, offset)

            monthCount += 1
            prevMonth = hourMonth
            prevTitle = hour.stats.starttime
            means = []
            times = []

    dist_subplot(fig, months, monthCount, prevTitle, times, means,
        sets, offset)

def dist_subplot(fig, months, monthCount, title, times, means, sets, offset):
    monthTotal = len(months)

    subplot = fig.add_subplot(int(str(sets*monthTotal) + "1"
        + str(monthCount + offset*monthTotal)))

    subplot.set_title("MHVs for " + title.strftime('%B %Y'))
    subplot.xaxis.set_major_locator(DayLocator([2,5,10,15,20,25,30]))
    subplot.xaxis.set_major_formatter(DateFormatter('%b %d %Y'))
    subplot.grid(True)

    times = matplotlib.dates.date2num(times)
    plot.plot(times, means, color='blue', marker='+', label='MHVs')

    pts = 0
    total = 0
    for mean in means:
        total += 1
        if not np.isnan(mean):
            pts += 1
    pts = str(pts) + " / " + str(total) + " pts (" + str(total/24) + " days)"
    plot.legend(loc='best', numpoints=1, frameon=False, title=pts)

    mean = np.nanmean(means)
    plot.plot([times[0], times[len(times)-1]],[mean, mean], lw=1,
        color='black', label='mean')

    stddev = np.nanmax(means) - np.nanmin(means)
    lower = mean - 1.0*stddev
    upper = mean + 1.0*stddev
    plot.fill_between(times, lower, upper, facecolor='red', alpha=0.1)
    lower = mean - 2.0*stddev
    upper = mean + 2.0*stddev
    plot.fill_between(times, lower, upper, facecolor='yellow', alpha=0.1)

def get_days(starttime, endtime):
    """Get days between (inclusive) starttime and endtime.

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
        # day = obspy.core.UTCDateTime(day.timestamp + 86400)
        day = UTCDateTime(day.timestamp + 86400)

    return days

def get_hours(starttime, endtime):
    """Get all of the hours between starttime and endtime.

    Parameters
    ----------
        starttime : obspy.core.utcdatetime.UTCDateTime
        endtime : obspy.core.utcdatetime.UTCDateTime

    Returns
    -------
        array_like
        List of times, one per hour, for all hours between and including
        ``starttime`` and ``endtime``.

    Raises
    ------
        TimeseriesFactoryException
            If starttime is after endtime
    """
    if starttime > endtime:
        raise AlgorithmException('The starttime must be before endtime.')

    hours = []
    date = starttime

    while date < endtime:
        hours.append(date)
        date += ONEHOUR

    hours.append(date)

    return hours

# def get_hours(day):
#     """
#         Get all of the hours in the given day.

#         Parameters
#         ----------
#         day : obspy.core.utcdatetime.UTCDateTime

#         Returns
#         -------
#         array_like
#             List of times, one per month, for all days between and including
#             ``starttime`` and ``endtime``.
#     """
#     hours = []
#     # date = np.datetime64(day)
#     date = day

#     for i in range(0, 24):
#         # hour = date + i * ONEHOUR
#         date = date + ONEHOUR
#         # hour = UTC.UTCDateTime(str(hour))

#         hours.append(date)

#     return hours

# def get_day_boundaries(date):
#     end = date + ONEDAY - ONEMINUTE

#     return { 'starttime': date, 'endtime': end }

# def get_hour_boundaries(date):
#     end = date + ONEHOUR - ONEMINUTE

#     return { 'starttime': date, 'endtime': end }

# def get_month_boundaries(month):
#     # Numpy doesn't know how to add a month...so work-around.
#     date = np.datetime64(month, timezone='UTC')
#     starttime = np.datetime_as_string(date, timezone='UTC')

#     date = pd.DatetimeIndex([date])
#     month = date.month + 1
#     year = date.year

#     if month > 12:
#         month = 1
#         year += 1

#     beginNextMonth = UTCDateTime(year, month, 1)
#     endOfMonth = np.datetime64(beginNextMonth, timezone='UTC') - ONEMINUTE
#     endtime = np.datetime_as_string(endOfMonth, timezone='UTC')

#     starttime = UTCDateTime(str(starttime))
#     endtime = UTCDateTime(str(endtime))
#     return { 'starttime': starttime, 'endtime': endtime }

def get_months(days):
    """
        Get months between (inclusive) starttime and endtime.

        Parameters
        ----------
        days : List of obspy.core.utcdatetime.UTCDateTime

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
            date = UTCDateTime(str(date))
            months.append(date)
            month = day.month

    return months

def kSubplot(fig, num, title, timeList, rLabel, mLabel, color='b', marker='s'):
    """
        Helper method for making subplots.

        Parameters
        ----------
        fig : magplotlib.figure.Figure
            Instance of plot.figure to attach the subplot to
        num : Integer
            Subplot rows, columns and position
        title : String
            Subplot title
        timeList :
            Array of times and stats to use for plot
        rLabel : String
            Range label
        mLabel : String
            Mean label
        color : String with valid marker color
            Marker color
        marker : String with valid marker string
            Marker style/shape
    """
    means = []
    times = []
    ranges = []

    subplot = fig.add_subplot(num)
    subplot.set_title(title)
    subplot.xaxis.set_major_formatter(DateFormatter('%b %y'))
    fig.autofmt_xdate(rotation=45)

    for time in timeList:
        times.append(time.stats.starttime)
        means.append(time.stats.statistics['average'])
        ranges.append(time.stats.statistics['maximum']
                - time.stats.statistics['minimum'])
    times = matplotlib.dates.date2num(times)

    plot.errorbar(times, means, ranges, color='cyan', label=rLabel)
    plot.plot(times, means, color=color, marker=marker, label=mLabel)

    pts = 0
    for mean in means:
        if not np.isnan(mean):
            pts += 1
    plot.legend(loc='best', numpoints=1, frameon=False, title=str(pts)+" pts")

    mean = np.nanmean(means)
    plot.plot([times[0], times[len(times)-1]],[mean, mean], lw=1, label='mean')

    stddev = np.nanmax(means) - np.nanmin(means)
    lower = mean - 0.5*stddev
    upper = mean + 0.5*stddev
    plot.fill_between(times, lower, upper, facecolor='red', alpha=0.05)
    lower = mean - 1.0*stddev
    upper = mean + 1.0*stddev
    plot.fill_between(times, lower, upper, facecolor='orange', alpha=0.1)
    lower = mean - 2.0*stddev
    upper = mean + 2.0*stddev
    plot.fill_between(times, lower, upper, facecolor='yellow', alpha=0.1)

def plot_all(monBefore, monAfter, dayBefore, dayAfter, hourBefore, hourAfter):
    """
        Plot montly, daily, and hourly statistics before and after cleaning.

        Parameters
        ----------
        monBefore : List <obspy.core.trac.Trace>
            List of monthly statistics before cleaning
        monAfter : List <obspy.core.trac.Trace>
            List of monthly statistics after cleaning
        dayBefore : List <obspy.core.trac.Trace>
            List of daily statistics before cleaning
        dayAfter : List <obspy.core.trac.Trace>
            List of daily statistics after cleaning
        hourBefore : List <obspy.core.trac.Trace>
            List of hourly statistics before cleaning
        hourAfter : List <obspy.core.trac.Trace>
            List of hourly statistics after cleaning
    """
    fig = plot.figure('Average nT Values')

    monthTitle = 'Monthly means (nT)'
    dayTitle = 'Daily means (nT)'
    hourTitle = 'Hourly means (nT)'

    monthLabel = 'Daily Mean Range'
    dayLabel = 'Hourly Mean Range'
    hourLabel = 'Minute Range'

    ### Set up all of the plots BEFORE the data has been cleaned. ###
    # Plot MONTHS before cleaning
    kSubplot(fig, 231, monthTitle, monBefore, monthLabel, 'Month', 'blue', 's')

    # Plot DAYS before cleaning
    kSubplot(fig, 232, dayTitle, dayBefore, dayLabel, 'Day', 'blue', '^')

    # Plot HOURS before cleaning
    kSubplot(fig, 233, hourTitle, hourBefore, hourLabel, 'MHVs', 'blue', '+')

    #### Set up all of the plots AFTER the data has been cleaned. ###
    # Plot MONTHS after cleaning
    kSubplot(fig, 234, monthTitle, monAfter, monthLabel, 'Month', 'green', 's')

    # Plot DAYS after cleaning
    kSubplot(fig, 235, dayTitle, dayAfter, dayLabel, 'Day', 'green', '^')

    # Plot HOURS after cleaning
    kSubplot(fig, 236, hourTitle, hourAfter, hourLabel, 'MHVs', 'green', '+')

    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.show()

def plot_days(dayBefore, dayAfter, beforeTitle='', afterTitle=''):
    """
        Plot daily statistics before and after cleaning.

        Parameters
        ----------
        dayBefore : List <obspy.core.trac.Trace>
            List of daily statistics before cleaning
        dayAfter : List <obspy.core.trac.Trace>
            List of daily statistics after cleaning
        beforeTitle: String
            Title description to append to "before" plot title
        afterTitle: String
            Title description to append to "after" plot title
    """
    fig = plot.figure('Average daily nT Values')

    dayTitle = ' - Daily means (nT)'
    dayLabel = 'Hourly Mean Range'

    beforeTitle = beforeTitle + dayTitle
    afterTitle = afterTitle + dayTitle

    kSubplot(fig, 211, beforeTitle, dayBefore, dayLabel, 'Day', 'blue', '^')
    kSubplot(fig, 212, afterTitle, dayAfter, dayLabel, 'Day', 'green', '^')

    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.show()

def plot_hours(hourBefore, hourAfter, beforeTitle='', afterTitle=''):
    """
        Plot hourly statistics before and after cleaning.

        Parameters
        ----------
        hourBefore : List <obspy.core.trac.Trace>
            List of hourly statistics before cleaning
        hourAfter : List <obspy.core.trac.Trace>
            List of hourly statistics after cleaning
        beforeTitle: String
            Title description to append to "before" plot title
        afterTitle: String
            Title description to append to "after" plot title
    """
    fig = plot.figure('MHVs (Mean Hourly nT Values)')

    hourTitle = ' - Hourly means (nT)'
    hourLabel = 'Minute Range'

    beforeTitle = beforeTitle + hourTitle
    afterTitle = afterTitle + hourTitle

    kSubplot(fig, 211, beforeTitle, hourBefore, hourLabel, 'MHVs', 'blue', '+')
    kSubplot(fig, 212, afterTitle, hourAfter, hourLabel, 'MHVs', 'green', '+')

    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.show()

def plot_months(monBefore, monAfter):
    """
        Plot monthly statistics before and after cleaning.

        Parameters
        ----------
        monBefore : List <obspy.core.trac.Trace>
            List of monthly statistics before cleaning
        monAfter : List <obspy.core.trac.Trace>
            List of monthly statistics after cleaning
    """
    fig = plot.figure('Average monthly nT Values')

    monthTitle = 'Monthly means (nT)'
    monthLabel = 'Daily mean range'

    kSubplot(fig, 211, monthTitle, monBefore, monthLabel, 'Month', 'blue', 's')
    kSubplot(fig, 212, monthTitle, monAfter, monthLabel, 'Month', 'green', 's')

    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.show()

def print_all(stats):
    """
        Print statistics for the entire trace to the terminal.

        Parameters
        ----------
        stats : trace.stats

    """
    statistics = stats.statistics
    print "Start Time            : " + str(stats.starttime)
    print "End Time              : " + str(stats.endtime)
    print "Total # of Minutes    : " + str(stats.npts) + " (" \
        + str(stats.npts/MINUTESPERDAY) + " * " + str(MINUTESPERDAY) + ")"
    print "Average of all Minutes: " + str(statistics['average']) + "nT"
    print "Std Dev of all Minutes: " \
        + str(statistics['standarddeviation']) + "nT"
    print "Range of all Minutes  : " \
        + str(statistics['maximum'] - statistics['minimum']) + "nT"

def print_days(days, format='wide'):
    """
        Print statistics for each day to terminal.
        ### Example output ###
            Day          : 2013-12-31T00:00:00.000000Z
            Daily Average: 20894.2173562
            Daily Std Dev: 9.39171243572
            Daily Range  : 44.319

        Parameters
        ----------
        days : List <obspy.core.trac.Trace>
            List of daily traces with statistics
        wide : String
            If 'wide', print more horizontal, else print more vertical.
    """
    for day in dailyStats:
        stats = day.stats
        statistics = stats.statistics
        if format == 'wide':
            print "  Day: " + str(stats.starttime) \
                  + "\tDay Avg: " + str(statistics['average']) \
                  + "\tDay Std Dev: " + str(statistics['standarddeviation']) \
                  + "\tDay Range : " + str(statistics['maximum'] \
                  - statistics['minimum'])
        else:
            print "    Day          : " + str(stats.starttime)
            print "    Daily Average: " + str(statistics['average'])
            print "    Daily Std Dev: " + str(statistics['standarddeviation'])
            print "    Daily Range  : " + str(statistics['maximum'] \
                - statistics['minimum']) + "\n"

def print_hours(hours, format='wide'):
    """
        Print statistics for each hour to terminal.
        ### Example output ###
              Hour          : 2014-01-02T17:00:00.000000Z
              Hourly Average: 20855.7571167
              Hourly Std Dev: 10.1907743067
              Hourly Range  : 36.883

        Parameters
        ----------
        hours : List <obspy.core.trac.Trace>
            List of hourly traces with statistics
        wide : String
            If 'wide', print more horizontal, else print more vertical.
    """
    for hour in hours:
        stats = hour.stats
        statistics = stats.statistics
        if format == 'wide':
            print "   Hour: " + str(stats.starttime) \
                + "\t Avg: " + str(statistics['average']) \
                + "\t   Std Dev: " + str(statistics['standarddeviation']) \
                + "\t Range: " + str(statistics['maximum'] \
                - statistics['minimum'])
        else:
            print "      Hour          : " + str(stats.starttime)
            print "      Hourly Average: " + str(statistics['average'])
            print "      Hourly Std Dev: " + str(statistics['standarddeviation'])
            print "      Hourly Range  : " + str(statistics['maximum'] \
                - statistics['minimum']) + "\n"

def print_months(monthlyStats):
    """
        Print statistics for each month to terminal.
        ### Example output ###
          Month          : 2013-12-31T00:00:00.000000Z
          Monthly Average: 20894.2173562
          Monthly Std Dev: 9.39171243572
          Monthly Range  : 44.319

        Parameters
        ----------
        monthlyStats : List <obspy.core.trac.Trace>
            List of monthly statistics
    """
    for month in monthlyStats:
        stats = month.stats
        statistics = stats.statistics
        print "  Month          : " + str(stats.starttime)
        print "  Monthly Average: " + str(statistics['average'])
        print "  Monthly Std Dev: " + str(statistics['standarddeviation'])
        print "  Monthly Range  : " + str(statistics['maximum'] \
            - statistics['minimum']) + "\n"

def get_traces(trace, interval='hours'):
    """Use array of times to slice up trace and collect statistics.

    Parameters
    ----------
        trace :
            a time-series trace of data
        interval: String
            Interval to use for trace boundaries.
            Trace should include complete intervals.
            'hours', 'days', 'months' are accepted
    Returns
    -------
        array-like list of traces with statistics
    """
    traces = []

    starttime = trace.stats.starttime
    endtime = trace.stats.endtime

    date = starttime

    while date < endtime:
        start = date
        if interval == 'hours':
            date = start + ONEHOUR
        elif interval == 'days':
            date = start + ONEDAY
        elif interval == 'months':
            year = date.year
            month = date.month + 1
            if month > 12:
                month = 1
                year += 1
            date = UTCDateTime(year, month, 1)
        end = date - ONEMINUTE

        localTrace = trace.slice(start, end)
        localTrace.stats.statistics = statistics(localTrace.data)

        traces.append(localTrace)

    return traces
    # traces = []

    # if interval == 'hours':
    #     return slice_hours(get_hours(trace.starttime, trace.endtime), trace)

    # for time in times:
    #     if interval == 'hours':
    #         boundaries = get_hour_boundaries(time)
    #     elif interval == 'days':
    #         boundaries = get_day_boundaries(time)
    #     elif interval == 'months':
    #         boundaries = get_month_boundaries(time)

    #     starttime = boundaries['starttime']
    #     endtime = boundaries['endtime']

    #     thisTrace = trace.slice(starttime, endtime)

    #     if (interval == 'days') and (thisTrace.stats.npts != MINUTESPERDAY):
    #         raise TimeseriesFactoryException(
    #                     'Entire calendar days of minute data required for K.')

    #     thisTrace.stats.statistics = statistics(thisTrace.data)
    #     traces.append(thisTrace)

    # return traces

def slice_hours(hours, trace):
    """
        Use array of "hours" to slice up "trace" and collect hourly statistics
        to be attached to dailyStats.

        Parameters
        ----------
        hours :
            array of days UTC
        trace :
            a time-series trace of data

        Returns
        -------
            array-like list of hourly traces with statistics
    """
    # hourTraces = []
    # for day in hours:
    #     for hour in day:
    #         end = np.datetime64(hour) + ONEHOUR - ONEMINUTE
    #         # TODO Look into using the raw time value instead of a string
    #         end = UTC.UTCDateTime(str(end))

    #         thisHour = trace.slice(hour, end)
    #         if thisHour.stats.npts != MINUTESPERHOUR:
    #             raise TimeseriesFactoryException(
    #                     '1 Hour should have 60 minutes.')

    #         thisHour.stats.statistics = statistics(thisHour.data)
    #         hourTraces.append(thisHour)
    # return hourTraces
    traces = []
    for hour in hours:
        end = hour + ONEHOUR - ONEMINUTE
        # TODO Look into using the raw time value instead of a string
        # end = UTC.UTCDateTime(str(end))

        trace = trace.slice(hour, end)
        if trace.stats.npts != MINUTESPERHOUR:
            raise TimeseriesFactoryException(
                    '1 Hour should have 60 minutes.')

        trace.stats.statistics = statistics(trace.data)
        traces.append(trace)
    return traces

def statistics(data):
    """
        Calculate average, standard deviation, minimum and maximum on given
        trace, add them to a 'statistics' object and attach them to the trace.

        Parameters
        ----------
        data : <numpy.ndarray>
            an array of time-series data

        Returns
        -------
            object with key/value pairs for statistics
    """
    mean = np.nanmean(data)
    # Skip some calculations if this entire array is NaN's.
    if not (np.isnan(mean)):
        # Ignoring any NaN's for these calculations.
        statistics = {
            'average': mean,
            'maximum': np.nanmax(data),
            'minimum': np.nanmin(data),
            'standarddeviation': np.nanstd(data)
        }
    else:
        statistics = {
            'average': np.nan,
            'maximum': np.nan,
            'minimum': np.nan,
            'standarddeviation': np.nan
        }
    return statistics
