"""
    Creates K-USGS Indices from H and D time-series data.
"""

import numpy as np
from obspy.core import Trace, Stats, Stream, UTCDateTime
import matplotlib
import matplotlib.pyplot as plot
from matplotlib.dates import DateFormatter, WeekdayLocator, DayLocator, MONDAY
from Algorithm import Algorithm
from geomagio import TimeseriesFactoryException
import copy

ONEMINUTE = 60
ONEHOUR = 60 * ONEMINUTE
ONEDAY = 24 * ONEHOUR

MINUTESPERHOUR = 60
MINUTESPERDAY = 24 * MINUTESPERHOUR

class KUSGSAlgorithm(Algorithm):
    """Algorithm for creating K-USGS indices.

    Parameters
    ----------
    """

    def __init__(self):
        Algorithm.__init__(self, inchannels=['H'], outchannels=['H','H','H','H'])

    def check_stream(self, timeseries, channels):
        """Checks a stream to make certain all the required channels exist.

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
        """Process the data to calculate K-USGS indices.

        Parameters
        ----------
            out_stream: obspy.core.Stream
                New stream object containing the converted coordinates.
        """
        out_stream = timeseries

        clean_MHVs(timeseries)

        return out_stream


def clean_MHVs(timeseries):
    """SR-Curve uses 24 Mean Hourly Values (MHVs) for 1 entire calendar
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

    hours = []
    rawHours = []
    days = []
    rangeLimit = 1.0        # 1 standard deviation TODO pass this in to the algorithm
    distributionLimit = 3.0 # 3 standard deviations TODO pass this in to the algorithm

    months = get_traces(trace, 'months')
    for month in months:
        avg = month.stats.statistics['average']
        stdD = month.stats.statistics['standarddeviation']
        maxRange = stdD * rangeLimit
        minimum = avg - distributionLimit * stdD
        maximum = avg + distributionLimit * stdD

        monthDays = get_traces(month, 'days')
        for day in monthDays:

            dayHours = get_traces(day, 'hours')
            for dayHour in dayHours:
                hour = clean_range(dayHour, maxRange)
                hour = clean_distribution(hour, minimum, maximum)
                hours.append(hour)
                rawHours.append(dayHour)
            days.append(day)

    # Uncomment any of the lines below to see data printed and/or plotted
    #   for evalutation purposes.
    # plot_ranges(rawHours, hours,
    #     'Before cleaning large minute ranges, all data',
    #     'After cleaning, all data')
    # plot_distribution(rawHours, hours, months)

    # print_times(hours, 'Hour', 'wide')
    # print_times(days, 'Day', 'wide')
    # print_times(months, 'Month', 'tall')
    print_all(trace.stats)

    # plot_days(dayBefore, days, 'Input data', 'After')

    # timeseries.plot() # This doesn't show anything
    # trace.plot()      # This also shows nothing...
    # plot_months(monthBefore, months)

    # plot_all(monthBefore, months, dayBefore, days, hourBefore, hours)

def clean_distribution(hour, minimum, maximum):
    """Elminiate any MHVs that are larger than maximum or smaller than minimum.

    Parameters
    ----------
        hour : Trace <obspy.core.trac.Trace>
            Hourly trace with statistics
        minimum : Float
        maximum : Float

    Returns
    -------

    """
    return hour

def clean_range(hour, maxRange):
    """Replace any hours within each month that have a range of minutes that
    is larger than maxRange with NaN. Also, eliminates any hours that have less
    than 30 minutes of valid (non-NaN) data.

    Parameters
    ----------
        hour : Trace <obspy.core.trac.Trace>
            Hourly trace with statistics
        maxRange : Float
            Number of standard deviations (from the monthly statistics) to
            use as the acceptable range of minutes within each hour in the
            month. Default is 1 Standard Deviation.

    Returns
    -------
        Trace with updated statistics
    """
    trace = []

    clearAvg = False

    maximum = hour.stats.statistics['maximum']
    minimum = hour.stats.statistics['minimum']
    hourRange = maximum - minimum

    if hourRange > maxRange:
        clearAvg = True

    else:
        nanMinuteCount = len(hour.data[hour.data == np.nan])

        # If half of the minute values are bad, the entire hour is bad.
        if nanMinuteCount >= 30:
            clearAvg = True

    if clearAvg:
        stats = Stats(hour.stats)
        stats.statistics = copy.deepcopy(stats.statistics)
        stats.statistics['average'] = np.nan
        return Trace(hour.data, stats)

    return hour

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

def kSubplot(fig, num, title, timeList, rLabel, mLabel, color='b', marker='s'):
    """Helper method for making subplots.

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
    """Plot montly, daily, and hourly statistics before and after cleaning.

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
    """Plot daily statistics before and after cleaning.

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

def plot_distribution(rawHours, hours, months):
    fig = plot.figure("Montly Statistics")

    dist_plot(fig, rawHours, months, 2, 0)

    dist_plot(fig, hours, months, 2, 1)

    # plot.subplots_adjust(hspace=0.23, wspace=0.01)
    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.tight_layout()
    plot.show()

def plot_months(monBefore, monAfter):
    """Plot monthly statistics before and after cleaning.

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

def plot_ranges(time1, time2, title1='', title2=''):
    """Plot hourly statistics before and after cleaning.

    Parameters
    ----------
        time1 : List <obspy.core.trac.Trace>
            List of hourly statistics before cleaning
        time2 : List <obspy.core.trac.Trace>
            List of hourly statistics after cleaning
        title1: String
            Title description to append to "before" plot title
        title2: String
            Title description to append to "after" plot title
    """
    fig = plot.figure('MHVs (Mean Hourly nT Values)')

    hourTitle = ' - Hourly means (nT)'
    rLabel1 = 'Minute Range Before'
    rLabel2 = 'Minute Range After'

    title1 = title1 + hourTitle
    title2 = title2 + hourTitle

    label1 = 'MHVs Before'
    label2 = 'MHVs After'

    plot_ranges_helper(fig, title1, title2, time1, time2,
                        rLabel1, rLabel2, label1, label2)

    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.tight_layout()
    plot.show()

def plot_ranges_helper(fig, title1, title2, list1, list2,
                        rLabel1, rLabel2, label1, label2):
    """Helper method for making subplots.

    Parameters
    ----------
        fig : magplotlib.figure.Figure
            Instance of plot.figure to attach the subplot to
        title1 : String
            1st half of subplot title
        title2 : String
            2nd half of subplot title
        list1 :
            Array of times and stats to use for before plot
        list2 :
            Array of times and stats to use for after plot
        rLabel1 : String
            Range label for before values
        rLabel2 : String
            Range label for after values
        label1 : String
            Mean label for before data
        label2 : String
            Mean label for after data
    """
    means1 = []
    means2 = []
    times1 = []
    times2 = []
    ranges1 = []
    ranges2 = []

    subplot = fig.add_subplot(111)
    subplot.set_title(title1 + " *VS* " + title2)
    subplot.xaxis.set_major_formatter(DateFormatter('%B %d, %Y'))
    subplot.xaxis.set_major_locator(DayLocator([5,15,25]))

    for time in list1:
        times1.append(time.stats.starttime)
        means1.append(time.stats.statistics['average'])
        ranges1.append(time.stats.statistics['maximum']
                - time.stats.statistics['minimum'])
    times1 = matplotlib.dates.date2num(times1)
    for time in list2:
        times2.append(time.stats.starttime)
        means2.append(time.stats.statistics['average'])
        ranges2.append(time.stats.statistics['maximum']
                - time.stats.statistics['minimum'])
    times2 = matplotlib.dates.date2num(times2)

    ptsTotal = 0
    for mean in means1:
        if not np.isnan(mean):
            ptsTotal += 1
    ptsRemain = 0
    for mean in means2:
        if not np.isnan(mean):
            ptsRemain += 1
    legendTitle = str(ptsRemain) + " of " +  str(ptsTotal) + " pts remaining"

    color1 = ['#ee95cf', '#e13636', 'red']
    plot.errorbar(times1, means1, ranges1, color=color1[0], label=rLabel1)
    mean1 = np.nanmean(means1)
    start1 = times1[0]
    end1 = times1[len(times1)-1]
    plot.plot([start1, end1], [mean1, mean1], label='Mean Before',
        color=color1[1], lw=2)
    # plot.plot(times1, means1, color=color1[2], marker='+', label=label1)

    color2 = ['#1523ea', '#000549', 'blue']
    plot.errorbar(times2, means2, ranges2, color=color2[0], label=rLabel2)
    mean2 = np.nanmean(means2)
    start2 = times2[0]
    end2 = times2[len(times2)-1]
    plot.plot([start2, end2], [mean2, mean2], label='Mean After',
        color=color2[1], lw=2)
    # plot.plot(times2, means2, color=color2[3], marker='+', label=label2)


    plot.legend(loc='best', numpoints=1, frameon=False, title=legendTitle)

    stddev = np.nanmax(means1) - np.nanmin(means1)
    lower = mean - 0.5*stddev
    upper = mean + 0.5*stddev
    plot.fill_between(times1, lower, upper, facecolor='red', alpha=0.04)
    lower = mean - 1.0*stddev
    upper = mean + 1.0*stddev
    plot.fill_between(times1, lower, upper, facecolor='orange', alpha=0.06)
    # lower = mean - 2.0*stddev
    # upper = mean + 2.0*stddev
    # plot.fill_between(times1, lower, upper, facecolor='yellow', alpha=0.08)

def print_all(stats):
    """Print statistics for the entire trace to the terminal.

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

def print_times(times, interval, format='wide'):
    """Print statistics for each time to terminal.

    Parameters
    ----------
        times : List <obspy.core.trac.Trace>
            List of time traces with statistics
        interval: String
            Display interval, typically: 'Hour', 'Day', 'Month'
        wide : String
            If 'wide', print more horizontal, else print more vertical.
    """
    for time in times:
        stats = time.stats
        statistics = stats.statistics
        if format == 'wide':
            print interval, ":", str(stats.starttime), \
                "\t Avg : ", str(statistics['average']), \
                "\t Std Dev: ", str(statistics['standarddeviation']), \
                "\t Range : ", str(statistics['maximum'] \
                - statistics['minimum'])
        else:
            print interval, "        : ", str(stats.starttime)
            print interval, " Average: ", str(statistics['average'])
            print interval, " Std Dev: ", str(statistics['standarddeviation'])
            print interval, " Range  : ", str(statistics['maximum'] \
                - statistics['minimum']), "\n"

def statistics(data):
    """Calculate average, standard deviation, minimum and maximum on given
    trace, add them to a 'statistics'.

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
