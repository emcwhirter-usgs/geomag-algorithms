"""Creates K-USGS Indices from H and D time-series data."""

import copy
import datetime
import matplotlib.pyplot as plot
import numpy as np

from Algorithm import Algorithm
from geomagio import TimeseriesFactoryException
from matplotlib.dates import \
        DateFormatter, WeekdayLocator, DayLocator, MONDAY, date2num
from obspy.core import Trace, Stats, Stream, UTCDateTime
from scipy import interpolate

ONEMINUTE = 60
ONEHOUR = 60 * ONEMINUTE
ONEDAY = 24 * ONEHOUR

MINUTESPERDAY = 24 * 60

class KUSGSAlgorithm(Algorithm):
    """Algorithm for creating K-USGS indices.

    Parameters
    ----------
    """

    def __init__(self, rangeLimit=1.0, distLimit=3.0):
        Algorithm.__init__(self, inchannels=['H', 'E'],
                        outchannels=['H','H','H','H'])

        self.rangeLimit = float(rangeLimit)
        self.distLimit = float(distLimit)

    def check_stream(self, timeseries, channels):
        """Checks a stream to make certain all the required channels exist.

        Parameters
        ----------
            timeseries: obspy.core.Stream
                Stream object containing input data.
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
        K-USGS is always calculated with entire months of data.

        Parameters
        ----------
            self.distLimit : Float
                Standard deviation limit to use for eliminating based on
                distribution.
            timeseries : obspy.core.Stream
                Stream object containing input data.
            self.rangeLimit : Float
                Standard deviation limit to use for eliminating based on ranges.
        """
        # TODO - Do all of the below with E in addition to H
        kVariationH = get_k_variation('H',
            timeseries, self.rangeLimit, self.distLimit)
        kVariationE = get_k_variation('E',
            timeseries, self.rangeLimit, self.distLimit)

        # TODO Translate
        translate(kVariationH, kVariationE)

        out_stream = timeseries

        return out_stream

def translate(kVariationH, kVariationD):
    """
    """
    print "translate"

def clean_distribution(hour, minimum, maximum, monthAverage):
    """Clean out MHVs at the edges of the monthly distribution, which is done
    by elminiating any MHVs that are larger than maximum or smaller than
    minimum. These eliminated values, along with any other NaNs in the data
    are replaced with the monthly average.

    Parameters
    ----------
        hour : Trace <obspy.core.trac.Trace>
            Hourly trace with statistics
        minimum : Float
            Minimum allowed hourly average
        maximum : Float
            Maximum allowed hourly average
        monthAverage : Float
            The monthly average for the month that the hour lies within

    Returns
    -------
        Trace with updated statistics.
    """
    clearAvg = False

    average = hour.stats.statistics['average']

    if (average > maximum) or (average < minimum) or (np.isnan(average)):
        clearAvg = True

    if clearAvg:
        stats = Stats(hour.stats)
        stats.statistics = copy.deepcopy(stats.statistics)
        stats.statistics['average'] = monthAverage
        return Trace(hour.data, stats)

    return hour

def clean_mhvs(channel, timeseries, rangeLimit, distributionLimit):
    """MHVs are excluded if they contain minute values that have a large range,
    defined by rangeLimit or they fall in the tails of the monthly MHV
    distribution, defined by distributionLimit. MHVs that are excluded or
    don't exist are replaced with the corresponding monthly mean.

    Parameters
    ----------
        distributionLimit : Float
            Standard deviation limit to use for eliminating based on
            distribution.
        rangeLimit : Float
            Standard deviation limit to use for eliminating based on ranges.
        timeseries : obspy.core.Stream
            Stream object containing input data.

    Returns
    -------
        List containing months with complete set of clean MHVs for the month
        attached as month.hours.
    """
    trace = timeseries.select(channel=channel)[0]
    trace.stats.statistics = statistics(trace.data)

    # This algorithm operates on entire calendar days of 1-Minute values.
    totalMinutes = trace.stats.npts
    if (totalMinutes % MINUTESPERDAY) != 0:
        raise TimeseriesFactoryException(
                'Entire calendar days of minute data required for K.')

    endtime = trace.stats.endtime
    starttime = trace.stats.starttime

    days = []
    hours = []
    rawHours = []

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
                # Clean out MHVs with ranges too large.
                hour = clean_range(dayHour, maxRange)
                # Clean out MHVs at the edges of the monthly distribution.
                hour = clean_distribution(hour, minimum, maximum, avg)
                hours.append(hour)
                rawHours.append(dayHour)      # Kept for printing/plotting only
            days.append(day)                  # Kept for printing/plotting only

        month.hours = hours
        hours = []

    """Uncomment any of the lines below to see data printed and/or plotted
       for evalutation purposes.
    """
    # Daily statistics - Average, Range and Standard Deviation #
    # print_stats(days, 'Day', 'wide')
    # plot_days(days)

    # Monthly statistics - Average, Range and Standard Deviation #
    # print_stats(months, 'Month', 'tall')
    # plot_months(months)

    # Total statistics - Average, Range and Standard Deviation #
    # print_all(trace.stats)
    # plot_all(months, days, rawHours)
    # print_stats(hours, 'Hour', 'wide')

    # Plot all data on 1 figure separated by month #
    # plot_distribution(rawHours, hours, months)

    return months

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
        Trace with updated statistics.
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

def clean_spline(x, y):
    """Remove any duplicate time points produced by np.linspace.

    Parameters
    ----------
        x : list
            Epoch times
        y : list
            nT values
    Returns
    -------
        Object with 'x' and 'y' defined as lists.
    """
    xClean = []
    yClean = []

    times = []
    for value in x:
        times.append(datetime.datetime.fromtimestamp(value))

    i = 0
    lastMinute = 0
    lastValue = -1
    for time in times:
        minute = time.minute

        # If it is the same time as previous, don't add it
        if minute == lastMinute:
            pass
        else:
            xClean.append(x[i])
            yClean.append(y[i])

        lastMinute = minute
        lastValue = y[i]
        i += 1

    return {'x': xClean, 'y': yClean}

def get_intercepts(lines):
    """Find intercepts of consecutive straight lines.

    Parameters
    ----------
        lines : List
            Array-like list of line segments defined by object with 'slope'
            and 'intercept' for y=mx+b

    Returns
    -------
        List of intercept objects with 'x' and 'y' defined.
    """
    xIntercepts = []
    yIntercepts = []

    for i in range(1, len(lines)-1):
        line0 = lines[i-1]
        line1 = lines[i]

        m0 = line0['slope']
        m1 = line1['slope']
        b0 = line0['intercept']
        b1 = line1['intercept']

        if ((m0 - m1) == 0):
            # Same slope, no intercept. Using the original point.
            x0 = line0['x']
            x1 = line1['x']
            # Don't add the point if is exactly the same as the last one.
            if ((b0 - b1) != 0) and ((x0 - x1) != 0):
                xIntercepts.append(line1['x'])
                yIntercepts.append(line1['y'])

        elif (m0 - m1) < 0.01:
            # Slopes are very close. Using the original point.
            # It looks like there are points not caught by 0.001
            xIntercepts.append(line1['x'])
            yIntercepts.append(line1['y'])

        else:
            x = (b1 - b0) / (m0 - m1)
            y = m0 * x + b0

            xIntercepts.append(x)
            yIntercepts.append(y)

    # plot_lines(lines, 0, 0)
    # plot_intercepts(xIntercepts, yIntercepts, 0, 0)
    return {'x-intercepts': xIntercepts, 'y-intercepts': yIntercepts}

def get_k_variation(channel, timeseries, rangeLimit, distLimit):
    """SR-Curve (Solar Regular Curve) uses 24 Mean Hourly Values (MHVs) for 1
    entire calendar day, plus the last 2 MHVs from the previous day and the
    first 2 MHVs from the following day. Thus the data is cleaned in daily
    chuncks using 28 MHVs. MHVs are excluded if they contain minute
    values that have a large range, or they fall in the tails of the
    monthly MHV distribution. MHVs that are excluded or don't exist
    are replaced with a daily or monthly mean.

    Parameters
    ----------
        channel : String
            H or D to specify which data channel to use.
        self.distLimit : Float
            Standard deviation limit to use for eliminating based on
            distribution.
        timeseries : obspy.core.Stream
            Stream object containing input data.
        self.rangeLimit : Float
            Standard deviation limit to use for eliminating based on ranges.

    Returns
    -------
        A trace of k-variation (SR-curve subtracted from every point).
    """
    # Clean up the data and make it continuous.
    months = clean_mhvs(channel, timeseries, rangeLimit, distLimit)

    # Get least square linear fit of sliding window over 3 MHVs at a time.
    lines = get_lines(months)

    # Get list of intercepts of all consecutive lines.
    intercepts = get_intercepts(lines)

    # Create a cubic spline from the intercepts.
    spline = get_spline(intercepts)

    # Subtract SR-curve (spline) from the input data.
    kTrace = remove_sr_curve(channel, timeseries, spline)

    return kTrace

def get_line(h0, h1, h2):
    """Find least squares fit of straight line of 3 points.

    Parameters
    ----------
        h0 : List <obspy.core.trac.Trace>
            First hour trace with statistics
        h1 : List <obspy.core.trac.Trace>
            Second hour trace with statistics
        h2 : List <obspy.core.trac.Trace>
            Third hour trace with statistics

    Returns
    -------
        Object with properties 'slope', 'intercept', 'x' and 'y' defined.
    """
    x0 = h0.stats.starttime.timestamp
    x1 = h1.stats.starttime.timestamp
    x2 = h2.stats.starttime.timestamp

    y0 = h0.stats.statistics['average']
    y1 = h1.stats.statistics['average']
    y2 = h2.stats.statistics['average']

    sumX = (x0 + x1 + x2)
    sumY = (y0 + y1 + y2)
    sumXY = (x0*y0 + x1*y1 + x2*y2)
    sumXX = (x0**2 + x1**2 + x2**2)

    slope = (3*sumXY - sumX * sumY) / (3*sumXX - sumX**2)
    intercept = (sumY*sumXX - sumX*sumXY) / (3*sumXX - sumX**2)

    return {'slope': slope, 'intercept': intercept, 'x': x1, 'y': y1}

def get_lines(months):
    """Create least squares fit of straight lines to sliding set of 3 MHVs.

    Parameters
    ----------
        months : List <obspy.core.trac.Trace>
            List containing months with MHVs for the month attached as
            month.hours.

    Returns
    -------
        List of line segments defined by 'slope' and 'intercept' for y=mx+b.
    """
    lines = []
    numMonths = len(months)

    if numMonths < 3:
        raise Exception('Must have data from 3 months in order to extract 2 '
            + 'additional points at the begin and end of the month for line '
            + 'segment calculations for desired month.')

    # If we pull 2 points instead of 4 from the month after, we
    # end up 2 hours short later in the process.
    if numMonths > 3:
        for i in range(1, len(months)-2):
            m0 = months[i-1]
            m1 = months[i]
            m2 = months[i+1]

            allHours = m0.hours[-2:] + m1.hours + m2.hours[0:4]

            for j in range(1, len(allHours)-2):
                h0 = allHours[j-1]
                h1 = allHours[j]
                h2 = allHours[j+1]

                lines.append(get_line(h0, h1, h2))
    else:
        m0 = months[0]
        m1 = months[1]
        m2 = months[2]

        allHours = m0.hours[-2:] + m1.hours + m2.hours[0:4]

        for j in range(1, len(allHours)-2):
            h0 = allHours[j-1]
            h1 = allHours[j]
            h2 = allHours[j+1]

            lines.append(get_line(h0, h1, h2))

    # plot_lines(lines, 0, 0)
    return lines

def get_spline(intercepts):
    """Create a cubic spline with list of consecutive best fit line intercepts.

    Parameters
    ----------
        intercepts : dictionary
            An object with 2 lists containing x and y-intercepts.
            x values are seconds since epoch.
            y values are nanoTeslas (nT).

    Returns
    -------
        Object with 'x' and 'y' defined as lists.
    """
    x = intercepts['x-intercepts']
    y = intercepts['y-intercepts']

    if len(x) != len(y):
        raise Exception('X and Y intercept lists must be the same length.')

    xLast = -1
    for value in x:
        if value < xLast:
            raise Exception('X (time) values must be consecutive (in order).')
        xLast = value

    f = interpolate.interp1d(x, y, kind='cubic')

    xnew = np.linspace(x[0], x[len(x)-1], 60*len(x))
    ynew = f(xnew)

    if len(xnew) % MINUTESPERDAY != 0:
        spline = clean_spline(xnew, ynew)
        xnew = spline['x']
        ynew = spline['y']

    # plot_spline(x, y, xnew, ynew)
    return {'x': xnew, 'y': ynew}

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
        array-like list of traces with statistics.
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

def plot_all(months, days, hours):
    """Plot montly, daily, and hourly statistics before and after cleaning.

    Parameters
    ----------
        months : List <obspy.core.trac.Trace>
            List of monthly statistics before cleaning
        days : List <obspy.core.trac.Trace>
            List of daily statistics before cleaning
        hours : List <obspy.core.trac.Trace>
            List of hourly statistics before cleaning
    """
    fig = plot.figure('Average nT Values')

    monthTitle = 'Monthly means (nT)'
    dayTitle = 'Daily means (nT)'
    hourTitle = 'Hourly means (nT)'

    monthLabel = 'Daily Mean Range'
    dayLabel = 'Hourly Mean Range'
    hourLabel = 'Minute Range'

    plot_subplot(fig, 311, monthTitle, months, monthLabel,
        'Month', 'blue', 's')

    plot_subplot(fig, 312, dayTitle, days, dayLabel,
        'Day', 'blue', '^')

    plot_subplot(fig, 313, hourTitle, hours, hourLabel,
        'MHVs', 'blue', '+')

    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.show()

def plot_days(days):
    """Plot daily statistics before and after cleaning.

    Parameters
    ----------
        days : List <obspy.core.trac.Trace>
            List of daily statistics before cleaning
    """
    fig = plot.figure('Average daily nT Values')

    dayTitle = ' - Daily means (nT)'
    dayLabel = 'Hourly Mean Range'

    beforeTitle = 'Daily means (nT) - entire data range'

    plot_subplot(fig, 111, beforeTitle, days, dayLabel,
        'Day', 'blue', '^')

    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.show()

def plot_distribution(rawHours, hours, months):
    """Plot monthly statistics before and after cleaning.

    Parameters
    ----------
        rawHours : List <obspy.core.trac.Trace>
            List of hourly statistics before cleaning
        hours : List <obspy.core.trac.Trace>
            List of hourly statistics after cleaning
        months : List <obspy.core.trac.Trace>
            List of monthly statistics to use to separate plots
    """
    fig = plot.figure("Montly Statistics")

    plot_dist_helper(fig, rawHours, hours, months)

    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.tight_layout()
    plot.show()

def plot_dist_helper(fig, beforeTimes, afterTimes, months):
    """Plot helper for distribution plotter. Separates subplots by month.

    Parameters
    ----------
        fig : magplotlib.figure.Figure
            Instance of plot.figure to attach the subplot to
        beforeTimes : List <obspy.core.trac.Trace>
            Array of times and stats before changes to use for plot
        afterTimes : List <obspy.core.trac.Trace>
            Array of times and stats after changes to use for plot
        months : List <obspy.core.trac.Trace>
            List of monthly statistics to use to separate plots
    """
    means1 = []
    means2 = []
    times1 = []
    times2 = []

    prevMonth = -1
    plotCount = 0

    count = 0
    lastMonth = beforeTimes[len(beforeTimes)-1].stats.starttime.month
    for hour in beforeTimes:
        means1.append(hour.stats.statistics['average'])
        times1.append(hour.stats.starttime)

        hourMonth = hour.stats.starttime.month
        if (prevMonth != hourMonth):
            if len(times1) > 1:
                monthHours = months[count].hours
                for monthHour in monthHours:
                    means2.append(monthHour.stats.statistics['average'])
                    times2.append(monthHour.stats.starttime)
                count += 1

                plot_dist_subplot(fig, months, plotCount, title,
                    times1, means1, times2, means2)

            plotCount += 1
            prevMonth = hourMonth
            title = hour.stats.starttime

            means1 = []
            times1 = []

            means2 = []
            times2 = []

    monthHours = months[count].hours
    for monthHour in monthHours:
        means2.append(monthHour.stats.statistics['average'])
        times2.append(monthHour.stats.starttime)

    plot_dist_subplot(fig, months, plotCount, title,
        times1, means1, times2, means2)

def plot_dist_subplot(fig, months, plotCount, title,
                    times1, means1, times2, means2):
    """
    Parameters
    ----------
        fig : magplotlib.figure.Figure
            Instance of plot.figure to attach the subplot to
        months : List <obspy.core.trac.Trace>
            List of monthly statistics to use to separate plots
        plotCount : Integer
            The number of the subplot to add to the figure
        title : String
            Title for the subplot
        times1 : List <obspy.core.trac.Trace>
            Array of times and stats before changes to use for plot
        means1 : List <obspy.core.trac.Trace>
            Array of means before changes to use for plot
        times2 : List <obspy.core.trac.Trace>
            Array of times and stats after changes to use for plot
        means2 : List <obspy.core.trac.Trace>
            Array of means before changes to use for plot
    """
    monthTotal = len(months)

    subplot = fig.add_subplot(int(str(monthTotal) + "1" + str(plotCount)))

    subplot.set_title("MHVs for " + title.strftime('%B %Y'))
    subplot.xaxis.set_major_locator(DayLocator([2,5,10,15,20,25,30]))
    subplot.xaxis.set_major_formatter(DateFormatter('%b %d %Y'))
    subplot.grid(True)

    times1 = date2num(times1)
    times2 = date2num(times2)
    plot.plot(times1, means1, color='#00bfff', marker='+', label='MHVs before')
    plot.plot(times2, means2, color='#104e8b', marker='^', label='MHVs after')

    pts = 0
    total = 0
    for mean in means2:
        total += 1
        if not np.isnan(mean):
            pts += 1
    pts = str(pts) + " / " + str(total) + " pts (" + str(total/24) + " days)"
    plot.legend(loc='best', numpoints=1, frameon=False, title=pts)

    mean = np.nanmean(means1)
    plot.plot([times1[0], times1[len(times1)-1]],[mean, mean], lw=1,
        color='black', label='mean')

    stddev = np.nanmax(means1) - np.nanmin(means1)
    lower = mean - 1.0*stddev
    upper = mean + 1.0*stddev
    plot.fill_between(times1, lower, upper, facecolor='red', alpha=0.07)
    lower = mean - 2.0*stddev
    upper = mean + 2.0*stddev
    plot.fill_between(times1, lower, upper, facecolor='yellow', alpha=0.08)

def plot_finalize(points):
    """Helper method for plotting

    Parameters
    ----------
        points : Integer
            Number of points displayed on the plot.
    """

    pts = str(points) + " pts"
    plot.legend(loc='best', title=pts)

    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()

    plot.tight_layout()
    plot.show()

def plot_initialize(fig, title, set_colors=False):
    """Helper method for plotting.

    Parameters
    ----------
        fig : magplotlib.figure.Figure
            Instance of plot.figure to attach the subplot to.

        title : String
            Title to appear on plot.
        set_colors : Boolean
            True to use colors defined here, False to use system defaults
    """
    subplot = fig.add_subplot(111)
    subplot.set_title(title)

    if set_colors:
        subplot.set_color_cycle(['red', 'orange', 'yellow'])

    subplot.xaxis.set_major_formatter(DateFormatter('%B %d, %Y'))
    subplot.xaxis.set_major_locator(DayLocator([5,15,25]))

def plot_intercepts(xIntercepts, yIntercepts, begin=0, cap=0):
    """Plot intercepts of line segments of least-squares-fit lines from the
    3 MHV sliding window.

    Parameters
    ----------
        begin : Integer
            Begin the plot at this point.
        cap : Integer
            Cap the plot at this many points.
        xIntercepts : List
            Epoch times.
        yIntercepts : List
            nT values.
    """

    if cap == 0:
        cap = len(xIntercepts)

    x = xIntercepts[begin:cap]
    i = 0
    for value in x:
        x[i] = datetime.datetime.fromtimestamp(value)
        i += 1
    y = yIntercepts[begin:cap]

    fig = plot.figure('Intercepts')
    plot_initialize(fig, 'Intercepts of line segments')

    plot.scatter(x, y, s=12, marker='D', color='green')
    plot.plot(x, y, label="Intercepts")

    plot_finalize(cap-begin)

def plot_lines(lines, begin=0, cap=0):
    """Plot line segments representing the least-squares-fit lines from the
    3 MHV sliding window.

    Parameters
    ----------
        begin : Integer
            Begin the plot at this point.
        cap : Integer
            Cap the plot at this many points.
        lines : Object
            Has 'x', 'y', 'slope', and 'intercept' defined for each point.
    """
    x = []
    y = []
    y2 = []

    count = 0
    if cap == 0:
        cap = len(lines)

    for line in lines:
        localX = line['x']
        localY = line['y']

        m = line['slope']
        b = line['intercept']
        localY2 = m * localX + b

        if count >= begin:
            localX = datetime.datetime.fromtimestamp(localX)
            x.append(localX)
            y.append(localY)
            y2.append(localY2)

        count += 1
        if count > cap:
            break

    fig = plot.figure('Line segments')
    plot_initialize(fig, 'Line segments from 3 consecutive MHVs', True)

    plot.scatter(x, y, s=8, marker='o', color='blue')
    plot.plot(x, y, label="x & y", color='blue')
    plot.plot(x, y2, label="x & y = mx + b", color='green')
    plot.scatter(x, y2, s=30, marker='h', color='green')

    count = 0
    for line in lines:
        localX = line['x']

        m = line['slope']
        b = line['intercept']
        localY2 = m * localX + b

        localX = datetime.datetime.fromtimestamp(localX)
        if count > 0 and count >= begin:
            xDiff = localX - prevX
            yDiff = localY2 - prevY
            plot.plot([prevX-xDiff, localX, localX+xDiff],
                    [prevY-yDiff, localY2, localY2+yDiff], lw=2)

        prevX = localX
        prevY = localY2

        count += 1
        if count > cap:
            break

    plot_finalize(count-begin-1)

def plot_months(months):
    """Plot monthly statistics before and after cleaning.

    Parameters
    ----------
        months : List <obspy.core.trac.Trace>
            List of monthly statistics before cleaning
    """
    fig = plot.figure('Average monthly nT Values')

    monthTitle = 'Monthly means (nT)'
    monthLabel = 'Daily mean range'

    plot_subplot(fig, 111, monthTitle, months, monthLabel, 'Month', 'blue', 's')

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
            Instance of plot.figure to attach the subplot to.
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
    times1 = date2num(times1)

    for time in list2:
        times2.append(time.stats.starttime)
        means2.append(time.stats.statistics['average'])
        ranges2.append(time.stats.statistics['maximum']
                - time.stats.statistics['minimum'])
    times2 = date2num(times2)

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

def plot_spline(x, y, xnew, ynew):
    """Plot cubic spline vs MHVs used to create it.

    Parameters
    ----------
        x : List
            MHV Epoch time values.
        y : List
            MHV nT values.
        xnew : List
            Spline Epoch time values.
        ynew : List
            Spline nT values.
    """

    times = []
    for time in x:
        times.append(datetime.datetime.utcfromtimestamp(time))
    times2 = []
    for time in xnew:
        times2.append(datetime.datetime.utcfromtimestamp(time))

    fig = plot.figure('MHV linear fit intercepts with cubic spline applied')

    plot.plot(times, y, color='blue', ls='none', marker='o', label='Intercepts')
    plot.plot(times2, ynew, color='green', linestyle='--', label='Cubic')
    plot.legend(loc='best', numpoints=1, frameon=False)

    mng = plot.get_current_fig_manager()
    mng.window.showMaximized()
    plot.tight_layout()
    plot.show()

def plot_subplot(fig, num, title, timeList, rLabel, mLabel,
                color='b', marker='s'):
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
    times = date2num(times)

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

def print_stats(times, interval, format='wide'):
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

def remove_sr_curve(channel, timeseries, spline):
    """Subtract the SR-curve from the original minute data.

    Parameters
    ----------
        channel : String
            H or D to specify which data channel is being used.
        spline : Dictionary
            Smoothed Solar-Regular (SR) curve (or non-K variation)
            Defined as an object/dictionary with a list of x's and a list of y's
        timeseries : obspy.core.Stream
            Stream object containing input data.

    Returns
    -------
        A trace of k-variation (SR-curve subtracted from every point).
    """
    trace = timeseries.select(channel=channel)[0]

    x = spline['x']
    y = spline['y']

    times = []
    for time in x:
        times.append(datetime.datetime.utcfromtimestamp(time))

    months = get_traces(trace, 'months')
    numMonths = len(months)

    if numMonths == 3:
        trace = months[1]
    elif numMonths != 1:
        raise Exception('SR curve will only be removed from 1 month of data.')

    if len(x) != len(y):
        raise Exception('X & Y must have the same number of points.')

    totalMinutes = trace.stats.npts

    if len(x) != totalMinutes:
        raise Exception('Spline and trace must have the same number of points.')

    i = 0
    for point in trace:
        trace.data[i] = point - y[i]
        i += 1

    return trace

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
