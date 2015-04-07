"""
    Creates K-USGS Indices from H and D time-series data.
"""

import numpy
import obspy.core
import obspy.core.utcdatetime as UTC
from Algorithm import Algorithm
from geomagio import TimeseriesFactoryException

MINUTESPERHOUR = 60
MINUTESPERDAY = 1440  # 24h * 60min

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

    # type = <type 'list'>
    days = get_days(trace.stats.starttime, trace.stats.endtime)

    # type = <type 'numpy.timedelta64'>
    oneDay = numpy.timedelta64(1, 'D')
    oneMinute = numpy.timedelta64(1, 'm')
    oneHour = numpy.timedelta64(1, 'h')

    dailyStats = []
    hours = []
    for day in days:
        end = numpy.datetime64(day) + oneDay - oneMinute
        end = UTC.UTCDateTime(str(end))

        thisDay = trace.slice(day, end)
        if thisDay.stats.npts != MINUTESPERDAY:
            raise TimeseriesFactoryException(
                    'Entire calendar days of minute data required for K.')

        statistics(thisDay)
        dailyStats.append(thisDay)

    for day in days:
        hours.append(get_hours(day))

    hourlyStats = []
    for day in hours:
        for hour in day:
            end = numpy.datetime64(hour) + oneHour - oneMinute
            # TODO Look into using the raw time value instead of a string
            end = UTC.UTCDateTime(str(end))

            thisHour = trace.slice(hour, end)
            if thisHour.stats.npts != MINUTESPERHOUR:
                raise TimeseriesFactoryException(
                        '1 Hour should have 60 minutes.')

            statistics(thisHour)
            hourlyStats.append(thisHour)

    print_days(dailyStats)
    # print_hours(hourlyStats)

    print_all(trace.stats)


def get_hours(day):
    """
        Get Mean Hourly Values (MHVs).
    """
    hours = []
    oneHour = numpy.timedelta64(1, 'h')
    date = numpy.datetime64(day)

    for i in range(0, 24):
        hour = date + i * oneHour
        hour = UTC.UTCDateTime(str(hour))

        hours.append(hour)

    return hours

def statistics(trace):
    """
        Calculate average, standard deviation, minimum and maximum on given
        trace, add them to a 'statistics' object and attach them to the trace.
    """
    H = trace.data

    average = numpy.nanmean(H)
    standardDeviation = numpy.nanstd(H)
    minimum = numpy.amin(H)
    maximum = numpy.amax(H)
    Range = maximum - minimum

    trace.stats.statistics = {
        'average': average,
        'standarddeviation': standardDeviation,
        'minimum': minimum,
        'maximum': maximum
    }

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
    print "Total # of Minutes    : " + str(stats.npts)
    print "Average of all Minutes: " + str(statistics['average']) + "nT"
    print "Std Dev of all Minutes: " + str(statistics['standarddeviation'])
    print "Range of all Minutes  : " + str(statistics['maximum'] - statistics['minimum']) + "nT"

def print_days(dailyStats):
    ### Example output ###
    #  Day          : 2013-12-31T00:00:00.000000Z
    #  Daily Average: 20894.2173562
    #  Daily Std Dev: 9.39171243572
    #  Daily Range  : 44.319
    for day in dailyStats:
        stats = day.stats
        statistics = stats.statistics
        print "  Day          : " + str(stats.starttime)
        print "  Daily Average: " + str(statistics['average'])
        print "  Daily Std Dev: " + str(statistics['standarddeviation'])
        print "  Daily Range  : " + str(statistics['maximum'] - statistics['minimum']) + "\n"

def print_hours(hourlyStats):
    ### Example output ###
    #    Hour          : 2014-01-02T17:00:00.000000Z
    #    Hourly Average: 20855.7571167
    #    Hourly Std Dev: 10.1907743067
    #    Hourly Range  : 36.883
    for hour in hourlyStats:
        stats = hour.stats
        statistics = stats.statistics
        print "    Hour          : " + str(stats.starttime)
        print "    Hourly Average: " + str(statistics['average'])
        print "    Hourly Std Dev: " + str(statistics['standarddeviation'])
        print "    Hourly Range  : " + str(statistics['maximum'] - statistics['minimum']) + "\n"
