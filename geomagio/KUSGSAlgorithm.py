"""
    Creates K-USGS Indices from H and D time-series data.
"""

import numpy
import obspy.core
import obspy.core.utcdatetime as UTC
from Algorithm import Algorithm
from geomagio import TimeseriesFactoryException

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
    print str(trace) + "\n"

    # type = <'numpy.ndarray'>
    H = trace.data

    numMinutes = H.size
    # This algorithm operates on entire calendar days of 1-Minute values.
    if (numMinutes % MINUTESPERDAY) != 0:
        raise TimeseriesFactoryException(
                'Entire calendar days of minute data required for K.')

    average = numpy.nanmean(H)
    stdDev = numpy.nanstd(H)
    minMinute = numpy.amin(H)
    maxMinute = numpy.amax(H)
    rangeMinutes = maxMinute - minMinute

    # hours = []
    # type = <class 'obspy.core.utcdatetime.UTCDateTime'>
    starttime = trace.stats.starttime
    endtime = trace.stats.endtime
    # type = <type 'list'>
    days = get_days(starttime, endtime)

    i = 0
    # type = <type 'numpy.timedelta64'>
    oneDay = numpy.timedelta64(1, 'D')
    oneMinute = numpy.timedelta64(1, 'm')
    for day in days:
        begin = numpy.datetime64(starttime) + i * oneDay
        begin = UTC.UTCDateTime(str(begin))

        end = numpy.datetime64(begin) + oneDay - oneMinute
        end = UTC.UTCDateTime(str(end))

        daily_stats(day, trace.slice(begin, end))
        i += 1
    # for day in days:
    #     hours.append(get_hours(day))

    # print "H Trace:"
    # print H
    print "Total # of Minutes    : " + str(numMinutes)
    print "Average of all Minutes: " + str(average) + "nT"
    print "Std Dev of all Minutes: " + str(stdDev)
    print "Range of all Minutes  : " + str(rangeMinutes) + "nT"

def daily_stats(day, trace):
    print "  Daily stats for " + str(day)
    print "    " + str(trace) + "\n"
    H = trace.data

    # numMinutes = H.length
    # average = numpy.nanmean(H)
    # stdDev = numpy.nanstd(H)
    # minMinute = numpy.amin(H)
    # maxMinute = numpy.amax(H)
    # rangeMinutes = maxMinute - minMinute

    # print "    Number of Minutes:        " + str(numMinutes)
    # print "    Average of day's Minutes: " + str(average) + "nT"
    # print "    Std Dev of day's Minutes: " + str(stdDev)
    # print "    Range of day's Minutes:   " + str(rangeMinutes) + "nT"

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

def get_hours(day):
    """
        Get Mean Hourly Values (MHVs).
    """
    hours = []
    delta = numpy.timedelta64(1, 'h')
    date = numpy.datetime64(day)

    for i in range(0, 24):
        hour = date + i*delta
        hours.append(hour)

    return hours
