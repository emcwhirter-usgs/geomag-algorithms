
from cStringIO import StringIO
from geomagio import TimeseriesFactoryException, ChannelConverter
import numpy
import IAGA2002Parser
import textwrap
from datetime import datetime


class IAGA2002Writer(object):
    """IAGA2002 writer.
    """

    def __init__(self, empty_value=IAGA2002Parser.NINES):
        self.empty_value = empty_value

    def write(self, out, timeseries, channels):
        stats = timeseries[0].stats
        out.write(self._format_headers(stats, channels))
        out.write(self._format_comments(stats))
        out.write(self._format_channels(channels, stats.station))
        out.write(self._format_data(timeseries, channels))
        pass

    def _format_headers(self, stats, channels):
        buf = []
        buf.append(self._format_header('Format', 'IAGA-2002'))
        buf.append(self._format_header('Source of Data', stats.agency_name))
        buf.append(self._format_header('Station Name', stats.station_name))
        buf.append(self._format_header('IAGA CODE', stats.station))
        buf.append(self._format_header('Geodetic Latitude',
                stats.geodetic_latitude))
        buf.append(self._format_header('Geodetic Longitude',
                stats.geodetic_longitude))
        buf.append(self._format_header('Elevation', stats.elevation))
        buf.append(self._format_header('Reported', ''.join(channels)))
        buf.append(self._format_header('Sensor Orientation',
                stats.sensor_orientation))
        buf.append(self._format_header('Digital Sampling',
                str(1 / stats.sensor_sampling_rate) + ' second'))
        buf.append(self._format_header('Data Interval Type',
                stats.data_interval_type))
        buf.append(self._format_header('Data Type', stats.data_type))
        return ''.join(buf)

    def _format_comments(self, stats):
        # build comments
        comments = []
        if 'declination_base' in stats:
            comments.append('DECBAS               {:<8d}'
                    '(Baseline declination value in tenths of minutes East'
                    ' (0-216,000)).'.format(stats.declination_base))
        if 'filter_comments' in stats:
            comments.extend(stats.filter_comments)
        if 'comments' in stats:
            comments.extend(stats.comments)
        if 'is_intermagnet' in stats and stats.is_intermagnet:
            comments.append('Final data will be available on the' +
                    ' INTERMAGNET DVD.')
            comments.append('Go to www.intermagnet.org for details on' +
                    ' obtaining this product.')
        if 'conditions_of_use' in stats:
            comments.append('CONDITIONS OF USE: ' + stats.conditions_of_use)
        # generate comment output
        buf = []
        for comment in comments:
            buf.append(self._format_comment(comment))
        return ''.join(buf)

    def _format_header(self, name, value):
        prefix = ' '
        suffix = ' |\n'
        return ''.join((prefix, name.ljust(23), value.ljust(44), suffix))

    def _format_comment(self, comment):
        buf = []
        prefix = ' # '
        suffix = ' |\n'
        lines = textwrap.wrap(comment, 65)
        for line in lines:
            buf.extend((prefix, line.ljust(65), suffix))
        return ''.join(buf)

    def _format_channels(self, channels, iaga_code):
        """Format channel header line.

        Parameters
        ----------
        channels : sequence
            list and order of channel values to output.
        iaga_code : str
            observatory code, which is prefixed to channel name in output.

        Returns
        -------
        str
            Channel header line as a string (including newline)
        """
        if len(iaga_code) != 3:
            raise TimeseriesFactoryException(
                    'iaga_code "{}" is not 3 characters'.format(iaga_code))
        if len(channels) != 4:
            raise TimeseriesFactoryException(
                    'more than 4 channels {}'.format(channels))
        buf = ['DATE       TIME         DOY']
        for channel in channels:
            if len(channel) != 1:
                raise TimeseriesFactoryException(
                        'channel "{}" is not 1 character'.format(channel))
            buf.append('     %s%s ' % (iaga_code, channel))
        buf.append('  |\n')
        return ''.join(buf)

    def _format_data(self, timeseries, channels):
        """Format all data lines.

        Parameters
        ----------
        timeseries : obspy.core.Stream
            stream containing traces with channel listed in channels
        channels : sequence
            list and order of channel values to output.
        """
        buf = []
        if timeseries.select(channel='D'):
            d = timeseries.select(channel='D')
            d[0].data = ChannelConverter.get_minutes_from_radians(d[0].data)
        traces = [timeseries.select(channel=c)[0] for c in channels]
        starttime = float(traces[0].stats.starttime)
        delta = traces[0].stats.delta
        for i in xrange(len(traces[0].data)):
            buf.append(self._format_values(
                datetime.utcfromtimestamp(starttime + i * delta),
                (t.data[i] for t in traces)))
        return ''.join(buf)

    def _format_values(self, time, values):
        """Format one line of data values.

        Parameters
        ----------
        time : datetime
            timestamp for values
        values : sequence
            list and order of channel values to output.
            if value is NaN, self.empty_value is output in its place.

        Returns
        -------
        unicode
            Formatted line containing values.
        """
        tt = time.timetuple()
        return '{0.tm_year:0>4d}-{0.tm_mon:0>2d}-{0.tm_mday:0>2d} ' \
                '{0.tm_hour:0>2d}:{0.tm_min:0>2d}:{0.tm_sec:0>2d}.{1:0>3d} ' \
                '{0.tm_yday:0>3d}   ' \
                '{2:10.2f}{3:10.2f}{4:10.2f}{5:10.2f}\n'.format(
                tt, int(time.microsecond / 1000),
                *[self.empty_value if numpy.isnan(val) else val
                        for val in values])

    @classmethod
    def format(self, timeseries, channels):
        """Get an IAGA2002 formatted string.

        Calls write() with a StringIO, and returns the output.

        Parameters
        ----------
        timeseries : obspy.core.Stream

        Returns
        -------
        unicode
          IAGA2002 formatted string.
        """
        out = StringIO()
        writer = IAGA2002Writer()
        writer.write(out, timeseries, channels)
        return out.getvalue()
