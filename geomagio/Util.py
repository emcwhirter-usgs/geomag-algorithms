import urllib2
import numpy
import obspy


class ObjectView(object):
    """
    Wrap a dictionary so its properties can be accessed as an object.

    Parameters
    ----------
    d : dictionary
        The dictionary to wrap.
    """
    def __init__(self, d):
        self.__dict__ = d

    def __str__(self):
        """
        Override string representation to output wrapped dictionary.
        """
        return str(self.__dict__)


def read_url(url):
    """Open and read url contents.

    Parameters
    ----------
    url : str
        A urllib2 compatible url, such as http:// or file://.

    Returns
    -------
    str
        contents returned by url.

    Raises
    ------
    urllib2.URLError
        if any occurs
    """
    response = urllib2.urlopen(url)
    content = None
    try:
        content = response.read()
    except urllib2.URLError, e:
        print e.reason
        raise
    finally:
        response.close()
    return content


def create_empty_trace(trace, channel):
    """
    Utility to create a trace containing the given numpy array.

    Parameters
    ----------
    stream: obspy.core.stream

    Returns
    -------
    obspy.core.Trace
        Trace a duplicated empty channel.
    """
    stats = obspy.core.Stats(trace.stats)
    stats.channel = channel
    count = len(trace.data)
    numpy_data = numpy.full((count), numpy.nan)
    return obspy.core.Trace(numpy_data, stats)
