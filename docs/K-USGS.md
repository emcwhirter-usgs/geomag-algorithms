
# Algorithm Theoretical Basis for "K-USGS"

E. Joshua Rigler &lt;[erigler@usgs.gov](mailto:erigler@usgs.gov)&gt;
Edward A. McWhirter Jr. &lt;[emcwhirter@usgs.gov](mailto:emcwhirter@usgs.gov)&gt;


## Summary

K-Indices are used as an approximate measure of magnetic activity at an
observatory. A scale is adopted for each observatory based on the typical
distribution of magnetic activity at that location. The scale is devided into
intervals that a translated to values 0 through 9. The K-USGS algorithm differs
from other K-Indices estimates a non-K variation curve for each day instead of
using the typical Solar Quiet (SQ) algorithm as an input.


## Background and Motivation

The USGS has been producing digital K-Indices since 1979. The "K-USGS" algorithm
was written to attempt to simulate the previously used hand-scaled process as
closely as possible. This method is described in "[An Evaluation of Digitally
Derived K-Indices][]", J. Geomag. Geoelectr., 39, 1987 by Lanny R. Wilson pages
97-109. The algorithm was tested and found acceptable in 1990 using IAGA WG-5
data. It tested again in 2010 using a subset of USGS data, and was found to
still be acceptable.

[An Evaluation of Digitally Derived K-Indices]: https://www.jstage.jst.go.jp/article/jgg1949/39/2/39_2_97/_article

The 3-hour K-Index was introduced by Bartels (1938) as a measure of irregular
and rapid storm-time magnetic activity. This same process was defined in detail
by Mayaud in 1957. It was designed to be insensitive to longer term components
of magnetic variation and to normalize the occurence frequency of individual K
values among many observatories, over many years. Thus, with this method, there
is a separate K-Index for each observatory. It has come to be used much more
generally as a measure of the magnetic activity at an observatory at any given
time, as opposed to just during magnetic storms.

In general the K-USGS algorithm consists of these steps:
1. Determine the number of minutes and means to use based on the selected time
   range.
1. Determine an acceptable range for the selected data. Replace any values
   outside of that range with DEAD values.
1. Replace all DEAD values with values interpolated from valid values near the
   DEAD values.
1. Create line segments and intercepts.
1. Create a Spline.
1. Create a Cubic Spline.

> Subset of USGS used for testing in 2010 included:
> CMO 1992; FRD 1985-1994,1997; GUA 1992-1994; SJG 1992-1994; TUC 1992-1993


## Math and Theory

### Inputs ###
- Station ID
- Start Time
- End Time
- H Time-Series Data
- D Time-Series Data

The algorithm for computing K-USGS is based on a CubicSpline. Before this can
be properly calculated, all data gaps, dead values and bad data must be
replaced with values.

<TODO: define equations>
<TODO: explain equations and why they are valid for this use>
<TODO: include image if needed>
<TODO: discuss edge cases>


## Practical Considerations

### Magnetic Intensity Units

It is understood that all raw data inputs are provided in units of nanoTesla
(nT). Of course this is not required for the equations to be valid, but it is
incumbent on the programmer to make sure all input data units are the same, and
that output units are defined accurately.

### Data Flags

It should go without saying that bad data in one coordinate system is bad data
in another. However, on occasion, operational USGS Geomagnetism Program code has
been discovered where coordinate transformations were applied before checking
data flags. This is not an issue if data flags are NaN (not-a-number values),
but more typical for Geomag data, these are values like 99999, which can lead to
seemingly valid, but erroneous values at times when the raw data were known to
be bad.

> Note: this library internally represents data gaps as NaN, and factories
> convert to this where possible.
