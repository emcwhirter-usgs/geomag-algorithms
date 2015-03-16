
# Algorithm Theoretical Basis for "K-USGS"

E. Joshua Rigler &lt;[erigler@usgs.gov](mailto:erigler@usgs.gov)&gt;
Edward A. McWhirter Jr. &lt;[emcwhirter@usgs.gov](mailto:emcwhirter@usgs.gov)&gt;


## Summary

K-Indices are used as an approximate measure of magnetic activity at an
observatory over a 3-hour window. A scale is adopted for each observatory based
on the typical distribution of magnetic activity at that location. The scale is
divided into intervals that a translated to values 0 through 9. The K-USGS
algorithm differs from other K-Indices algorithms in one regard. A non-K
variation curve is calculated for each day instead of using the typical Solar
Quiet (SQ) algorithm as an input, which means that K-USGS can be calculated with
just the raw time-series data as inputs.


## Background and Motivation

The USGS has been producing digital K-Indices since 1979. The K-USGS algorithm
was written to attempt to simulate the previously used hand-scaled process as
closely as possible. This method is described in
"[An Evaluation of Digitally Derived K-Indices][]", J. Geomag. Geoelectr., 39,
1987 by Lanny R. Wilson pages 97-109. The algorithm was tested and found
acceptable in 1990 using IAGA WG-5 data. It was tested again in 2010 using a
subset of USGS data, and was found to still be acceptable for producing
K-Indices.

[An Evaluation of Digitally Derived K-Indices]: https://www.jstage.jst.go.jp/article/jgg1949/39/2/39_2_97/_article

The 3-hour K-Index was introduced by Bartels (1938) as a measure of irregular
and rapid storm-time magnetic activity. This same process was defined in detail
by Mayaud in 1957. It was designed to be insensitive to longer term components
of magnetic variation and to normalize the occurence frequency of individual K
values among many observatories, over many years. Thus, with this method, there
is a separate K-Index for each observatory. It has come to be used much more
generally as a measure of the magnetic activity at an observatory at any given
time, as opposed to just during magnetic storms.

In general, the K-USGS algorithm consists of these steps:
 1. Determine the number of minutes and means to use based on the selected time
    range.
 2. Determine an acceptable range for the selected data. Replace any values
    outside of that range with DEAD values.
 3. Replace all DEAD values with values interpolated from valid values near the
    DEAD values.
 4. Create line segments and intercepts.
 5. Create a Spline.
 6. Create a Cubic Spline.

> Subset of USGS used for testing in 2010 included:
> CMO 1992; FRD 1985-1994,1997; GUA 1992-1994; SJG 1992-1994; TUC 1992-1993


## Math and Theory

### Inputs ###
- Station ID
- Start Time (ST)
- End Time (ET)
- H Time-Series Data
- D Time-Series Data

### Static Defaults and Configurable Values ###
* `MeansPerSegment  3`                          // Configurable 3 or 4
* `MinutesPerMean   60`                         // Configurable from 12 to 180
* `MinutesPerMean3HrList { 12, 20, 30, 36, 60 }`
* `MinutesPerMean24List { 12, 20, 24, 30, 36, 48, 60, 72, 90 }`
* `StationK9Value   2500`                       // Configurable from 50 to 3500
* `StdDevFor1MinuteValues 2.5`                  // Configurable 0.50 to 7.00
* `StdDevForMeanValues    1.8`                  // Configurable 0.50 to 7.00

### Computation ###
The algorithm for computing K-USGS is based on a CubicSpline. Before this can
be properly calculated, all data gaps, dead values and bad data must be
replaced with appropriate values. This is done by first calculating some
statistics about the selected data.

#### Initial Statistics ####

Determine the number of minutes to use.

* Equation 1: `DataDays = 1 + ET - ST`
* Equation 2: `DataMinutes = DataDays * 1440 minutes per day`

An extra buffer of data is added to the front and back of the data so that end
the end points of the selected data can be calculated.

* Equation 3: `DataBuffer = 2 * 60 minutes per mean`
* Equation 4: `DataWidth = DataMinutes + 2 * DataBuffer`

An extra point is added to the means to act as an anchor to the end of the data.

* Equation 5: `NumberOfMeans = 1 + DataWidth / 60 minutes per mean`

#### Counting Statistics ####

Next, loop over the means while skipping DEAD values. Each mean has 60 1-Minute
Values associated with it. Find:

Among all of the means.

* `MeanCount`
* `MeanSum`
* `HighestMeanValue`
* `LowestMeanValue`
* `SumSqrsMeans1MinValues`       // Sum of the squares of valid Mean Values
* `HighestRange`
* `LowestRange`

Within each mean.

* `HighMinute`
* `LowMinute`
* `MinuteCount`
* `MinuteSum`

If more than half of the minute values for a mean are DEAD, the mean is DEAD.
If `HighMinute - LowMinute > 5 * StationK9Value / 500`, the mean is DEAD.

#### Computed Statistics ####

Use the counting statistics to establish some computed values.

* `RangeOfMeansSum`              // Sum of Ranges of valid Means
* `SumSqrsRangeOfMeansSum`       // Sum of the squares of RangeOfMeansSum
* `DataMean`                     // X, Mean of the location of each point
* `MeanOfSqrs1MinValues`         // Mean of squares of 1 Minute Values
* `ValueMean`                    // Y, Mean of the values at each point
* `MeanOfSqrsMeanValues`         // Mean of squares of "Mean" values
* `RangeMean`                    // Mean of Ranges of "Means"
* `MeanOfSqrsMeanRange`          // Mean of squares of Ranges of "Means"

* Equation 6: `DataStdDev = sqrt(MeanOfSqrs1MinuteValues - DataMean^2)`
* Equation 7: `DataMax = DataMean + 2.5 * DataStdDev`
* Equation 8: `DataMin = DataMean - 2.5 * DataStdDevn`

* Equation 9: `ValueStdDev = sqrt(MeanOfSqrsMeanValues - ValueMean^2)`
* Equation 10: `ValueMax = ValueMean + 1.8 * ValueStdDev`
* Equation 11: `ValueMin = ValueMean - 1.8 * ValueStdDev`

* Equation 12: `RangeStdDev = sqrt(MeanOfSqrsMeanRange - RangeMean^2)`
* Equation 13: `RangeMax = RangeMean + 1.8 * RangeStdDev`

* Equation 14: `DataStandardError = DataStdDev / DataMeanCount`
* Equation 15: `DataMedian = LowestMeanValue + (HighestMeanValue - LowestMeanValue) / 2.0`

* Equation 16: `ValueStdError = ValueStdDev / MeanCount`
* Equation 17: `ValueMedian = LowMinute + (HighMinute - LowMinute) / 2.0`

* Equation 18: `RangeStdError = RangeStdDev / MeanCount`
* Equation 19: `RangeMedian = RangeLowValue + (RangeHighValue - RangeLowValue) / 2.0`

* Equation 20: `DataInterval = (HighestMeanValue - LowestMeanValue) / 10.0`
* Equation 21: `ValueInterval = (HighMinute - LowMinute) / 10.0`
* Equation 22: `RangeInterval = (HighestRange - LowestRange) / 10.0`

X (Time) values are fractional hours centered on the minute data used for the
mean. So, if MeansPerMinute is 60, the 1st range is 0.5hr.

* Equation 23: `Interval = MinutesPerMean / 60 min/hr = 60.0 / 60.0 = 1.0`
* Equation 24: `Offset = Interval / 2.0 = 1.0 / 2.0 = 0.5`

Once all of these values are found, loop over all of the means and the minutes
in each mean, still skipping DEAD values, to reject values based on found
statistics. If the MinuteValue is between HighestMeanValue and LowestMeanValue,
then:

* Equation 25: `ValueIndex = Floor( (MinuteValue-LowestMeanValue) / DateInterval )`

Otherwise the value is outside of the acceptable range and should be rejected.
Rejected values are replaced with DEAD values.

Also reject the mean and replace it with a DEAD value if more than half of the
minutes are DEAD, if the range of minutes is greater than RangeMax or if the
MeanValue is not between ValueMin and ValueMax.

Next, loop over the data again and replace all DEAD values with an
interpolated value based on the closest points in the data moving forward. If
there are no non-DEAD values forward in time, then look backwards instead.

<TODO: compare this to code again to ensure nothing was lost>
<TODO: include images if needed>
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
