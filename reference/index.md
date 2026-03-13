# Package index

## Main functions

Main functions for estimating kinematic parameters from midline data.
These are more or less in the order you would use them. Functions ending
in `_df` operate on a data frame and often add multiple columns to the
data frame. The data frame should have two identifier columns:

- `frame`: Identifies instants in time
- `point`: Identifies locations along the body, so that they can be
  ordered from head to tail.

To use identifier columns with different names, the functions take
options `.frame` and `.point`. To change the names of the output
columns, the functions use an option `.out`. `.out` can be a vector
containing the names of the columns, or a named list.

- [`arclength()`](https://tytell.github.io/fishmechr/reference/arclength.md)
  : Calculate arc length along a 2D curve
- [`interpolate_points_df()`](https://tytell.github.io/fishmechr/reference/interpolate_points_df.md)
  : Interpolates and smooths a 2D curve at new arc length
- [`interpolate_width()`](https://tytell.github.io/fishmechr/reference/interpolate_width.md)
  : Interpolates and scales fish body width
- [`get_midline_center_df()`](https://tytell.github.io/fishmechr/reference/get_midline_center_df.md)
  : Gets the center of a midline for many midlines in a data frame
- [`get_primary_swimming_axis_df()`](https://tytell.github.io/fishmechr/reference/get_primary_swimming_axis_df.md)
  : Gets the primary swimming axis for many midlines
- [`curvature()`](https://tytell.github.io/fishmechr/reference/curvature.md)
  : Estimates curvature for a single curve
- [`hilbert_phase()`](https://tytell.github.io/fishmechr/reference/hilbert_phase.md)
  : Compute phase of an oscillation using the Hilbert transform
- [`peak_phase()`](https://tytell.github.io/fishmechr/reference/peak_phase.md)
  : Compute phase of an oscillation by locating peaks and zero
  crossings.
- [`get_frequency()`](https://tytell.github.io/fishmechr/reference/get_frequency.md)
  : Estimates the cycle frequency based on time and phase
- [`get_wavelength()`](https://tytell.github.io/fishmechr/reference/get_wavelength.md)
  : Computes the body wavelength based on the phase at each point and
  the arc length
- [`get_body_cycle_numbers_df()`](https://tytell.github.io/fishmechr/reference/get_body_cycle_numbers_df.md)
  : Gets oscillation cycle numbers for a midline data set

## Helper functions

Other useful functions for estimating kinematic parameters.

- [`apply_filter()`](https://tytell.github.io/fishmechr/reference/apply_filter.md)
  : Apply a filter constructed with build_filter
- [`build_filter()`](https://tytell.github.io/fishmechr/reference/build_filter.md)
  : Constructs a smoothing filter
- [`deriv()`](https://tytell.github.io/fishmechr/reference/deriv.md) :
  Estimate first or second derivatives for dy/dx.
- [`find_gaps_df()`](https://tytell.github.io/fishmechr/reference/find_gaps_df.md)
  : Find gaps in a data series
- [`fishwidth`](https://tytell.github.io/fishmechr/reference/fishwidth.md)
  : Fish body width profiles
- [`get_cycle_numbers()`](https://tytell.github.io/fishmechr/reference/get_cycle_numbers.md)
  : Gets cycle numbers from a phase variable
- [`get_primary_swimming_axis()`](https://tytell.github.io/fishmechr/reference/get_primary_swimming_axis.md)
  : Gets the main swimming axis from a midline
- [`get_volume()`](https://tytell.github.io/fishmechr/reference/get_volume.md)
  : Gets the volume of segments of a cylindrical body with elliptical
  cross section
- [`interpolate_peak_location()`](https://tytell.github.io/fishmechr/reference/interpolate_peak_location.md)
  : Interpolate the location of a peak based on three points
- [`interpolate_points_frame()`](https://tytell.github.io/fishmechr/reference/interpolate_points_frame.md)
  : Interpolates x and y points on a curve to different arc lengths
- [`lampreydata`](https://tytell.github.io/fishmechr/reference/lampreydata.md)
  : Lamprey midline data
- [`pivot_kinematics_longer()`](https://tytell.github.io/fishmechr/reference/pivot_kinematics_longer.md)
  : Pivots a kinematics dataset into long format
- [`skip_na()`](https://tytell.github.io/fishmechr/reference/skip_na.md)
  : Skip NAs when running a function on a vector
- [`smooth_point()`](https://tytell.github.io/fishmechr/reference/smooth_point.md)
  : Applies a smoothing spline to a data series, potentially with gaps
- [`smooth_points_df()`](https://tytell.github.io/fishmechr/reference/smooth_points_df.md)
  : Smooths locations of points over time
- [`xmucosusdata`](https://tytell.github.io/fishmechr/reference/xmucosusdata.md)
  : Prickleback tracking data
- [`zebrafish_shape`](https://tytell.github.io/fishmechr/reference/zebrafish_shape.md)
  : Zebrafish body shape
- [`zfish_goodframes`](https://tytell.github.io/fishmechr/reference/zfish_goodframes.md)
  : Zebrafish good frame ranges
- [`zfishdata`](https://tytell.github.io/fishmechr/reference/zfishdata.md)
  : Zebrafish keypoint tracking data
