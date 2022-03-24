# PiecewiseBilinear

!syntax description /Functions/PiecewiseBilinear

## Description

The `PiecewiseBilinear` function reads a csv file and interpolates values based on the
data in the file.  The interpolation is based on x-y pairs.  If `axis` is given, time is
used as the y index.  Either `xaxis` or `yaxis` or both may be given.  Time is used as
the other index if one of them is not given.  If `radius` is given, `xaxis` and `yaxis`
are used to orient a cylindrical coordinate system, and the x-y pair used in the query
will be the radial coordinate and time.

The csv file `data_file` format expected is:

- first line holds the `x` values.
- each subsequent line holds the `y` value then the list of `z` values for this `y` and
  all values of `x`.

The csv file `data_file` may be substituted by specifying the `x`, `y` and `z` parameters.

## Example Input Syntax

In this example, the csv file below `fred.csv` is read with different parameters to achieve
different function definitions. Explanations are provided in the input file snippet.

!listing test/tests/utils/2d_linear_interpolation/fred.csv

!listing test/tests/utils/2d_linear_interpolation/2d_linear_interpolation_test.i block=Functions

!syntax parameters /Functions/PiecewiseBilinear

!syntax inputs /Functions/PiecewiseBilinear

!syntax children /Functions/PiecewiseBilinear
