# CSVReader

The CSVReader reads [CSV](https://en.wikipedia.org/wiki/Comma-separated_values) data from a file and
converts each column into a VectorPostprocessor vector. This object uses the
[DelimitedFileReader](MooseUtils.md#delimitedfilereader) utility to perform the reading of the file.

## Example Input Syntax

This example shows how the following CSV file containing dates can be loaded into
a simulation using a `CSVReader`.

!listing test/tests/vectorpostprocessors/csv_reader/example.csv

!listing test/tests/vectorpostprocessors/csv_reader/read.i block=VectorPostprocessors

!syntax parameters /VectorPostprocessors/CSVReader

!syntax inputs /VectorPostprocessors/CSVReader

!syntax children /VectorPostprocessors/CSVReader
