# Small Example Data

A small example dataset for testing and examples consisting of 121
observations in a (25m x 25m) square domain.

## Usage

``` r
data("small_example_data")
```

## Format

\## `small_example_data` A data frame with 121 rows and 3 columns:

- x:

  x coordinate

- y:

  y coordinate

- size:

  Size

## Source

Simulated dataset. Code to generate it can be found in
`data_raw/small_example_data.R`.

## Details

The dataset was generated using the example raster data and an
exponential decay size function.

The full code to generate this dataset is available in the package's
`data_raw` directory.
