# calculates acceptance for thinning mechanism during simulation

calculates acceptance for thinning mechanism during simulation

## Usage

``` r
thin_st_fast(data, params)
```

## Arguments

- data:

  NumericMatrix with columns (time, x, y). Assumed sorted by time
  ascending.

- params:

  NumericVector length 3: (alpha3, beta3, gamma3

## Value

LogicalVector length n of whether to keep each point (true) or thin it
(false).
