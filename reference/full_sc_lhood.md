# calculates full self-correcting log-likelihood

calculates full self-correcting log-likelihood

## Usage

``` r
full_sc_lhood(xgrid, ygrid, tgrid, tobs, data, params, bounds)
```

## Arguments

- xgrid:

  a vector of grid values for x.

- ygrid:

  a vector of grid values for y.

- tgrid:

  a vector of grid values for t.

- tobs:

  a vector of observed values for t.

- data:

  a matrix of times and locations.

- params:

  a vector of parameters.

- bounds:

  a vector of bounds for time, x, and y.

## Value

evaluation of full log-likelihood.
