# Evaluate reference self-correcting log-likelihood

Reference implementation used for parity checks and validation. For
production estimation, prefer
[`full_sc_lhood_fast()`](https://lanedrew.github.io/ldmppr/reference/full_sc_lhood_fast.md).

## Usage

``` r
full_sc_lhood(xgrid, ygrid, tgrid, tobs, data, params, bounds)
```

## Arguments

- xgrid:

  NumericVector of x-grid values.

- ygrid:

  NumericVector of y-grid values.

- tgrid:

  NumericVector of integration-time grid values.

- tobs:

  NumericVector of observed event times.

- data:

  NumericMatrix with columns (time, x, y), sorted by nondecreasing time.

- params:

  NumericVector of model parameters (alpha1, beta1, gamma1, alpha2,
  beta2, alpha3, beta3, gamma3).

- bounds:

  NumericVector of integration bounds (bt, bx, by).

## Value

Full self-correcting log-likelihood value.
