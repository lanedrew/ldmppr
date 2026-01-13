# Simulate the temporal component of the self-correcting model

Simulate the temporal component of the self-correcting model

## Usage

``` r
sim_temporal_sc(Tmin = 0, Tmax = 1, params = as.numeric(c(0, 0, 0)))
```

## Arguments

- Tmin:

  minimum time value.

- Tmax:

  maximum time value.

- params:

  a vector of parameters (alpha_1, beta_1, gamma_1).

## Value

a vector of thinned and unthinned temporal samples.
