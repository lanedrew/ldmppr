# Simulate the spatial component of the self-correcting model

Simulate the spatial component of the self-correcting model

## Usage

``` r
sim_spatial_sc(M_n, params, nsim_t, xy_bounds)
```

## Arguments

- M_n:

  a vector of (x,y)-coordinates for largest point.

- params:

  a vector of parameters (alpha_2, beta_2).

- nsim_t:

  number of points to simulate.

- xy_bounds:

  vector of lower and upper bounds for the domain (2 for x, 2 for y).

## Value

a matrix of point locations in the (x,y)-plane.
