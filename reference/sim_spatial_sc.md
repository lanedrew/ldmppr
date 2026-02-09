# Simulate the spatial component of the self-correcting model (faster)

Simulate the spatial component of the self-correcting model (faster)

## Usage

``` r
sim_spatial_sc(M_n, params, nsim_t, xy_bounds)
```

## Arguments

- M_n:

  a vector of (x,y)-coordinates for anchor/first point.

- params:

  a vector of parameters (alpha_2, beta_2).

- nsim_t:

  number of points to simulate.

- xy_bounds:

  vector: (ax, bx, ay, by).

## Value

a matrix nsim_t x 2 of point locations (x,y).
