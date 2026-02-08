# Simulate from the self-correcting model

Allows the user to simulate a realization from the self-correcting model
given a set of parameters and a point to condition on.

## Usage

``` r
simulate_sc(
  t_min = 0,
  t_max = 1,
  sc_params = NULL,
  anchor_point = NULL,
  xy_bounds = NULL
)
```

## Arguments

- t_min:

  minimum value for time.

- t_max:

  maximum value for time.

- sc_params:

  a vector of parameter values corresponding to
  \\(\alpha_1,\beta_1,\gamma_1,\alpha_2,\beta_2,\alpha_3,\beta_3,\gamma_3)\\
  (i.e., alpha_1, beta_1, gamma_1, alpha_2, beta_2, alpha_3, beta_3,
  gamma_3).

- anchor_point:

  a vector of (x,y) coordinates of point to condition on.

- xy_bounds:

  a vector of domain bounds (2 for x, 2 for y).

## Value

a list containing the thinned and unthinned simulation realizations.

## Examples

``` r
# Specify the generating parameters of the self-correcting process
generating_parameters <- c(2, 8, .02, 2.5, 3, 1, 2.5, .2)

# Specify an anchor point
M_n <- c(10, 14)

# Simulate the self-correcting process
generated_locs <- simulate_sc(
  t_min = 0,
  t_max = 1,
  sc_params = generating_parameters,
  anchor_point = M_n,
  xy_bounds = c(0, 25, 0, 25)
)
```
