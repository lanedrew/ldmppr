# Fast spatio-temporal interaction for the self-correcting model

Computes \$g_i = exp(-alpha3 \* sum\_{j\<i} 1\[ \|\|x_i-x_j\|\| \<=
beta3 AND (t_i - t_j) \>= gamma3 \])\$ for i = 1..n, with g_0 = exp(0) =
1.

## Usage

``` r
interaction_st_fast(data, params)
```

## Arguments

- data:

  NumericMatrix with columns (time, x, y). Assumed sorted by time
  ascending.

- params:

  NumericVector length 3: (alpha3, beta3, gamma3)

## Value

NumericVector length n of exp(-alpha3 \* counts)
