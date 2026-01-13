# ldmppr ![](reference/figures/ldmppr_logo_hex5.png)

`ldmppr` is an `R` package for working with location dependent marked
point processes. The package includes a suite of tools for model
estimation, model fit assessment, visualization, and simulation for
marked point processes with dependence between the marks and locations
and regularity in the pattern.

### Workflow Overview

1.  Estimate the parameters of a self-correcting point process given a
    reference dataset.
2.  Train a mark model using simulated or real-world data.
3.  Check the fit of the model using various non-parametric summaries
    for point processes and global envelope tests.
4.  Simulate and visualize datasets from the fitted model.

For additional details on implementing the package workflow, run
[`vignette("ldmppr_howto")`](https://lanedrew.github.io/ldmppr/articles/ldmppr_howto.md)
in `R` after installing the package.

### Installation

You can install the development version of `ldmppr` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lanedrew/ldmppr", build_vignettes = TRUE)
```

You can install the stable version of `ldmppr` from CRAN:

``` r
install.packages("ldmppr")
```

For details on how to install the `terra` package that `ldmppr` depends
on, please visit the [`terra` installation
page](https://github.com/rspatial/terra).
