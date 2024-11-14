
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sesp <a href="https://stscl.github.io/sesp/"><img src="man/figures/logo.png" align="right" height="139" alt="sesp website" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-cyan.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN](https://www.r-pkg.org/badges/version/sesp)](https://CRAN.R-project.org/package=sesp)
[![R-universe](https://stscl.r-universe.dev/badges/sesp?color=cyan)](https://stscl.r-universe.dev/sesp)
<!-- badges: end -->

**Spatially Explicit Stratified Power Model**

## Installation

<!-- - Install from [CRAN](https://CRAN.R-project.org/package=sesp) with: -->
<!-- ``` r -->
<!-- install.packages("sesp") -->
<!-- ``` -->

- Install development binary version from
  [R-universe](https://stscl.r-universe.dev/sesp) with:

``` r
install.packages('sesp',
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dep = TRUE)
```

- Install development source version from
  [GitHub](https://github.com/stscl/sesp) with:

``` r
if (!requireNamespace("devtools")) {
    install.packages("devtools")
}
devtools::install_github("stscl/sesp",
                         build_vignettes = TRUE,
                         dep = TRUE)
```
