
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ICIKendallTau

<!-- badges: start -->

[![ICIKendallTau status
badge](https://moseleybioinformaticslab.r-universe.dev/badges/ICIKendallTau)](https://moseleybioinformaticslab.r-universe.dev)
<!-- badges: end -->

You can see the pkgdown site
[here](https://moseleybioinformaticslab.github.io/ICIKendallTau).

## Installation

You can install the current version of ICIKendallTau via GitHub:

``` r
remotes::install_github("MoseleyBioinformaticsLab/ICIKendallTau")
```

You can also install Windows or Mac binaries using our r-universe:

``` r
options(repos = c(
    moseleybioinformaticslab = 'https://moseleybioinformaticslab.r-universe.dev',
    CRAN = "https://cloud.r-project.org"))
install.packages("ICIKendallTau")
```

## Problem

-   How to handle missing data (i.e. `NA`’s) in calculating a
    correlation between two variables.
-   Current calculations of correlation are based on having all pairs of
    observations for two variables.
    -   However, whether an observation is present or missing is
        semi-quantitative information for many analytical measurements
        with sensitivity limits.
    -   i.e. in many cases, missing observations are not
        “missing-at-random”, but “missing-not-at-random” due to falling
        below the detection limit.
    -   In these cases, NA is informative.
    -   Therefore, in **most** analytical measurements (gene expression,
        proteomics, metabolomics), missing measurements should be
        included, and contribute to the correlation.

If you want to read more on **how** we solve this problem, see the
package vignette.

## Package Functions

The functions that implement this include:

-   `ici_kt`: the C++ workhorse, actually calculating a correlation
    between an X and Y.
    -   The option `perspective` will control how the `NA` values
        influence ties.
    -   When comparing samples, you likely want to use
        `perspective = "global"`.
-   `ici_kendallt`: Handles comparisons for a large matrix.
    -   Rows are samples, columns are features.
    -   Implicitly parallel, but have to call:
        -   library(furrr)
        -   plan(multiprocess)
    -   Otherwise will only use a single core.

## Examples

The most common case is a large matrix of independent samples (columns)
and measured features in each of the samples (i.e. gene expression).

Here we will make some artificial data to show how the correlation
changes as we introduce missingness.

``` r
set.seed(1234)
library(ICIKendallTau)

s1 = sort(rnorm(1000, mean = 100, sd = 10))
s2 = s1 + 10 

matrix_1 = cbind(s1, s2)

r_1 = ici_kendalltau(t(matrix_1))
r_1$cor
#>    s1 s2
#> s1  1  1
#> s2  1  1
```

Now we introduce some missing values at the low end of each one. We will
just do the simplest thing and introduce `NA` values in the bottom set.

``` r
s3 = s1
s3[sample(100, 50)] = NA

s4 = s2
s4[sample(100, 50)] = NA

matrix_2 = cbind(s3, s4)
r_2 = ici_kendalltau(t(matrix_2))
r_2$cor
#>           s3        s4
#> s3 1.0000000 0.9944616
#> s4 0.9944616 1.0000000
```

## Is It Fast?

The C++ code implementation (thanks {Rcpp}!) is based on the SciPy
implementation, which uses two merge sorts of the ranks of each vector,
and then looks for differences between them. This is the fastest method
we know of, and has a complexity of O(nlogn). The naive way of computing
it, which explicitly examines all of the pairs, has a complexity of n^2.
Our implementation was compared to the {pcaPP::cov.fk} function, and the
use of {Rcpp} and our inefficient copying of vectors makes ours 3X
slower than that one. Which honestly isn’t too bad.

``` r
library(microbenchmark)
x = rnorm(1000)
y = rnorm(1000)

x2 = rnorm(40000)
y2 = rnorm(40000)

library(microbenchmark)

microbenchmark(
  cor(x, y, method = "kendall"),
  ici_kt(x, y, "global"),
  ici_kt(x2, y2, "global"),
  times = 5
)
#> Unit: microseconds
#>                           expr       min        lq       mean
#>  cor(x, y, method = "kendall") 19619.802 20277.542 20434.6350
#>         ici_kt(x, y, "global")   310.463   346.523   405.0156
#>       ici_kt(x2, y2, "global") 19563.587 19949.583 20738.8768
#>    median        uq       max neval cld
#>  20297.02 20778.512 21200.297     5   b
#>    350.09   371.053   646.949     5  a 
#>  20195.84 20957.648 23027.727     5   b
```

In the case of 40,000 features, the average time on a modern CPU is 13
milliseconds.

Of course, if you want to use it to calculate Kendall-tau-b without
incorporating missingness, it can do that just fine as well.

``` r
k_tau = ici_kt(x, y, "global")
all.equal(k_tau[[1]] ,cor(x, y, method = "kendall"))
#> [1] TRUE
```

## Code of Conduct

Please note that the ICIKendallTau project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
