
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

- How to handle missing data (i.e. `NA`’s) in calculating a correlation
  between two variables.
- Current calculations of correlation are based on having all pairs of
  observations for two variables.
  - However, whether an observation is present or missing is
    semi-quantitative information for many analytical measurements with
    sensitivity limits.
  - i.e. in many cases, missing observations are not
    “missing-at-random”, but “missing-not-at-random” due to falling
    below the detection limit.
  - In these cases, NA is informative.
  - Therefore, in **most** analytical measurements (gene expression,
    proteomics, metabolomics), missing measurements should be included,
    and contribute to the correlation.

If you want to read more on **how** we solve this problem, see the
package vignette.

## Package Functions

The functions that implement this include:

- `ici_kt`: the C++ workhorse, actually calculating a correlation
  between an X and Y.
  - The option `perspective` will control how the `NA` values influence
    ties.
  - When comparing samples, you likely want to use
    `perspective = "global"`.
- `ici_kendallt`: Handles comparisons for a large matrix.
  - Rows are features, columns are samples.
  - Implicitly parallel, but have to call:
    - `library(furrr)`
    - `plan(multiprocess)`
  - Otherwise will only use a single core.

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

r_1 = ici_kendalltau(matrix_1)
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
r_2 = ici_kendalltau(matrix_2)
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
#>                           expr       min        lq       mean    median        uq       max
#>  cor(x, y, method = "kendall") 11506.697 11670.094 12169.6628 12006.418 12482.883 13182.222
#>         ici_kt(x, y, "global")   243.866   250.125   294.6542   275.104   320.058   384.118
#>       ici_kt(x2, y2, "global") 13467.011 13739.312 14658.5050 14945.446 14987.140 16153.616
#>  neval
#>      5
#>      5
#>      5
```

In the case of 40,000 features, the average time on a modern CPU is 14
milliseconds.

Of course, if you want to use it to calculate Kendall-tau-b without
incorporating missingness, it can do that just fine as well.

``` r
k_tau = ici_kt(x, y, "global")
all.equal(k_tau[[1]] ,cor(x, y, method = "kendall"))
#> [1] TRUE
```

We also provide the `kt_fast` function, if you want something that
treats `NA` values similarly to `stats::cor`.

``` r
k_tau_fast = kt_fast(x, y)
k_tau_fast
#>          tau       pvalue 
#> -0.003411411  0.871672260
```

## Parallelism

If you have {future} and the {furrr} packages installed, then it is also
possible to split up the a set of matrix comparisons across compute
resources for any multiprocessing engine registered with {future}.

``` r
library(furrr)
future::plan(multicore, workers = 4)
r_3 = ici_kendalltau(matrix_2)
```

## Many Many Comparisons

In the case of hundreds of thousands of comparisons to be done, the
result matrices can become very, very large, and require lots of memory
for storage. They are also inefficient, as both the lower and upper
triangular components are stored. An alternative storage format is as a
`data.frame`, where there is a single row for each comparison performed.
This is actually how the results are stored internally, and then they
are converted to a matrix form if requested (the default).s To keep the
`data.frame` output, add the argument `return_matrix=FALSE` to the call
of `ici_kendalltau`.

``` r
r_4 = ici_kendalltau(matrix_2, return_matrix = FALSE)
r_4
#> $cor
#>   s1 s2 core       raw pvalue   taumax       cor
#> 1 s3 s4    1 0.9924359      0 0.997963 0.9944616
#> 2 s3 s3    0 1.0000000      0 1.000000 1.0000000
#> 3 s4 s4    0 1.0000000      0 1.000000 1.0000000
#> 
#> $run_time
#> [1] 0.01606894
```

## Code of Conduct

Please note that the ICIKendallTau project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
