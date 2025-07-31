
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ICIKendallTau

<!-- badges: start -->

[![ICIKendallTau status
badge](https://moseleybioinformaticslab.r-universe.dev/badges/ICIKendallTau)](https://moseleybioinformaticslab.r-universe.dev)
<!-- badges: end -->

You can see the pkgdown site
[here](https://moseleybioinformaticslab.github.io/ICIKendallTau/).

## Citation

This package has an [associated preprint on
bioRxiv](https://doi.org/10.1101/2022.02.24.481854), please cite it if
you use the package in your own work.

> Flight RM, Bhatt PS, Moseley HNB (2025). “Information-Content-Informed Kendall-tau Correlation Methodology: Interpreting Missing Values as Useful Information.” _bioRxiv_. doi:10.1101/2022.02.24.481854 

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
package
[vignette](https://moseleybioinformaticslab.github.io/ICIKendallTau/articles/ici-kendalltau.html).

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

We’ve also included a function for testing if the missingness in your
data comes from left-censorship, `test_left_censorship`. We walk through
creating example data and testing it in the vignette [Testing for Left
Censorship](https://moseleybioinformaticslab.github.io/ICIKendallTau/articles/testing-for-left-censorship).
In addition to testing, you can also visualize the missing data pattern
by feature rank using the `rank_order_data` function, and use
`visdat::vis_miss()` on the original and reordered missing data.

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

microbenchmark(
  cor(x, y, method = "kendall"),
  ici_kt(x, y, "global"),
  ici_kt(x2, y2, "global"),
  times = 5
)
#> Unit: microseconds
#>                           expr       min        lq      mean    median        uq       max neval
#>  cor(x, y, method = "kendall") 12020.945 12227.957 12778.679 12241.762 13459.533 13943.197     5
#>         ici_kt(x, y, "global")   331.001   346.866   456.299   361.061   401.216   841.351     5
#>       ici_kt(x2, y2, "global") 18462.550 18468.493 20342.581 19215.913 22400.049 23165.899     5
```

In the case of 40,000 features, the average time on a modern CPU is 14
milliseconds.

Of course, if you want to use it to calculate Kendall-tau-b without
incorporating missingness, it can do that just fine as well.

``` r
k_tau = ici_kt(x, y, "global")
all.equal(k_tau[[1]], cor(x, y, method = "kendall"))
#> [1] TRUE
```

We also provide the `kt_fast` function, if you want something that
treats `NA` values similarly to `stats::cor`.

``` r
k_tau_fast = kt_fast(x, y)
k_tau_fast
#> $tau
#>              x            y
#> x  1.000000000 -0.003411411
#> y -0.003411411  1.000000000
#> 
#> $pvalue
#>           x         y
#> x 0.0000000 0.8716723
#> y 0.8716723 0.0000000
#> 
#> $run_time
#> [1] 0.02020574
```

## P-Values

ICI-Kt functions only calculates the tau-b variant that handles ties.
P-value calculations use the asymptotic approximation in all cases, and
thus may vary slightly from the p-values returned by R’s `cor.test` and
Python’s `scipy.stats.kendalltau` depending on the number of values in
*x* and *y*.

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
are converted to a matrix form if requested (the default). To keep the
`data.frame` output, add the argument `return_matrix=FALSE` to the call
of `ici_kendalltau`.

``` r
r_4 = ici_kendalltau(matrix_2, return_matrix = FALSE)
r_4
#> $cor
#>   s1 s2 core       raw pvalue   taumax completeness       cor
#> 1 s3 s4    1 0.9924359      0 0.997963        0.921 0.9944616
#> 2 s3 s3    0 1.0000000      0 1.000000        0.950 1.0000000
#> 3 s4 s4    0 1.0000000      0 1.000000        0.950 1.0000000
#> 
#> $run_time
#> [1] 0.01733327
```

## Other Correlations

`ici_kendalltau` and `ici_kt` calculate the p-value of the correlation
as part of the overall calculation. `stats::cor` does not, and
`stats::cor.test` can only calculate the p-value for a single comparison
of two vectors. It is sometimes advantageous to obtain p-values for a
large number of correlations. We provide `cor_fast`, which works
analogously to `kt_fast`, with the ability to choose `pearson` or
`spearman` as the method. Note that if a matrix is provided, the columns
must be named.

``` r
r_5 = cor_fast(x, y, method = "pearson")
r_5
#> $rho
#>            x          y
#> x 1.00000000 0.00720612
#> y 0.00720612 1.00000000
#> 
#> $pvalue
#>           x         y
#> x 0.0000000 0.8199608
#> y 0.8199608 0.0000000
#> 
#> $run_time
#> [1] 0.02877402
```

``` r
m_3 = cbind(x, y, x)
colnames(m_3) = c("s1", "s2", "s3")
r_6 = cor_fast(m_3)
r_6
#> $rho
#>            s1         s2         s3
#> s1 1.00000000 0.00720612 1.00000000
#> s2 0.00720612 1.00000000 0.00720612
#> s3 1.00000000 0.00720612 1.00000000
#> 
#> $pvalue
#>           s1        s2        s3
#> s1 0.0000000 0.8199608 0.0000000
#> s2 0.8199608 0.0000000 0.8199608
#> s3 0.0000000 0.8199608 0.0000000
#> 
#> $run_time
#> [1] 0.02604246
```

## Code of Conduct

Please note that the ICIKendallTau project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
