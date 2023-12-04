# ICIKendallTau 0.3.0

- Handling parallel execution differently to avoid really large matrix issues on each core.
- Introduces the `return_matrix` parameter to `ici_kendalltau` that allows to return the results in the form of a `data.frame` instead of lists of matrices.

# ICIKendallTau 0.2.10

- Provide the `kt_fast` function that handles missing or `NA` values similarly to the `stats::cor` function, but uses our `ici_kt` fast function underneath.

# ICIKendalltau 0.2.8

- Made all error outputs the same length as the default output but containing `NA`.
- Check for the case when one of the variables in `ici_kt` has all identical arguments, warns the user and returns `NA`.

# ICIKendallTau 0.2.1

- Added new argument `include_only` to `ici_kendalltau` that allows the user to define which of the pairwise correlations to actually do.

# ICIKendallTau 0.1.17

- switched theme.
- updated installation instructions to use r-universe.

# ICIKendallTau 0.1.16

- updating documentation and examples

# ICIKendallTau 0.1.4

- Switched `pairwise_completeness` to also use the `global_na` parameter.
- Oh yes, there is a function `pairwise_completeness` to enable scaling by the completeness between any two samples.

# ICIKendallTau 0.1.2

-   Updated the API to use a single variable, `global_na` that defines all the values that should be set to `NA` for the correlation calculation.
-   This is a big API change, so we bumped the version up to 0.1.

# ICIKendallTau 0.0.6

-   Fixing a bug where if there are 55,000 elements in the vector, some of the match overflows from 32 bit to 64 bit, and results make no sense.

# ICIKendallTau 0.0.5

-   Fixed a bug where instead of returning a two element vector, it only returned a zero length value. This probably only happened when one of the entries contained all NA values, or you tried to pass "vectors" with less than two entries.

# ICIKendallTau 0.0.4

-   Fixing things that came up with R CMD check in the documentation. Should be all good now (I hope).

# ICIKendallTau 0.0.3

-   Added tests!
-   Removed a bunch of code that wasn't necessary because it was using the incorrect formula's, and we have the correct version in the `ici_kt_pairs` function for reference purposes.

# ICIKendallTau 0.0.2

-   Fixed an issue where a warning would be issued if {furrr} was not installed.

# ICIKendallTau 0.0.1

-   First release. It's been tested and used rather thoroughly by myself, but I'd say it's still rather alpha. Even though I use it all the time myself.
-   Added a `NEWS.md` file to track changes to the package.
