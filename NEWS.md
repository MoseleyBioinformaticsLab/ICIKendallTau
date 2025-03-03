# ICIKendallTau 1.2.10

- Added the number of non-missing values used to calculate correlation in data.frame output of `cor_fast`.

# ICIKendallTau 1.2.9

- Fixing an issue where **any** `NA` observations would be removed when using `cor_fast` and `pairwise.complete.obs`, instead of just those for the pairs (closes #25).
- Adding `include_only` function argument to `cor_fast` to allow specification of which correlations to include (closes #24).

# ICIKendallTau 1.2.7

- Fixing a test that fails on r-universe.

# ICIKendallTau 1.2.5

- Fixed a bug where the `check_timing` wasn't taking the right arguments, and wasn't actually calculating what it should.
- Added the `completeness` metric to the output of `ici_kt` (and therefore to `ici_kendalltau`) (closes #23).
- Added more tests (closes #22 for now).

# ICIKendallTau 1.2.2

- Fixed a bug where passing a two vector list is **supposed** to restrict the comparisons to just those provided, but instead did all possible pairwise comparisons between the two things and all others available.

# ICIKendallTau 1.2.0

- Refactored much of `ici_kendalltau`, making the code more consistent and easier to extend, as well as providing more informative error messages. Thanks to @njtierney for suggestions.
- Also refactored `kt_fast` to be more consistent and use more functions internally. **Note**: now returns matrix or data.frame regardless of whether passing simply two vectors or a matrix input.
- Added `alternative` and `continuity` arguments to both `ici_kendalltau` and `kt_fast`.
- Added `cor_fast` to allow running many iterations of `cor.test` on large matrix inputs if desired, with parallel processing to speed things up.

# ICIKendallTau 1.1.3

- makes `rank_order_data` take a sample class argument to enable splitting out by class.

# ICIKendallTau 1.1.2

- adds the `rank_order_data` function to help with visualizing the pattern of missingness with respect to the rank of the rows.
- adds the `yeast_missing` dataset as a more realistic dataset for running the `test_left_censorship` and `rank_order_data`.

# ICIKendallTau 1.1.0

- adds the function `test_left_censorship` to verify if `ici_kendalltau` is appropriate to use on the data or not.

# ICIKendallTau 1.0.0

- Calculates correlation between columns of the matrix, **not** the rows.

# ICIKendallTau 0.3.20

- `kt_fast` now uses the same data.frame format for output as `ici_kendalltau`, but returns a matrix by default. The data.frame is useful when large amounts of comparisons are run.

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
