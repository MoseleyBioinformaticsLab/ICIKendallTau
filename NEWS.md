# ICIKendallTau 0.0.6

* Fixing a bug where if there are 55,000 elements in the vector, some of the match overflows from 32 bit to 64 bit, and results make no sense.

# ICIKendallTau 0.0.5

* Fixed a bug where instead of returning a two element vector, it only returned a zero length value.
This probably only happened when one of the entries contained all NA values, or you tried to pass "vectors" with less than two entries.

# ICIKendallTau 0.0.4

* Fixing things that came up with R CMD check in the documentation.
Should be all good now (I hope).

# ICIKendallTau 0.0.3

* Added tests!
* Removed a bunch of code that wasn't necessary because it was using the incorrect formula's, and we have the correct version in the `ici_kt_pairs` function for reference purposes.

# ICIKendallTau 0.0.2

* Fixed an issue where a warning would be issued if {furrr} was not installed.

# ICIKendallTau 0.0.1

* First release. It's been tested and used rather thoroughly by myself, but I'd say it's still rather alpha. Even though I use it all the time myself.
* Added a `NEWS.md` file to track changes to the package.
