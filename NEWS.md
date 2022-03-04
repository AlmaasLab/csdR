# csdR 1.1.3

Added citation to the `csdR` article which is now printed.

# csdR 1.1.2

Made it explicit in the documentation that missing values are not allows and
wrote a test for this case.

# csdR 1.1.1

Made some minor modifications to the `README` and the vignette.

# csdR 0.99.6

* Bugfix: Open MP was not working correctly because of missing compiler flags. For this reason, the `Makevars` file has been created.
* Calculation of column ranks now uses `matrixStats::colRanks` instead of an apply statement with `base::rank`.


# csdR 0.99.0

* Added a `NEWS.md` file to track changes to the package. This is the first public version of the package.
