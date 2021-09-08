# csdR 0.99.6

* Bugfix: Open MP was not working correctly because of missing compiler flags. For this reason, the `Makevars` file has been created.
* Calculation of column ranks now uses `matrixStats::colRanks` instead of an apply statement with `base::rank`.


# csdR 0.99.0

* Added a `NEWS.md` file to track changes to the package. This is the first public version of the package.
