# csdR

This R-package implements the CSD algorithm presented by [Voigt et al. 2017](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005739) in an efficient manner.
* Requrirements: Requirements: R (version 4.1.0 or higher) with packages `WGCNA`, `optparse`, `glue`, `magrittr` and `Rcpp` (and of course a C++ compiler) installed. Additionally, having an optimized Blas library such as openBlas is highly recommended for performance reasons (see [this link](https://www.r-bloggers.com/2010/06/faster-r-through-better-blas/) for more info).

## Installation

In order to install the release version from Bioconductor, do type in the R terminal:

```
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("csdR")
```

If you want to install the development version instead:

```
# install.packages("devtools")
devtools::install_github("AlmaasLab/csdR", ref = "main")
```

## Usage
Please see the [package vignette](https://almaaslab.github.io/csdR/articles/csdR.html). Additionally check out the article on `csdR`:

Pettersen, J.P., Almaas, E. csdR, an R package for differential co-expression analysis. BMC Bioinformatics 23, 79 (2022). https://doi.org/10.1186/s12859-022-04605-1


## Issues and feedback
Please use the repository's [issue tracker](https://github.com/AlmaasLab/csdR/issues) if you cannot make the package work or if you have suggestions for improvements.