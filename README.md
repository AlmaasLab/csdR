This repository implements the CSD algorithm presented by [Voigt et al. 2017](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005739) in an efficient manner in R.
* Requrirements: Requirements: R (version 3.5.0 or higher) with packages `WGCNA`, `optparse`, `glue`, `magrittr` and `Rcpp` (and of course a C++ compiler) installed. Additionally, having an optimized Blas library such as openBlas is highly recommended for performance reasons (see [this link](https://www.r-bloggers.com/2010/06/faster-r-through-better-blas/) for more info).
## Usage
Currently, there are two user options. Let us assume that you want to compare the expression profiles in the two files `healthy.txt` and `sick.txt`, where the samples are in columns and genes are in rows. We want to run 2000 iterations, set the random seed to 123, use 5 cores and save the results to a file named `results.txt`. The file `welford.cpp`  must be present for the program to work, but is compiled automatically by Rcpp, and does hence not need to be touched by the user.
### Command line option: 
```
./find_rho_and_var -a healthy.txt -b sick.txt -p 5 -o results.txt -s 123
```
However, this outputs a huge file the disk and still leaves to us to find the top links.
### Script approach
This is more flexible and allows all steps to be taken at once, eliminating the need for writing a large file. See the  [example script](csd_example.R) for this procedure.
