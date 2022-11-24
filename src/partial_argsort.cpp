#include <Rcpp.h>
using namespace Rcpp;

//' @title Extract indecies corresponding to the largest elements
//' @description Extracts the indecies of the \eqn{n} largest elements of the input.
//' This procedure is equivalent to \code{order(x, decreasing = TRUE)[1:n_elements]},
//' but is much faster and avoids the overhead of sorting discarded elements.
//' This function is useful for extracting the rows in a data frame having the
//' largest values in one of the columns.
//' @param x Numeric vector, the vector containing the numbers to sort.
//' @param n_elements Integer scalar, the number of indecies to return.
//' @return Numeric vector, the indecies of the largest elements (in sorted order) in
//' \code{x}.
//' @examples
//' x <- c(10L,5L,-2L,12L,15L)
//' max_indecies <- partial_argsort(x,3L)
//' max_indecies
//' x[max_indecies]
//' order(x)[1:3]
//' mtcars[partial_argsort(mtcars$hp,5L),]
//' @export 
// [[Rcpp::export]]
NumericVector partial_argsort(NumericVector x, int n_elements){
  int sort_upper_bound;
  // Handles the special case when the requested number of elements is higher than
  // the actual size of the array
  if (n_elements > x.length()){
    sort_upper_bound = x.length();
  } else {
    sort_upper_bound = n_elements;
  }
  std::vector < R_xlen_t > idx (x.length());
  std::iota(idx.begin(),idx.end(),0);
  std::nth_element(idx.begin(),idx.begin() + sort_upper_bound, idx.end(),
                   [&x](R_xlen_t i1, R_xlen_t i2) {return x[i1] > x[i2];});
  std::sort(idx.begin(), idx.begin() + sort_upper_bound,
            [&x](R_xlen_t i1, R_xlen_t i2) {return x[i1] > x[i2];});
  NumericVector res (n_elements);
  for(R_xlen_t i = 0; i < n_elements; i++){
    res[i] = idx[i] + 1;
    }
  // If the requested number of elements is higher than the actual size of the array,
  // pad the reminder of the answer with NAs
  for(R_xlen_t i = sort_upper_bound; i < n_elements; i++){
    res[i] = NA_REAL;
  }
  return res;
}
