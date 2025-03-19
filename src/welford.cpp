#define RCPP_NO_BOUNDS_CHECK
#include <Rcpp.h>
using namespace Rcpp;

/**
 * Note: This function updates means and variances IN PLACE as opposed to
 * the normal behaviour of R functions
 */
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
SEXP welford_update(NumericVector mu, NumericVector var,NumericVector cor_matrix,
                    int iteration, int n_threads) {
  R_xlen_t n_elements = mu.length();
#if defined(_OPENMP)
  #pragma omp parallel for simd num_threads(n_threads)
#endif
  for(R_xlen_t i = 0; i < n_elements; i++){
      double cor = cor_matrix[i];
      double old_mean = mu[i];
      double old_var = var[i];
      double new_mean = old_mean + (cor - old_mean) / iteration;
      double new_var = old_var + (cor - old_mean) * (cor - new_mean);
      mu[i] = new_mean;
      var[i] = new_var;
  }
  return R_NilValue;
}
