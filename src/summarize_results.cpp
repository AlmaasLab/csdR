#define RCPP_NO_BOUNDS_CHECK
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
DataFrame summarizeResults(List res_1, List res_2, int n_threads) {
  NumericMatrix rho_1_mat = res_1["rho"];
  NumericMatrix rho_2_mat = res_2["rho"];
  NumericMatrix var_1_mat = res_1["var"];
  NumericMatrix var_2_mat = res_2["var"];
  int n_genes = rho_1_mat.ncol();
  R_xlen_t n_pairs = n_genes*(n_genes - 1) / 2;
  CharacterVector gene_names = colnames(rho_1_mat);
  CharacterVector Gene1_vec (n_pairs);
  CharacterVector Gene2_vec (n_pairs);
  NumericVector rho1_vec (n_pairs);
  NumericVector rho2_vec (n_pairs);
  NumericVector var1_vec (n_pairs);
  NumericVector var2_vec (n_pairs);
  NumericVector cVal (n_pairs);
  NumericVector sVal (n_pairs);
  NumericVector dVal (n_pairs);
/*
  What is this evil trick?
  Originally it said
  for(int i = 0; i < n_genes; i++){
    for(int j = i  + 1; j < n_genes; j++){
      ...
  }
    }
  However, we parallelize the loop by fusing it:
  i, j >= 0
  j < i
  k = (i-1)*i/2 + j = (i² - i) / 2 +j
  Then, we prove i = floor(sqrt(2*(k+1))+0.5)
  i-0.5 < sqrt(2*(k+1)) < i+0.5
  i² - i + 0.25 < 2k + 2 < i² + i + 0.25 
  i² - i + 0.25 < (i² - i) + 2j + 2 < i² + i + 0.25
  0.25 < 2j + 2 < 2i + 0.25
  which is true because 2i >= 2(j+1) = 2j + 2 such that
  2i + 0.25 > 2j + 2

  */
#if defined(_OPENMP)
  #pragma omp parallel for simd num_threads(n_threads)
#endif
  for(R_xlen_t k = 0; k < n_pairs; k++){
      int i = floor(sqrt(2*(k+1))+.5);
      int j = k - (i-1)*i/2;
      Gene1_vec[k] = gene_names[j];
      Gene2_vec[k] = gene_names[i];
      double rho_1 = rho_1_mat(i,j);
      rho1_vec[k] = rho_1;
      double rho_2 = rho_2_mat(i,j);
      rho2_vec[k] = rho_2;
      double var_1 = var_1_mat(i,j);
      var1_vec[k] = var_1;
      double var_2 = var_2_mat(i,j);
      var2_vec[k] = var_2;
      double std_estimate = sqrt(var_1 + var_2);
      cVal[k] = fabs(rho_1 + rho_2) / std_estimate;
      sVal[k] = fabs(fabs(rho_1) - fabs(rho_2)) / std_estimate;
      dVal[k] = fabs(fabs(rho_1) + fabs(rho_2) - fabs(rho_1 + rho_2)) / std_estimate;
      }
  DataFrame res = DataFrame::create(Named("Gene1") = Gene1_vec, Named("Gene2") = Gene2_vec,
                                    Named("rho1") = rho1_vec, Named("rho2") = rho2_vec,
                                    Named("var1") = var1_vec, Named("var2") = var2_vec,
                                    Named("cVal") = cVal, Named("sVal") = sVal, Named("dVal") = dVal);
  return res;
}
