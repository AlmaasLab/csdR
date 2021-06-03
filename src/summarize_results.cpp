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
#if defined(_OPENMP)
#pragma omp parallel for simd num_threads(n_threads)
#endif
  for(int i = 0; i < n_genes; i++){
    for(int j = i  + 1; j < n_genes; j++){
      R_xlen_t this_index = ((n_genes - 1) + (n_genes - i)) * i / 2
      + j - (i + 1);
      Gene1_vec[this_index] = gene_names[i];
      Gene2_vec[this_index] = gene_names[j];
      double rho_1 = rho_1_mat(i,j);
      rho1_vec[this_index] = rho_1;
      double rho_2 = rho_2_mat(i,j);
      rho2_vec[this_index] = rho_2;
      double var_1 = var_1_mat(i,j);
      var1_vec[this_index] = var_1;
      double var_2 = var_2_mat(i,j);
      var2_vec[this_index] = var_2;
      double std_estimate = sqrt(var_1 + var_2);
      cVal[this_index] = fabs(rho_1 + rho_2) / std_estimate;
      sVal[this_index] = fabs(fabs(rho_1) - fabs(rho_2)) / std_estimate;
      dVal[this_index] = fabs(fabs(rho_1) + fabs(rho_2) - fabs(rho_1 + rho_2)) / std_estimate;
    }
  }
  DataFrame res = DataFrame::create(Named("Gene1") = Gene1_vec, Named("Gene2") = Gene2_vec,
                                    Named("rho1") = rho1_vec, Named("rho2") = rho2_vec,
                                    Named("var1") = var1_vec, Named("var2") = var2_vec,
                                    Named("cVal") = cVal, Named("sVal") = sVal, Named("dVal") = dVal);
  return res;
}
