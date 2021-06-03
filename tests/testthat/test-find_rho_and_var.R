set.seed(432)
n_it <- 10000
N_SAMPLES <- 100
gene_1 <- rnorm(N_SAMPLES,mean = 0,sd = 1)
gene_2 <- 1 / sqrt(2)*gene_1 + 1 / sqrt(2)*rnorm(N_SAMPLES,mean = 0,sd = 1)
combined_matrix <- cbind(gene_1,gene_2)
gene_names <- c("gene_1","gene_2")
gene_names_dimnames <- list(gene_names,gene_names)
expected_res <- replicate(n_it, {
  bootstrap_sample <- sample(1:N_SAMPLES,size = N_SAMPLES,replace = TRUE)
  stats::cor(x=gene_1[bootstrap_sample], y =gene_2[bootstrap_sample], method = "spearman")
}
)
expected_rho <- mean(expected_res)
expected_var <- var(expected_res)
expected_result <- list(rho = matrix(c(1,expected_rho,expected_rho,1),byrow = FALSE,ncol = 2,dimnames = gene_names_dimnames),
                        var = matrix(c(0, expected_var,expected_var,0),byrow = FALSE,ncol = 2,dimnames = gene_names_dimnames))
test_that("run_cor_bootstrap works as intended",
          {
            cor_bootstrap_res <- run_cor_bootstrap(combined_matrix,n_it = 1000,nThreads = 2L, verbose = FALSE)
            expect_equal(object = cor_bootstrap_res, expected = expected_result,tolerance = 0.05)
          }
          )

gene_1_sick <- rnorm(N_SAMPLES,mean = 0,sd = 2)
gene_2_sick <- 1 / sqrt(2)*gene_1 - 1 / sqrt(2)*rnorm(N_SAMPLES,mean = 0,sd = 2)
combined_matrix_sick <- cbind(gene_1_sick,gene_2_sick)
colnames(combined_matrix_sick) <- gene_names
expected_res_sick <- replicate(n_it, {
  bootstrap_sample <- sample(1:N_SAMPLES,size = N_SAMPLES,replace = TRUE)
  stats::cor(x=gene_1_sick[bootstrap_sample], y = gene_2_sick[bootstrap_sample], method = "spearman")
}
)
expected_rho_sick <- mean(expected_res_sick)
expected_var_sick <- var(expected_res_sick)
expected_sd_estimate <- sqrt(expected_var+expected_var_sick)
expected_c <- abs(expected_rho + expected_rho_sick) / expected_sd_estimate
expected_s <- abs(abs(expected_rho) - abs(expected_rho_sick)) / expected_sd_estimate
expected_d <- abs(abs(expected_rho) + abs(expected_rho_sick) -
                       abs(expected_rho + expected_rho_sick)) / expected_sd_estimate
expected_csd_res <- structure(list(Gene1 = "gene_1",Gene2 = "gene_2",rho1 = expected_rho,
                                   rho2 = expected_rho_sick,
                                   var1 = expected_var, var2 = expected_var_sick,
                                   cVal = expected_c, sVal = expected_s, dVal = expected_d),
                              class = c("csd_res","data.frame"), row.names = 1L
                              )

test_that("overall CSD algorithm works with a small dataset",{
  expect_equal(object = run_csd(x_1 = combined_matrix, x_2 = combined_matrix_sick, n_it = n_it, nThreads = 2, verbose = FALSE),
              expected = expected_csd_res,tolerance = 0.05
               )
}
)

# col_subset <- 1:10
test_that("overall CSD algorithm gives a reasonable result for more realistic data",
          {
            csd_res <- run_csd(x_1 = csdR::normal_expression, x_2 = csdR::sick_expression,n_it = 20,
                               nThreads = 2,verbose = FALSE
                                )
            expect_equal(object = nrow(csd_res), as.integer(1000*999 / 2))
          }
)
