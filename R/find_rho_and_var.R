log_progress <- function(msg) {
    message(glue("{date()} => {msg}"))
}

fast_spearman <- function(x, nThreads = 1L) {
    ranks <- colRanks(x,ties.method = "average",preserveShape = TRUE)
    cor(ranks, nThreads = nThreads)
}

#' @title Run bootstrapping of Spearman correlations within a dataset
#' @description This function provides the more low-level functionality
#' of bootstrapping the Spearman correlations of the columns within a dataset.
#' Only use this function if you want
#' a low-level interface, else \code{\link{run_csd}}
#' provides a more streamlined approach if you want to do a CSD analysis.
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom glue glue
#' @importFrom WGCNA cor
#' @importFrom matrixStats colRanks
#' @param x Numeric matrix, the gene expression matrix to analyse.
#' Genes are in columns, samples are in rows.
#' @inheritParams run_csd
#' @return A list with two fields \describe{
#' \item{rho}{Numeric matrix constaining the bootstrapped
#'  mean of the Spearman correlation between each column}
#' \item{var}{Numeric matrix constaining the bootstrapped
#'  variance of the Spearman correlation between each column}
#' }
#' @examples
#' data("normal_expression")
#' cor_res <- run_cor_bootstrap(
#'     x = normal_expression,
#'     n_it = 100, nThreads = 2L
#' )
#' @export
run_cor_bootstrap <- function(x, n_it = 20L, nThreads = 1L,
                            verbose = TRUE, iterations_gap = 1L) {
    blas_set_num_threads(nThreads)
    n_res <- ncol(x)
    gene_names <- colnames(x)
    rho_matrix <- array(0, c(n_res, n_res),
        dimnames = list(gene_names, gene_names)
    )
    var_matrix <- array(0, c(n_res, n_res),
        dimnames = list(gene_names, gene_names)
    )
    for (i in seq_len(n_it)) {
        if (verbose && i %% iterations_gap == 0) {
            log_progress(glue("Running bootstrap iteration {i} of {n_it}..."))
        }
        bootstrap_ind <- sample(seq_len(nrow(x)), replace = TRUE)
        cor_matrix <- fast_spearman(x[bootstrap_ind, ], nThreads)
        welford_update(
            mu = rho_matrix, var = var_matrix,
            cor_matrix = cor_matrix,
            iteration = i, n_threads = nThreads
        )
    }
    var_matrix <- var_matrix / (n_it - 1)
    list(rho = rho_matrix, var = var_matrix)
}

validate_csd_input <- function(x_1, x_2) {
    if (anyNA(x_1)) {
        stop("Argument x_1 has missing values")
    }
    if (anyNA(x_2)) {
        stop("Argument x_2 has missing values")
    }
    if (is.null(colnames(x_1)) || is.null(colnames(x_2))) {
        stop("The input matrices must be labelled with gene names")
    }
    genes_to_analyze <- intersect(colnames(x_1), colnames(x_2))
    shared_x_1 <- intersect(colnames(x_1), genes_to_analyze)
    if (length(shared_x_1) != ncol(x_1)) {
        warning(glue("{ncol(x_1)-length(shared_x_1)} genes
                    for the first condition were not found in the second"))
    }
    shared_x_2 <- intersect(colnames(x_2), genes_to_analyze)
    if (length(shared_x_2) != ncol(x_2)) {
        warning(glue("{ncol(x_2)-length(shared_x_2)} genes
                    for the second condition were not found in the first"))
    }
    list(x_1 = x_1[, genes_to_analyze], x_2 = x_2[, genes_to_analyze],
        genes_to_analyze = genes_to_analyze)
}

#' @title Run CSD analysis
#' @description This function implements the CSD algorithm based on the
#' one presented by Voigt et al. 2017. All pairs of genes are
#' first compared within each condition by the Spearman correlation
#' and the correlation and its variance are
#' estimated by bootstrapping.
#' Finally, the results for the two conditions are
#' compared and C-, S- and D-values are computed
#' and returned.
#' @param x_1 Numeric matrix, the gene expression matrix for
#' the first condition.
#' Genes are in columns, samples are in rows.
#' The columns must be named with the name of the genes.
#' Missing values are not allowed.
#' @param x_2 Numeric matrix, the gene expression matrix for
#' the second condition.
#' @param n_it Integer, number of bootstrap iterations
#' @param nThreads Integer, number of threads to use for computations
#' @param verbose Logical, should progress be printed?
#' @param iterations_gap If output is verbose - Number of iterations between
#' each status message
#' (Default=1 - Displayed only if \code{verbose=TRUE})
#' @return A \code{data.frame} with
#' the additional class attribute \code{csd_res} with the
#'  results of the CSD analysis.
#'  This frame has a row for each pair of genes and has the
#'  following columns: \describe{
#'  \item{\code{Gene1}}{Character, the name of the first gene}
#'  \item{\code{Gene2}}{Character, the name of the second gene}
#'  \item{\code{rho1}}{Mean correlation of the two genes in the first condition}
#'  \item{\code{rho2}}{Mean correlation of the two genes in the second condition}
#'  \item{\code{var1}}{The estimated variance of \code{rho1}
#'   determined by bootstrapping}
#'  \item{\code{var2}}{The estimated variance of \code{rho2}
#'   determined by bootstrapping}
#'  \item{\code{cVal}}{Numeric, the conserved score.
#'   A high value indicates that the co-expression of the two genes
#'   have the same sign in both conditions}
#'   \item{\code{sVal}}{Numeric, the specific score.
#'   A high value indicates that the co-expression of the two genes
#'   have a high degree of co-expression in
#'   one condition, but not the other.}
#'   \item{\code{dVal}}{Numeric, the differentiated score.
#'   A high value indicates that the co-expression of the two genes have
#'   a high degree of co-expression in
#'   both condition, but the sign of co-expression is different.}
#' }
#' @details The gene names in \code{x_1} and \code{x_2}
#' do not need to be in the same order,
#' but must be in the same namespace.
#' Only genes present in both datasets will be considered for the analysis.
#' The parallelism gained by `nThreads` applies to the computations
#' within a single iteration. The iterations are run is serial in order
#' to reduce the memory footprint.
#' @references Voigt A, Nowick K and Almaas E
#' 'A composite network of conserved and tissue specific
#' gene interactions reveals possible genetic interactions in glioma'
#' In: \emph{PLOS Computational Biology} 13(9): e1005739.
#' (doi: \url{https://doi.org/10.1371/journal.pcbi.1005739})
#' @examples
#' data("sick_expression")
#' data("normal_expression")
#' cor_res <- run_csd(
#'     x_1 = sick_expression, x_2 = normal_expression,
#'     n_it = 100, nThreads = 2L
#' )
#' c_max <- max(cor_res$cVal)
#' @export
run_csd <- function(x_1, x_2, n_it = 20L, nThreads = 1L,
                    verbose = TRUE, iterations_gap = 1L) {
    validated_matrices <- validate_csd_input(x_1, x_2)
    x_1 <- validated_matrices$x_1
    x_2 <- validated_matrices$x_2
    genes_to_analyze <- validated_matrices$genes_to_analyze
    if (verbose) {
        log_progress(glue("Running CSD with
\t                  {nrow(x_1)} samples from condition 1
\t                  {nrow(x_2)} samples from condition 2
\t                  Number of genes: {length(genes_to_analyze)}
\t                  Number of bootstrap iterations: {n_it}
\t                  Number of threads: {nThreads}"))
    }
    if (verbose) {
        log_progress("Running correlation bootstrapping on first condition...")
    }
    res_1 <- run_cor_bootstrap(x_1,
        n_it = n_it,
        nThreads = nThreads, verbose = verbose,
        iterations_gap = iterations_gap
    )
    if (verbose) {
        log_progress("Running correlation bootstrapping on second condition...")
    }
    res_2 <- run_cor_bootstrap(x_2,
        n_it = n_it,
        nThreads = nThreads, verbose = verbose,
        iterations_gap = iterations_gap
    )
    if (verbose) {
        log_progress("Summarizing results")
    }
    csd_df <- summarizeResults(
        res_1 = res_1, res_2 = res_2,
        n_threads = nThreads
    )
    class(csd_df) <- c("csd_res", class(csd_df))
    if (verbose) {
        log_progress("Returning from CSD procedure...")
    }
    csd_df
}
