
#'
#' Sample expression matrices for CSD
#'
#' Sample expression matrices of thyroid gland tissue for
#' thyroid cancer patients and healthy individuals. These datasets
#' were pre-processed by Gulla et al. (2019).
#' Due to size requirements, only 1000 randomly selected genes
#' are provided in the dataset. Number of samples are 399 and 504
#' in the healthy controls and the sick samples, respectively.
#'
#'
#' @format Numeric matrices of normalized gene expression.
#' Genes are in columns, whereas samples are in rows.
#' @docType data
#' @aliases normal_expression sick_expression
#' @name data-expression
#' @usage 
#' data(normal_expression)
#' data(sick_expression)
#' @source For the expression matrix for healthy individuals,
#' GenotypeTissue Expression (GTEx) V7.
#' For the thyroid cancer patients,
#' the data are obtained for the Thyroid Cancer
#' project (THCA) from The Cancer Genome Atlas (TCGA).
#' @references 
#' Gulla, Almaas, Eivind, & Voigt, André (2019). 
#' An integrated systems biology approach to investigate
#' transcriptomic data of thyroid carcinoma.
#' NTNU.
#' \url{http://hdl.handle.net/11250/2621725}
#' @keywords datasets
NULL
