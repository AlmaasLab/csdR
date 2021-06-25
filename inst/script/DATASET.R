## code to prepare dataset goes here
library(readr)
library(magrittr)
set.seed(35397290)
normal_data <- read.table("normal_an10_full_correct.txt", header = TRUE, row.names = 1) %>%
    as.matrix() %>%
    t()
sick_data <- read.table("sick_an10_full.txt", header = TRUE, row.names = 1) %>%
    as.matrix() %>%
    t()
shared_gene_names <- intersect(colnames(normal_data), colnames(sick_data))
genes_to_pick <- sample(shared_gene_names, size = 1000, replace = FALSE)
normal_expression <- normal_data[, genes_to_pick]
sick_expression <- sick_data[, genes_to_pick]
usethis::use_data(normal_expression, sick_expression, overwrite = TRUE)
