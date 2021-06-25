# Instructions for downloading and pre-processing the datasets belonging to the package `csdR`
## Background
The purpose of providing the two datasets `sick_expression` and
`normal_expression` is to have some toy dataset for examples and tinkering.
The two datasets are gene expression datasets from patients with thyroid cancer
and healthy controls, respectively.
They were originally pre-processed and analyzed as a part of
Marie Gulla's [master thesis](http://hdl.handle.net/11250/2621725).
However, due to the sheer size of the data and the computing power
it takes to run examples on them,
in this package we only choose a subset of 1000 randomly selected genes.
Due to their arbitrary and limited covering of genes, it is
**not** recommended to use the datasets for real differential gene expression analysis. 

## Download
The TMP values for the healthy controls can be downloaded from GTEx V7 at (https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz).

For the thyroid cancer cases, the expression files (FPKM) in HTSeq-format can be found at the GDC cancer portal (https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-THCA%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22HTSeq%20-%20FPKM%22%5D%7D%7D%5D%7D&searchTableTab=files). However, the files referring to healthy thyroid tissue were omitted.


## Pre-processing
Please consult Chapter 3 of the [master thesis](http://hdl.handle.net/11250/2621725). The datasets were produced by Analysis 2 with soft filtering.

## Trimming the size of the datasets
In order to prepare the datasets for this package,
trimming was done choosing 1000 genes at random. The script for doing this is
shown in the file `DATASET.R`.
