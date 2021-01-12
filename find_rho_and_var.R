#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(Rcpp))
sourceCpp("welford.cpp")
# make_option(c("-h", "--help"), action="store", default=FALSE, 
#                      help="Show this help message and exit")
options_list  <- list(
		      make_option(c("-a","--condition1"), action="store",type="character", help= "Expression file in tab deliminated format for first condition"),
		      make_option(c("-b","--condition2"), action="store",type="character",help = "Expression file in tab deliminated format   for second condition"),
		      make_option(c("-p","--threads"), action="store",type="integer", default = 1L, help = "Number of threads for the correlation calculations"),
		      make_option(c("-n","--iterations"),action="store",type="integer",default = 20L, help = "Number of bootstrap iterations to run"),
		      make_option(c("-o","--output"), action="store",type="character", default = "csd_res.txt", help = "Name of output file"),
		      make_option(c("-q","--quiet"),action="store_false",default=TRUE,type="logical",dest="verbose",help = "Run without printing concurrent progress"),
		      make_option(c("-s","--seed"),action="store",type="integer",help="Random seed for reproducable results")
		      )

log_progress  <- function(msg){
	message(glue("{date()} => {msg}"))
}

fast_spearman  <- function(x,nThreads=1L){
	ranks <- apply(X=x,MARGIN=2L,FUN=base::rank)
	WGCNA::cor(ranks,nThreads=nThreads)
}

run_cor_bootstrap  <- function(x,n_it=20L,nThreads=1L,verbose = TRUE){
	n_res  <- ncol(x)
	gene_names <- colnames(x)
	rho_matrix  <- array(0,c(n_res,n_res),
	                     dimnames = list(gene_names,gene_names))
	var_matrix  <- array(0,c(n_res,n_res),
	                     dimnames = list(gene_names,gene_names))
	for(i in seq_len(n_it)){
		if(verbose){
		log_progress(glue("Running bootstrap iteration {i} of {n_it}..."))
		}
		bootstrap_ind  <- sample(seq_len(nrow(x)),replace=TRUE)
		cor_matrix <- fast_spearman(x[bootstrap_ind,],nThreads)
		welford_update(mu = rho_matrix,var = var_matrix,cor_matrix = cor_matrix,
		               iteration = i,n_threads = nThreads)
		
	}
	# The procedure might be numerically unstable, so we make sure the estimated variances are always positive
	var_matrix  <- var_matrix / (n_it - 1)
	list(rho=rho_matrix,var=var_matrix)
}

run_csd  <- function(x_1,x_2,n_it=20L,nThreads=1L,verbose=TRUE){
	if(is.null(colnames(x_1)) || is.null(colnames(x_2))){
		   stop("The input files must be labelled with gene names")
		      }
	genes_to_analyze  <- intersect(colnames(x_1),colnames(x_2))
	shared_x_1  <- intersect(colnames(x_1),genes_to_analyze)
	if(length(shared_x_1) != ncol(x_1)){
		   warning(glue("{ncol(x_1)-length(shared_x_1)} genes for the first conditions were not found in the second"))
	}
	shared_x_2  <- intersect(colnames(x_2),genes_to_analyze)
	if(length(shared_x_2) != ncol(x_2)){
		   warning(glue("{ncol(x_2)-length(shared_x_2)} genes for the second conditions were not found in the first"))
	}
	x_1  <- x_1[,genes_to_analyze]
	x_2  <- x_2[,genes_to_analyze]
	if(verbose){
	log_progress(glue("Running CSD with\n{nrow(x_1)} samples from condition 1\n{nrow(x_2)} samples from condition 2\nNumber of genes: {length(genes_to_analyze)}\nNumber of bootstrap iterations: {n_it}\nNumber of threads: {nThreads}"))
}
	if(verbose){
		log_progress("Running correlation bootstrapping on first condition...")
	}
	res_1  <- run_cor_bootstrap(x_1,n_it=n_it,nThreads=nThreads,verbose = verbose)
	if(verbose){
		log_progress("Running correlation bootstrapping on second condition...")
	}
	res_2  <- run_cor_bootstrap(x_2,n_it=n_it,nThreads=nThreads, verbose = verbose)
	if(verbose){
		log_progress("Summarizing results")
	}
	upper_tri_matrix <- upper.tri(res_1$rho)
	Gene1 <- rownames(res_1$rho)[row(res_1$rho)[upper_tri_matrix]]
	Gene2 <- colnames(res_1$rho)[col(res_1$rho)[upper_tri_matrix]]
	rho1 <- res_1$rho[upper_tri_matrix]
	rho2 <- res_2$rho[upper_tri_matrix]
	var1 <- res_1$var[upper_tri_matrix]
	var2 <- res_2$var[upper_tri_matrix]
	std_estimate  <- sqrt(var1 + var2)
	cVal  <- abs(rho1 + rho2) / std_estimate
	sVal <- abs(abs(rho1) - abs(rho2)) / std_estimate
	dVal <- abs(abs(rho1) + abs(rho2) - abs(rho1 + rho2)) / std_estimate
	csd_df  <- data.frame(Gene1,
			      Gene2,
			      rho1,
			      rho2,
			      var1,
			      var2,
			      cVal,
			      sVal,
			      dVal)
	class(csd_df)  <- c("csd_res", class(csd_df))
	csd_df
}

if(sys.nframe() == 0L){
opt  <- parse_args(OptionParser(option_list=options_list))
if(!is.null(opt$seed)){
	set.seed(opt$seed)
}
if(opt$verbose){
	log_progress("Parsing file from first condition...")
}
x_1  <- read.table(opt$condition1,header = TRUE, sep = "", row.names = 1) %>% as.matrix() %>% t()
if(opt$verbose){
	log_progress("Parsing file from second condition...")
}
x_2  <- read.table(opt$condition2,header = TRUE, sep = "", row.names = 1) %>% as.matrix %>% t()
csd_df  <- run_csd(x_1,x_2,n_it = opt$iterations, nThreads = opt$threads,verbose=opt$verbose)
if(opt$verbose){
	log_progress("Writing results to file...")
}
write.table(x = csd_df, file = opt$output, sep = '\t', row.names = FALSE, quote = FALSE)
if(opt$verbose){
	log_progress("DONE")
}
}

