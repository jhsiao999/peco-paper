#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}

ind <- args[1]
ngenes <- as.numeric(args[2])
ncores <- as.numeric(args[3])
grids <- as.numeric(args[4])


library(Biobase)
library(peco)


# getting Leng data...
dir <- "/project2/gilad/joycehsiao/fucci-seq"


# --- prediction for samples from one individual cell ine
library(peco)
# ind <- "NA19098"
# ngenes <- 5

dir <-"/project2/gilad/joycehsiao/fucci-seq"
#source(file.path(dir,"code/working/run_methods.R"))

data_train <- readRDS(file.path(dir,paste0("analysis/revision/output/data_training_test/ind_",
                                           ind, "_data_training.rds")))
data_test <- readRDS(file.path(dir,paste0("analysis/revision/output/data_training_test/ind_",
                                          ind, "_data_test.rds")))

fits_all <- readRDS(file.path(dir,"output/npreg-trendfilter-quantile.Rmd/fit.quant.rds"))
genes_all <- names(fits_all)[order(sapply(fits_all,"[[",3), decreasing=T)]

# ngenes=1
which_genes <- genes_all[1:ngenes]
# Y_train <- with(data_train, log2cpm_quant_train[which(rownames(log2cpm_quant_train) %in% which_genes), ])
# Y_train <- matrix(Y_train, nrow=1, ncol=length(Y_train))
# colnames(Y_train) <- colnames(data_train$log2cpm_quant_train)
# rownames(Y_train) <- which_genes
fit_train <- cycle_npreg_insample(
  Y = with(data_train, log2cpm_quant_train[which(rownames(log2cpm_quant_train) %in% which_genes), ]),
  theta = with(data_train, theta_train),
  polyorder=2,
  ncores=ncores,
  method.trend="trendfilter")

# str(data_train$theta_train)
# str(data_train$log2cpm.quant_train)

fit_test <- cycle_npreg_outsample(
  Y_test=with(data_test, log2cpm_quant_test[which(rownames(log2cpm_quant_test) %in% which_genes), ]),
  sigma_est=with(fit_train, sigma_est),
  funs_est=with(fit_train, funs_est),
  method.grid = "uniform",
  method.trend="trendfilter",
  polyorder=2,
  ncores=ncores, grids=grids)

out <- list(fit_train=fit_train,
            fit_test=fit_test)

saveRDS(out, file.path(dir,
                       paste0("analysis/revision/output/ourdata_ind_peco/ourdata_peco_",
                              ind, "_top", sprintf("%03d", ngenes),"genes_grids_",grids,".rds")))


#
# a=readRDS(file.path(dir,
#                        paste0("analysis/revision/output/ourdata_ind_peco/ourdata_peco_",
#                               "NA18855", "_top", 500,"genes.rds")))
