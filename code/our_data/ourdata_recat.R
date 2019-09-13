#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.
#
# To run this script on the command line
# sbatch

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}

seed <- args[1]

library(doParallel)
library(mclust)
library(cluster)
library(TSP)

# get our data -----
dir <- "/project2/gilad/joycehsiao/fucci-seq"
log2cpm.quant <- readRDS(file.path(dir,
                         "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))

# run reCAT on our data ------
dir_recat <- "/project2/gilad/joycehsiao/reCAT/R"
setwd(dir_recat)

input <- t(log2cpm.quant)

set.seed(seed)
input <- input[sample(1:nrow(input)),]

source("get_test_exp.R")
test_exp <- get_test_exp(t(input))

source("get_ordIndex.R")
res_ord <- get_ordIndex(test_exp, 10)

ordIndex <- res_ord$ordIndex

save(test_exp, ordIndex,
     file = file.path(dir, paste0("analysis/revision/output/ourdata_recat_",seed,".rda")))




