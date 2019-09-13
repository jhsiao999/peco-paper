#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


start_ngene <- as.numeric(args[1])
end_ngene <- as.numeric(args[2])


library(peco)
library(Biobase)


# getting Leng data...
dir <- "/project2/gilad/joycehsiao/fucci-seq"

eset <- readRDS(file.path(dir,
                          "data/cellcycle-genes-previous-studies/rnaseq/leng/HumanLengESC.rds"))
pdata_leng <- pData(eset)
counts_leng <- exprs(eset)

counts_leng_fucci <- counts_leng[,pdata_leng$cell_state != "H1"]

pdata_leng_fucci <- pdata_leng[pdata_leng$cell_state != "H1",]
pdata_leng_fucci$cell_state <- droplevels(pdata_leng_fucci$cell_state)

# counts_f ucci
log2cpm_leng_fucci <- t(log2(10^6 * t(counts_leng_fucci + 1)/colSums(counts_leng_fucci)))
keep_genes <- which(rowMeans(log2cpm_leng_fucci) > 1)
log2cpm_leng_fucci_filt <- log2cpm_leng_fucci[keep_genes,]

# quantile normalize data
log2cpm_leng_fucci_quant <- do.call(rbind, lapply(1:nrow(log2cpm_leng_fucci_filt), function(g) {
  yy <- log2cpm_leng_fucci_filt[g,]
  is.zero <- which(yy == 0)
  qq.map <- qqnorm(yy, plot.it = F)
  yy.qq <- qq.map$x
  yy.qq[is.zero] <- sample(qq.map$x[is.zero])
  return(yy.qq)
}) )
rownames(log2cpm_leng_fucci_quant) <- rownames(log2cpm_leng_fucci_filt)
colnames(log2cpm_leng_fucci_quant) <- colnames(log2cpm_leng_fucci_filt)



# # run reCAT on our data ------
# dir_recat <- "/project2/gilad/joycehsiao/reCAT/R"
# setwd(dir_recat)
#
# input_data <- t(log2cpm_leng_fucci_quant)
#
# source("get_test_exp.R")
# # gene by cell
# test_exp <- get_test_exp(t(input_data))
#
# source("get_ordIndex.R")
# res_ord <- get_ordIndex(test_exp, 10)
#
# ordIndex <- res_ord$ordIndex
# ensembleResultLst <- res_ord$ensembleResultLst
#
# source("get_score.R")
# score_result <- get_score(t(test_exp))
# score_result$bayes_score
# score_result$mean_score
#
#
# source("plot.R")
# # plot_bayes(score_result$bayes_score, ordIndex)
# myord <- c(50:247, 1:50)
# # plot_bayes(score_result$bayes_score, ordIndex[myord])
#
# message("get_hmm.R", "\n")
# source("get_hmm.R")
# rdata = as.matrix(data.frame(start = c(40, 130, 220),
#                              end = c(50, 140, 230)))
#
# message("get_hmm_order", "\n")
# hmm_result <- get_hmm_order(bayes_score = score_result$bayes_score,
#                             mean_score = score_result$mean_score,
#                             ordIndex = ordIndex, cls_num = 3, myord = myord, rdata = rdata)
#
# # message("get clusters", "\n")
# # source("get_cluster_result.R")
# # cls_result = get_cluster_result(test_exp = test_exp,
# #                                 ensembleResultLst = ensembleResultLst,
# #                                 resultLst = resultLst, cls_num = 20)
# save(test_exp, score_result, ordIndex, ensembleResultLst,
#      myord, rdata, hmm_result, #cls_result,
#      file = file.path(dir, "analysis/revision/output/leng_fucci_recat.rda"))

load(file = file.path(dir, "analysis/revision/output/leng_fucci_recat.rda"))
# fitting trends -----------------------
# log2cpm_nest_filt_seurat <- log2cpm_nest_filt[,match(rownames(out_seurat),
#                                              colnames(log2cpm_nest_filt))]

#source(file.path(dir,"code/utility.R"))
#fit.trend <- mclapply(1:nrow(log2cpm_nest_filt_seurat), function(g) {
fit.trend <- mclapply(start_ngene:end_ngene, function(g) {
#  fit.trend <- mclapply(1:10, function(g) {
  fit_trendfilter_generic(log2cpm_leng_fucci_quant[g, ordIndex],
                          polyorder = 2)
}, mc.cores=10)
names(fit.trend) <- rownames(log2cpm_leng_fucci_quant)[start_ngene:end_ngene]

saveRDS(fit.trend,
        file.path(dir,
        paste0("analysis/revision/output/leng_fucci_recat/leng_fucci_recat_",
               start_ngene, "_", end_ngene, ".rds")))



