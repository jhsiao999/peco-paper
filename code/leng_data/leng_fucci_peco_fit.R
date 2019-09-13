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


library(Biobase)
library(peco)



# getting Leng data...
dir <- "/project2/gilad/joycehsiao/fucci-seq"

eset <- readRDS(file.path(dir,
            "data/cellcycle-genes-previous-studies/rnaseq/leng/HumanLengESC.rds"))
pdata_leng <- pData(eset)
counts_leng <- exprs(eset)

counts_leng_fucci <- counts_leng[,pdata_leng$cell_state != "H1"]
counts_leng_h1 <- counts_leng[,pdata_leng$cell_state == "H1"]

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



fit_predict <- readRDS(file.path(dir,
                                 "analysis/revision/output/leng_fucci_peco/top_10.rds"))$fit_predict
phase_peco <- fit_predict$cell_times_est

#all.equal(names(phase_peco), colnames(log2cpm_leng_fucci_quant))


# fit trendfilter to each gene
fit.trend <- mclapply(start_ngene:end_ngene, function(g) {
  fit_trendfilter_generic(log2cpm_leng_fucci_quant[g, order(phase_peco)],
                          polyorder = 2)
}, mc.cores=10)
names(fit.trend) <- rownames(log2cpm_leng_fucci_quant)[start_ngene:end_ngene]

saveRDS(fit.trend,
        file.path(dir, paste0("analysis/revision/output/leng_fucci_peco_fit/leng_peco_",
                              start_ngene, "_", end_ngene, ".rds")))






