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


# get Cyclone time ------
# library(scran)

# # get pre-trained marker genes
# hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
#
# # get ENSG ID
# library(biomaRt)
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# geneinfo <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
#                filters = "hgnc_symbol",
#                values = rownames(log2cpm_leng_fucci_quant), bmHeader = T, mart = mart)
# colnames(geneinfo) <- c("symbol", "ensg")
# write.table(geneinfo,
#       file = file.path(dir,
#          "data/cellcycle-genes-previous-studies/rnaseq/leng/leng_geneinfo.txt"))

# geneinfo <- read.table(file.path(dir,
#           "data/cellcycle-genes-previous-studies/rnaseq/nestorowa/nestorawa_geneinfo.txt"))
#
# indata <- which(geneinfo$ensg %in% unlist(hs.pairs))
# geneinfo_indata <- geneinfo[indata,]
#
# Y_cyclone <- log2cpm_leng_fucci_quant[rownames(log2cpm_leng_fucci_quant) %in% as.character(geneinfo_indata$symbol),]
# rownames(Y_cyclone) <- geneinfo_indata$ensg[match(rownames(Y_cyclone),
#                                                   geneinfo_indata$symbol)]
#
# message("run cyclone", "\n")
# out_cyclone <- cyclone(Y_cyclone, pairs = hs.pairs,
#                        gene.names=rownames(Y_cyclone),
#                        iter=1000, min.iter=100, min.pairs=50,
#                        BPPARAM=SerialParam(), verbose=T, subset.row=NULL)
# names(out_cyclone$phases) <- colnames(Y_cyclone)
# all.equal(names(out_cyclone$phases), colnames(log2cpm_leng_fucci_quant))
#
# # Infer phase using Seurat
# pc_cyclone <- prcomp(out_cyclone$normalized.scores, scale. = T, center = T)
#
# library(peco)
# phase_cyclone <- intensity2circle(pc_cyclone$x, plot.it = F, method = "trig")
# names(phase_cyclone) <- colnames(log2cpm_leng_fucci_quant)
#
# message("save cyclone phase", "\n")
# saveRDS(phase_cyclone,
#         file = file.path(dir, "analysis/revision/output/leng_fucci_phase_cyclone.rds"))

phase_cyclone <- readRDS(file = file.path(dir,
                                  "analysis/revision/output/leng_fucci_phase_cyclone.rds"))
#all.equal(names(phase_cyclone), colnames(log2cpm_nest_quant))

library(peco)
#source(file.path(dir,"code/utility.R"))
#fit.trend <- mclapply(1:nrow(log2cpm.quant_matched), function(g) {
#fit.trend <- mclapply(1:10, function(g) {
fit.trend <- mclapply(start_ngene:end_ngene, function(g) {
  fit_trendfilter_generic(log2cpm_leng_fucci_quant[g,order(phase_cyclone)],
                          polyorder = 2)
}, mc.cores=10)
names(fit.trend) <- rownames(log2cpm_leng_fucci_quant)[start_ngene:end_ngene]

saveRDS(fit.trend,
        file.path(dir, paste0("analysis/revision/output/leng_fucci_cyclone/leng_fucci_cyclone_",
                              start_ngene, "_", end_ngene, ".rds")))


