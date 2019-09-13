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

# # get Seurat-based phase
# library(Seurat)
# cc.genes <- readLines(con = "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt")
#
# s.genes <- cc.genes[1:43]
# g2m.genes <- cc.genes[44:97]
#
# obj <- CreateSeuratObject(counts = counts_leng_fucci)
# obj <- NormalizeData(obj)
# obj <- FindVariableFeatures(obj, selection.method = "vst")
# obj <- ScaleData(obj, features = rownames(obj))
#
# # obj <- RunPCA(obj, features = VariableFeatures(obj), ndims.print = 1,
# #               nfeatures.print = 1)
# #DimHeatmap(obj, dims = c(8, 10))
# #
# obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#
# out_seurat <- obj[[]]
# all.equal(rownames(out_seurat), colnames(counts_leng_fucci))
# out_seurat <- out_seurat[match(colnames(counts_leng_fucci), rownames(out_seurat)),]
#
# # Infer phase using Seurat
# pc_seurat <- prcomp(cbind(out_seurat$S.Score,
#                           out_seurat$G2M.Score), scale. = T, center = T)
#
#
# phase_seurat <- intensity2circle(pc_seurat$x, plot.it = F, method = "trig")
# #phase_seurat_rotated <- rotation(pdata$theta, phase_seurat)
# names(phase_seurat) <- colnames(counts_leng_fucci)
#
#
# out_seurat <- cbind(out_seurat,
#                     phase_peco=phase_seurat[match(rownames(out_seurat),
#                                                names(phase_seurat))])
# saveRDS(out_seurat,
#         file = "analysis/revision/output/leng_phase_seurat.rds")


out_seurat <- readRDS(file.path(dir,
                        "analysis/revision/output/leng_phase_seurat.rds"))
phase_seurat <- out_seurat$phase_peco


# fitting trends -----------------------
# log2cpm_nest_filt_seurat <- log2cpm_nest_filt[,match(rownames(out_seurat),
#                                              colnames(log2cpm_nest_filt))]

#source(file.path(dir,"code/utility.R"))
#fit.trend <- mclapply(1:nrow(log2cpm_nest_filt_seurat), function(g) {
fit.trend <- mclapply(start_ngene:end_ngene, function(g) {
#  fit.trend <- mclapply(1:10, function(g) {
  fit_trendfilter_generic(log2cpm_leng_fucci_quant[g,order(phase_seurat)],
                          polyorder = 2)
}, mc.cores=10)
names(fit.trend) <- rownames(log2cpm_leng_fucci_quant)[start_ngene:end_ngene]

saveRDS(fit.trend,
        file.path(dir,
        paste0("analysis/revision/output/leng_fucci_seurat/leng_seurat_",
               start_ngene, "_", end_ngene, ".rds")))



