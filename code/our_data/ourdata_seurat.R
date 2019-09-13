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


library(parallel)
library(genlasso)
library(peco)

dir <- "/project2/gilad/joycehsiao/fucci-seq"

log2cpm.quant <- readRDS(file.path(dir, "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))

# # get Seurat-based phase
# library(Seurat)
# cc.genes <- readLines(con = "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt")
#
# s.genes <- cc.genes[1:43]
# g2m.genes <- cc.genes[44:97]
#
# obj <- CreateSeuratObject(counts = counts)
# obj <- NormalizeData(obj)
# obj <- FindVariableFeatures(obj, selection.method = "vst")
# obj <- ScaleData(obj, features = rownames(obj))
#
# obj <- RunPCA(obj, features = VariableFeatures(obj), ndims.print = 1:10,
#               nfeatures.print = 10)
# #DimHeatmap(obj, dims = c(8, 10))
#
# obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#
# out_seurat <- obj[[]]
# all.equal(rownames(out_seurat), colnames(log2cpm.quant))
# out_seurat <- out_seurat[match(colnames(log2cpm.quant), rownames(out_seurat)),]
# pdata <- pdata[match(rownames(out_seurat), rownames(pdata)),]
#
# # Infer phase using Seurat
# pc_seurat <- prcomp(cbind(out_seurat$S.Score,
#                           out_seurat$G2M.Score), scale. = T, center = T)
#
# library(peco)
# phase_seurat <- intensity2circle(pc_seurat$x, plot.it = F, method = "trig")
# phase_seurat_rotated <- rotation(pdata$theta, phase_seurat)
# names(phase_seurat_rotated) <- colnames(counts)
#
# match(phase_seurat_rotated)
# out_seurat <- cbind(out_seurat,
#                     phase_peco=phase_seurat_rotated[match(rownames(out_seurat),
#                                                names(phase_seurat_rotated))])
# saveRDS(out_seurat,
#         file = "analysis/revision/output/ourdata_phase_seurat.rds")


out_seurat <- readRDS(file.path(dir,
                                "analysis/revision/output/ourdata_phase_seurat.rds"))
phase_seurat_rotated <- out_seurat$phase_peco
names(phase_seurat_rotated) <- rownames(out_seurat)
phase_seurat_rotated <- phase_seurat_rotated[order(phase_seurat_rotated)]
log2cpm.quant_seurat <- log2cpm.quant[,match(names(phase_seurat_rotated),
                                             colnames(log2cpm.quant))]

library(peco)
#source(file.path(dir,"code/utility.R"))
fit.trend <- mclapply(1:nrow(log2cpm.quant_seurat), function(g) {
#fit.trend <- mclapply(1:10, function(g) {
  fit_trendfilter_generic(log2cpm.quant_seurat[g,],
                              polyorder = 2)
}, mc.cores=10)
names(fit.trend) <- rownames(log2cpm.quant_seurat)

saveRDS(fit.trend,
        file.path(dir,"analysis/output/revision/ourdata_seurat.rds"))


