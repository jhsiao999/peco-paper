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


start_ngene <- as.numeric(args[1])
end_ngene <- as.numeric(args[2])

library(Biobase)
library(parallel)
library(genlasso)
library(scran)

dir <- "/project2/gilad/joycehsiao/fucci-seq"

log2cpm.quant <- readRDS(file.path(dir, "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))

eset <- readRDS(file=file.path(dir,"data/eset-final.rds"))
fdata <- fData(eset)
fdata <- fdata[grep("ENSG", rownames(eset)), ]


# # get Cyclone time ------
# get pre-trained marker genes
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

Y_cyclone <- log2cpm.quant
rownames(Y_cyclone) <- rownames(fdata)

# message("run cyclone", "\n")
# out_cyclone <- cyclone(Y_cyclone, pairs = hs.pairs,
#                        gene.names=rownames(Y_cyclone),
#                        iter=1000, min.iter=100, min.pairs=50,
#                        BPPARAM=SerialParam(), verbose=T, subset.row=NULL)
# names(out_cyclone$phases) <- colnames(Y_cyclone)
# #all.equal(names(out_cyclone$phases), rownames(pdata))
# all.equal(names(out_cyclone$phases), colnames(log2cpm.quant))
#
# # Infer phase using Seurat
# pc_cyclone <- prcomp(out_cyclone$normalized.scores, scale. = T, center = T)
#
# library(peco)
# phase_cyclone <- intensity2circle(pc_cyclone$x, plot.it = T, method = "trig")
# phase_cyclone_rotated <- rotation(pdata$theta, phase_cyclone)
#
# names(phase_cyclone_rotated) <- colnames(counts)
#
# out_cyclone$phase_cyclone_rotated <- phase_cyclone_rotated
# # message("save cyclone phase", "\n")
# saveRDS(out_cyclone,
#         file = "analysis/revision/output/ourdata_phase_cyclone.rds")


#-----
out_cyclone <- readRDS(file.path(dir, "analysis/revision/output/ourdata_phase_cyclone.rds"))
phases <- out_cyclone$phase_cyclone_rotated
names(phases) <- names(out_cyclone$phases)
#all.equal(names(phases), colnames(log2cpm.quant))
phases <- phases[order(phases)]
log2cpm.quant_matched <- log2cpm.quant[,match(names(phases),
                                             colnames(log2cpm.quant))]

library(peco)
#source(file.path(dir,"code/utility.R"))
#fit.trend <- mclapply(1:nrow(log2cpm.quant_matched), function(g) {
#fit.trend <- mclapply(1:10, function(g) {
fit.trend <- mclapply(start_ngene:end_ngene, function(g) {
  fit_trendfilter_generic(log2cpm.quant_matched[g,],
                              polyorder = 2)
}, mc.cores=10)
names(fit.trend) <- rownames(log2cpm.quant_matched)[start_ngene:end_ngene]

saveRDS(fit.trend,
        file.path(dir, paste0("analysis/revision/output/ourdata_cyclone/ourdata_cyclone_",
                              start_ngene, "_", end_ngene, ".rds")))


