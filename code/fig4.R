# Figure 4

dir <- "/project2/gilad/joycehsiao/fucci-seq"

library(yarrr)
library(ggplot2)

# ------ getting data
library(Biobase)
eset <- readRDS(file.path(dir,
                          "data/cellcycle-genes-previous-studies/rnaseq/leng/HumanLengESC.rds"))
pdata_leng <- pData(eset)
counts_leng <- exprs(eset)
fdata_leng <- fData(eset)

counts_leng_fucci <- counts_leng[,pdata_leng$cell_state != "H1"]
pdata_leng_fucci <- pdata_leng[pdata_leng$cell_state != "H1",]
pdata_leng_fucci$cell_state <- droplevels(pdata_leng_fucci$cell_state)
pdata_leng_fucci$cell_state <- factor(pdata_leng_fucci$cell_state,
                                      levels = c("G1", "S", "G2"),
                                      labels = c("G1", "S", "G2M"))
# counts_f ucci
log2cpm_leng_fucci <- t(log2(10^6 * t(counts_leng_fucci + 1)/colSums(counts_leng_fucci)))
keep_genes <- which(rowMeans(log2cpm_leng_fucci) > 1)
log2cpm_leng_fucci_filt <- log2cpm_leng_fucci[keep_genes,]

# quantile normalize data
log2cpm_leng_fucci_quant <- do.call(rbind, lapply(1:nrow(log2cpm_leng_fucci_filt), function(g) {
  yy <- log2cpm_leng_fucci_filt[g,]
  is.zero <- which(yy == 0)
  qq.map <- qqnorm(yy, plot.it = FALSE)
  yy.qq <- qq.map$x
  yy.qq[is.zero] <- sample(qq.map$x[is.zero])
  return(yy.qq)
}) )
rownames(log2cpm_leng_fucci_quant) <- rownames(log2cpm_leng_fucci_filt)
colnames(log2cpm_leng_fucci_quant) <- colnames(log2cpm_leng_fucci_filt)


# Use top 4 genes for HIST1H4E was deprecated in genome build hg38
labs <- c("CDK1", "UBE2C", "TOP2A", "HIST1H4C")
genes <- c("ENSG00000170312","ENSG00000175063",
           "ENSG00000131747", "ENSG00000197061")


# Figure 4A
out_peco <- readRDS(file=file.path(dir,
                                   paste0("analysis/revision/output/leng_fucci_peco/top_5.rds")))
subdata_plot <- do.call(rbind, lapply(1:4, function(g) {
  gindex <- which(rownames(out_peco$fit_predict$Y) == labs[g])
  gexp <- out_peco$fit_predict$Y[gindex,]
  data.frame(gexp=gexp, gene=labs[g], cell_state=pdata_leng_fucci$cell_state)
}))
ggplot(subdata_plot, aes(x=cell_state, y = gexp, fill=cell_state)) +
  geom_boxplot() + facet_wrap(~gene, ncol=4) +
  ylab("Quantile-normalized \n gene expression levels") + xlab("") +
  scale_fill_manual(values=as.character(yarrr::piratepal("espresso")[3:1]),
                    name="Phase") + theme_bw()


# ------ peco
# code/leng_data/leng_fucci_peco.R predicts cell phase using top X genes from the training data
# code/leng_data/lleng_fucc_peco_fit.R estimate cyclic trend for all genes in the Leng fucci data

out_peco <- readRDS(file=file.path(dir,
                                   paste0("analysis/revision/output/leng_fucci_peco/top_5.rds")))

fits_peco <- lapply(1:4, function(g) {
  gindex <- which(rownames(out_peco$fit_predict$Y) == labs[g])
  fit_g <- data.frame(
    gexp=out_peco$fit_predict$Y[gindex,],
    phase=shift_origin(out_peco$fit_predict$cell_times_est, origin = pi/2))
  fit_g <- fit_g[order(fit_g$phase),]

  fit_trend <- fit_trendfilter_generic(fit_g$gexp, polyorder = 2)
  fit_g$trend.yy <- fit_trend$trend.yy
  fun_g <- approxfun(x=as.numeric(fit_g$phase),
                     y=as.numeric(fit_g$trend.yy), rule=2)
  fit_out <- list(fit_g=fit_g,
                  fun_g = fun_g)
  return(fit_out)
})
names(fits_peco) <- labs

all.equal(rownames(fits_peco[[1]]$fit_g), rownames(pdata_leng_fucci))
cell_state_col <- data.frame(cell_state = pdata_leng_fucci$cell_state[match(rownames(fits_peco[[1]]$fit_g),
                                                                            rownames(pdata_leng_fucci))],
                             sample_id = rownames(fits_peco[[1]]$fit_g))
cell_state_col$cols <- as.character(yarrr::piratepal("espresso")[3:1])[cell_state_col$cell_state]

par(mfrow=c(1,4), mar=c(2,2,2,1))
for (g in 1:4) {
  if (g>1) { ylab <- ""} else { ylab <- "Quantile-normalized expression value"}
  #    if (i==1) { xlab <- "Predicted phase"} else { xlab <- "Inferred phase"}
  res_g <- fits_peco[[g]]
  plot(x= res_g$fit_g$phase,
       y = res_g$fit_g$gexp, pch = 16, cex=.5, ylab=ylab,
       xlab="Predicted phase", col = cell_state_col$cols,
       main = names(fits_peco)[g], axes=F, ylim=c(-3,3))
  abline(h=0, lwd=.5)
  lines(x = seq(0, 2*pi, length.out=200),
        y = res_g$fun_g(seq(0, 2*pi, length.out=200)),
        col=wesanderson::wes_palette("FantasticFox1")[1], lty=1, lwd=2)
  axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
                labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                         expression(2*pi)))
}



# Oscope
library(openxlsx)
oscope29genes <- unlist(read.xlsx(file.path(dir,
                        "data/cellcycle-genes-previous-studies/nmeth.3549-S2.xlsx"), colNames = F))

# run Oscope on FUCCI cells ------------------------------
# load pre-computed results
load(file = file.path(dir, "analysis/revision/output/leng_fucci_oscope_29genes.rda"))
# library(Oscope)
# Sizes <- MedianNorm(counts_leng_fucci)
# DataNorm <- GetNormalizedMat(counts_leng_fucci, Sizes)
# DataNorm_sub <- DataNorm[which(rownames(DataNorm) %in% oscope29genes),]
# DataInput <- NormForSine(DataNorm_sub)
# SineRes <- OscopeSine(DataInput, parallel=F)
# KMResUse <- list(cluster1=oscope29genes)
# ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput,
#                     NCThre = 100, parallel=T)
# save(DataInput,
#      KMResUse, ENIRes, pdata_leng_fucci,
#      file = file.path(dir, "analysis/revision/output/leng_fucci_oscope_29genes.rda"))

phase_pred <- seq(0, 2*pi, length.out= length(ENIRes$cluster1))
names(phase_pred) <- colnames(DataInput)[ENIRes[["cluster1"]]]
fits_leng_oscope <- lapply(1:4, function(g) {
  gindex <- which(rownames(log2cpm_leng_fucci_quant) == labs[g])
  fit_g <- data.frame(
    gexp=log2cpm_leng_fucci_quant[gindex, match(names(phase_pred), colnames(log2cpm_leng_fucci_quant))],
    phase=shift_origin((-phase_pred+2*pi), origin = pi/2))
  fit_g <- fit_g[order(fit_g$phase),]

  fit_trend <- fit_trendfilter_generic(fit_g$gexp)
  fit_g$trend.yy <- fit_trend$trend.yy
  fun_g <- approxfun(x=as.numeric(fit_g$phase),
                     y=as.numeric(fit_g$trend.yy), rule=2)
  fit_out <- list(fit_g=fit_g,
                  pve = fit_trend$pve,
                  fun_g = fun_g)
  return(fit_out)
})
names(fits_leng_oscope) <- labs

all.equal(rownames(fits_leng_oscope[[1]]$fit_g), rownames(pdata_leng_fucci))
cell_state_oscope_col <- data.frame(cell_state = pdata_leng_fucci$cell_state[match(rownames(fits_leng_oscope[[1]]$fit_g),
                                                                                   rownames(pdata_leng_fucci))],
                                    sample_id = rownames(fits_leng_oscope[[1]]$fit_g) )
cell_state_oscope_col$cols <- as.character(yarrr::piratepal("espresso")[3:1])[cell_state_oscope_col$cell_state]

par(mfrow=c(1,4), mar=c(2,2,2,1))
for (g in 1:4) {
  if (g>1) { ylab <- ""} else { ylab <- "Normalized log2CPM"}
  #    if (i==1) { xlab <- "Predicted phase"} else { xlab <- "Inferred phase"}
  res_g <- fits_leng_oscope[[g]]
  plot(x= res_g$fit_g$phase,
       y = res_g$fit_g$gexp, pch = 16, cex=.5, ylab=ylab,
       xlab="Predicted phase", col=cell_state_oscope_col$cols,
       main = names(fits_leng_oscope)[g], axes=F, ylim=c(-3,3))
  abline(h=0, lwd=.5)
  lines(x = seq(0, 2*pi, length.out=200),
        y = res_g$fun_g(seq(0, 2*pi, length.out=200)),
        col = wesanderson::wes_palette("Darjeeling2")[2], lty=1, lwd=2)
  axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
                labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                         expression(2*pi)))
}







# --- recat
# load pre-computed results
# see leng_data/leng_fucci_recat.R
out_recat <- load(file.path(dir,"analysis/revision/output/leng_fucci_recat.rda"))
sample_ordered <- rownames(test_exp)[ordIndex]
phase_pred <- seq(0, 2*pi, length.out= length(ordIndex))
names(phase_pred) <- sample_ordered

fits_recat <- lapply(1:4, function(g) {
  gindex <- which(rownames(log2cpm_leng_fucci_quant) == labs[g])
  fit_g <- data.frame(
    gexp=log2cpm_leng_fucci_quant[gindex, match(sample_ordered, colnames(log2cpm_leng_fucci_quant))],
    phase=(-phase_pred)+2*pi )
  fit_g <- fit_g[order(fit_g$phase),]

  fit_trend <- fit_trendfilter_generic(fit_g$gexp, polyorder = 2)
  fit_g$trend.yy <- fit_trend$trend.yy
  #    fit_g$pve <- fit$pve
  fun_g <- approxfun(x=as.numeric(fit_g$phase),
                     y=as.numeric(fit_g$trend.yy), rule=2)
  fit_out <- list(fit_g=fit_g,
                  #                  pve = fit$pve,
                  fun_g = fun_g)
  return(fit_out)
})
names(fits_recat) <- labs

all.equal(rownames(fits_recat[[1]]$fit_g), rownames(pdata_leng_fucci))
cell_state_recat_col <- data.frame(cell_state = pdata_leng_fucci$cell_state[match(rownames(fits_recat[[1]]$fit_g),
                                                                                  rownames(pdata_leng_fucci))],
                                   sample_id = rownames(fits_recat[[1]]$fit_g) )
cell_state_recat_col$cols <- as.character(yarrr::piratepal("espresso")[3:1])[cell_state_recat_col$cell_state]

par(mfrow=c(1,4), mar=c(2,2,2,1))
for (g in 1:4) {
  if (g>1) { ylab <- ""} else { ylab <- "Quantile-normalized expression value"}
  res_g <- fits_recat[[g]]
  plot(x= res_g$fit_g$phase,
       y = res_g$fit_g$gexp, pch = 16, cex=.5, ylab=ylab,
       xlab="Predicted phase", col = cell_state_recat_col$cols,
       main = names(fits_recat)[g], axes=F, ylim=c(-3,3))
  abline(h=0, lwd=.5)
  lines(x = seq(0, 2*pi, length.out=200),
        y = res_g$fun_g(seq(0, 2*pi, length.out=200)),
        col=wesanderson::wes_palette("FantasticFox1")[5], lty=1, lwd=2)
  axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
                labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                         expression(2*pi)))
}




# ------ seurat
library(Seurat)
cc.genes <- readLines(con = file.path(dir,
                                      "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt"))

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

obj <- CreateSeuratObject(counts = counts_leng_fucci)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst")
obj <- ScaleData(obj, features = rownames(obj))

obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
out_seurat <- obj[[]]
all.equal(rownames(out_seurat), colnames(counts_leng_fucci))
out_seurat <- out_seurat[match(colnames(counts_leng_fucci), rownames(out_seurat)),]

out_seurat <- readRDS(file = file.path(dir, "analysis/revision/output/leng_fucci_seurat.rds"))
out_seurat$gates <- factor(out_seurat$gates,
                           levels = c("G1", "S", "G2"),
                           labels = c("G1", "S", "G2M"))
out_seurat$Phase <- factor(out_seurat$Phase,
                           levels = c("G1", "S", "G2M"))
tab_seurat <- table(out_seurat$gates, out_seurat$Phase)
sum(tab_seurat[1,c(2,3)])/rowSums(tab_seurat)[1]
sum(tab_seurat[2,c(1,3)])/rowSums(tab_seurat)[2]
sum(tab_seurat[3,c(1,2)])/rowSums(tab_seurat)[3]



#--- cyclone
library(scran)

# get pre-trained marker genes
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

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

# geneinfo <- read.table(write.table(geneinfo,
#       file = file.path(dir,
#          "data/cellcycle-genes-previous-studies/rnaseq/leng/leng_geneinfo.txt"))

geneinfo <- read.table(file.path(dir,
                                 "data/cellcycle-genes-previous-studies/rnaseq/leng/leng_geneinfo.txt"))

indata <- which(geneinfo$ensg %in% unlist(hs.pairs))
geneinfo_indata <- geneinfo[indata,]

Y_cyclone <- log2cpm_leng_fucci_quant[rownames(log2cpm_leng_fucci_quant) %in% as.character(geneinfo_indata$symbol),]
rownames(Y_cyclone) <- geneinfo_indata$ensg[match(rownames(Y_cyclone),
                                                  geneinfo_indata$symbol)]

out_cyclone <- cyclone(Y_cyclone, pairs = hs.pairs,
                       gene.names=rownames(Y_cyclone),
                       iter=1000, min.iter=100, min.pairs=50,
                       BPPARAM=SerialParam(), verbose=T, subset.row=NULL)
names(out_cyclone$phases) <- colnames(Y_cyclone)
all.equal(names(out_cyclone$phases), colnames(log2cpm_leng_fucci_quant))
all.equal(names(out_cyclone$phases), rownames(pdata_leng_fucci))

pdata_leng_fucci$cell_state <- factor(pdata_leng_fucci$cell_state,
                                      levels = c("G1", "S", "G2"),
                                      labels = c("G1", "S", "G2M"))
out_cyclone$phases <- factor(out_cyclone$phases,
                             levels = c("G1", "S", "G2M"))
tab_cyclone <- table(pdata_leng_fucci$cell_state, out_cyclone$phases)
sum(tab_cyclone[1,c(2,3)])/rowSums(tab_cyclone)[1]
sum(tab_cyclone[2,c(1,3)])/rowSums(tab_cyclone)[2]
sum(tab_cyclone[3,c(1,2)])/rowSums(tab_cyclone)[3]




