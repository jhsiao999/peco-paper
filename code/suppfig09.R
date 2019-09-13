# Supplemental Figure 9

library(Biobase)
library(circular)
library(peco)
library(cluster)
library(ggplot2)

dir <- "/project2/gilad/joycehsiao/fucci-seq"
eset <- readRDS(file.path(dir, "data/eset-final.rds"))

# Fig 2A. Derive FUCCI-based cell cycle phase
pca <- prcomp(cbind(pData(eset)$rfp.median.log10sum.adjust,
                    pData(eset)$gfp.median.log10sum.adjust))
(pca$sdev^2)/sum(pca$sdev^2)

plot(pca$x[,1], pca$x[,2], pch=16, cex=.5, xlim=c(-4, 4), ylim=c(-4,4),
     xlab="PC1 (67%)", ylab="PC2 (33%)",
     main = "fucci intensities PC1 vs PC2", col="gray50", axes=F)
axis(1);axis(2)
abline(h=0,v=0, col="gray50", lty=2)
par(new=TRUE)
theta <- coord2rad(pca$x)
plot(circular(theta), stack=T, shrink=1.3, cex=.5, bins=200)
theta_final <- shift_origin(as.numeric(theta), 3*pi/4)

# get PAM-based clusters
theta_final <- shift_origin(as.numeric(theta), 3*pi/4)
pam_res <- pam(cbind(pData(eset)$rfp.median.log10sum.adjust,
                     pData(eset)$gfp.median.log10sum.adjust), k=3)
clust <- data.frame(clust=pam_res$clustering,
                    sample_id=rownames(pData(df)))

# get all significant cyclic genes assigend as cell cycle genes
# get significant cyclic genes that are also classifed as cell cycle genes
# permutation-based PVE
perm.lowmiss <- readRDS(file.path(dir,
                                  "output/npreg-trendfilter-quantile.Rmd/fit.trend.perm.lowmiss.rds"))
pve.perm.lowmiss <- sapply(perm.lowmiss, "[[", "trend.pve")
fit.quant <- readRDS(file.path(dir, "output/npreg-trendfilter-quantile.Rmd/fit.quant.rds"))
pval_peco <- data.frame(pve=sapply(fit.quant, "[[", "trend.pve"),
                        pval_perm = sapply(fit.quant, function(x) (1+sum(pve.perm.lowmiss > x$trend.pve))/(1+ length(pve.perm.lowmiss)) ),
                        row.names = fdata$name[match(names(fit.quant),
                                                     rownames(fdata))])
pval_peco_ordered <- pval_peco[order(pval_peco$pve, decreasing = T),]


macosko <- readRDS(file.path(dir,
                             "data/cellcycle-genes-previous-studies/cycle_genes/macosko-2015.rds"))
metrics <- pval_peco_ordered$pval_perm
cutoff <- .001
cyclegenes <- toupper(macosko$hgnc)
allgenes <- toupper(rownames(pval_peco_ordered))

which_sig_cycle <- allgenes[metrics < cutoff][allgenes[metrics < cutoff] %in% cyclegenes]
which_sig_cycle_genes <- data.frame(names=which_sig_cycle,
                                    ensg=rownames(fData(eset))[match(which_sig_cycle,
                                                                     fData(eset)$name)])

all.equal(as.character(clust$sample_id), colnames(data_quant_ord))
#subdata_plot <- do.call(rbind, lapply(1:nrow(which_sig_cycle_genes), function(g) {

# top 30 of 54 significant cyclic genes also known as cell cycle genes
subdata_plot <- do.call(rbind, lapply(1:30, function(g) {
  gindex <- which(rownames(data_quant_ord) == which_sig_cycle_genes$ensg[g])
  gexp <- data_quant_ord[gindex,]
  data.frame(gexp=gexp, gindex=g, gene=which_sig_cycle_genes$names[g],
             cell_state=as.character(clust$clust)[match(colnames(data_quant_ord),
                                                        as.character(clust$sample_id))])
}))
subdata_plot$cell_state <- factor(subdata_plot$cell_state,
                                  labels=c("G1", "S", "G2M"))

subdata_plot$gene <- factor(subdata_plot$gene,
                            levels=which_sig_cycle_genes$names,
                            labels=which_sig_cycle_genes$names)
ggplot(subdata_plot, aes(x=cell_state, y = gexp, fill=cell_state)) +
  geom_boxplot() + facet_wrap(~gene, nrow=6) +
  ylab("Quantile-normalized \n gene expression levels") + xlab("") +
  scale_fill_manual(values=cols,
                    name="Phase") + theme_bw()



# bottom 24 of 54 significant cyclic genes also known as cell cycle genes
subdata_plot <- do.call(rbind, lapply(31:54, function(g) {
  gindex <- which(rownames(data_quant_ord) == which_sig_cycle_genes$ensg[g])
  gexp <- data_quant_ord[gindex,]
  data.frame(gexp=gexp, gindex=g, gene=which_sig_cycle_genes$names[g],
             cell_state=as.character(clust$clust)[match(colnames(data_quant_ord),
                                                        as.character(clust$sample_id))])
}))
subdata_plot$cell_state <- factor(subdata_plot$cell_state,
                                  labels=c("G1", "S", "G2M"))

subdata_plot$gene <- factor(subdata_plot$gene,
                            levels=which_sig_cycle_genes$names,
                            labels=which_sig_cycle_genes$names)
ggplot(subdata_plot, aes(x=cell_state, y = gexp, fill=cell_state)) +
  geom_boxplot() + facet_wrap(~gene, nrow=6) +
  ylab("Quantile-normalized \n gene expression levels") + xlab("") +
  scale_fill_manual(values=cols,
                    name="Phase") + theme_bw()
