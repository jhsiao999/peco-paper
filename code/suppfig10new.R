# Supplemental Figure 10

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
                    sample_id=rownames(pData(eset)))


# check all cell cycle phase markers in Whitfield et al.
data_quant <- readRDS(file.path(dir, "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))
sample_ord <- rownames(pData(eset))[order(theta_final)]
data_quant_ord <- data_quant[,match(sample_ord,colnames(data_quant))]


macosko <- readRDS(file.path(dir,
                             "data/cellcycle-genes-previous-studies/cycle_genes/macosko-2015.rds"))

subdata_plot = do.call(rbind, lapply(1:nrow(macosko), function(g) {
  gindex = which(rownames(data_quant_ord) == macosko$ensembl[g])
  if (length(gindex) != 0 ) {
    gexp = data_quant_ord[gindex,]
    tmp = data.frame(gexp=gexp,
                     ensembl=macosko$ensembl[g], hgnc = macosko$hgnc[g],
                     cell_state=as.character(clust$clust)[match(colnames(data_quant_ord),
                                                                as.character(clust$sample_id))])
    tmp2 = data.frame(max_state = which.max(sapply(1:3, function(i) mean(tmp$gexp[tmp$cell_state==i]))),
                      ensembl=macosko$ensembl[g], hgnc = macosko$hgnc[g])
    return(tmp2)
  }
}))
subdata_plot$max_state = factor(subdata_plot$max_state,
                                labels=c("G1", "S", "G2M"))
macosko_sub = macosko[match(subdata_plot$ensembl, macosko$ensembl),]
table(subdata_plot$max_state, macosko_sub$phase)

# top100 = rownames(data_quant_ord)[order(rowVars(data_quant_ord),
#                                         decreasing = TRUE)][1:500]

# table(subdata_plot[as.character(subdata_plot$ensembl) %in% top100,]$max_state,
#       macosko_sub[as.character(macosko_sub$ensembl) %in% top100,]$phase)

table(subdata_plot$max_state,
      macosko_sub$phase)

# compute mis-classfication rate
# G1
(39+13+27)/(53+39+13+27+22)
# S
(14+24+25)/(8+14+14+24+25)
# G2M
(26+51)/(26+51+99+91+53)

