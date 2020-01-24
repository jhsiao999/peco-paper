# Figure 2

library(SingleCellExperiment)
library(circular)
library(peco)
library(cluster)
library(ggplot2)

dir <- "/project2/gilad/joycehsiao/fucci-seq"
sce <- readRDS(file.path(dir, "data/sce-final.rds"))

# Fig 2A. Derive FUCCI-based cell cycle phase
pca <- prcomp(cbind(colData(sce)$rfp.median.log10sum.adjust,
                    colData(sce)$gfp.median.log10sum.adjust))
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

# Fig 2B. FUCCI phase summarizes fluorescent intensities

# get PAM-based clusters
theta_final <- shift_origin(as.numeric(theta), 3*pi/4)
pam_res <- pam(cbind(colData(sce)$rfp.median.log10sum.adjust,
                     colData(sce)$gfp.median.log10sum.adjust), k=3)
clust <- data.frame(clust=pam_res$clustering,
                    sample_id=rownames(colData(sce)))
plot(theta_final, clust$clust)
b1 <- mean(max(range(theta_final[clust$clust==2])), min(range(theta_final[clust$clust==3])))
b2 <- mean(max(range(theta_final[clust$clust==3])), max(range(theta_final[clust$clust==1])))
b3 <- mean(min(range(theta_final[clust$clust==2])), min(range(theta_final[clust$clust==1])))
abline(v=c(b1,b2,b3), lty=2)
# clust 1 = G1
# clust 2 = S
# clust 3 = G2/M

par(mfrow=c(1,1))
plot(x=shift_origin(as.numeric(theta), 3*pi/4),
     y=colData(sce)$gfp.median.log10sum.adjust, col="forestgreen",
     pch=c(16,1,4)[clust$clust],
     ylim=c(-1.5, 1.5), cex=.5,
     xlab="FUCCI phase", ylab="FUCCI scores",
     main="FUCCI scores", axes=F)
axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
              labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                       expression(2*pi)))
abline(h=0, lty=1, col="black", lwd=.7)
points(x=shift_origin(as.numeric(theta), 3*pi/4),
       y=colData(sce)$rfp.median.log10sum.adjust, col="firebrick",
       ylim=c(-1.5, 1.5),
       pch=c(16,1,4)[clust$clust],
       cex=.5)
abline(v=c(b1,b2,b3), lty=2)


# Fig 2C. FUCCI phase summarizes fluorescent intensities
# get top 5 cyclic genes
fits_all <- readRDS("data/fit.quant.rds")
pve_all <- sapply(fits_all, "[[", 3)
pve_all_ord <- pve_all[order(pve_all, decreasing = T)]
genes <- names(pve_all_ord)[1:5]
labs <- c("CDK1", "UBE2C", "TOP2A", "HIST1H4E", "HIST1H4C")

data_quant <- readRDS("data/log2cpm.quant.rds")
sample_ord <- rownames(colData(df))[order(theta_final)]
data_quant_ord <- data_quant[,match(sample_ord,colnames(data_quant))]
fits_tmp <- lapply(1:5, function(i) {
  ii <- which(rownames(data_quant_ord)==genes[i])
  fit_g <- fit_trendfilter_generic(data_quant_ord[ii,])
  fun_g <- approxfun(x=as.numeric(theta_final[order(theta_final)]),
                     y=as.numeric(fit_g$trend.yy), rule=2)
  return(fun_g)
})
names(fits_tmp) <- genes

par(mfrow=c(1,5), mar=c(2,2,2,1))
for (i in 1:5) {
  plot(x=theta_final[order(theta_final)],
       y=data_quant_ord[rownames(data_quant_ord)==genes[i],], col="gray50",
       xlab="FUCCI phase",
       ylab="", axes=F, cex=.3, pch=16,
       main = labs[i])
  lines(x = seq(0, 2*pi, length.out=200),
        y = fits_tmp[[which(names(fits_tmp) == genes[i])]](seq(0, 2*pi, length.out=200)),
        col=wesanderson::wes_palette("FantasticFox1")[1], lty=1, lwd=2)
  axis(2);
  axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                               expression(3*pi/2), expression(2*pi)))
  abline(h=0, col="black", lty=1, lwd=.7)
  abline(v=c(b1,b2,b3), lty=2)
}