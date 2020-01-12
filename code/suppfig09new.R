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
                    sample_id=rownames(pData(eset)))
plot(theta_final, clust$clust)

ggplot(data.frame(gfp=pData(eset)$gfp.median.log10sum.adjust,
                  rfp=pData(eset)$rfp.median.log10sum.adjust,
                  clust=factor(clust$clust)), aes(x=gfp, y=rfp, shape=clust)) +
  geom_point() + theme_classic() +
  xlab("EGFP score") +
  ylab("mCherry score") +
  scale_shape_manual(values=c(16,1,4),
                     labels = c("G1", "S", "G2M"),
                     name = "") +
  ggtitle("")


par(mfrow = c(1,3))
for (i in 1:3) {
  hist(pData(eset)$dapi.median.log10sum.adjust[clust$clust==i],
       xlab = "DAPI score", main = "", xlim = c(-1.6, .7),
       lty = "blank", col = "gray40", breaks = 60)
  abline(v = mean(pData(eset)$dapi.median.log10sum.adjust[clust$clust==i]) +
           c(1, 0, -1) * sd(pData(eset)$dapi.median.log10sum.adjust[clust$clust==i]),
         col = "orange", lty = c(2, 1, 2))
}
