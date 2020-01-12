# Supplemental Figure 8

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
hist(pData(eset)$dapi.median.log10sum.adjust[clust$clust==1],
     xlab = "DAPI scores", main = "", xlim = c(-1.6, .7),
     lty = "blank", col = "gray40", breaks = 60)
abline(v = mean(pData(eset)$dapi.median.log10sum.adjust[clust$clust==1]) +
         c(1, 0, -1) * sd(pData(eset)$dapi.median.log10sum.adjust[clust$clust==1]),
       col = "orange", lty = c(2, 1, 2))
hist(pData(eset)$dapi.median.log10sum.adjust[clust$clust==2],
     xlab = "DAPI scores", main = "", xlim = c(-1.6, .7),
     lty = "blank", col = "gray40", breaks = 60)
abline(v = mean(pData(eset)$dapi.median.log10sum.adjust[clust$clust==2]) +
         c(1, 0, -1) * sd(pData(eset)$dapi.median.log10sum.adjust[clust$clust==2]),
       col = "orange", lty = c(2, 1, 2))
hist(pData(eset)$dapi.median.log10sum.adjust[clust$clust==3],
     xlab = "DAPI scores", main = "", xlim = c(-1.6, .7),
     lty = "blank", col = "gray40", breaks = 60)
abline(v = mean(pData(eset)$dapi.median.log10sum.adjust[clust$clust==3]) +
         c(1, 0, -1) * sd(pData(eset)$dapi.median.log10sum.adjust[clust$clust==3]),
       col = "orange", lty = c(2, 1, 2))


par(mfrow = c(2,2))
plot(x = pData(eset)$dapi.median.log10sum.adjust, ylim = c(-1.6, 1.6),
     y = pData(eset)$gfp.median.log10sum.adjust, pch = 16, col = "gray30", cex = .7,
     xlab = "DAPI score", ylab = "EGFP score")
plot(x = pData(eset)$dapi.median.log10sum.adjust, ylim = c(-1.6, 1.6),
     y = pData(eset)$rfp.median.log10sum.adjust, pch = 16, col = "gray30", cex = .7,
     xlab = "DAPI score", ylab = "mCherry score")
hist(pData(eset)$dapi.median.log10sum.adjust, breaks = 40, main = "",
     xlab = "DAPI score")
cor.test(pData(eset)$dapi.median.log10sum.adjust,
         pData(eset)$gfp.median.log10sum.adjust, method = "spearman")
cor.test(pData(eset)$dapi.median.log10sum.adjust,
         pData(eset)$rfp.median.log10sum.adjust, method = "spearman")



par(mfrow = c(2,2))
plot(x = 10^pData(eset)$dapi.median.log10sum, ylim = c(0, 2500), xlim = c(0, 2600),
     y = 10^pData(eset)$gfp.median.log10sum, pch = 16, col = "gray30", cex = .7,
     xlab = "DAPI intensity", ylab = "EGFP intensity")
abline(0, 1, col = "red")
plot(x = 10^pData(eset)$dapi.median.log10sum, ylim = c(0, 2500), xlim = c(0, 2600),
     y = 10^pData(eset)$rfp.median.log10sum, pch = 16, col = "gray30", cex = .7,
     xlab = "DAPI intensity", ylab = "mCherry intensity")
abline(0, 1, col = "red")
hist(10^pData(eset)$dapi.median.log10sum, breaks = 40, main = "",
     xlab = "DAPI intensity", xlim = c(0, 2600))
cor.test(pData(eset)$dapi.median.log10sum.adjust,
         pData(eset)$gfp.median.log10sum.adjust, method = "spearman")
cor.test(pData(eset)$dapi.median.log10sum.adjust,
         pData(eset)$rfp.median.log10sum.adjust, method = "spearman")

par(mfrow = c(2,2))
plot(x = 10^pData(eset)$dapi.median.log10sum, #ylim = c(-1.6, 1.6),
     y = 10^pData(eset)$gfp.median.log10sum, pch = 16, col = "gray30", cex = .7,
     xlab = "DAPI intensity sum", ylab = "EGFP intensity sum")
plot(x = 10^pData(eset)$dapi.median.log10sum.adjust, ylim = c(-1.6, 1.6),
     y = 10^pData(eset)$rfp.median.log10sum.adjust, pch = 16, col = "gray30", cex = .7,
     xlab = "DAPI score", ylab = "mCherry score")
hist(10^pData(eset)$dapi.median.log10sum, breaks = 40, main = "",
     xlab = "DAPI sum")
cor.test(pData(eset)$dapi.median.log10sum.adjust,
         pData(eset)$gfp.median.log10sum.adjust, method = "spearman")
cor.test(pData(eset)$dapi.median.log10sum.adjust,
         pData(eset)$rfp.median.log10sum.adjust, method = "spearman")
