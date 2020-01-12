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

# basically difference is just a rotation
# but a bit nonlinear transfomration
# what's the impact of this difference?
theta_raw = coord2rad(scale(cbind(pData(eset)$rfp.median.log10sum.adjust,
                                  pData(eset)$gfp.median.log10sum.adjust)))

plot(as.numeric(theta_raw), as.numeric(theta))
plot(2*pi-as.numeric(theta_raw), theta_raw)
plot(shift_origin(2*pi-as.numeric(theta_raw), pi+.08), theta)
abline(0,1, col = "red")
plot(shift_origin(shift_origin(2*pi-as.numeric(theta_raw), pi+.08),
                  2.715*pi/4), theta_final)

# shifting theta_raw (based on intensities) to be the same as theta_final
theta_raw_final = shift_origin(shift_origin(2*pi-as.numeric(theta_raw), pi+.08),
             2.715*pi/4)



fits_all <- readRDS(file.path(dir, "output/npreg-trendfilter-quantile.Rmd/fit.quant.rds"))
pve_all <- sapply(fits_all, "[[", 3)
pve_all_ord <- pve_all[order(pve_all, decreasing = T)]
genes <- names(pve_all_ord)[1:5]
labs <- c("CDK1", "UBE2C", "TOP2A", "HIST1H4E", "HIST1H4C")

data_quant <- readRDS(file.path(dir, "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))
sample_ord <- rownames(pData(eset))[order(theta_raw_final)]
data_quant_ord <- data_quant[,match(sample_ord,colnames(data_quant))]
fits_tmp <- lapply(1:5, function(i) {
  ii <- which(rownames(data_quant_ord)==genes[i])
  fit_g <- fit_trendfilter_generic(data_quant_ord[ii,])
  fun_g <- approxfun(x=as.numeric(theta_raw_final[order(theta_raw_final)]),
                     y=as.numeric(fit_g$trend.yy), rule=2)
  return(fun_g)
})
names(fits_tmp) <- genes


par(mfrow=c(1,5), mar=c(2,2,2,1))
for (i in 1:5) {
  plot(x=theta_final[order(theta_raw_final)],
       y=data_quant_ord[rownames(data_quant_ord)==genes[i],], col="gray50",
       xlab="Fucci phase",
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


