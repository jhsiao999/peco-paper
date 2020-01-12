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

fits_tmp <- lapply(seq_len(nrow(which_sig_cycle_genes)), function(i) {
  ii <- which(rownames(data_quant_ord)==which_sig_cycle_genes$ensg[i])
  fit_g <- fit_trendfilter_generic(data_quant_ord[ii,])
  fun_g <- approxfun(x=as.numeric(theta_final[order(theta_final)]),
                     y=as.numeric(fit_g$trend.yy), rule=2)
  return(fun_g)
})
names(fits_tmp) <- which_sig_cycle_genes$names

par(mfrow=c(4,5), mar=c(2,2,2,1))
for (i in 1:20) {
  plot(x=theta_final[order(theta_final)],
       y=data_quant_ord[rownames(data_quant_ord)==which_sig_cycle_genes$ensg[i],],
       col="gray50",
       xlab="Fucci phase",
       ylab="", axes=F, cex=.3, pch=16,
       main = paste(which_sig_cycle_genes$names[i], paste0("(", i ,")")))
  lines(x = seq(0, 2*pi, length.out=200),
        y = fits_tmp[[which(names(fits_tmp) == which_sig_cycle_genes$names[i])]](seq(0, 2*pi, length.out=200)),
        col=wesanderson::wes_palette("FantasticFox1")[1], lty=1, lwd=2)
  axis(2);
  axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                               expression(3*pi/2), expression(2*pi)))
  abline(h=0, col="black", lty=1, lwd=.7)
  abline(v=c(b1,b2,b3), lty=2)
}

par(mfrow=c(4,5), mar=c(2,2,2,1))
for (i in 21:40) {
  plot(x=theta_final[order(theta_final)],
       y=data_quant_ord[rownames(data_quant_ord)==which_sig_cycle_genes$ensg[i],],
       col="gray50",
       xlab="Fucci phase",
       ylab="", axes=F, cex=.3, pch=16,
       main = paste(which_sig_cycle_genes$names[i], paste0("(", i ,")")))
  lines(x = seq(0, 2*pi, length.out=200),
        y = fits_tmp[[which(names(fits_tmp) == which_sig_cycle_genes$names[i])]](seq(0, 2*pi, length.out=200)),
        col=wesanderson::wes_palette("FantasticFox1")[1], lty=1, lwd=2)
  axis(2);
  axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                               expression(3*pi/2), expression(2*pi)))
  abline(h=0, col="black", lty=1, lwd=.7)
  abline(v=c(b1,b2,b3), lty=2)
}


par(mfrow=c(4,5), mar=c(2,2,2,1))
for (i in 41:54) {
  plot(x=theta_final[order(theta_final)],
       y=data_quant_ord[rownames(data_quant_ord)==which_sig_cycle_genes$ensg[i],],
       col="gray50",
       xlab="Fucci phase",
       ylab="", axes=F, cex=.3, pch=16,
       main = paste(which_sig_cycle_genes$names[i], paste0("(", i ,")")))
  lines(x = seq(0, 2*pi, length.out=200),
        y = fits_tmp[[which(names(fits_tmp) == which_sig_cycle_genes$names[i])]](seq(0, 2*pi, length.out=200)),
        col=wesanderson::wes_palette("FantasticFox1")[1], lty=1, lwd=2)
  axis(2);
  axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                               expression(3*pi/2), expression(2*pi)))
  abline(h=0, col="black", lty=1, lwd=.7)
  abline(v=c(b1,b2,b3), lty=2)
}


