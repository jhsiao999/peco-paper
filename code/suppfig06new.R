# Supplemental Figure 6

library(Biobase)
library(circular)
library(peco)
library(cluster)
library(ggplot2)

dir <- "/project2/gilad/joycehsiao/fucci-seq"
eset <- readRDS(file.path(dir, "data/eset-final.rds"))

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

