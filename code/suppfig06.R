# Supplemental Figure 6

library(SingleCellExperiment)
library(circular)
library(peco)
library(cluster)
library(ggplot2)

sce <- readRDS("data/sce-final.rds")
pdata = data.frame(colData(sce))

par(mfrow = c(2,2))
plot(x = 10^pdata$dapi.median.log10sum, ylim = c(0, 2500), xlim = c(0, 2600),
     y = 10^pdata$gfp.median.log10sum, pch = 16, col = "gray30", cex = .7,
     xlab = "DAPI intensity", ylab = "EGFP intensity")
abline(0, 1, col = "red")
plot(x = 10^pdata$dapi.median.log10sum, ylim = c(0, 2500), xlim = c(0, 2600),
     y = 10^pdata$rfp.median.log10sum, pch = 16, col = "gray30", cex = .7,
     xlab = "DAPI intensity", ylab = "mCherry intensity")
abline(0, 1, col = "red")
hist(10^pdata$dapi.median.log10sum, breaks = 40, main = "",
     xlab = "DAPI intensity", xlim = c(0, 2600))
cor.test(pdata$dapi.median.log10sum.adjust,
         pdata$gfp.median.log10sum.adjust, method = "spearman")
cor.test(pdata$dapi.median.log10sum.adjust,
         pdata$rfp.median.log10sum.adjust, method = "spearman")

