# Supplemental Figure 6
#   FUCCI and DAPI intensity before C1 plate effect correction of the 888 samples in the data.

library(SingleCellExperiment)
library(circular)
library(peco)
library(cluster)
library(ggplot2)

sce <- readRDS("data/sce-final.rds")
pdata = data.frame(colData(sce))


# split samples into 40 bins based on DAPI intensity distribution
# then plot GFP/RFP distribution in the 40 bins
pdata %>%
  mutate(dapi_bins=cut(10^pdata$dapi.median.log10sum, breaks=40)) %>%
  group_by(dapi_bins) %>%
  ggplot(., aes(x = dapi_bins, y = 10^pdata$gfp.median.log10sum)) +
    geom_boxplot() +
    ylab("EGFP intensity") +
    scale_x_discrete("DAPI intensity (40 equally-spaced bins)",
        breaks = levels(cut(10^pdata$dapi.median.log10sum, breaks=40))[c(1,5,10,15,20,25,30,35,40)],
        labels = levels(cut(10^pdata$dapi.median.log10sum, breaks=40))[c(1,5,10,15,20,25,30,35,40)]) +
    theme(axis.text.x = element_text(size = 10, angle = 90,
                                     vjust=.5, hjust=1))

# split samples into 40 bins based on DAPI intensity distribution
# then plot GFP/RFP distribution in the 40 bins
pdata %>%
  mutate(dapi_bins=cut(10^pdata$dapi.median.log10sum, breaks=40)) %>%
  group_by(dapi_bins) %>%
  ggplot(., aes(x = dapi_bins, y = 10^pdata$rfp.median.log10sum)) +
  geom_boxplot() +
  ylab("mCherry intensity") +
  scale_x_discrete("DAPI intensity (40 equally-spaced bins)",
    breaks = levels(cut(10^pdata$dapi.median.log10sum, breaks=40))[c(1,5,10,15,20,25,30,35,40)],
    labels = levels(cut(10^pdata$dapi.median.log10sum, breaks=40))[c(1,5,10,15,20,25,30,35,40)]) +
  theme(axis.text.x = element_text(size = 10, angle = 90,
                                   vjust=.5, hjust=1))

# histogram displaying distribution of DAPI intensities
hist(10^pdata$dapi.median.log10sum, breaks = 40, main = "",
     xlab = "DAPI intensity", xlim = c(0, 2600))


# compute spearman correlation between DAPI intensities and GFP/RFP intensities
cor.test(pdata$dapi.median.log10sum.adjust,
         pdata$gfp.median.log10sum.adjust, method = "spearman")
cor.test(pdata$dapi.median.log10sum.adjust,
         pdata$rfp.median.log10sum.adjust, method = "spearman")
cor.test(pdata$dapi.median.log10sum.adjust,
         pdata$gfp.median.log10sum.adjust, method = "pearson")
cor.test(pdata$dapi.median.log10sum.adjust,
         pdata$rfp.median.log10sum.adjust, method = "pearson")

cor.test(pdata$gfp.median.log10sum.adjust,
         pdata$rfp.median.log10sum.adjust, method = "pearson")
