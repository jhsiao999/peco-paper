# Supplemental Figure 4
# Quality metrics in individual cell lines

dir <- "/project2/gilad/joycehsiao/fucci-seq"

eset <- readRDS(file.path(dir, "data/eset-final.rds"))

counts <- exprs(eset)
counts <- counts[grep("ENSG", rownames(counts)), ]
pdata <- pData(eset)

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

library(ggplot2)
library(cowplot)
pdata$dropout <- colMeans(log2cpm.all >0)
pdata$experiment <- factor(pdata$experiment)
pdata$chip_id <- factor(pdata$chip_id)
ggplot(pdata, aes(x=chip_id, y=log10(molecules),
                  group=chip_id, col=chip_id)) +
  geom_boxplot() + coord_flip() +
  ylab("log10 molecule count") + xlab("Individual") +
  scale_color_brewer(palette="Dark2") +
  ggtitle("Library size and individual") + theme(legend.position="none")
ggplot(pdata, aes(x=chip_id, y=dropout,
                group=chip_id, col=chip_id)) +
  geom_boxplot() + coord_flip() +
  ylab("Sample detection rate") + xlab("Individual") +
  scale_color_brewer(palette="Dark2") +
  ggtitle("Detection rate and individual") + theme(legend.position="none")
ggplot(pdata, aes(x=chip_id, fill=chip_id)) +
  geom_bar() + coord_flip() +
  ylab("Singletons") + xlab("Individual") +
  scale_fill_brewer("Individual",palette="Dark2") +
  ggtitle("Singletons and individual") + theme(legend.position="none")


