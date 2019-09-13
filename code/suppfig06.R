# Supplemental Figure 6

dir <- "/project2/gilad/joycehsiao/fucci-seq"

eset <- readRDS(file.path(dir, "data/eset-final.rds"))

counts <- exprs(eset)
counts <- counts[grep("ENSG", rownames(counts)), ]
pdata <- pData(eset)


# Intensity batch effect ---------------------------------------------------------
lm.rfp <- lm(rfp.median.log10sum~factor(chip_id)+factor(experiment) + factor(image_label),
             data = pdata)
lm.gfp <- lm(gfp.median.log10sum~factor(chip_id)+factor(experiment) + factor(image_label),
             data = pdata)
lm.dapi <- lm(dapi.median.log10sum~factor(chip_id)+factor(experiment) + factor(image_label),
              data = pdata)

library(car)
aov.lm.rfp <- Anova(lm.rfp, type = "III")
aov.lm.gfp <- Anova(lm.gfp, type = "III")
aov.lm.dapi <- Anova(lm.dapi, type = "III")

aov.lm.rfp
aov.lm.gfp

library(ggplot2)
library(RColorBrewer)
library(cowplot)
pdata$experiment <- factor(pdata$experiment)
plot_grid(
  ggplot(pdata,
         aes(x=experiment, y=gfp.median.log10sum,
             fill=experiment)) +
    geom_violin(alpha=.5) + geom_boxplot(width=.1, fill="gray90") +
    ylab("FUCC score") + xlab("Individual") +
    ggtitle("GFP before batch-correction") + labs(fill="experiment") +
    theme(axis.text.x  = element_text(angle=90, vjust=1, hjust=1,size=10),
          axis.text.y  = element_text(angle=0, vjust=0.2, size=10)) +
    theme(legend.position="none"),
  ggplot(pdata,
         aes(x=experiment, y=rfp.median.log10sum, fill=experiment)) +
    geom_violin(alpha=.5) + geom_boxplot(width=.1, fill="gray90") +
    ylab("FUCC score") + xlab("Individual") +
    ggtitle("RFP before batch-correction") + labs(fill="experiment") +
    theme(axis.text.x  = element_text(angle=90, vjust=1, hjust=1,size=10),
          axis.text.y  = element_text(angle=0, vjust=0.2, size=10)) +
    theme(legend.position="none"),
  ggplot(pdata,
         aes(x=experiment, y=gfp.median.log10sum.adjust, fill=experiment)) +
    geom_violin(alpha=.5) + geom_boxplot(width=.1, fill="gray90") +
    ylab("FUCC score") + xlab("Individual") +
    ggtitle("GFP after batch-correction") + labs(fill="experiment") +
    theme(axis.text.x  = element_text(angle=90, vjust=1, hjust=1,size=10),
          axis.text.y  = element_text(angle=0, vjust=0.2, size=10)) +
    theme(legend.position="none"),
  ggplot(pdata,
         aes(x=experiment, y=rfp.median.log10sum.adjust, fill=experiment)) +
    geom_violin(alpha=.5) + geom_boxplot(width=.1, fill="gray90") +
    ylab("FUCC score") + xlab("Individual") +
    ggtitle("RFP after batch-correction") + labs(fill="experiment") +
    theme(axis.text.x  = element_text(angle=90, vjust=1, hjust=1,size=10),
          axis.text.y  = element_text(angle=0, vjust=0.2, size=10)) +
    theme(legend.position="none"),
  labels=LETTERS[1:4])

