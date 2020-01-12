# Figure 3

library(Biobase)
library(circular)
library(peco)

dir <- "/project2/gilad/joycehsiao/fucci-seq"
eset <- readRDS(file=file.path(dir, "data/eset-final.rds"))
fdata <- fData(eset)[grep("ENSG", rownames(fData(eset))),]

# Figure 3B. Prediction error in the test data

# reCAT results
# see code/ourdata_recat.R
load(file.path(dir, "analysis/revision/output/ourdata_recat.rda"))
samples_ordered <- rownames(test_exp)[ordIndex]
phase_pred <- seq(0, 2*pi, length.out= length(ordIndex))
names(phase_pred) <- samples_ordered
diff_recat <- data.frame(phase_ref=pData(eset)$theta[match(samples_ordered, colnames(eset))])
diff_recat$ind <- pData(eset)$chip_id[match(names(phase_pred), colnames(eset))]
diff_recat$phase_pred_rot <- rotation(ref_var= diff_recat$phase_ref, shift_var = phase_pred)
rownames(diff_recat) <- samples_ordered
diff_recat$diff_time <- circ_dist(diff_recat$phase_ref, diff_recat$phase_pred_rot)
diff_recat$method <- "reCAT"
diff_recat <- diff_recat[,c(1,3,4,2,5)]


# Oscope results
# see code/ourdata_oscope.R
load(file = file.path(dir, "analysis/revision/output/ourdata_oscope_366genes.rda"))

samples_ordered <- colnames(DataInput)[ENIRes[["cluster2"]]]
phase_oscope <- seq(0, 2*pi, length.out= length(ENIRes[[2]]))

diff_oscope <- data.frame(phase_ref=pData(eset)$theta[match(samples_ordered, colnames(eset))])
diff_oscope$phase_pred_rot <- rotation(ref_var= diff_oscope$phase_ref, shift_var = phase_pred)
diff_oscope$diff_time <- circ_dist(diff_oscope$phase_ref, diff_oscope$phase_pred_rot)
diff_oscope$ind <- pData(eset)$chip_id[match(names(phase_pred), colnames(eset))]
diff_oscope$method <- "Oscope"


# cyclone results
# see code/ourdata_cyclone.R
inds <- c("NA18511", "NA18855", "NA18870", "NA19098", "NA19101", "NA19160")
diff_cyclone <- do.call(rbind, lapply(1:length(inds), function(i) {
  out_cyclone <- readRDS(file.path(dir,
                                   paste0("analysis/revision/output/ourdata_ind_cyclone/ourdata_cyclone_",
                                          inds[i],".rds")))
  phase_peco <- pData(eset)$theta[match(rownames(out_cyclone), rownames(pData(eset)))]
  phase_cyclone_rot <- rotation(ref_var= phase_peco, shift_var = out_cyclone$phase_pred)
  return(data.frame(phase_ref = phase_peco,
                    phase_pred_rot = phase_cyclone_rot,
                    diff_time = circ_dist(phase_peco, phase_cyclone_rot),
                    ind = inds[i],
                    method = "Cyclone"))
}) )


# Seurat results
# see code/ourdata_seurat.R
source(file.path(dir, "peco/R/run_seurat.R"))
seurat.genes <- readLines(con = file.path(dir,
                                          "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt"))
seurat.genes <- list(s.genes=seurat.genes[1:43],
                     g2m.genes=seurat.genes[44:97])

diff_seurat <- do.call(rbind, lapply(1:length(inds), function(i) {
  #  print(i)
  data <- readRDS(file.path(dir,paste0("analysis/revision/output/data_training_test/ind_",
                                       inds[i], "_data_test.rds")))
  Y_fit <- data$log2cpm_test
  rownames(Y_fit) <- fdata$name[match(rownames(Y_fit), rownames(fdata))]
  fit <- run_seurat(Y_fit,
                    s.genes=seurat.genes$s.genes,
                    g2m.genes=seurat.genes$g2m.genes)
  fit <- as.list(fit)
  seurat.pca <- prcomp(cbind(fit$G2M, fit$S), scale=TRUE)
  phase_pred <- intensity2circle(seurat.pca$x, plot.it = F, method = "trig")
  names(phase_pred) <- colnames(Y_fit)
  phase_ref <- shift_origin(pData(eset)$theta[match(names(phase_pred),
                    rownames(pData(eset)))], origin = 3*pi/4)
  phase_pred_rot <- rotation(ref_var=phase_ref, shift_var=phase_pred)
  return(data.frame(phase_ref = phase_ref,
                    phase_pred_rot = phase_pred_rot,
                    diff_time = circ_dist(phase_ref, phase_pred_rot),
                    class_seurat = as.character(fit$assignments),
                    ind = inds[i],
                    method = "Seurat", stringsAsFactors = FALSE))
}) )



# peco
# see code/ourdata_peco.R
inds <- c("NA19098","NA18511","NA18870","NA19101","NA18855","NA19160")
ngene <- 5
diff_peco <- do.call(rbind, lapply(1:length(inds), function(i) {
  ind <- inds[i]
  fl_name <- list.files(file.path(dir, "data/results/finalizing"),
                        pattern=paste0("ind_",ind,"_results_overallcyclical.top",ngene,".rds"),
                        full.names = TRUE)
  df <- readRDS(fl_name)
  df_pred <- df$fit.test$fit.supervised
  return(data.frame(phase_ref = df_pred$ref_time,
                    phase_pred_rot = df_pred$pred_time_shift,
                    diff_time = df_pred$diff_time,
                    ind = ind,
                    method = "peco"))
}))



# gather prediction results from the different methods
diff_cyclone <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_cyclone.rds"))
diff_recat <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_recat.rds"))
diff_peco <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_peco.rds"))
diff_seurat <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_seurat.rds"))
diff_oscope <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_oscope.rds"))

diff_all <- rbind(diff_recat,
                  diff_cyclone,
                  diff_peco,
                  diff_seurat,
                  diff_oscope)

# compare prediction results between methods
pval_difftime <- do.call(rbind, lapply(unique(diff_all$ind), function(ind) {
  data.frame(ind=ind,
             p_cyclone=wilcox.test(diff_peco[diff_peco$ind==ind,]$diff_time/2/pi,
                                   diff_cyclone[diff_cyclone$ind==ind,]$diff_time/2/pi)$p.value,
             p_seurat=wilcox.test(diff_peco[diff_peco$ind==ind,]$diff_time/2/pi,
                                  diff_seurat[diff_peco$ind==ind,]$diff_time/2/pi)$p.value,
             p_oscope=wilcox.test(diff_peco[diff_peco$ind==ind,]$diff_time/2/pi,
                                  diff_oscope[diff_peco$ind==ind,]$diff_time/2/pi)$p.value,
             p_recat=wilcox.test(diff_peco[diff_peco$ind==ind,]$diff_time/2/pi,
                                 diff_recat[diff_peco$ind==ind,]$diff_time/2/pi)$p.value,
             ratio_cyclone=mean(diff_peco[diff_peco$ind==ind,]$diff_time)/mean(diff_cyclone[diff_cyclone$ind==ind,]$diff_time),
             ratio_seurat=mean(diff_peco[diff_peco$ind==ind,]$diff_time)/mean(diff_seurat[diff_seurat$ind==ind,]$diff_time),
             ratio_oscope=mean(diff_peco[diff_peco$ind==ind,]$diff_time)/mean(diff_oscope[diff_oscope$ind==ind,]$diff_time),
             ratio_recat=mean(diff_peco[diff_peco$ind==ind,]$diff_time)/mean(diff_recat[diff_recat$ind==ind,]$diff_time))
}))
round(pval_difftime[,-1],4)

mean(pval_difftime$ratio_cyclone)
mean(pval_difftime$ratio_seurat)
mean(pval_difftime$ratio_oscope)
mean(pval_difftime$ratio_recat)
# NA19098: cyclone < .05; else < .00001
# NA18511: cyclone < .05; else

# all 6 methods
library(tidyverse)
library(ggplot2)
library(wesanderson)
cols_method <- c(wes_palette("FantasticFox1")[1:3],
                 wes_palette("Darjeeling1")[2],
                 wes_palette("FantasticFox1")[5])
diff_all$method <- factor(diff_all$method, levels=c("peco", "Cyclone", "Seurat", "Oscope", "reCAT"))
diff_all %>% group_by(ind, method) %>%
  summarise(diff_mean = mean(diff_time/2/pi),
            diff_se = sd(diff_time/2/pi)/sqrt(length(diff_time/2/pi))) %>%
  ggplot(., aes(x=method,y=diff_mean, col=method, group=method)) +
  geom_point(aes(x=method, y=diff_mean)) +
  #    geom_vline(xintercept=seq(1, 5, 1), col="gray90", lty=1) +
  geom_hline(yintercept=seq(.1, .25, .01), col="gray90", alpha=.6, lty=1) +
  geom_errorbar(aes(ymin=diff_mean-diff_se,
                    ymax=diff_mean+diff_se, col=method), width=.2, alpha=.5) +
  #    geom_hline(yintercept=.25, col="gray30", lty=1) +
  scale_colour_manual(values = c(as.vector(cols_method))) +
  ylab("Prediction error (% circle)") + xlab("Method") +
  labs(color="Test data (Cell line)") + ylim(0,.25) +
  facet_wrap(~ind) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1,size=10))


# only 3 continuous methods
library(tidyverse)
library(ggplot2)
library(wesanderson)
cols_method <- c(wes_palette("FantasticFox1")[1],
                 wes_palette("Darjeeling1")[2],
                 wes_palette("FantasticFox1")[5])
diff_all$method <- factor(diff_all$method, levels=c("peco", "Oscope", "reCAT"))
diff_all %>% filter(method %in% c("peco", "Oscope", "reCAT")) %>%
  group_by(ind, method) %>%
  summarise(diff_mean = mean(diff_time/2/pi),
            diff_se = sd(diff_time/2/pi)/sqrt(length(diff_time/2/pi))) %>%
  ggplot(., aes(x=method,y=diff_mean, col=ind, group=ind)) +
  geom_point(aes(x=method, y=diff_mean), size = .5) +
  geom_hline(yintercept=seq(.1, .25, .01), col="gray90", alpha=.6, lty=1) +
  geom_errorbar(aes(ymin=diff_mean-diff_se,
                    ymax=diff_mean+diff_se, col=ind), width=.2, alpha=.5) +
#  scale_colour_manual(values = c(as.vector(cols_method))) +
  ylab("Prediction error (% circle)") + xlab("Method") +
  labs(color="Test data (Cell line)") + ylim(0,.3) +
  geom_hline(yintercept = .25, col = "red") + theme_classic() +
#  facet_wrap(~method) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1,size=10))


# Figure 3C. Estimate cyclic trend of gene expression for cells from NA18511 samples
# See code/suppfig15.R for reCAT and Oscope
pca = prcomp(cbind(pData(eset)$rfp.median.log10sum.adjust,
                    pData(eset)$gfp.median.log10sum.adjust))
theta = coord2rad(pca$x)
theta_final = shift_origin(as.numeric(theta), 3*pi/4)
names(theta_final) = rownames(pData(eset))
library(cluster)
pam_res <- pam(cbind(pData(eset)$rfp.median.log10sum.adjust,
                     pData(eset)$gfp.median.log10sum.adjust), k=3)
clust <- data.frame(clust=pam_res$clustering,
                    sample_id=rownames(pData(eset)))
plot(theta_final, clust$clust)
b1 <- mean(max(range(theta_final[clust$clust==2])), min(range(theta_final[clust$clust==3])))
b2 <- mean(max(range(theta_final[clust$clust==3])), max(range(theta_final[clust$clust==1])))
b3 <- mean(min(range(theta_final[clust$clust==2])), min(range(theta_final[clust$clust==1])))
abline(v=c(b1,b2,b3), lty=2)

# seurat discrete classes vs FUCCI phase
library(Seurat)
eset <- readRDS(file=file.path(dir,"data/eset-final.rds"))
fdata <- fData(eset)
fdata <- fdata[grep("ENSG", rownames(eset)), ]
cc.genes <- readLines(con = file.path(dir,
      "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt"))
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
counts = exprs(eset)[grep("ENSG", rownames(eset)), ]
rownames(counts) = fdata$name[match(rownames(counts),rownames(fdata))]
obj <- CreateSeuratObject(counts = counts)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst")
obj <- ScaleData(obj, features = rownames(obj))
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
out_seurat <- obj[[]]

out_seurat2 = out_seurat[match(names(theta_final), rownames(out_seurat)),]
all.equal(rownames(out_seurat2), names(theta_final))

xr1

# Cyclone discrete classes vs FUCCI phase
out_cyclone <- readRDS(file.path(dir, "analysis/revision/output/ourdata_phase_cyclone.rds"))
out_cyclone2 = out_cyclone$phases[match(names(theta_final),
                                        names(out_cyclone$phases))]
all.equal(names(out_cyclone2), names(theta_final))

par(mfrow = c(1,3))
hist(theta_final[out_cyclone2=="G1"],
     breaks = 80, ylim = c(0,40), axes = F, ann = F, lty = "blank",
     col = as.character(yarrr::piratepal("espresso")[3]))
axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
              labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                       expression(2*pi)))
abline(v=c(b1,b2,b3), lty=2)
hist(theta_final[out_cyclone2=="S"],
     breaks = 80, ylim = c(0,40), axes = F, ann = F, lty = "blank",
     col = as.character(yarrr::piratepal("espresso")[2]))
axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
              labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                       expression(2*pi)))
abline(v=c(b1,b2,b3), lty=2)
hist(theta_final[out_cyclone2=="G2M"],
     breaks = 90, ylim = c(0,40), axes = F, ann = F, lty = "blank",
     col = as.character(yarrr::piratepal("espresso")[1]))
axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
              labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                       expression(2*pi)))
abline(v=c(b1,b2,b3), lty=2)

