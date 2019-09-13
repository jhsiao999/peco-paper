# Supplemental Figure S14

library(Biobase)
library(peco)

dir <- "/project2/gilad/joycehsiao/fucci-seq"

eset <- readRDS(file=file.path(dir, "data/eset-final.rds"))
pdata <- pData(eset)
fdata <- fData(eset)

log2cpm.quant <- readRDS(file.path(dir, "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))

inds <- c("NA18511", "NA18855", "NA18870", "NA19098", "NA19101", "NA19160")


# recat
load(file.path(dir, paste0("analysis/revision/output/ourdata_recat.rda")))
sample_ordered <- rownames(test_exp)[ordIndex]
phase_pred <- seq(0, 2*pi, length.out= length(ordIndex))
names(phase_pred) <- sample_ordered

diff_recat <- data.frame(phase_ref=pData(eset)$theta[match(samples_ordered, colnames(eset))])
diff_recat$ind <- pData(eset)$chip_id[match(names(phase_pred), colnames(eset))]

diff_recat$phase_pred_rot <- rotation(ref_var= diff_recat$phase_ref, shift_var = phase_pred)
rownames(diff_recat) <- samples_ordered

diff_recat$diff_time <- circ_dist(diff_recat$phase_ref, diff_recat$phase_pred_rot)
diff_recat$method <- "reCAT"
diff_recat <- diff_recat[,c(1,3,4,2,5)]
#saveRDS(diff_recat, file=file.path(dir, "analysis/revision/output/fit_diff_recat.rds"))


# cyclone
diff_cyclone <- do.call(rbind, lapply(1:length(inds), function(i) {
  out_cyclone <- readRDS(file.path(dir,
                                   paste0("analysis/revision/output/ourdata_ind_cyclone/ourdata_cyclone_",
                                          inds[i],".rds")))
  phase_peco <- pdata$theta[match(rownames(out_cyclone), rownames(pdata))]
  phase_cyclone_rot <- rotation(ref_var= phase_peco, shift_var = out_cyclone$phase_pred)
  return(data.frame(phase_ref = phase_peco,
                    phase_pred_rot = phase_cyclone_rot,
                    diff_time = circ_dist(phase_peco, phase_cyclone_rot),
                    ind = inds[i],
                    method = "Cyclone"))
}) )


# Seurat
source(file.path(dir, "peco/R/run_seurat.R"))
library(Biobase)
library(peco)
eset <- readRDS(file=file.path(dir, "data/eset-final.rds"))
fdata <- fData(eset)[grep("ENSG", rownames(fData(eset))),]
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
  phase_ref <- pdata$theta[match(names(phase_pred), rownames(pdata))]
  phase_pred_rot <- rotation(ref_var=phase_ref, shift_var=phase_pred)
  return(data.frame(phase_ref = phase_ref,
                    phase_pred_rot = phase_pred_rot,
                    diff_time = circ_dist(phase_ref, phase_pred_rot),
                    ind = inds[i],
                    method = "Seurat"))
}) )


# Oscope results
eset <- readRDS(file=file.path(dir, "data/eset-final.rds"))
load(file = file.path(dir, "analysis/revision/output/ourdata_oscope_366genes.rda"))

samples_ordered <- colnames(DataInput)[ENIRes[["cluster2"]]]
phase_oscope <- seq(0, 2*pi, length.out= length(ENIRes[[2]]))

diff_oscope <- data.frame(phase_ref=pData(eset)$theta[match(samples_ordered, colnames(eset))])
diff_oscope$phase_pred_rot <- rotation(ref_var= diff_oscope$phase_ref, shift_var = phase_pred)
diff_oscope$diff_time <- circ_dist(diff_oscope$phase_ref, diff_oscope$phase_pred_rot)
diff_oscope$ind <- pData(eset)$chip_id[match(names(phase_pred), colnames(eset))]
diff_oscope$method <- "Oscope"




# compre predictions between methods


diff_cyclone <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_cyclone.rds"))
diff_recat <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_recat.rds"))
diff_peco <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_peco.rds"))
diff_seurat <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_seurat.rds"))
diff_oscope <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_oscope.rds"))

diff_peco %>% group_by(ind) %>% summarise(mn=mean(diff_time)/2/pi,
                                          sd=sd(diff_time/2/pi)/sqrt(sum(phase_pred_rot>0)))

diff_cyclone %>% group_by(ind) %>% summarise(mn=mean(diff_time)/2/pi,
                                             sd=sd(diff_time/2/pi)/sqrt(sum(phase_pred_rot>0)))

diff_seurat %>% group_by(ind) %>% summarise(mn=mean(diff_time)/2/pi,
                                            sd=sd(diff_time/2/pi)/sqrt(sum(phase_pred_rot>0)))

diff_oscope %>% group_by(ind) %>% summarise(mn=mean(diff_time)/2/pi,
                                            sd=sd(diff_time/2/pi)/sqrt(sum(phase_pred_rot>0)))

diff_recat %>% group_by(ind) %>% summarise(mn=mean(diff_time)/2/pi,
                                           sd=sd(diff_time/2/pi)/sqrt(sum(phase_pred_rot>0)))

diff_all <- rbind(diff_recat,
                  diff_cyclone,
                  diff_peco,
                  diff_seurat,
                  diff_oscope)

pval_difftime <- do.call(rbind, lapply(inds, function(ind) {
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
cbind(inds, round(pval_difftime[,-1],4))

mean(pval_difftime$ratio_cyclone)
mean(pval_difftime$ratio_seurat)
mean(pval_difftime$ratio_oscope)
mean(pval_difftime$ratio_recat)

# peco results
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



diff_cyclone <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_cyclone.rds"))
diff_recat <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_recat.rds"))
diff_peco <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_peco.rds"))
diff_seurat <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_seurat.rds"))
diff_oscope <- readRDS(file=file.path(dir, "analysis/revision/output/fit_diff_oscope.rds"))

diff_peco %>% group_by(ind) %>% summarise(mn=mean(diff_time)/2/pi,
                                          sd=sd(diff_time/2/pi)/sqrt(sum(phase_pred_rot>0)))

diff_cyclone %>% group_by(ind) %>% summarise(mn=mean(diff_time)/2/pi,
                                             sd=sd(diff_time/2/pi)/sqrt(sum(phase_pred_rot>0)))

diff_seurat %>% group_by(ind) %>% summarise(mn=mean(diff_time)/2/pi,
                                            sd=sd(diff_time/2/pi)/sqrt(sum(phase_pred_rot>0)))

diff_oscope %>% group_by(ind) %>% summarise(mn=mean(diff_time)/2/pi,
                                            sd=sd(diff_time/2/pi)/sqrt(sum(phase_pred_rot>0)))

diff_recat %>% group_by(ind) %>% summarise(mn=mean(diff_time)/2/pi,
                                           sd=sd(diff_time/2/pi)/sqrt(sum(phase_pred_rot>0)))
