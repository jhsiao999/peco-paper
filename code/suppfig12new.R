# Supplemental Fig S10

# import data
dir <- "/project2/gilad/joycehsiao/fucci-seq"
eset <- readRDS(file=file.path(dir, "data/eset-final.rds"))
log2cpm_quant <- readRDS(file.path(dir,
                       "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))
log2cpm_quant <- log2cpm_quant[, match(colnames(eset), colnames(log2cpm_quant))]
theta <- pData(eset)$theta




# partition data to training and test datasets
makedata_supervised <- function(eset, log2cpm_quant,
                                theta) {

  library(Biobase)
  library(peco)
  pdata <- pData(eset)
  fdata <- fData(eset)
  counts <- exprs(eset)[grep("ENSG", rownames(eset)), ]
  log2cpm <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

  for (ind in unique(pdata$chip_id)) {
    ii_test <- c(1:nrow(pdata))[which(pdata$chip_id == ind)]
    ii_train <- c(1:nrow(pdata))[which(pdata$chip_id != ind)]

    pdata_test <- pdata[ii_test,]
    pdata_train <- pdata[ii_train,]

    log2cpm_quant_test <- log2cpm_quant[,ii_test]
    log2cpm_quant_train <- log2cpm_quant[,ii_train]
    theta <- pdata$theta
    names(theta) <- rownames(pdata)

    log2cpm_test <- log2cpm[,ii_test]
    log2cpm_train <- log2cpm[,ii_train]

    counts_test <- counts[,ii_test]
    counts_train <- counts[,ii_train]

    theta_test <- theta[ii_test]
    theta_train <- theta[ii_train]

    #sig.genes <- readRDS("output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.476.rds")
    data_training <- list(theta_train=theta_train,
                          log2cpm_quant_train=log2cpm_quant_train,
                          log2cpm_train=log2cpm_train,
                          counts_train=counts_train,
                          pdata_train=pdata_train,
                          fdata=fdata)

    data_test <- list(theta_test=theta_test,
                      log2cpm_quant_test=log2cpm_quant_test,
                      log2cpm_test=log2cpm_test,
                      counts_test = counts_test,
                      pdata_test=pdata_test,
                      fdata=fdata)

    saveRDS(data_training,
            file=file.path(dir, paste0("analysis/revision/output/data_training_test/ind_",ind,"_data_training.rds")))
    saveRDS(data_test,
            file=file.path(dir, paste0("analysis/revision/output/data_training_test/ind_",ind,"_data_test.rds")))
  }
}
makedata_supervised(eset, log2cpm_quant, theta)


# make thinned data
makedata_thinned <- function(eset, thinlog2_rate) {

  dir <- "/project2/gilad/joycehsiao/fucci-seq"

  library(Biobase)
  library(peco)
  pdata <- pData(eset)
  fdata <- fData(eset)
  counts <- exprs(eset)[grep("ENSG", rownames(eset)), ]

  # ---- thinning!
  library(seqgendiff)
  nsamp <- ncol(counts)
  inputmat <- counts
  thinlog2 <- rexp(nsamp, rate = thinlog2_rate)
  outmat <- thin_lib(inputmat, thinlog2 = thinlog2)

  counts_thin <- outmat$mat
  dimnames(counts_thin) <- dimnames(counts)
  libsize <- colSums(counts_thin)

  keep_samples <- which(colMeans(counts_thin==0)>.01)
  keep_genes <- which(rowMeans(counts_thin == 0 ) < .5)
  counts_thin <- counts_thin[keep_genes, keep_samples]
  libsize <- libsize[keep_samples]

  log2cpm_thin <- t(log2(1+(10^6)*t(counts_thin)/libsize))
  dimnames(log2cpm_thin) <- dimnames(counts_thin)

  log2cpm_thin_quant <- do.call(rbind,
                                lapply(1:nrow(log2cpm_thin), function(i) {
                                  qqnorm(log2cpm_thin[i,], plot.it = F )$x }) )
  dimnames(log2cpm_thin_quant) <- dimnames(counts_thin)

  saveRDS(list(counts_thin=counts_thin,
               log2cpm_thin= log2cpm_thin,
               log2cpm_thin_quant=log2cpm_thin_quant),
          file=file.path(dir, paste0("analysis/revision/output/data_training_test/data_thinlog2_",
                                     sprintf("%03d", 100*thinlog2_rate),".rds")))

  pdata_thin <- pdata[keep_samples,]

  for (i in 1:length(unique(pdata$chip_id))) {
    ind <- unique(pdata_thin$chip_id)[i]
    ii_test <- c(1:nrow(pdata_thin))[which(pdata_thin$chip_id == ind)]
    counts_test <- counts[,ii_test]
    counts_thin_test <- counts_thin[,ii_test]
    log2cpm_thin_test <- log2cpm_thin[,ii_test]
    log2cpm_thin_quant_test <- log2cpm_thin_quant[,ii_test]

    data_thin_test <- list(counts_test=counts_test,
                           counts_thin_test=counts_thin_test,
                           log2cpm_thin_test=log2cpm_thin_test,
                           log2cpm_thin_quant_test=log2cpm_thin_quant_test)

    saveRDS(data_thin_test,
            file=file.path(dir, paste0("analysis/revision/output/data_training_test/ind_",ind,"_data_test_thinlog2_",
                                       sprintf("%03d", 100*thinlog2_rate),".rds")))
  }
}

makedata_thinned(eset, thinlog2_rate = .8)
makedata_thinned(eset, thinlog2_rate = .33)

data_thinlog2_033 <- readRDS(file.path(dir,
                                       "analysis/revision/output/data_training_test/data_thinlog2_033.rds"))

data_thinlog2_080 <- readRDS(file.path(dir,
                                       "analysis/revision/output/data_training_test/data_thinlog2_080.rds"))

mean(colSums(data_thinlog2_033$counts_thin))/mean(colSums(exprs(eset)))

mean(colSums(data_thinlog2_080$counts_thin))/mean(colSums(exprs(eset)))

mean(colSums(data_thinlog2_033$counts_thin/colSums(exprs(eset))))
mean(colSums(exprs(eset))/(1/.23))
mean(colSums(data_thinlog2_033$counts_thin))
1/.23

mean(colSums(data_thinlog2_080$counts_thin/colSums(exprs(eset))))
mean(colSums(exprs(eset))/(1/.45))
mean(colSums(data_thinlog2_080$counts_thin))
1/.45


# number of genes detected in each data
top5genes <- genes
do.call(rbind, lapply(1:length(inds), function(i) {
  ind <- inds[i]
  samps <- rownames(pdata)[which(pdata$chip_id == ind)]
  exp_tmp <- data_thinlog2_080$counts_thin[,colnames(data_thinlog2_080$counts_thin) %in% samps]
  data.frame(ind = ind,
             g_detect = min(rowMeans(exp_tmp > 0)),
             g_detect_top50 = mean(rowMeans(exp_tmp[which(rownames(exp_tmp) %in% top5genes),] > 0)))
}))


do.call(rbind, lapply(1:length(inds), function(i) {
  ind <- inds[i]
  samps <- rownames(pdata)[which(pdata$chip_id == ind)]
  exp_tmp <- data_thinlog2_033$counts_thin[,colnames(data_thinlog2_033$counts_thin) %in% samps]
  data.frame(ind = ind,
             g_detect = min(rowMeans(exp_tmp > 0)))
}))

eset <- readRDS(file=file.path(dir, "data/eset-final.rds"))
counts <- exprs(eset)
counts_filt <- counts[grep("ENSG", rownames(counts)),]
do.call(rbind, lapply(1:length(inds), function(i) {
  ind <- inds[i]
  samps <- rownames(pdata)[which(pdata$chip_id == ind)]
  exp_tmp <- counts_filt[,colnames(counts_filt) %in% samps]
  data.frame(ind = ind,
             g_detect = min(rowMeans(exp_tmp > 0)))
}))


# properties of thinned data
# checking...
thing_080 <- readRDS(file.path(dir,
                               "analysis/revision/output/data_training_test/data_thinlog2_080.rds"))
thing_033 <- readRDS(file.path(dir,
                               "analysis/revision/output/data_training_test/data_thinlog2_033.rds"))

subsam <- do.call(rbind, list(data.frame(libsize=colSums(exprs(eset)), thinlog2="1"),
                              data.frame(libsize=colSums(thing_080$counts_thin), thinlog2="080"),
                              data.frame(libsize=colSums(thing_033$counts_thin), thinlog2="033")))
subsam$thinlog2 <- factor(subsam$thinlog2,
                          levels=c("033", "080", "1"),
                          labels=c("4.4", "2.2", "1 (none)"))
ggplot(subsam, aes(x=thinlog2, y=libsize, group=thinlog2, fill=thinlog2)) +
  geom_violin() +
  geom_boxplot(width=.2, col="black") +
  labs(fill="Thinning factor") +
  ylab("Sample molecule count") +
  xlab("Thinning factor")
subsam %>% group_by(thinlog2) %>% summarize(mn = mean(libsize), sd = sd(libsize))

# --- results using peco
eset <- readRDS(file=file.path(dir, "data/eset-final.rds"))



# Supp Fig S10A
# not thinned
grids = "100"
out_notthin <- do.call(rbind, lapply(unique(pData(eset)$chip_id), function(ind) {
  #  ind = unique(pData(eset)$chip_id)[3]
  out_ind <- do.call(rbind, lapply(2:50, function(g) {
    res <- readRDS(file.path(dir,
                             paste0("analysis/revision/output/ourdata_ind_peco/ourdata_peco_",
                                    ind, "_top",sprintf("%03d", g), "genes_grids_", grids, ".rds")))
    data <- readRDS(file.path(dir,
                              paste0("analysis/revision/output/data_training_test/ind_",ind,"_data_test.rds")))
    diff_time <- circ_dist(data$theta_test,
                           rotation(data$theta_test, res$fit_test$cell_times_est))
    data.frame(phase_pred_rot=rotation(data$theta_test, res$fit_test$cell_times_est),
               phase_ref=data$theta_test,
               diff_time=diff_time,
               ind = ind,
               grids = grids,
               ngenes = g,
               thinlog2 = "1",
               method = "peco")
  }))
  return(out_ind)
}))
library(ggplot2)
library(dplyr)
out_notthin %>% group_by(ind, ngenes) %>%
  summarise(diff_mean = mean(diff_time/2/pi),
            diff_se = sd(diff_time/2/pi)/sqrt(length(diff_time/2/pi))) %>%
  ggplot(., aes(x=factor(ngenes), y=diff_mean, group = ind)) +
  #  geom_vline(xintercept=seq(5, 50,5)-1, col="gray90", lty=1) +
  geom_hline(yintercept=seq(.1, .2, .01), col="gray90", lty=1) +
  geom_line(aes(col=ind), lwd=.7) +
  scale_fill_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin=diff_mean-diff_se,
                    ymax=diff_mean+diff_se, col=ind), width=.2, alpha=.5) +
  stat_summary(fun.y=mean,geom="line",lwd=.5, group=1) +
  ylim(0,.3) + geom_hline(yintercept=.25, col="red") +
  labs(color="Test data (Cell line)") +
  xlab("Number of cyclic genes used in peco prediction") +
  ylab("Prediction error (% circle)") +
  scale_x_discrete(breaks=c(2:50),
                   labels=c(rep("",3),5, rep("",4), 10, rep("",4), 15,
                            rep("",4), 20, rep("",4), 25,
                            rep("",4), 30, rep("",4), 35,
                            rep("",4), 40, rep("",4), 45, rep("",4), 50))


# not thinned up to 500
grids = "100"
out_notthin_500 <- do.call(rbind, lapply(unique(pData(eset)$chip_id), function(ind) {
  #  ind = unique(pData(eset)$chip_id)[3]
  out_ind <- do.call(rbind, lapply(c(2:50, seq(60, 500, by=10)), function(g) {
    res <- readRDS(file.path(dir,
                             paste0("analysis/revision/output/ourdata_ind_peco/ourdata_peco_",
                                    ind, "_top",sprintf("%03d", g), "genes_grids_", grids, ".rds")))
    data <- readRDS(file.path(dir,
                              paste0("analysis/revision/output/data_training_test/ind_",ind,"_data_test.rds")))
    diff_time <- circ_dist(data$theta_test,
                           rotation(data$theta_test, res$fit_test$cell_times_est))
    data.frame(phase_pred_rot=rotation(data$theta_test, res$fit_test$cell_times_est),
               phase_ref=data$theta_test,
               diff_time=diff_time,
               ind = ind,
               grids = grids,
               ngenes = g,
               thinlog2 = "1",
               method = "peco")
  }))
  return(out_ind)
}))
library(ggplot2)
library(dplyr)
out_notthin_500 %>% group_by(ind, ngenes) %>%
  summarise(diff_mean = mean(diff_time/2/pi),
            diff_se = sd(diff_time/2/pi)/sqrt(length(diff_time/2/pi))) %>%
  ggplot(., aes(x=ngenes, y=diff_mean, group = ind)) +
  #  geom_vline(xintercept=seq(5, 50,5)-1, col="gray90", lty=1) +
  geom_hline(yintercept=seq(.1, .2, .01), col="gray90", lty=1) +
  geom_line(aes(col=ind), lwd=.7) +
  scale_fill_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin=diff_mean-diff_se,
                    ymax=diff_mean+diff_se, col=ind), width=.2, alpha=.5) +
  stat_summary(fun.y=mean,geom="line",lwd=.5, group=1) +
  ylim(0,.3) + geom_hline(yintercept=.25, col="red") +
  labs(color="Test data (Cell line)") +
  xlab("Number of cyclic genes used in peco prediction") +
  ylab("Prediction error (% circle)") +
  scale_x_continuous(breaks=c(2, 50, seq(100, 500, by=50)),
                     labels=c(2, 50, seq(100, 500, by=50)))


# Supp Fig S11C
grids = "100"
thinlog2 = "033"
out_thinlog2_033 <- do.call(rbind, lapply(unique(pData(eset)$chip_id), function(ind) {
  #  ind = unique(pData(eset)$chip_id)[3]
  out_ind <- do.call(rbind, lapply(2:50, function(g) {
    res <- readRDS(file.path(dir,
                             paste0("analysis/revision/output/ourdata_ind_peco_thin/ourdata_peco_",
                                    ind, "_top",sprintf("%03d", g), "genes_thinlog2_", thinlog2,
                                    "_grids_", grids, ".rds")))
    data <- readRDS(file.path(dir,
                              paste0("analysis/revision/output/data_training_test/ind_",ind,"_data_test.rds")))
    diff_time <- circ_dist(data$theta_test,
                           rotation(data$theta_test, res$fit_test$cell_times_est))
    data.frame(phase_pred_rot=rotation(data$theta_test, res$fit_test$cell_times_est),
               phase_ref=data$theta_test,
               diff_time=diff_time,
               ind = ind,
               grids = grids,
               ngenes = g,
               thinlog2 = thinlog2)
  }))
  return(out_ind)
}))
library(ggplot2)
library(dplyr)
out_thinlog2_033 %>% group_by(ind, ngenes) %>%
  summarise(diff_mean = mean(diff_time/2/pi),
            diff_se = sd(diff_time/2/pi)/sqrt(length(diff_time/2/pi))) %>%
  ggplot(., aes(x=factor(ngenes), y=diff_mean, group = ind)) +
  #  geom_vline(xintercept=seq(5, 50,5)-1, col="gray90", lty=1) +
  geom_hline(yintercept=seq(.1, .2, .01), col="gray90", lty=1) +
  geom_line(aes(col=ind), lwd=.7) + #ggtitle("thinlog2 = .33") +
  scale_fill_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin=diff_mean-diff_se,
                    ymax=diff_mean+diff_se, col=ind), width=.2, alpha=.5) +
  stat_summary(fun.y=mean,geom="line",lwd=.5, group=1) +
  ylim(0,.3) + geom_hline(yintercept=.25, col="red") +
  labs(color="Test data (Cell line)") +
  xlab("Number of cyclic genes used in peco prediction") +
  ylab("Prediction error (% circle)") +
  scale_x_discrete(breaks=c(2:50),
                   labels=c(rep("",3),5, rep("",4), 10, rep("",4), 15,
                            rep("",4), 20, rep("",4), 25,
                            rep("",4), 30, rep("",4), 35,
                            rep("",4), 40, rep("",4), 45, rep("",4), 50))



# Supp Fig S11D
grids = "100"
thinlog2 = "080"
out_thinlog2_080 <- do.call(rbind, lapply(unique(pData(eset)$chip_id), function(ind) {
  #  ind = unique(pData(eset)$chip_id)[3]
  out_ind <- do.call(rbind, lapply(2:50, function(g) {
    res <- readRDS(file.path(dir,
                             paste0("analysis/revision/output/ourdata_ind_peco_thin/ourdata_peco_",
                                    ind, "_top",sprintf("%03d", g), "genes_thinlog2_", thinlog2,
                                    "_grids_", grids, ".rds")))
    data <- readRDS(file.path(dir,
                              paste0("analysis/revision/output/data_training_test/ind_",ind,"_data_test.rds")))
    diff_time <- circ_dist(data$theta_test,
                           rotation(data$theta_test, res$fit_test$cell_times_est))
    data.frame(phase_pred_rot=rotation(data$theta_test, res$fit_test$cell_times_est),
               phase_ref=data$theta_test,
               diff_time=diff_time,
               ind = ind,
               grids = grids,
               ngenes = g,
               thinlog2 = thinlog2)
  }))
  return(out_ind)
}))
library(ggplot2)
library(dplyr)
out_thinlog2_080 %>% group_by(ind, ngenes) %>%
  #  filter(ngenes <= 20) %>%
  summarise(diff_mean = mean(diff_time/2/pi),
            diff_se = sd(diff_time/2/pi)/sqrt(length(diff_time/2/pi))) %>%
  ggplot(., aes(x=factor(ngenes), y=diff_mean, group = ind)) +
  #  geom_vline(xintercept=seq(5, 50,5)-1, col="gray90", lty=1) +
  geom_hline(yintercept=seq(.1, .2, .01), col="gray90", lty=1) +
  geom_line(aes(col=ind), lwd=.7) + # ggtitle("thinlog2 = .80") +
  scale_fill_brewer(palette="Dark2") +
  geom_errorbar(aes(ymin=diff_mean-diff_se,
                    ymax=diff_mean+diff_se, col=ind), width=.2, alpha=.5) +
  stat_summary(fun.y=mean,geom="line",lwd=.5, group=1) +
  ylim(0,.3) + geom_hline(yintercept=.25, col="red") +
  labs(color="Test data (Cell line)") +
  xlab("Number of cyclic genes used in peco prediction") +
  ylab("Prediction error (% circle)") +
  scale_x_discrete(breaks=c(2:50),
                   labels=c(rep("",3),5, rep("",4), 10, rep("",4), 15,
                            rep("",4), 20, rep("",4), 25,
                            rep("",4), 30, rep("",4), 35,
                            rep("",4), 40, rep("",4), 45, rep("",4), 50))


mean(out_thinlog2_080$diff_time/2/pi)
mean(out_thinlog2_033$diff_time/2/pi)




# Supp Fig S10B
inds <- c("NA19098","NA18511","NA18870","NA19101","NA18855","NA19160")
ngene <- 5
eval_fit <- lapply(1:length(inds), function(i) {
  ind <- inds[i]
  fl_name <- file.path(dir,"data/results/finalizing",
                        paste0("ind_",ind,"_results_overallcyclical.top",ngene,".rds"))
  df <- readRDS(fl_name)
  return(df$fit.test)
})
names(eval_fit) <- inds


# sample size and prediction error
ns_error <- do.call(rbind, lapply(1:length(eval_fit), function(i) {
  data.frame(ind=names(eval_fit)[i],
             error=mean(eval_fit[[i]]$fit.supervised$diff_time),
             ns=ncol(eval_fit[[i]]$fit.supervised$Y_reordered),
             stringsAsFactors = F)
}))
cols <- RColorBrewer::brewer.pal(9,"Dark2")
par(mfrow=c(1,1))
plot(y=(ns_error$error/2/pi),x=ns_error$ns, pch=16,
     ylab="Prediction error",
     xlab="Sample size in test sample", axes=F,
     xlim=c(80, 220), col=cols[1:6])
axis(1); axis(2)
text(labels=ns_error$ind[1],
     y=(ns_error$error[1]/2/pi),
     x=ns_error$ns[1], pos=4, offset=1, col=cols[1])
text(labels=ns_error$ind[2],
     y=(ns_error$error[2]/2/pi),
     x=ns_error$ns[2], pos=3, offset=1, col=cols[2])
text(labels=ns_error$ind[3],
     y=(ns_error$error[3]/2/pi),
     x=ns_error$ns[3], pos=3, offset=1, col=cols[3])
text(labels=ns_error$ind[4],
     y=(ns_error$error[4]/2/pi),
     x=ns_error$ns[4], pos=1, offset=1, col=cols[4])
text(labels=ns_error$ind[5],
     y=(ns_error$error[5]/2/pi),
     x=ns_error$ns[5], pos=3, offset=1, col=cols[5])
text(labels=ns_error$ind[6],
     y=(ns_error$error[6]/2/pi),
     x=ns_error$ns[6], pos=1, offset=1, col=cols[6])
cor.test(ns_error$error, ns_error$ns, method = "spearman")
title("Supp Fig S10B")
