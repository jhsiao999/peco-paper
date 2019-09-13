#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}

ngenes <- as.numeric(args[1])


library(Biobase)
library(peco)



# getting Leng data...
dir <- "/project2/gilad/joycehsiao/fucci-seq"

eset <- readRDS(file.path(dir,
            "data/cellcycle-genes-previous-studies/rnaseq/leng/HumanLengESC.rds"))
pdata_leng <- pData(eset)
counts_leng <- exprs(eset)

counts_leng_fucci <- counts_leng[,pdata_leng$cell_state != "H1"]
counts_leng_h1 <- counts_leng[,pdata_leng$cell_state == "H1"]

pdata_leng_fucci <- pdata_leng[pdata_leng$cell_state != "H1",]
pdata_leng_fucci$cell_state <- droplevels(pdata_leng_fucci$cell_state)
# counts_f ucci
log2cpm_leng_fucci <- t(log2(10^6 * t(counts_leng_fucci + 1)/colSums(counts_leng_fucci)))
keep_genes <- which(rowMeans(log2cpm_leng_fucci) > 1)
log2cpm_leng_fucci_filt <- log2cpm_leng_fucci[keep_genes,]

# quantile normalize data
log2cpm_leng_fucci_quant <- do.call(rbind, lapply(1:nrow(log2cpm_leng_fucci_filt), function(g) {
  yy <- log2cpm_leng_fucci_filt[g,]
  is.zero <- which(yy == 0)
  qq.map <- qqnorm(yy, plot.it = F)
  yy.qq <- qq.map$x
  yy.qq[is.zero] <- sample(qq.map$x[is.zero])
  return(yy.qq)
}) )
rownames(log2cpm_leng_fucci_quant) <- rownames(log2cpm_leng_fucci_filt)
colnames(log2cpm_leng_fucci_quant) <- colnames(log2cpm_leng_fucci_filt)



# getting our data! -----------------------
# expression set of top 100 cyclic genes
eset <- readRDS(file=file.path(dir,"data/eset-final.rds"))
pdata <- pData(eset)
fdata <- fData(eset)

# select endogeneous genes
counts <- exprs(eset)[grep("ENSG", rownames(eset)), ]
fdata <- fdata[grep("ENSG", rownames(eset)), ]

# use gene symbols as indicators
rownames(counts) <- fdata$name
log2cpm.quant <- readRDS(file = file.path(dir,
                      "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))
log2cpm.quant <- log2cpm.quant[grep("ENSG", rownames(log2cpm.quant)), ]

# All 11,040 genes included in our dataset
# ordered by cyclic trend (strong to weak)
fits_all <- readRDS("/project2/gilad/joycehsiao/fucci-seq/output/npreg-trendfilter-quantile.Rmd/fit.quant.rds")
pve <- sapply(fits_all, "[[", 3)
cyclic_gene <- data.frame(ensg=names(fits_all),
                          pve=sapply(fits_all, "[[", 3),
                          stringsAsFactors = F)
cyclic_gene <- cyclic_gene[order(cyclic_gene$pve, decreasing = T),]
genes_cyclic_list <- cyclic_gene

# learning cyclic functions
Y_train <- log2cpm.quant[match(rownames(cyclic_gene),rownames(log2cpm.quant)),]
theta_train <- pdata$theta
names(theta_train) <- rownames(pdata)
theta_train <- theta_train[order(theta_train)]
Y_train <- Y_train[,match(names(theta_train), colnames(Y_train))]
#all.equal(names(theta_train), colnames(Y_train))
rownames(Y_train) <- fdata$name[match(rownames(Y_train), rownames(fdata))]


# applying peco -----------------------
# training
topX_genes <- rownames(Y_train)[1:ngenes]
genes_training <- which(topX_genes %in% rownames(log2cpm_leng_fucci_quant))
Y_train_sub <- Y_train[genes_training,]
fit_train <- cycle_npreg_insample(Y = Y_train_sub,
                                  theta = theta_train,
                                  polyorder=2,
                                  ncores=1,
                                  method.trend="trendfilter")
#save(Y_train, fit_train, file = "output/benchmarking.Rmd/peco_predtop10.rda")

# predicting
genes_predict <- which(rownames(log2cpm_leng_fucci_quant) %in% topX_genes)
Y_predict <- log2cpm_leng_fucci_quant[genes_predict,]
fit_predict <- cycle_npreg_outsample(Y_test=Y_predict,
                                     sigma_est=fit_train$sigma_est,
                                     funs_est=fit_train$funs_est,
                                     method.grid = "uniform",
                                     method.trend="trendfilter",
                                     polyorder=2, ncores=2, get_trend_estimates = T)

out <- list(fit_train=fit_train,
            fit_predict=fit_predict)
#}
#names(fits) <- paste0("fold.", 1:length(fold_indices))
saveRDS(out,
        file=file.path(dir,
        paste0("analysis/revision/output/leng_fucci_peco/top_",ngenes,".rds")))





