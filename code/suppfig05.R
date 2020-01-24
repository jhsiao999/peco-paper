# Supplemental Figure S5

library("cowplot")
library("dplyr")
library("edgeR")
library("ggplot2")
library("heatmap3")
library("reshape2")
library("SingleCellExperiment")
library(heatmap3)
library(RColorBrewer)

#' @title Perform Principal Components Analysis (PCA).
#
#' @param x gene-by-sample matrix
#' @param retx center, scale - see ?prcomp
#
#' @return a list with the following elements:
#   PCs - sample-by-PC matrix of principal components
#   explained - proportion of variance explained by each PC
#
#' @references Zhang et al. 2009 (http://www.ncbi.nlm.nih.gov/pubmed/19763933)
#'
#' @export
run_pca <- function(x, retx = TRUE, center = TRUE, scale = TRUE) {
  library("testit")

  pca <- prcomp(t(x), retx = retx, center = center, scale. = scale)
  variances <- pca$sdev^2
  explained <- variances / sum(variances)
  assert("Variance explained is calculated correctly.",
         explained[1:2] - summary(pca)$importance[2, 1:2] < 0.0001)
  return(list(PCs = pca$x, explained = explained))
}

#' @title Plot PCA results.
#
#' @param x numeric matrix of PCs
#' @param pcx PC to plot on x-axis (default: 1)
#' @param pcy PC to plot on y-axis (default: 2)
#' @param explained numeric vector of fractions of variance explained by each PC
#' @param metadata data frame or matrix that contains the metadata used to annotate
#             the plot
#' @param color, shape, size: column name of metadata used to pass column to ggplot
#                       aesthetic
#' @param factors character vector which contains the column names of metadata that
#            need to be explicitly converted to a factor
#' @param ...  Additional arguments passed to geom_point
#'
#' @export
plot_pca <- function(x, pcx = 1, pcy = 2, explained = NULL, metadata = NULL,
                     color = NULL, shape = NULL, factors = NULL,
                     ...) {
  library("ggplot2")
  library("testit")

  # Prepare data frame to pass to ggplot
  if (!is.null(metadata)) {
    assert("PC and metadata have same number of rows.",
           nrow(x) == nrow(metadata))
    plot_data <- cbind(x, metadata)
    plot_data <- as.data.frame(plot_data)
    # Convert numeric factors to class "factor"
    for (f in factors) {
      plot_data[, f] <- as.factor(plot_data[, f])
    }
  } else {
    plot_data <- as.data.frame(x)
  }
  # Prepare axis labels
  if (!is.null(explained)) {
    assert("Number of PCs differs between x and explained.",
           length(explained) == ncol(x))
    xaxis <- sprintf("PC%d (%.2f%%)", pcx, round(explained[pcx] * 100, 2))
    yaxis <- sprintf("PC%d (%.2f%%)", pcy, round(explained[pcy] * 100, 2))
  } else {
    xaxis <- paste0("PC", pcx)
    yaxis <- paste0("PC", pcy)
  }
  # Plot
  p <- ggplot(plot_data, aes_string(x = paste0("PC", pcx),
                                    y = paste0("PC", pcy),
                                    color = color,
                                    shape = shape)) +
    geom_point(...) +
    labs(x = xaxis, y = yaxis)
  p
}

## a function of pca vs technical factors
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}


## PCA
sce_filtered <- readRDS("data/sce-filtered.rds")
counts = assay(sce_filtered)
pdata <- data.frame(colData(sce_filtered))
pdata$experiment <- as.factor(pdata$experiment)

# Compute log2 CPM based on the library size before filtering.
log2cpm <- edgeR::cpm(counts, log = TRUE)

pca_log2cpm <- run_pca(log2cpm)


# PCA variation -------------------------------------------------------------------
# selection of technical factor

covariates <- pdata %>% dplyr::select(experiment, well, chip_id,
                                            concentration, raw:unmapped,
                                            starts_with("detect"),  molecules)
# look at the first 6 PCs
pca_log2cpm <- prcomp(t(log2cpm), scale. = TRUE, center = TRUE)
pcs <- pca_log2cpm$x[, 1:6]


# R-square between PCs and the covariates
r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, cov], pcs[, pc])
  }
}

# plot heatmap
heatmap3(r2, cexRow=1, cexCol=1, margins=c(4,15), scale="none",
         ylab="", main = "", symm=F,
         Colv=NA, showColDendro = F,
         labRow = c("C1: C1 batch",
                    "Well: C1 capture site (well) ID",
                    "Individual",
                    "cDNA: cDNA concentration",
                    "Reads raw: Raw reads",
                    "Reads UMI: Reads with valid UMI",
                    "Reads mapped: Reads with valid UMI mapped to genome",
                    "Reads unmapped: Reads with valid UMI not mapped to genome",
                    "ERCC: ERCC gene proportion detected",
                    "ENSG: Endogeneous gene proportion detected",
                    "Molecules: Total sample molecule count"
         ))

# r-squared between PC1 and individual
f1 <- lm(pcs[,1] ~ covariates$chip_id)
f2 <- lm(pcs[,1] ~ covariates$molecules)
summary(f1)

summary(f1)$adj.r.squared
summary(f2)$adj.r.squared


# PC proportions
100*((pca_log2cpm$sdev^2)/sum(pca_log2cpm$sdev^2))[1:6]


# Correlation between Technical factors
cor_tech <- cor(as.matrix(covariates[,4:11]),use="pairwise.complete.obs")
heatmap3(cor_tech, symm = TRUE, margins=c(7,15),
         col=brewer.pal (9, "Blues" ), cexRow=1, cexCol=1, scale="none",
         labRow = c("cDNA: cDNA concentration",
                    "Reads raw: Raw reads",
                    "Reads UMI: Reads with valid UMI",
                    "Reads mapped: Reads with valid UMI mapped to genome",
                    "Reads unmapped: Reads with valid UMI not mapped to genome",
                    "ERCC: ERCC gene proportion detected",
                    "ENSG: Endogeneous gene proportion detected",
                    "Molecules: Total sample molecule count"),
         labCol = c("cDNA",
                    "Reads raw",
                    "Reads UMI",
                    "Reads mapped",
                    "Reads unmapped",
                    "ERCC",
                    "ENSG",
                    "Molecules"))

summary(lm(covariates$detect_hs ~ covariates$experiment))

summary(lm(covariates$detect_hs ~ covariates$mapped))

summary(lm(covariates$detect_hs ~ covariates$chip_id))

summary(lm(covariates$detect_hs ~ covariates$molecules))

summary(lm(covariates$detect_hs ~ covariates$concentration))

