# Supplemental Figure S14
# Comparison of FUCCI phase with phase assignment by Seurat and Cyclone

library(SingleCellExperiment)
library(peco)
sce <- readRDS("data/sce-final.rds")
pdata <- data.frame(colData(sce))

# derive cell cycle phase
# use PCA to rotate the axes
pca <- prcomp(cbind(pdata$rfp.median.log10sum.adjust,
                    pdata$gfp.median.log10sum.adjust))
theta <- coord2rad(pca$x)
theta_final <- shift_origin(as.numeric(theta), 3*pi/4)
names(theta_final) = rownames(pdata)

# get PAM-based clusters
pam_res <- pam(cbind(pdata$rfp.median.log10sum.adjust,
                     pdata$gfp.median.log10sum.adjust), k=3)
clust <- data.frame(clust=pam_res$clustering,
                    sample_id=rownames(pdata))
plot(theta_final, clust$clust)

# find boundaries
b1 <- mean(max(range(theta_final[clust$clust==2])), min(range(theta_final[clust$clust==3])))
b2 <- mean(max(range(theta_final[clust$clust==3])), max(range(theta_final[clust$clust==1])))
b3 <- mean(min(range(theta_final[clust$clust==2])), min(range(theta_final[clust$clust==1])))
abline(v=c(b1,b2,b3), lty=2)


# Supp Fig 14A
# Seurat and PAM-based classification obtained from FUCCI scores
out_seurat <- readRDS("data/ourdata_phase_seurat.rds")
phase_seurat_rot <- rotation(pdata$theta,
                             out_seurat$phase_peco[match(rownames(pdata),
                                                         rownames(out_seurat))])
out_seurat$Phase <- factor(out_seurat$Phase,
                           levels = c("G1", "S", "G2M"))
out_seurat_ordered <- out_seurat[match(rownames(pdata), rownames(out_seurat)),]

all.equal(as.character(clust$sample_id), rownames(out_seurat_ordered))
table(out_seurat_ordered$Phase, clust$clust)

# Supp Fig S14B
# Cyclone and PAM-based classification obtained from FUCCI scores
out_cyclone <- readRDS("data/ourdata_phase_cyclone.rds")
phase_cyclone_rot <- rotation(pdata$theta,
                              out_cyclone$phase_cyclone_rotated[match(rownames(pdata),
                                                                      rownames(out_cyclone))])
out_cyclone$phases <- factor(out_cyclone$phases,
                             levels = c("G1", "S", "G2M"))
out_cyclone_ordered_phases <- out_cyclone$phases[match(rownames(pdata), names(out_cyclone$phases))]

all.equal(as.character(clust$sample_id), names(out_cyclone_ordered_phases))
table(out_cyclone_ordered_phases, clust$clust)


# Supp Fig S14C
# Seurat-based classification and FUCCI phase
# seurat discrete classes vs FUCCI phase
out_seurat2 = out_seurat$Phase[match(names(theta_final), rownames(out_seurat))]
all.equal(rownames(out_seurat2), names(theta_final))

phases = c("G1", "S", "G2M")
cols = as.character(yarrr::piratepal("espresso")[c(3,2,1)])
par(mfrow = c(1,3))
for (i in seq_along(phases)) {
  hist(theta_final[out_seurat2==phases[i]],
       breaks = 80, ylim = c(0,40), axes = F, ann = F, lty = "blank",
       col = cols[i])
  title(main = paste(sum(out_seurat2==phases[i]), "cells"))
  axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
                labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                         expression(2*pi)))
  abline(v=c(b1,b2,b3), lty=2)
}


# Supp Fig S14D
# Cyclone-based classification and FUCCI phase
out_cyclone2 = out_cyclone$phases[match(names(theta_final),
                                        names(out_cyclone$phases))]
all.equal(names(out_cyclone2), names(theta_final))

phases = c("G1", "S", "G2M")
cols = as.character(yarrr::piratepal("espresso")[c(3,2,1)])
par(mfrow = c(1,3))
for (i in seq_along(phases)) {
  hist(theta_final[out_cyclone2==phases[i]],
       breaks = 80, ylim = c(0,40), axes = F, ann = F, lty = "blank",
       col = cols[i])
  title(main = paste(sum(out_cyclone2==phases[i]), "cells"))
  axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
                labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                         expression(2*pi)))
  abline(v=c(b1,b2,b3), lty=2)
}




