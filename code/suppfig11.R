# Supplemental Figure S11

dir <- "/project2/gilad/joycehsiao/fucci-seq"
eset <- readRDS(file=file.path(dir, "data/eset-final.rds"))

# Supp Fig S11A
# Cyclone and PAM-based classification obtained from FUCCI scores
out_cyclone <- readRDS(file.path(dir, "analysis/revision/output/ourdata_phase_cyclone.rds"))
phase_cyclone_rot <- rotation(pdata$theta,
                              out_cyclone$phase_cyclone_rotated[match(rownames(pdata),
                                                                      rownames(out_cyclone))])
out_cyclone$phases <- factor(out_cyclone$phases,
                             levels = c("G1", "S", "G2M"))
out_cyclone_ordered_phases <- out_cyclone$phases[match(rownames(pData(df)), names(out_cyclone$phases))]

all.equal(as.character(clust$sample_id), names(out_cyclone_ordered_phases))
table(out_cyclone_ordered_phases, clust$clust)

par(mfrow=c(1,1))
all.equal(names(out_cyclone_ordered_phases), rownames(pData(df)) )
plot(x=shift_origin(as.numeric(theta), 3*pi/4),
     y=pData(df)$gfp.median.log10sum.adjust, col=cols[out_cyclone_ordered_phases],
     ylim=c(-1.5, 1.5), cex=.5, pch=16,
     xlab="FUCCI phase", ylab="FUCCI intensities",
     main="Cyclone phases and FUCCI intensities", axes=F)
axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
              labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                       expression(2*pi)))
abline(h=0, lty=1, col="black", lwd=.7)
points(x=shift_origin(as.numeric(theta), 3*pi/4),
       y=pData(df)$rfp.median.log10sum.adjust, col=cols[out_cyclone_ordered_phases],
       ylim=c(-1.5, 1.5), cex=.5, pch=2)




# Supp Fig 11B
# Seurat and PAM-based classification obtained from FUCCI scores
out_seurat <- readRDS(file.path(dir, "analysis/revision/output/ourdata_phase_seurat.rds"))
phase_seurat_rot <- rotation(pData(eset)$theta,
                             out_seurat$phase_peco[match(rownames(pData(eset)),
                                                         rownames(out_seurat))])
out_seurat$Phase <- factor(out_seurat$Phase,
                           levels = c("G1", "S", "G2M"))
out_seurat_ordered <- out_seurat[match(rownames(pData(eset)), rownames(out_seurat)),]

all.equal(as.character(clust$sample_id), rownames(out_seurat_ordered))
table(out_seurat_ordered$Phase, clust$clust)

par(mfrow=c(1,1))
all.equal(rownames(out_seurat_ordered), rownames(pData(eset)) )
plot(x=shift_origin(as.numeric(theta), 3*pi/4),
     y=pData(eset)$gfp.median.log10sum.adjust, col=cols[out_seurat_ordered$Phase],
     ylim=c(-1.5, 1.5), cex=.5, pch=16,
     xlab="FUCCI phase", ylab="FUCCI intensities",
     main="Seurat phases FUCCI intensities", axes=F)
axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
              labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                       expression(2*pi)))
abline(h=0, lty=1, col="black", lwd=.7)
points(x=shift_origin(as.numeric(theta), 3*pi/4),
       y=pData(df)$rfp.median.log10sum.adjust, col=cols[out_seurat_ordered$Phase],
       ylim=c(-1.5, 1.5), cex=.5, pch=2)


# Supp Fig S11C, S11D
# Oscope and FUCCI score
load(file = file.path(dir, "analysis/revision/output/ourdata_oscope_366genes.rda"))
samples_ordered <- colnames(DataInput)[ENIRes[["cluster2"]]]
phase_oscope <- seq(0, 2*pi, length.out= length(ENIRes[[2]]))
names(phase_oscope) <- samples_ordered
ii_match <- match(samples_ordered, rownames(pData(eset)) )

theta_final <- shift_origin(pData(eset)$theta, origin=3*pi/4)
names(theta_final) <- rownames(pData(eset))

map_oscope <- data.frame(theta_final=theta_final,
                         phase_oscope=phase_oscope[match(names(theta_final), names(phase_oscope))],
                         sample_id=names(theta_final))
map_oscope$phase_oscope_rot <- rotation(ref_var=map_oscope$theta_final,
                                        shift_var=map_oscope$phase_oscope)

# Supp Fig S11C
par(mfrow=c(1,1))
plot(x=map_oscope$theta_final,
     y=map_oscope$phase_oscope_rot%%(2*pi),
     xlab="FUCCI phase", ylab="Oscope predicted phase",
     pch=16, cex=.6, axes=F, col="gray50")
abline(0,1, lty=1)
axis(1, at=c(0,pi/2, pi, 3*pi/2, 2*pi),
     labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
              expression(2*pi)))
axis(2, at=c(0,pi/2, pi, 3*pi/2, 2*pi),
     labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
              expression(2*pi)))
abline(a=pi/2, b=1, lty=2, col="black", lwd=.7)
abline(a=-pi/2, b=1, lty=2, col="black", lwd=.7)
title("Supp Fig S11C")


# Supp Fig S11D
par(mfrow=c(1,1))
plot(x=phase_oscope,
     y=pData(eset)$gfp.median.log10sum.adjust[ii_match], col="forestgreen",
     ylim=c(-1.5, 1.5), cex=.5, pch=16,
     xlab="Oscope phase", ylab="FUCCI intensities",
     main="Oscope phase and FUCCI intensities", axes=F)
axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
              labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                       expression(2*pi)))
abline(h=0, lty=1, col="black", lwd=.7)
points(x=phase_oscope,
       y=pData(eset)$rfp.median.log10sum.adjust[ii_match], col="firebrick",
       ylim=c(-1.5, 1.5), cex=.5, pch=16)
gfp <- pData(eset)$gfp.median.log10sum.adjust[ii_match]
gfp_ordered <- gfp[order(phase_oscope)]
rfp <- pData(eset)$rfp.median.log10sum.adjust[ii_match]
rfp_ordered <- rfp[order(phase_oscope)]
fit_gfp_oscope <- fit_trendfilter_generic(gfp_ordered)
fit_rfp_oscope <- fit_trendfilter_generic(rfp_ordered)
title("Supp Fig S11D")




# Supp Fig S11E, S11F
# reCAT phase and FUCCI score
load(file.path(dir, paste0("analysis/revision/output/ourdata_recat.rda")))
samples_ordered <- rownames(test_exp)[ordIndex]
phase_recat <- seq(0, 2*pi, length.out= length(ordIndex))
names(phase_recat) <- samples_ordered
ii_match <- match(samples_ordered, rownames(pData(eset)) )

theta_final <- shift_origin(pData(eset)$theta, origin=3*pi/4)
names(theta_final) <- rownames(pData(eset))

map_recat <- data.frame(theta_final=theta_final,
                        phase_recat=phase_recat[match(names(theta_final), names(phase_recat))],
                        sample_id=names(theta_final))
map_recat$phase_recat_rot <- rotation(ref_var=map_recat$theta_final,
                                      shift_var=map_recat$phase_recat)

par(mfrow=c(1,1))
plot(x=map_recat$theta_final,
     y=map_recat$phase_recat_rot %% (2*pi),
     xlab="FUCCI phase", ylab="reCAT predicted phase",
     pch=16, cex=.6, axes=F, col="gray50")
abline(0,1, lty=1)
abline(a=pi/2, b=1, lty=2, col="black", lwd=.7)
abline(a=-pi/2, b=1, lty=2, col="black", lwd=.7)
axis(1, at=c(0,pi/2, pi, 3*pi/2, 2*pi),
     labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
              expression(2*pi)))
axis(2, at=c(0,pi/2, pi, 3*pi/2, 2*pi),
     labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
              expression(2*pi)))
title("Supp Fig S11E")


par(mfrow=c(1,1))
plot(x=phase_recat,
     y=pData(df)$gfp.median.log10sum.adjust[ii_match], col="forestgreen",
     ylim=c(-1.5, 1.5), cex=.5, pch=16,
     xlab="reCAT phase", ylab="FUCCI intensities",
     main="reCAT phase and FUCCI intensities", axes=F)
axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
              labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                       expression(2*pi)))
abline(h=0, lty=1, col="black", lwd=.7)
points(x=phase_recat,
       y=pData(df)$rfp.median.log10sum.adjust[ii_match], col="firebrick",
       ylim=c(-1.5, 1.5), cex=.5, pch=16)

gfp <- pData(df)$gfp.median.log10sum.adjust[ii_match]
gfp_ordered <- gfp[order(phase_recat)]
rfp <- pData(df)$rfp.median.log10sum.adjust[ii_match]
rfp_ordered <- rfp[order(phase_recat)]
fit_gfp_recat <- fit_trendfilter_generic(gfp_ordered)
fit_rfp_recat <- fit_trendfilter_generic(rfp_ordered)
title("Supp Fig S11F")

