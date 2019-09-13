#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.
#
# To run this script on the command line
# sbatch

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


library(Oscope)
library(Biobase)
library(peco)

dir <- "/project2/gilad/joycehsiao/fucci-seq"

eset <- readRDS(file.path(dir,
        "data/cellcycle-genes-previous-studies/rnaseq/leng/HumanLengESC.rds"))
pdata_leng <- pData(eset)
counts_leng <- exprs(eset)

counts_leng_fucci <- counts_leng[,pdata_leng$cell_state != "H1"]
pdata_leng_fucci <- pdata_leng[pdata_leng$cell_state != "H1",]
pdata_leng_fucci$cell_state <- droplevels(pdata_leng_fucci$cell_state)

counts_leng_h1 <- counts_leng[,pdata_leng$cell_state == "H1"]


# run Oscope on FUCCI cells ------------------------------
Sizes <- MedianNorm(counts_leng_fucci)
DataNorm <- GetNormalizedMat(counts_leng_fucci, Sizes)
MV <- CalcMV(Data = as.matrix(counts_leng_fucci), Sizes = Sizes, Plot=F)
#str(MV$GeneToUse)
DataNorm_high_var <- DataNorm[MV$GeneToUse,]
DataInput <- NormForSine(DataNorm_high_var)
SineRes <- OscopeSine(DataInput, parallel=T)
KMRes <- OscopeKM(SineRes, maxK = NULL)
KMResUse <- KMRes
ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput,
		    NCThre = 100, parallel=T)
save(DataInput,
     KMResUse, ENIRes, pdata_leng_fucci,
     file = file.path(dir, "analysis/revision/output/leng_fucci_oscope_2903genes.rda"))




# Oscope example ------
#data(OscopeExampleData)
#Sizes <- MedianNorm(OscopeExampleData)
#DataNorm <- GetNormalizedMat(OscopeExampleData, Sizes)
#MV <- CalcMV(Data = OscopeExampleData, Sizes = Sizes)
#DataSubset <- DataNorm[MV$GeneToUse,]
#DataInput <- NormForSine(DataNorm)
#SineRes <- OscopeSine(DataInput, parallel=T)
#KMRes <- OscopeKM(SineRes, maxK = 10)
#ToRM <- FlagCluster(SineRes, KMRes, DataInput)
#KMResUse <- KMRes[-ToRM$FlagID]
#ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput, NCThre = 100, parallel=T)
#save(DataNorm, KMResUse, ENIRes, file = "oscope.rda")






