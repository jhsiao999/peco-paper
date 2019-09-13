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

dir <- "/project2/gilad/joycehsiao/fucci-seq"

# get our data -----
eset <- readRDS(file=file.path(dir,"data/eset-final.rds"))

library(Biobase)
pdata <- pData(eset)
fdata <- fData(eset)

# select endogeneous genes
counts <- exprs(eset)[grep("ENSG", rownames(eset)), ]
fdata <- fdata[grep("ENSG", rownames(eset)), ]

# use gene symbols as indicators
rownames(counts) <- fdata$name

# take cells from one individual
# select_ind <- which(pData(eset)$chip_id == "NA18511")

# run Oscope on 366 genes (across all), MeanCutLow=10 ------
Sizes <- MedianNorm(counts)
DataNorm <- GetNormalizedMat(counts, Sizes)
MV <- CalcMV(Data = counts, Sizes = Sizes, MeanCutLow = 1)
DataNorm_high_var <- DataNorm[MV$GeneToUse,]
DataInput <- NormForSine(DataNorm_high_var)

# select samples from one individual
#DataInput_oneind <- DataInput[,select_ind]
#SineRes <- OscopeSine(DataInput_oneind, parallel=F)

# select samples from all individuals
SineRes <- OscopeSine(DataInput, parallel=T)
KMRes <- OscopeKM(SineRes, maxK = NULL)

# to flag clusters...
#ToRM <- FlagCluster(SineRes, KMRes, DataInput)
#KMResUse <- KMRes[-ToRM$FlagID]
KMResUse <- KMRes
ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput,
		    NCThre = 100, parallel=T)

save(DataInput,
     KMResUse, ENIRes,
     file = file.path(dir, "analysis/revision/output/ourdata_oscope_2895genes.rda"))


# # run Oscope on 76 genes ------
# Sizes <- MedianNorm(counts)
# DataNorm <- GetNormalizedMat(counts, Sizes)
# MV <- CalcMV(Data = counts, Sizes = Sizes, MeanCutLow = 30)
# DataNorm_high_var <- DataNorm[MV$GeneToUse,]
# DataInput <- NormForSine(DataNorm_high_var)
#
# # select samples from one individual
# #DataInput_oneind <- DataInput[,select_ind]
# #SineRes <- OscopeSine(DataInput_oneind, parallel=F)
#
# # select samples from all individuals
# SineRes <- OscopeSine(DataInput, parallel=F)
# KMRes <- OscopeKM(SineRes, maxK = NULL)
#
# # to flag clusters...
# #ToRM <- FlagCluster(SineRes, KMRes, DataInput)
# #KMResUse <- KMRes[-ToRM$FlagID]
# KMResUse <- KMRes
# ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput,
# 		    NCThre = 100, parallel=T)
#
# save(DataInput,
#      KMResUse, ENIRes, file = "output/ourdata_oscope_76genes.rda")


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






