#!/bin/bash

# Submit batch job running npreg.R

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joyce.hsiao1@gmail.com
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2GB

set -e

#source activate fucci-seq

seed=$1

Rscript ourdata_recat.R ${seed}

