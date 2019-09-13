#!/bin/bash

# Submit batch job running job_nestorowa_peco.R

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joyce.hsiao1@gmail.com
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2GB

set -e

ind=$1
ngenes=$2
ncores=$3
grids=$4

#source activate fucci-seq

Rscript ourdata_ind_peco.R ${ind} ${ngenes} ${ncores} ${grids}
