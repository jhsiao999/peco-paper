#!/bin/bash

# Submit batch job running job_nestorowa_peco.R

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=2GB

set -e

ngenes=$1

#source activate fucci-seq

start_ngene=$1
end_ngene=$2

Rscript leng_fucci_seurat.R ${start_ngene} ${end_ngene}



