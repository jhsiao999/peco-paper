#!/bin/bash

# Submit batch job running npreg.R

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=mstephens
#SBATCH --account=pi-mstephens
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2GB

set -e

#source activate fucci-seq

start_ngene=$1
end_ngene=$2

Rscript ourdata_cyclone.R ${start_ngene} ${end_ngene}
