#!/bin/bash

# Submit batch job running npreg.R

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=broadwl
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=2GB

set -e

#source activate fucci-seq

start_ngene=$1
end_ngene=$2

Rscript leng_fucci_cyclone.R ${start_ngene} ${end_ngene}
