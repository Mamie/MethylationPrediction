#!/bin/bash

#SBATCH --job-name=HM27subset
#SBATCH --output=log/HM27subset.%A.out
#SBATCH --error=log/HM27subset.%A.err
#SBATCH --time=0-1:00:00
#SBATCH -p normal
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=64000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R

Rscript task_subset_HM27_centered.R
