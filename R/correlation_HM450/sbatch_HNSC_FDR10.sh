#!/bin/bash

#SBATCH --job-name=corFDR10
#SBATCH --output=log/corFDR10.%A.out
#SBATCH --error=log/corFDR10.%A.err
#SBATCH --time=0-1:00:00
#SBATCH -p owners,normal
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=64000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R

Rscript task_FDR10.R /oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/HNSC/processed
