#!/bin/bash

#SBATCH --job-name=sbatchHNSC
#SBATCH --output=log/sbatchHNSC.%A.out
#SBATCH --error=log/sbatchHNSC.%A.err
#SBATCH --time=0-1:00:00
#SBATCH -p normal
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --cores=10
#SBATCH --mem=64000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R

Rscript task_subset.R /oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/HNSC/processed
