#!/bin/bash

#SBATCH --job-name=randIdx
#SBATCH --output=log/randIdx.%A.out
#SBATCH --error=log/randIdx.%A.err
#SBATCH --time=0-12:00:00
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

Rscript compareClustering.R
