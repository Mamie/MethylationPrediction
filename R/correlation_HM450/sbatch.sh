#!/bin/bash

#SBATCH --job-name=computeCorrelation
#SBATCH --output=log/computeCorrelation.%A.out
#SBATCH --error=log/computeCorrelation.%A.err
#SBATCH --time=2-0:00:00
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

Rscript task.R
