#!/bin/bash

#SBATCH --job-name=sbatchHM450
#SBATCH --output=log/sbatchHM450.%A.out
#SBATCH --error=log/sbatchHM450.%A.err
#SBATCH --time=0-16:00:00
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

Rscript task.R
