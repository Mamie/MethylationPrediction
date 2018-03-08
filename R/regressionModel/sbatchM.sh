#!/bin/bash

#SBATCH --job-name=sbatchM
#SBATCH --output=log/sbatchM.%A.out
#SBATCH --error=log/sbatchM.%A.err
#SBATCH --time=0-16:00:00
#SBATCH -p bigmem
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --cores=10
#SBATCH --mem=500000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R

Rscript task.R
