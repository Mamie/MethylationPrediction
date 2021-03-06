#!/bin/bash

#SBATCH --job-name=HM450centered
#SBATCH --output=log/HM450centered.%A.out
#SBATCH --error=log/HM450centered.%A.err
#SBATCH --time=0-2:00:00
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

Rscript task_centered.R
