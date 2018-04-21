#!/bin/bash

#SBATCH --job-name=vbsrHM450centered
#SBATCH --output=log/vbsrHM450centered.%A.out
#SBATCH --error=log/vbsrHM450centered.%A.err
#SBATCH --time=0-2:00:00
#SBATCH -p andrewg
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --cores=10
#SBATCH --mem=64000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R

Rscript task_vbsr_HM450M.R
