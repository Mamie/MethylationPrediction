#!/bin/bash

#SBATCH --job-name=preprocessHM270
#SBATCH --output=preprocessHM27.%A.out
#SBATCH --error=preprocessHM27.%A.err
#SBATCH --time=2-0:00:00
#SBATCH --qos=normal
#SBATCH -p owners,normal
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=64000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R
Rscript preprocessing_HM27K.R
