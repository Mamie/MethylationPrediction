#!/bin/bash

#SBATCH --job-name=HM450HNSC
#SBATCH --output=log/HM450HNSC.%A.out
#SBATCH --error=log/HM450HNSC.%A.err
#SBATCH --time=0-1:00:00
#SBATCH --qos=normal
#SBATCH -p andrewg
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=64000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R
Rscript subsetHNSC.R /oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/HNSC/processed/HM450.RData 100

