#!/bin/bash

#SBATCH --job-name=HM450HNSC
#SBATCH --output=HM450HNSC.%A.out
#SBATCH --error=HM450HNSC.%A.err
#SBATCH --time=0-1:00:00
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
Rscript preprocessing_HM450K.R /oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/HNSC/HM450/HNSC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt /oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/HNSC/processed

