#!/bin/bash

#SBATCH --job-name=mRNAseqHNSC
#SBATCH --output=mRNAseqHNSC.%A.out
#SBATCH --error=mRNAseqHNSC.%A.err
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
Rscript preprocessing_mRNAseq.R /oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/HNSC/mRNAseq/HNSC.uncv2.mRNAseq_RSEM_normalized_log2.txt /oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/HNSC/processed

