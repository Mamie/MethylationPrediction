#!/bin/bash

#SBATCH --job-name=filterFDR
#SBATCH --out=log/filterFDR.%A_%a.out
#SBATCH --error=log/filterFDR.%A_%a.err
#SBATCH --time=0-01:00:00
#SBATCH -p bigmem
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=100000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R

scriptDir="/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/R/correlation_HM450/"
script=${scriptDir}/filterFDR.R
lookupFile=${scriptDir}/filterFDR-list.tsv

taskID=${SLURM_ARRAY_TASK_ID}

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2

cancer=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $1 }')
methylation=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $2 }')
rnaseq=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $3 }')
resultsPath=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $4 }')
fdr=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $5 }')

echo $cancer >&2
echo Rscript $script $methylation $rnaseq $resultsPath $fdr >&2
Rscript $script $methylation $rnaseq $resultsPath $fdr >&2

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2
