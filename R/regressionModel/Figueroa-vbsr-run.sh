#!/bin/bash

#SBATCH --job-name=Figueroavbsr
#SBATCH --output=log/Figueroavbsr.%A_%a.out
#SBATCH --error=log/Figueroavbsr.%A_%a.err
#SBATCH --time=0-01:00:00
#SBATCH -p normal
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=32000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R

scriptDir="/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/R/regressionModel/"
resultDir="/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/results/"
script=${scriptDir}/Figueroa-vbsr-run.R
lookupFile=${scriptDir}/Figueroa-vbsr-list.tsv

taskID=${SLURM_ARRAY_TASK_ID}

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2

model=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $1 }')
methylation=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $2 }')
rnaseq=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $3 }')
resultDir=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $4 }')

echo Rscript $script $methylation $rnaseq $resultDir >&2
Rscript $script $methylation $rnaseq $resultDir >&2

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2
