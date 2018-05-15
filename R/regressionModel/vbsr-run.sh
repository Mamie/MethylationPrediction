#!/bin/bash

#SBATCH --job-name=vbsr
#SBATCH --output=log/vbsr.%A_%a.out
#SBATCH --error=log/vbsr.%A_%a.err
#SBATCH --time=0-01:00:00
#SBATCH -p andrewg
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=15000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R

scriptDir="/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/R/regressionModel/"
resultDir="/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/results/"
script=${scriptDir}/vbsr-run.R
lookupFile=${scriptDir}/vbsr-list.tsv

taskID=${SLURM_ARRAY_TASK_ID}

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2

cancer=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $1 }')
resultDir=${resultDir}/${cancer}/
methylation=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $2 }')
rnaseq=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $3 }')
convertToM=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $4 }')
CGIprobes=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $5 }')
samplecode=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $6 }')

echo Rscript $script $resultDir $methylation $rnaseq $convertToM $CGIprobes $samplecode >&2
Rscript $script $resultDir $methylation $rnaseq $convertToM $CGIprobes $samplecode >&2

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2
