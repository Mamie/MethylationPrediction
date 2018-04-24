#!/bin/bash

#SBATCH --job-name=HM450
#SBATCH --error=log/HM450.%A_%a.err
#SBATCH --out=log/HM450.%A_%a.out
#SBATCH --time=0-01:00:00
#SBATCH -p andrewg
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R

scriptDir="/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/R/preprocess/"
script=${scriptDir}/HM450K.R
lookupFile=${scriptDir}/HM450K-list.tsv

taskID=${SLURM_ARRAY_TASK_ID}

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2

cancer=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $1 }')
methylation=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $2 }')
outdir=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $3 }')

echo $cancer >&2
echo Rscript $script $methylation $outdir >&2
Rscript $script $methylation $outdir >&2

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2
