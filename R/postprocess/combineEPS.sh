#!/bin/bash

#SBATCH --job-name=combineEPS
#SBATCH --error=log/combineEPS.%A_%a.err
#SBATCH --out=log/combineEPS.%A_%a.out
#SBATCH --time=0-00:10:00
#SBATCH -p normal
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

module load R

scriptDir="/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/R/postprocess/"
lookupFile=${scriptDir}/combineEPS-list.tsv

taskID=${SLURM_ARRAY_TASK_ID}

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2

cancer=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $1 }')
folder=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $2 }')
outfile=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $3 }')

echo $cancer >&2
echo bash /oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/R/postprocess/gsCombineEPS.sh $folder $outfile >&2
bash /oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/R/postprocess/gsCombineEPS.sh $folder $outfile >&2

echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2
