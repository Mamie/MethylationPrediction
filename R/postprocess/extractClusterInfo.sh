#!/bin/bash

#SBATCH --job-name=writeCSVs
#SBATCH --error=log/writeCSVs.%A_%a.err
#SBATCH --out=log/writeCSVs.%A_%a.out
#SBATCH --time=0-00:10:00
#SBATCH -p andrewg
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=END,FAIL
#################
# Usage: $ sbatch --array=1-2,4%1 sbatch.sh

set -beEu -o pipefail

#module load R

scriptDir="/Users/wangmeng/Documents/Research/Gentles/repos/MethylationPrediction/R/postprocess" #"/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/R/postprocess/"
script=${scriptDir}/extractClusterInfo.R
lookupFile=${scriptDir}/extractClusterInfo-list.tsv

for taskID in 1 2 3 4 5 6 7 8  #${SLURM_ARRAY_TASK_ID}
do
    SLURM_JOBID=1
    SLURM_ARRAY_TASK_ID=$taskID

    echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2

    cancer=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $1 }')
    data=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $2 }')
    methylationProbes=$(cat $lookupFile | awk -F"\t" -v taskID=$taskID '(NR == taskID) { print $3 }')

    echo $cancer 
    echo Rscript $script $data $methylationProbes 
    Rscript $script $data $methylationProbes

    echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}" >&2
done
