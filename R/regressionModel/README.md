`vbsr-list.tsv` contains a list of arguments to models (8 models in total). Each line corresponds to a task in the job array. Information on job array can be found [here](https://slurm.schedmd.com/job_array.html) 

`vbsr-run.sh` is batch script file for the job array. To submit the jobs, use command 

```
sbatch --array=1-8 vbsr-run.sh   # submit task 1 - 8
sbatch --array=1,3-4 vbsr-run.sh # submit task 1, 3, 4
``` 
