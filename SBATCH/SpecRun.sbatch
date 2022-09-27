#!/bin/bash -l

#SBATCH --partition=prod
#SBATCH --nodes 1 
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 12:00:00
#SBATCH --mem=16000
#SBATCH --array=1-60

for i in {1..60}
do
   srun ~/SPEC/xspec Input_scan${SLURM_ARRAY_TASK_ID}_${i}.sp
done 

