#!/bin/bash

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J prof2lambda1
#! Which project should be charged:
#SBATCH -A DENGEN-SL2-CPU
#SBATCH -D /home/ja850/rds/hpc-work/DENV_reinfections/  # working directory
#SBATCH --output=prof2lambda1_%A.log
#SBATCH --error=prof2lambda1_%A.err
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=64G
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time=12:00:00

# Print job info
echo "Running on host: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK"

# run the script
Rscript cluster_profile2.R
