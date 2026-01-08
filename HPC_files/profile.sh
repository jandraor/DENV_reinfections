#!/bin/bash

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J proflik
#! Which project should be charged:
#SBATCH -A DENGEN-SL2-CPU
#SBATCH -D /home/ja850/rds/hpc-work/DENV_reinfections/  # working directory
#SBATCH --output=proflik_%A.log
#SBATCH --error=proflik_%A.err
#SBATCH -p sapphire
#SBATCH --nodes=1
#SBATCH --cpus-per-task=110
#SBATCH --mem=64G
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time=24:00:00

# Print job info
echo "Job started at: $(date)"
echo "Running on host: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK"

# run the script
Rscript cluster_profile.R "$1"

echo "========================================"
echo "Job finished at: $(date)"
echo "========================================"
