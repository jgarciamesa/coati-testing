#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 16           # number of cores
#SBATCH --mem=32G       # memory limit
#SBATCH -t 6-00:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public       # QOS
#SBATCH -o logs/slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e logs/slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --export=NONE   # Purge the job-submitting shell environment

module purge
module add r-4.3.0-aocc-3.1.0

make -k -j ${SLURM_CPUS_PER_TASK} $1
