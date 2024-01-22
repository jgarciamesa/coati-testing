#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 1            # number of cores
#SBATCH -t 0-01:00:00   # time in d-hh:mm:ss
#SBATCH --mem=16G       # memory limit
#SBATCH -p htc          # partition
#SBATCH -q public       # QOS
#SBATCH -o logs/slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e logs/slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=FAIL # Send an e-mail when a job fails
#SBATCH --export=NONE   # Purge the job-submitting shell environment

module purge
module add r-4.3.0-aocc-3.1.0

stem=$1

make -f align.mk "benchmark_fasta_aligned/coati-tri-mg/${stem}.coati-tri-mg.fasta"
make -f align.mk "benchmark_fasta_aligned/coati-tri-ecm/${stem}.coati-tri-ecm.fasta"
make -f align.mk "benchmark_fasta_aligned/coati-dna/${stem}.coati-dna.fasta"
make -f align.mk "benchmark_fasta_aligned/coati-mar-mg/${stem}.coati-mar-mg.fasta"
make -f align.mk "benchmark_fasta_aligned/coati-mar-ecm/${stem}.coati-mar-ecm.fasta"
make -f align.mk "benchmark_fasta_aligned/rev-coati-tri-mg/${stem}.rev-coati-tri-mg.fasta"
make -f align.mk "benchmark_fasta_aligned/mafft/${stem}.mafft.fasta"
make -f align.mk "benchmark_fasta_aligned/macse/${stem}.macse.fasta"
make -f align.mk "benchmark_fasta_aligned/prank/${stem}.prank.fasta"
make -f align.mk "benchmark_fasta_aligned/clustalo/${stem}.clustalo.fasta"
