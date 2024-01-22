#!/bin/bash

find benchmark_fasta_nogaps -type f -name 'TEST*.fasta' -exec \
    sh -c 'STEM=$(basename {} .nogaps.fasta); echo sbatch -J "align_${STEM}" slurm/run_align_test.bash ${STEM}' \;
