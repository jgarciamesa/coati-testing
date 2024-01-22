#!/bin/bash

find benchmark_fasta_nogaps -type f -name 'TEST*.fasta' -exec \
    sh -c 'STEM=$(basename {} .nogaps.fasta); echo sbatch -J "align_${STEM}" sbatch/sbatch_align_tests.bash ${STEM}' \;

