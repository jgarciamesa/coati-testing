# Description

################################################################################
# Modify values in this section												   #
################################################################################

# species options: gorilla, mouse, dmelanogaster
species=gorilla

# n - number of sequences. Default: 4000.
#  MAX values: gorilla 16809, mouse 15922, dmelanogaster 3157
n=4000

# len - max length of sequences. Default: 3500 nucs/seq.
len=7000

# number of parallel jobs
j=1

################################################################################
################################################################################

# Download geneid list
#make download_geneid SPECIES=${species}

# Download sequences
tput setaf 11; echo "Download sequences"
tput setaf 15
make download_genes SPECIES=${species} N=${n} -j${j}
tput setaf 11; echo "Sequences downloaded"
tput setaf 15

# Filter sequences
tput setaf 11; echo "Filter sequences"
tput setaf 15
make filter SPECIES=${species} N=${n} LEN=${len}
tput setaf 11; echo "Sequences filtered"
tput setaf 15

# Initial alignment
tput setaf 11; echo "Initial alignment"
tput setaf 15
make initial_alignment SPECIES=${species} N=${n} -j${j} -i
tput setaf 11; echo "Initial alignment done"
tput setaf 15

# Identify which alignments have gaps
tput setaf 11; echo "Identify which alignments have gaps   "
tput setaf 15
make data/${species}/gaps.csv SPECIES=${species} N=${n}

# create gap CIGAR strings
tput setaf 11; echo "Create biologically-like gap maps     "
tput setaf 15
rm -f data/${species}/gaps_cigar.csv
make data/${species}/gaps_cigar.csv SPECIES=${species} N=${n} -j${j}

# create list of no gap sequences to be converted into reference alignments
tput setaf 11; echo "Create list of sequences with no gaps"
tput setaf 15
make data/${species}/nogaps.csv SPECIES=${species} N=${n}

# create reference alignments
tput setaf 11; echo "Simulate reference alignments (true data set)"
tput setaf 15
make reference SPECIES=${species} N=${n} -j${j}

# remove gaps from reference
tput setaf 11; echo "Remove gaps from true alignments, creating benchmark data set"
tput setaf 15
make no_gaps_reference SPECIES=${species} N=${n}
#for file in data/${species}/no_gaps_ref/*; do mv $file data/${species}/no_gaps_ref/$(basename $file .ref); done

# align reference alignments
tput setaf 11; echo "Align benchmark data set sequences"
tput setaf 15
make align_ref SPECIES=${species} N=${n} -j${j}

# compute statistics
tput setaf 11; echo "Compute summary statistics"
tput setaf 15
models=(mecm macse mafft mcoati clustalo prank coati dna ecm)
make data/${species}/dseq_summary.csv SPECIES=${species} N=${n}
Rscript --vanilla scripts/number_alignments.R ${species}
Rscript --vanilla scripts/kaks.R ${species} ${models[*]}
