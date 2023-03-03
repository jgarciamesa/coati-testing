# Description

################################################################################
# Modify values in this section												   #
################################################################################

# n - number of sequences. Default: 4000.
#  MAX value: gorilla 16809
n=16800

# len - max length of sequences. Default: 3500 nucs/seq.
len=6000

# number of parallel jobs
j=4

################################################################################
################################################################################

# Download geneid list
#make download_geneid

# Download sequences
tput setaf 11; echo "Download sequences"
tput setaf 15
make download_genes N=${n} -j${j}
tput setaf 11; echo "Sequences downloaded"
tput setaf 15

# Filter sequences
tput setaf 11; echo "Filter sequences"
tput setaf 15
make filter N=${n} LEN=${len}
tput setaf 11; echo "Sequences filtered"
tput setaf 15

# Initial alignment
tput setaf 11; echo "Initial alignment"
tput setaf 15
make initial_alignment N=${n} -j${j} -i
tput setaf 11; echo "Initial alignment done                "
tput setaf 15

# Identify which alignments have gaps
tput setaf 11; echo "Identify which alignments have gaps   "
tput setaf 15
make data/gaps.csv N=${n}

# create gap CIGAR strings
tput setaf 11; echo "Create biologically-like gap maps     "
tput setaf 15
rm -f data/gaps_cigar.csv
make data/gaps_cigar.csv N=${n} -j${j}

# create list of no gap sequences to be converted into reference alignments
tput setaf 11; echo "Create list of sequences with no gaps"
tput setaf 15
make data/nogaps.csv N=${n}

# create reference alignments
tput setaf 11; echo "Simulate reference alignments (true data set)"
tput setaf 15
make reference N=${n} -j${j}

# remove gaps from reference
tput setaf 11; echo "Remove gaps from true alignments, creating benchmark data set"
tput setaf 15
make no_gaps_reference N=${n}

# align reference alignments
tput setaf 11; echo "Align benchmark data set sequences"
tput setaf 15
make align_ref N=${n} -j${j} -i

# compute statistics
tput setaf 11; echo "Compute summary statistics           "
tput setaf 15
make results/results_summary.csv
make supplementary_materials.pdf
