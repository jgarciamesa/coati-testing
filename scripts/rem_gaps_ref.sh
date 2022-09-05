mkdir -p ./data/no_gaps_ref

if [ $# -lt 1 ]
then
	echo "At least one argument is required"
fi

species=$1

mkdir -p data/${species}/no_gaps_ref

for file in $(ls data/${species}/ref_alignments/*)
do
	cat ${file} | tr -d '-' > data/${species}/no_gaps_ref/$(basename ${file})
done
# create a file with the name of the pairwise aligned ENSEMBL sequences
#  that DON'T contain gaps

