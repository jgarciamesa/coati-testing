# create a file with the name of the pairwise aligned ENSEMBL sequences
#  that DON'T contain gaps

for file in $(ls data/ref_alignments/*)
do
	cat ${file} | tr -d '-' > data/no_gaps_ref/$(basename ${file})
done

