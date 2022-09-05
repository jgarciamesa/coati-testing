# create a file with the name of the pairwise aligned ENSEMBL sequences
#  that DON'T contain gaps

if [ $# -lt 1 ]
then
	echo "At least one argument is required"
fi

species=$1
N=$2
output="data/${species}/nogaps.csv"

rm -f ${output}

for file in $(head -n${N} raw_data/${species}_geneId.tsv | cut -f1)
do
	c=$(grep -c $(basename ${file}) data/${species}/gaps.csv)
	if [ ${c} -eq 0 ]
	then
		echo $(basename ${file}) >> ${output}
	fi
done

if [ ! -f ${output} ]
then
	echo "No alignment without gaps was found."
fi
