# create a file with the name of the pairwise aligned ENSEMBL sequences
#  that contain gaps

if [ $# -lt 2 ]
then
	echo "At least two arguments are required"
fi

species="$1"
output="data/${species}/gaps.csv"

shift 1

mkdir -p "data/${species}/"
#rm -f "data/${species}/gaps_${species}_tmp.csv"
rm -f ${output}

for arg in $@
do
	for file in aln/${arg}/*.fasta
	do
		c=$(grep -c '-' ${file})
		if [ ${c} -gt 0 ]
		then
			echo ${file} >> ${output}
		fi
	done
done


