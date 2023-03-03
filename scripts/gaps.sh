# create a file with the name of the pairwise aligned ENSEMBL sequences
#  that contain gaps

if [ $# -lt 1 ]
then
	echo "At least one argument is required"
fi

output="data/gaps.csv"

mkdir -p "data/"
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


