# create a file with the name of the pairwise aligned ENSEMBL sequences
#  that DON'T contain gaps

if [ $# -lt 1 ]
then
	echo "At least one argument is required"
fi

species=$1
output="data/${species}/nogaps.csv"

rm -f ${output}

for file in $(cat data/${species}/filtered.csv | cut -f1)
do
	c=$(grep -c ${file} data/${species}/gaps.csv)
	if [ ${c} -eq 0 ]
	then
		echo ${file} >> ${output}
	fi
done

if [ ! -f ${output} ]
then
	echo "No alignment without gaps was found."
fi
