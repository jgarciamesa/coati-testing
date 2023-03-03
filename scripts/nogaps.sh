# create a file with the name of the pairwise aligned ENSEMBL sequences
#  that DON'T contain gaps

output="data/nogaps.csv"

rm -f ${output}

for file in $(cat data/filtered.csv | cut -f1)
do
	c=$(grep -c ${file} data/gaps.csv)
	if [ ${c} -eq 0 ]
	then
		echo ${file} >> ${output}
	fi
done

if [ ! -f ${output} ]
then
	echo "No alignment without gaps was found."
fi
