# filter sequences by length - discard long sequences

SPECIES=$1
N=$2
SIZE=$3
OUTPUT=$4

rm -f ${OUTPUT}

for file in $(ls raw_data/${SPECIES}/* | head -n${N} )
do
	s=$(stat --printf="%s" ${file})
	if [ ${s} -le ${SIZE} ]
	then
		echo $(basename ${file}) >> ${OUTPUT}
	fi
done

if [ ! -f ${OUTPUT} ]
then
	echo "No files smaller than ${SIZE} found."
fi