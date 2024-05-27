#!/bin/bash

for i in $(find ../rawdata -name "*TSS*1.fq.gz" | sort); do
	date
	echo "Processing $i"

	filename=$(basename "$i")
	id=${filename%_1.fq.gz}

	# Use sed to escape spaces in directory paths
	dir=$(dirname "$i" | sed "s/ /\\ /g")

	CMD="../../data/gzp -I $dir/${id}_1.fq.gz -i $dir/${id}_2.fq.gz \
		-O ${id}_1.fq.gz -o ${id}_2.fq.gz \
		-p CTTCCGATCT -P GAGTACATGG  -u 28 -t 20 -r true"

	echo "Running command:"
	echo "$CMD"

	eval "$CMD"
	echo 
done
