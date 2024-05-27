#!/bin/bash

source ~/mambaforge/bin/activate scRCAT-seq2

for i in $(find ./ -maxdepth 1 -type f -name "*TES*1.fq.gz" | sort); do
	date
	echo "STAR MAPPING on $i"

	sample_filename=$(basename $i)
	id=${sample_filename%_1.fq.gz}
	dirname=$(dirname $i)
	echo "Sample ID: $id"

	CMD="STAR --genomeDir /home/xyh/index/STARSOLO/gencode_v44/gencode_v44_starsolo \
		--readFilesIn ${dirname}/${id}_2.fq.gz $i --readFilesCommand zcat \
		--soloType CB_UMI_Simple --soloCBstart 11 --soloCBlen 16 --soloUMIstart 27 --soloUMIlen 12 \
		--soloFeatures Gene GeneFull SJ \
		--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
		--runThreadN 40 \
		--outSAMtype BAM SortedByCoordinate \
		--soloBarcodeMate 1 \
		--clip5pNbases 38 0 \
		--sjdbGTFfile /home/xyh/index/STARSOLO/gencode_v44/gencode.v44.annotation.gtf \
		--outFileNamePrefix ./$id/ \
		--outSAMattrRGline ID:$id SM:TSS LB:Long \
		--soloCBwhitelist /data2/xyh/zsy/D45/first_BIG/process/737K-arc-v1.txt"

	echo "Running command:"
	echo "$CMD"

	eval "$CMD"
	echo
done
