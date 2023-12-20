#!/bin/bash

seqtk seq -r source_files/test_dataset/all-TSS.R2.fastq > temp-TSS.R2_rev.fastq
python find_smart3_and_trim_umi.py -i temp-TSS.R2_rev.fastq -o temp_TSS_umi_barcode_trim.fastq

python find_smart3_and_trim_TSS_cdna.py -i source_files/test_dataset/all-TSS.R1.fastq -b temp_TSS_cdna_trim_withOligo.fastq
python find_smart3_and_trim_TSS_cdna_noOligo.py -i source_files/test_dataset/all-TSS.R1.fastq -b temp_TSS_cdna_trim_noOligo.fastq
perl cmpfastq_pe.pl temp_TSS_cdna_trim_withOligo.fastq temp_TSS_umi_barcode_trim.fastq
mv temp_TSS_umi_barcode_trim.fastq-common.out temp_TSS_umi_barcode_trim_withOligo_R1.fastq-common.out
perl cmpfastq_pe.pl temp_TSS_cdna_trim_noOligo.fastq temp_TSS_umi_barcode_trim.fastq
mv temp_TSS_umi_barcode_trim.fastq-common.out temp_TSS_umi_barcode_trim_noOligo_R1.fastq-common.out

python find_smart3_and_trim_TSS_cdna.py -i temp-TSS.R2_rev.fastq -b temp_TSS_cdna_trim_withOligo_R2.fastq
perl cmpfastq_pe.pl temp_TSS_cdna_trim_withOligo_R2.fastq temp_TSS_umi_barcode_trim.fastq
mv temp_TSS_umi_barcode_trim.fastq-common.out temp_TSS_umi_barcode_trim_R2.fastq-common.out
awk '{ if (NR%4==1) gsub(" ","R2 "); print }' temp_TSS_umi_barcode_trim_R2.fastq-common.out > suffix_TSS_umi_barcode_trim_R2.fastq-common.out
awk '{ if (NR%4==1) gsub(" ","R2 "); print }' temp_TSS_cdna_trim_withOligo_R2.fastq-common.out  > suffix_TSS_cdna_trim_withOligo_R2.fastq-common.out

cat temp_TSS_umi_barcode_trim_withOligo_R1.fastq-common.out suffix_TSS_umi_barcode_trim_R2.fastq-common.out > all_TSS_umi_barcode_trim_withOligo.fastq-common.out
cat temp_TSS_cdna_trim_withOligo.fastq-common.out suffix_TSS_cdna_trim_withOligo_R2.fastq-common.out > all_TSS_cdna_trim_withOligo.fastq-common.out
seqtk seq -r all_TSS_cdna_trim_withOligo.fastq-common.out > temp_TSS_cdna_trim_withOligo.rc.fastq
seqtk seq -r temp_TSS_cdna_trim_noOligo.fastq-common.out > temp_TSS_cdna_trim_noOligo.rc.fastq
awk '{ if (NR%4==1) gsub(" ","_exonCap "); print }' temp_TSS_umi_barcode_trim_noOligo_R1.fastq-common.out > renamed_TSS_umi_barcode_trim_noOligo.fastq
awk '{ if (NR%4==1) gsub(" ","_exonCap "); print }' temp_TSS_cdna_trim_noOligo.rc.fastq > renamed_TSS_cdna_trim_noOligo.rc.fastq
awk '{ if (NR%4==1) gsub(" ","_TSS "); print }' all_TSS_umi_barcode_trim_withOligo.fastq-common.out > renamed_TSS_umi_barcode_trim_withOligo.fastq
awk '{ if (NR%4==1) gsub(" ","_TSS "); print }' temp_TSS_cdna_trim_withOligo.rc.fastq > renamed_TSS_cdna_trim_withOligo.rc.fastq

python find_smart3_and_trim_umi.py -i source_files/test_dataset/all-TES.R2.fastq -o temp_TES_trimed.R2.fastq
perl cmpfastq_pe.pl temp_TES_trimed.R2.fastq source_files/test_dataset/all-TES.R1.fastq
python3 find_withA10_and_noA10_v2.py -i all-TES.R1.fastq-common.out -w temp_TES_trimed.R1.withA10.fastq -n temp_TES_trimed.R1.noA10.fastq
perl cmpfastq_pe.pl temp_TES_trimed.R1.withA10.fastq temp_TES_trimed.R2.fastq-common.out
mv temp_TES_trimed.R2.fastq-common.out-common.out temp_TES_trimed.R2.fastq-common.out_withA10
perl cmpfastq_pe.pl temp_TES_trimed.R1.noA10.fastq temp_TES_trimed.R2.fastq-common.out
mv temp_TES_trimed.R2.fastq-common.out-common.out temp_TES_trimed.R2.fastq-common.out_noA10

awk '{ if (NR%4==1) gsub(" ","_TES "); print }' temp_TES_trimed.R1.withA10.fastq-common.out > renamed_temp_TES_trimed.R1.withA10.fastq
awk '{ if (NR%4==1) gsub(" ","_TES "); print }' temp_TES_trimed.R2.fastq-common.out_withA10 > renamed_temp_TES_trimed.R2.fastq-common.out_withA10
awk '{ if (NR%4==1) gsub(" ","_exonTail "); print }' temp_TES_trimed.R1.noA10.fastq-common.out > renamed_temp_TES_trimed.R1.noA10.fastq
awk '{ if (NR%4==1) gsub(" ","_exonTail "); print }' temp_TES_trimed.R2.fastq-common.out_noA10 > renamed_temp_TES_trimed.R2.fastq-common.out_noA10

cat renamed_temp_TES_trimed.R1.withA10.fastq renamed_temp_TES_trimed.R1.noA10.fastq renamed_TSS_cdna_trim_noOligo.rc.fastq renamed_TSS_cdna_trim_withOligo.rc.fastq > renamed_all.R1.fastq
cat renamed_temp_TES_trimed.R2.fastq-common.out_withA10 renamed_temp_TES_trimed.R2.fastq-common.out_noA10 renamed_TSS_umi_barcode_trim_noOligo.fastq renamed_TSS_umi_barcode_trim_withOligo.fastq > renamed_all.R2.fastq

gzip renamed_all.R2.fastq
gzip renamed_all.R1.fastq
seqkit seq --name --only-id renamed_all.R1.fastq.gz renamed_all.R2.fastq.gz | sort | uniq -d >id.txt
seqkit grep --pattern-file id.txt renamed_all.R1.fastq.gz -o sorted_R1.fastq.gz
seqkit grep --pattern-file id.txt renamed_all.R2.fastq.gz -o sorted_R2.fastq.gz

rm temp*
rm *out
rm renamed*


if test -e zUMIs_output ; then
rm -r zUMIs_output
fi
bash /home/zjw/biosoft/zUMIs-main/zUMIs.sh -y run_human_20210423.yaml
if ! test -e ./zUMIs_output/ESkept_barcodes_binned.txt ; then
cp ./zUMIs_output/ESkept_barcodes.txt ./zUMIs_output/ESkept_barcodes_binned.txt
bash /home/zjw/biosoft/zUMIs-main/zUMIs.sh -y run_human_20210423.yaml
fi



python ss3_isofrom.py -i ES.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam -c hg38_SIRVercc.conf -e isorestructure -o ./ -p 20 -s hg38 -P -Q













