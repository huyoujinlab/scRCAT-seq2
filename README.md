# Description

This pipeline is designed for processing scRCAT-seq2 raw fastq data, enabling the identification of isoforms at high-throughput levels.

# System Requirements

This pipeline has been successfully tested on CentOS and Ubuntu operating systems. The following software and their specified versions are required:

- Python (version 3.7.12)
- STAR (version 2.7.11b)
- samtools(version 1.17)

## Installation

We recommend installing the required software and packages using Conda.

### Download Anaconda

Download the latest version of Anaconda suitable for your system.

### Install Anaconda

Follow the instructions to install Anaconda on your system.

we recommond usage mamba(https://github.com/mamba-org/mamba) to install these package quicly.

### Build Environment

Create and activate a Conda environment:

```bash
### create environment

conda create -n scRCAT-seq2 python=3.7.12
conda activate scRCAT-seq2

### Install Software and Packages

Install the necessary software and packages using Conda:

conda install -c bioconda star=2.7.11b
conda install -c bioconda samtools=1.17

### Please note that the installation time will vary based on your system and network speed.
```

# Demo

This section demonstrates the pipeline's application for identifying isoforms in hESC scRCAT-seq2 data.

## library structure
```
TSS:
@LH00524:15:223WTVLT4:4:1101:18924:1098 2:N:0:CTCTCTAC+TACTCCTT
CGCGTGCCGTATAGCGTTAACTTACTTCAGTCGTTCCGCTACTCTGCGTTGATACCACTGCTTGCAATGAAGTCGCAGGGTTGGGGGGGAATTCAGATAAAACGAATAGCTCGTAACCAAACATGCACAGCGGTCAAACAGTATGTCCCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9II

TES:
@LH00524:15:223WTVLT4:8:1101:6487:1168 2:N:0:CTCTCTAC+AGGCTTAG
CAACCCTGCGACTTCATTGCAAGCAGTGGTATCAACGCAGAGTCCACTCCGCTGCGGGGGCAGTTAACGCTATACGGCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCAGAACTCAAGACGTTAAACGTTCTTGGCGCAAATA
+
IIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIII9IIIIIIIIIII-IIIII9--I-9III----99--99-IIII9----9-9999
```
Both TSS and TES data have the same barcode and UMI sequence structure in R2 read, but the barcode and UMI sequences reverse complement, see more details on our manuscript to get more details on the barcode and UMI sequence structure (suppmentary figure 2 and 3).

## Usage

Run the program gzp to identify and extract valid reads containing cDNA or cell barcode, UMI sequences based on a fixed 10-base tag from the raw paired-end sequencing fastq files, with the command demonstrated in the identify_valid_reads.sh script.

```bash
usage: gzp --help
gzp 1.0
Searches for sequence patterns in paired-end FASTQ files

USAGE:
    gzp [OPTIONS] --input1 <INPUT1> --input2 <INPUT2> --output1 <OUTPUT1> --output2 <OUTPUT2> --pattern1 <PATTERN1> --umi <UMI>

OPTIONS:
    -a, --remove-polya <REMOVE_POLYA>
            If true, remove the polyA sequence from the cDNA sequence [default: false] [possible
            values: true, false]
    -d, --only-cut-adpater <ONLY_CUT_ADPATER>
            If true, only trim the sequence ahead of adapter [default: false] [possible values:
            true, false]
    -h, --help
            Print help information
    -i, --input1 <INPUT1>
            Input FASTQ file for R1, optionally gzipped. R1 contains the cDNA sequence
    -I, --input2 <INPUT2>
            Input FASTQ file for R2, optionally gzipped. R2 contains the Cell Barcode and UMI
            sequence
    -m, --max-mismatches <MAX_MISMATCHES>
            Maximum number of mismatches allowed when searching for patterns [default: 1]
    -o, --output1 <OUTPUT1>
            Output FASTQ file for R1, containing the cDNA sequence from the R1 file
    -O, --output2 <OUTPUT2>
            Output FASTQ file for R2, containing the cDNA sequence from the R2 file
    -p, --pattern1 <PATTERN1>
            DNA sequence to search for in R2, used to locate the UMI sequence
    -P, --pattern3 <PATTERN3>
            DNA sequence to search for in R1, used for trimming the TSO and other sequences in R1
            sequence and find the cDNA sequence in R2
    -q, --pattern2 <PATTERN2>
            DNA sequence to search for in R2, used as a helper sequence to locate the CB and UMI
            [default: ]
    -r, --reverse <REVERSE>
            If false, search for the CB and UMI pattern in the forward direction; if true, search in
            the reverse direction [default: false] [possible values: true, false]
    -t, --threads <THREADS>
            Number of threads to use for parallel processing [default: 8]
    -u, --umi <UMI>
            Length of UMI & Cell Barcode sequence between pattern1 and pattern2 (if exists)
    -V, --version
            Print version information
        --verbose <VERBOSE>
            [default: false] [possible values: true, false]

# for TSS data, which contains 
TSS="./gzp -i R1_reads -I R2_reads \
    -o R1_out -O R2_out \
    -p CTTCCGATCT -P GAGTACATGG  -u 28 -t 20 -r true"
TES="./gzp -i R1_reads -I R2_reads \
        -o R1_out -O R2_out \
        -p AACGCAGAGT -q AGTTAACGCT -a true -u 18 -t 20"
eval $TSS 
eval $TES
```

The STAR software was employed to align the identified valid reads to the genome with parameters as demonstrated, see `https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md` to find more help.

```bash
CMD="STAR --genomeDir /home/xyh/index/STARSOLO/gencode_v44/gencode_v44_starsolo 
    --readFilesIn${dirname}/${id}_2.fq.gz $i  --readFilesCommand zcat 
    --soloType CB_UMI_Simple --soloCBstart 11 --soloCBlen 16 --soloUMIstart 27 --soloUMIlen 12
    --soloFeatures Gene GeneFull SJ
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN
    --runThreadN 40
    --outSAMtype BAM SortedByCoordinate
    --soloBarcodeMate 1
    --clip5pNbases 38 0
    --sjdbGTFfile /home/xyh/index/STARSOLO/gencode_v44/gencode.v44.annotation.gtf
    --outFileNamePrefix ./$id/
    --outSAMattrRGline ID:$id SM:TSS LB:Long
    --soloCBwhitelist /data2/xyh/zsy/D45/first_BIG/process/737K-arc-v1.txt"
eval $CMD
```

The program scRCAT-seq was executed with input from the BAM file aligned by STAR, combining reads with identical UMI and cell barcode into an isoform, and annotating the reconstructed isoforms based on the reference transcriptome.The command example for running the process is as follows:

```bash
usage: scRCAT-seq --help
scRCAT-seq 1.0
Your Name

USAGE:
    scRCAT-seq [OPTIONS] --bam <BAM> --gtf <GTF> --out <OUT>

OPTIONS:
    -b, --bam <BAM>          Input bam files
    -g, --gtf <GTF>          Input gtf files
    -h, --help               Print help information
    -o, --out <OUT>          Output file path, output is csv format
    -V, --version            Print version information

./scRCAT-seq -b ./BIG_second_sub_rep2_1_0.25.bam -g /home/xyh/index/STARSOLO/gencode_v44/gencode.v44.annotation.gtf.gz -o ./BIG_second_sub_rep2_1_0.25_summary.csv
```

## Output

The output files of scRCAT-seq2 is the *.csv file, which contains the isoform information and the corresponding gene annotation.