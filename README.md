
# Description

This pipeline is designed for processing scRCAT-seq2 raw fastq data, enabling the identification of transcripts at the single-cell level.

# System Requirements

This pipeline has been successfully tested on CentOS and Ubuntu operating systems. The following software and their specified versions are required:

- Cell Ranger (version 6.1.1)
- Python (version 3.7.12)
- R (version 4.1.3)
- zUMIs (version 2.9.4f)
- STAR (version 2.7.9a)
- deeptools (version 3.5.1)
- minimap2 (version 2.17)
- gffcompare (version 0.11.2)
- Seurat (version 4.1.0)
- Signac (version 1.6.0)
- Monocle (version 2.27.0)
- chromVAR (version 1.16.0)
- rMATS (version 4.1.2)
- rmats2sashimiplot (version 2.0.4)
- bedtoolsR (version 2.30.0-2)
- bedtools (version 2.30.0)

## Installation

We recommend installing the required software and packages using Conda.

### Download Anaconda

Download the latest version of Anaconda suitable for your system.

### Install Anaconda

Follow the instructions to install Anaconda on your system.

### Build Environment

Create and activate a Conda environment:

```bash
conda create -n scRCAT_seq2 python=3.7
conda activate scRCAT_seq2
```

### Install Software and Packages

Install the necessary software and packages using Conda:

```
conda install -c bioconda star=2.7.9a
conda install -c bioconda deeptools=3.5.1
conda install -c bioconda minimap2=2.17
conda install -c bioconda gffcompare=0.11.2
conda install -c bioconda bedtools=2.30.0

# Install R and R packages
conda install -c conda-forge r-base=4.1.3
conda install -c bioconda r-seurat=4.1.0
conda install -c bioconda r-signac=1.6.0
conda install -c bioconda bioconductor-monocle=2.27.0
conda install -c bioconda bioconductor-chromvar=1.16.0
conda install -c bioconda r-rmats=4.1.2
conda install -c conda-forge bedtools
```

Please note that the installation time will vary based on your system and network speed.

# Demo

This section demonstrates the pipeline's application for identifying TSSs and TESs in human embryonic stem cells.

## Usage

Run the provided shell script to map fastq files to the genome and obtain UMI counts for each transcript:

```
sh fastq_extraction_reconstruct.sh
```

Next, use the `isoform_construct.R` script to construct the transcript structure:

```
Rscript isoform_construct.R
```

## Output

The output files from this demo pipeline will be stored in the `output/` directory. The constructed bed file, named `all_exon_12.bed`, will be found in the `assigned_isoforms` folder.
