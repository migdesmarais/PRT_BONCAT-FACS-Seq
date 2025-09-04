# Read Quality Control

This document describes the steps and best practices for performing quality control (QC) on sequencing reads in the `PRT_BONCAT-FACS-Seq` project.

## Overview

Quality control of sequencing reads is a critical step to ensure reliable downstream analyses. QC processes remove low-quality bases, adapter contamination, and artifacts, resulting in high-confidence data.

## Recommended Tools

- **FastQC**: For initial QC assessment of raw reads.
- **Trimmomatic/Cutadapt/fastp**: For trimming adapters and low-quality bases.
- **MultiQC**: For aggregating QC reports across samples.

## Workflow

### 1. FastQC Analysis

Run FastQC on all raw FASTQ files:

```
conda activate fastqc
fastqc -version

mkdir -p /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/raw_data/prefastqc

fastqc -o /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/raw_data/prefastqc /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/raw_data/250813_DTSA1131_1132_1133_NovaX25B/*fastq.gz
```

Review output for:
- Per base sequence quality
- Adapter content
- Overrepresented sequences
```
scp -r mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/raw_data/prefastqc ~/Downloads/
```

### 2. Read Trimming

Trim adapters and low-quality bases using your preferred tool. Example with Trimmomatic:

```
conda activate trimmomatic_env

mkdir -p /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/trimmed_reads

for file in $(cat samples.txt); do
  trimmomatic PE -phred33 -threads 12 \
    ${file}_R1_001.fastq.gz ${file}_R2_001.fastq.gz \
    /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/trimmed_reads/${file}_paired_R1.fastq.gz \
    /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/trimmed_reads/${file}_unpaired_R1.fastq.gz \
    /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/trimmed_reads/${file}_paired_R2.fastq.gz \
    /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/trimmed_reads/${file}_unpaired_R2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```

### 3. Post-trimming QC

Run FastQC again on trimmed reads:

```
conda activate fastqc

mkdir -p /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/trimmed_reads/postfastqc

fastqc /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/trimmed_reads/*paired.fastq.gz -o /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/trimmed_reads/postfastqc
```

### 4. MultiQC Aggregation

Aggregate all FastQC reports for easy comparison:

```
conda create -n qc -y -c conda-forge -c bioconda python=3.10 fastqc multiqc
conda activate qc
multiqc --version
fastqc --version

multiqc /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/trimmed_reads/postfastqc
```

## Additional Recommendations

- Always inspect QC reports before and after trimming.
- Adjust trimming parameters as needed based on initial QC.
- Document all QC steps and parameters for reproducibility.

## References

- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [MultiQC](https://multiqc.info/)
