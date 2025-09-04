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

```bash
fastqc *.fastq.gz -o qc_reports/
```

Review output for:
- Per base sequence quality
- Adapter content
- Overrepresented sequences

### 2. Read Trimming

Trim adapters and low-quality bases using your preferred tool. Example with Trimmomatic:

```bash
trimmomatic PE \
  input_R1.fastq.gz input_R2.fastq.gz \
  output_R1_paired.fastq.gz output_R1_unpaired.fastq.gz \
  output_R2_paired.fastq.gz output_R2_unpaired.fastq.gz \
  ILLUMINACLIP:adapters.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 3. Post-trimming QC

Run FastQC again on trimmed reads:

```bash
fastqc output_R1_paired.fastq.gz output_R2_paired.fastq.gz -o qc_reports/
```

### 4. MultiQC Aggregation

Aggregate all FastQC reports for easy comparison:

```bash
multiqc qc_reports/
```

## Additional Recommendations

- Always inspect QC reports before and after trimming.
- Adjust trimming parameters as needed based on initial QC.
- Document all QC steps and parameters for reproducibility.

## References

- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [MultiQC](https://multiqc.info/)