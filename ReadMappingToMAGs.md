# Workflow: Mapping Trimmed BONCAT-FACS-Seq Reads to High-Quality in situ MAGs

This workflow describes how to map quality-controlled sequencing reads from BONCAT-FACS-Seq sorted cells to high-quality in situ metagenome-assembled genomes (MAGs), to determine which MAGs are transcriptionally and translationally active.

---

## 1. Prepare MAG Reference Sequences

**Ensure each MAG FASTA has unique headers for downstream mapping and quantification.**

```bash
mkdir -p renamed_derep_mags

for MAG in /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/PRT_MAGs/derep_mags/final_derep/*.fa; do
  ID=$(basename "$MAG" .fa)
  awk -v prefix="$ID" '/^>/{print ">" prefix "_" substr($0,2)} !/^>/' "$MAG" > renamed_derep_mags/"$ID.fa"
done
```

**Concatenate all renamed MAGs into a single reference fasta:**

```bash
cat renamed_dered_mags/*.fa > all_MAGs_unique.fa
```

---

## 2. Build Bowtie2 Index

Create an index for rapid mapping.

```bash
conda activate bowtie2_env2
bowtie2-build all_MAGs_unique.fa all_MAGs_index
```

---

## 3. Map Cleaned Reads to MAGs

Use Bowtie2 and Samtools to map paired-end reads and sort/filter mapped reads.

```bash
conda activate bowtie2_env2

#!/bin/bash
set -euo pipefail

BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS_DIR="$BASE/trimmed_reads"
OUT="$BASE/magmap_out"
THREADS=12

# Correct locations (as in your screenshots)
INDEX_DIR="$BASE/PRT_MAGs/all_MAGs"
IDX="$INDEX_DIR/all_MAGs_index"            # bowtie2 prefix (no extension)
REF="$INDEX_DIR/all_MAGs_unique.fa"        # concatenated MAGs FASTA

mkdir -p "$OUT"/{logs,bam,counts}

# sanity checks
ls -lh "${IDX}."*.bt2*   >/dev/null
[ -f "$REF" ] || { echo "REF not found: $REF"; exit 1; }
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# map ONE sample first (quick test)
R1=$(ls "$READS_DIR"/*_paired_R1.fastq.gz | head -n1)
R2="${R1/_paired_R1/_paired_R2}"
[ -f "$R1" ] && [ -f "$R2" ] || { echo "Missing mates in $READS_DIR"; exit 1; }
SAMPLE="$(basename "$R1" | sed 's/_paired_R1\.fastq\.gz$//')"

echo ">>> Test mapping $SAMPLE"
bowtie2 --very-sensitive -p "$THREADS" --no-unal --no-mixed --no-discordant -k 1 \
  -x "$IDX" -1 "$R1" -2 "$R2" 2> "$OUT/logs/${SAMPLE}_bowtie2.log" \
| samtools view -h -b -q 30 -F 4 -F 256 -F 2048 \
| samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.q30.primary.bam"

samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"
samtools calmd -bAr "$OUT/bam/${SAMPLE}.q30.primary.bam" "$REF" > "$OUT/bam/${SAMPLE}.tmp.bam"
mv -f "$OUT/bam/${SAMPLE}.tmp.bam" "$OUT/bam/${SAMPLE}.q30.primary.bam"
samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"
samtools idxstats "$OUT/bam/${SAMPLE}.q30.primary.bam" > "$OUT/counts/${SAMPLE}_idxstats.tsv"

# if that worked, do the loop
for R1 in "$READS_DIR"/*_paired_R1.fastq.gz; do
  [ -e "$R1" ] || { echo "No R1 files in $READS_DIR"; break; }
  R2="${R1/_paired_R1/_paired_R2}"
  [ -f "$R2" ] || { echo "Missing mate for $R1"; exit 1; }
  SAMPLE="$(basename "$R1" | sed 's/_paired_R1\.fastq\.gz$//')"
  echo ">>> Mapping $SAMPLE"

  bowtie2 --very-sensitive -p "$THREADS" --no-unal --no-mixed --no-discordant -k 1 \
    -x "$IDX" -1 "$R1" -2 "$R2" 2> "$OUT/logs/${SAMPLE}_bowtie2.log" \
  | samtools view -h -b -q 30 -F 4 -F 256 -F 2048 \
  | samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.q30.primary.bam"

  samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"
  samtools calmd -bAr "$OUT/bam/${SAMPLE}.q30.primary.bam" "$REF" > "$OUT/bam/${SAMPLE}.tmp.bam"
  mv -f "$OUT/bam/${SAMPLE}.tmp.bam" "$OUT/bam/${SAMPLE}.q30.primary.bam"
  samtools index "$OUT/bam/${SAMPLE}.q30.primary.bam"
  samtools idxstats "$OUT/bam/${SAMPLE}.q30.primary.bam" > "$OUT/counts/${SAMPLE}_idxstats.tsv"
done

```

---

## 4. Verify MD/NM Tags for CoverM

CoverM requires NM (edit distance) and MD (mismatch string) tags in BAMs.

```bash
# Test one BAM
samtools view -h $OUT/bam/OB129_S51_L003.q30.primary.bam | grep -m1 -E "NM:i|MD:Z" || echo "No tags found"

# If missing, add them to all BAMs
for bam in $OUT/bam/*.bam; do
  samtools calmd -bAr "$bam" $FASTAS > "${bam%.bam}.tmp.bam"
  mv "${bam%.bam}.tmp.bam" "$bam"
done
```

---

## 5. Quantify MAG Activity with CoverM

Estimate MAG coverage, RPKM, and relative abundance from mapped reads.

```bash
coverm genome \
  --bam-files /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/magmap_out/bam/*.bam \
  --genome-fasta-directory /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/dereplicated_genomes/renamed_mags \
  --methods covered_bases rpkm relative_abundance \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 75 \
  --output-file /scratch/mdesmarais/OB_BONCAT-FACS-SEQ/magmap_out/mag_coverage_summary.tsv \
  --threads 12
```

---

## Notes & Recommendations

- MAGs should be high quality (e.g., >90% completeness, <5% contamination).
- Use additional filtering or normalization as needed for your downstream analysis.
- Document file paths and sample IDs for reproducibility.
- Adjust Bowtie2 and CoverM parameters for your project or experimental design.

---

## References

- [Bowtie2 Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [Samtools Documentation](http://www.htslib.org/doc/samtools.html)
- [CoverM Documentation](https://github.com/wwood/CoverM)
