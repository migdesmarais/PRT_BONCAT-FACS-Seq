# FILTERED MAGS
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
cat renamed_derep_mags/*.fa > all_MAGs_unique.fa
```

---

## 2. Build Bowtie2 Index

Create an index for rapid mapping.

```bash
conda activate bowtie2_env2
bowtie2-build all_MAGs_unique.fa all_MAGs_index
```

---

## 3. Map Cleaned AND Combined (2 sequencing runs) Reads to MAGs

Use Bowtie2 and Samtools to map paired-end reads and sort/filter mapped reads.

----------------> NEEDS TO BE CHANGE TO Q20

```bash
conda activate bowtie2_env2

#!/bin/bash
set -euo pipefail

BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS_DIR="$BASE/trimmed_reads_FINAL"
OUT="$BASE/magmap_out"
THREADS=12

# Correct locations (as in your screenshots)
INDEX_DIR="$BASE/PRT_MAGs/derep_mags"
IDX="$INDEX_DIR/all_MAGs_index"            # bowtie2 prefix (no extension)
REF="$INDEX_DIR/all_MAGs_unique.fa"        # concatenated MAGs FASTA

mkdir -p "$OUT"/{logs,bam,counts}

# sanity checks
ls -lh "${IDX}."*.bt2*   >/dev/null
[ -f "$REF" ] || { echo "REF not found: $REF"; exit 1; }
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# if that worked, do the loop Q30
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

# if that worked, do the loop Q20
BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS_DIR="$BASE/trimmed_reads_FINAL"
OUT="$BASE/magmap_out_all_q20"
THREADS=12

# Correct locations (as in your screenshots)
INDEX_DIR="$BASE/MAGs/ALL_MAGs"
IDX="$INDEX_DIR/all_MAGs_index"            # bowtie2 prefix (no extension)
REF="$INDEX_DIR/all_MAGs_unique.fa"        # concatenated MAGs FASTA

mkdir -p "$OUT"/{logs,bam,counts}

for R1 in "$READS_DIR"/*_paired_R1.fastq.gz; do
  [ -e "$R1" ] || { echo "No R1 files in $READS_DIR"; break; }
  R2="${R1/_paired_R1/_paired_R2}"
  [ -f "$R2" ] || { echo "Missing mate for $R1"; exit 1; }

  SAMPLE="$(basename "$R1" | sed 's/_paired_R1\.fastq\.gz$//')"
  echo ">>> Mapping $SAMPLE"

  bowtie2 --very-sensitive -p "$THREADS" --no-unal --no-mixed --no-discordant -k 1 \
    -x "$IDX" -1 "$R1" -2 "$R2" 2> "$OUT/logs/${SAMPLE}_bowtie2.q20.log" \
  | samtools view -h -b -q 20 -F 4 -F 256 -F 2048 \
  | samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.q20.primary.bam"

  samtools index "$OUT/bam/${SAMPLE}.q20.primary.bam"

  samtools calmd -bAr "$OUT/bam/${SAMPLE}.q20.primary.bam" "$REF" > "$OUT/bam/${SAMPLE}.tmp.bam"
  mv -f "$OUT/bam/${SAMPLE}.tmp.bam" "$OUT/bam/${SAMPLE}.q20.primary.bam"

  samtools index "$OUT/bam/${SAMPLE}.q20.primary.bam"
  samtools idxstats "$OUT/bam/${SAMPLE}.q20.primary.bam" > "$OUT/counts/${SAMPLE}_idxstats.q20.tsv"
done
```

Summary of mapping rate (per sample).

```bash
OUT=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out_q20
DEST="$OUT/mapping_from_bowtie2.tsv"

echo -e "sample\ttotal_reads\tmapped_reads\toverall_alignment_rate" > "$DEST"

for f in "$OUT"/logs/*_bowtie2.log; do
  sample=$(basename "$f" _bowtie2.q20.log)

  # total reads
  total=$(awk '/ reads; of these:/{print $1; exit}' "$f")

  # mapped reads (concordant pairs only, matches --no-mixed --no-discordant)
  exact1=$(awk '/aligned concordantly exactly 1 time/{print $1; exit}' "$f")
  gt1=$(awk '/aligned concordantly >1 times/{print $1; exit}' "$f")
  mapped=$(( (${exact1:-0}) + (${gt1:-0}) ))

  # overall alignment rate string (e.g. "0.35% overall alignment rate")
  rate=$(awk '/overall alignment rate/{print $1; exit}' "$f")   # keeps the %

  printf "%s\t%s\t%s\t%s\n" "$sample" "${total:-0}" "${mapped:-0}" "${rate:-0%}" >> "$DEST"
done

# pretty view
column -t -s$'\t' "$DEST" | less -S
```

```bash
OUT=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out_q20
DEST="$OUT/mapping_from_bowtie2.q20.tsv"

echo -e "sample\ttotal_reads\tmapped_reads\toverall_alignment_rate" > "$DEST"

shopt -s nullglob
for f in "$OUT"/logs/*_bowtie2.q20.log; do
  sample=$(basename "$f" _bowtie2.q20.log)

  # total reads
  total=$(awk '/ reads; of these:/{print $1; exit}' "$f")

  # mapped reads (concordant pairs only, matches --no-mixed --no-discordant)
  exact1=$(awk '/aligned concordantly exactly 1 time/{print $1; exit}' "$f")
  gt1=$(awk '/aligned concordantly >1 times/{print $1; exit}' "$f")
  mapped=$(( (${exact1:-0}) + (${gt1:-0}) ))

  # overall alignment rate (keeps %)
  rate=$(awk '/overall alignment rate/{print $1; exit}' "$f")

  printf "%s\t%s\t%s\t%s\n" "$sample" "${total:-0}" "${mapped:-0}" "${rate:-0%}" >> "$DEST"
done

# pretty view
column -t -s$'\t' "$DEST" | less -S
```

## 4. Verify MD/NM Tags for CoverM

CoverM requires NM (edit distance) and MD (mismatch string) tags in BAMs.

```bash
# Test one BAM
OUT=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out
samtools view -h $OUT/bam/540AT2ALL.q30.primary.bam | grep -m1 -E "NM:i|MD:Z" || echo "No tags found"
```
---

## 5. Quantify MAG Activity with CoverM

Estimate MAG coverage, RPKM, and relative abundance from mapped reads.

Strict setting 95 75
```bash
coverm genome \
  --bam-files /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out/bam/*.bam \
  --genome-fasta-directory /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/MAGs/derep_mags/renamed_derep_mags \
  --genome-fasta-extension fa \
  --methods covered_bases covered_fraction mean rpkm relative_abundance \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 75 \
  --output-file /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out/mag_coverage_summary_strict.tsv \
  --threads 12
```

Strict setting 95 75
```bash
coverm genome \
  --bam-files /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out_q20/bam/*.bam \
  --genome-fasta-directory /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/MAGs/derep_mags/renamed_derep_mags \
  --genome-fasta-extension fa \
  --methods covered_bases covered_fraction mean rpkm relative_abundance \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 75 \
  --output-file /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out_q20/mag_coverage_summary_strict.tsv \
  --threads 12
```

## 6 RUN CHECKM on MAGs

```
#!/usr/bin/env bash
set -euo pipefail

conda activate checkm_env

BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
MAGS_DIR="$BASE/MAGs/derep_mags/renamed_derep_mags"     # <-- EDIT if needed
OUTDIR="$BASE/MAGs/checkm_derep"
THREADS=24
EXT=fa                                          # <-- fa or fna, etc.

mkdir -p "$OUTDIR"

# Run CheckM
checkm lineage_wf -x "$EXT" -t "$THREADS" --pplacer_threads "$THREADS" \
  "$MAGS_DIR" "$OUTDIR" 2>&1 | tee "$OUTDIR/checkm_lineage_wf.log"

# Produce an easy-to-read QC table (optional but useful)
checkm qa -o 2 -f "$OUTDIR/checkm_QA.tsv" \
  "$OUTDIR/lineage.ms" "$OUTDIR" 2>&1 | tee "$OUTDIR/checkm_qa.log"

echo "[done] CheckM output: $OUTDIR"
echo "[done] QA table: $OUTDIR/checkm_QA.tsv"
```






















# ALL MAGS
# Workflow: Mapping Trimmed BONCAT-FACS-Seq Reads to non-filtered in situ MAGs

This workflow describes how to map quality-controlled sequencing reads from BONCAT-FACS-Seq sorted cells to high-quality in situ metagenome-assembled genomes (MAGs), to determine which MAGs are transcriptionally and translationally active.

---

## CheckM on them and assign taxonomy with graftM, make table
```
conda activate checkm_env
checkm lineage_wf -x fa /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/PRT_MAGs/mags_all_completeness_contam/renamed_all_mags /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/PRT_MAGs/mags_all_completeness_contam/renamed_all_mags/checkm -t 20

????
checkm qa -o 2 -f $RUN/checkm/qa.tsv $RUN/checkm/lineage.ms $RUN/checkm
```

# 6) GTDB-Tk (taxonomy for bins)
```
conda activate gtdbtk_env
export GTDBTK_DATA_PATH=/data_store/gtdbtk_db/release226   # your DB

EXT=fa                          # change if your bins are .fasta or .fna

gtdbtk classify_wf \
  --genome_dir "/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/PRT_MAGs/mags_all_completeness_contam/renamed_all_mags/_no_allMAGs_unique" \
  --extension "$EXT" \
  --out_dir "/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/PRT_MAGs/mags_all_completeness_contam/renamed_all_mags/_no_allMAGs_unique/gtdbtk_classification" \
  --cpus 24 \
  --skip_ani_screen
```

## 1. Prepare MAG Reference Sequences

**Ensure each MAG FASTA has unique headers for downstream mapping and quantification.**

```bash
mkdir -p renamed_all_mags

for MAG in /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/PRT_MAGs/mags_all_completeness_contam/*.fa; do
  ID=$(basename "$MAG" .fa)
  awk -v prefix="$ID" '/^>/{print ">" prefix "_" substr($0,2)} !/^>/' "$MAG" > renamed_all_mags/"$ID.fa"
done
```

**Concatenate all renamed MAGs into a single reference fasta:**

```bash
cat *.fa > all_MAGs_unique.fa
```

---

## 2. Build Bowtie2 Index

Create an index for rapid mapping.

```bash
conda activate bowtie2_env2
bowtie2-build all_MAGs_unique.fa all_MAGs_index
```
---

## 3. Map Cleaned Reads to MAGs - Q30

Use Bowtie2 and Samtools to map paired-end reads and sort/filter mapped reads.

```bash
#!/bin/bash
set -euo pipefail
# conda activate bowtie2_env2

# --------- paths & params ---------
BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS_DIR="$BASE/trimmed_reads"
OUT_BASE="$BASE/magmap_out_all"
RUN_TAG=q30
OUT="$OUT_BASE/$RUN_TAG"
THREADS=12
MAPQ=30

INDEX_DIR="$BASE/PRT_MAGs/mags_all_completeness_contam/renamed_all_mags"
IDX="$INDEX_DIR/all_MAGs_index"      # bowtie2 index prefix (no extension)
REF="$INDEX_DIR/all_MAGs_unique.fa"  # concatenated MAGs FASTA
# ----------------------------------

mkdir -p "$OUT"/{logs,bam,counts}

# sanity checks
ls -lh "${IDX}."*.bt2* >/dev/null
[ -f "$REF" ] || { echo "REF not found: $REF"; exit 1; }
[ -f "${REF}.fai" ] || samtools faidx "$REF"

for R1 in "$READS_DIR"/*_paired_R1.fastq.gz; do
  [ -e "$R1" ] || { echo "No R1 files in $READS_DIR"; break; }
  R2="${R1/_paired_R1/_paired_R2}"
  [ -f "$R2" ] || { echo "Missing mate for $R1"; exit 1; }
  SAMPLE="$(basename "$R1" | sed 's/_paired_R1\.fastq\.gz$//')"
  echo ">>> Mapping $SAMPLE  (RUN_TAG=$RUN_TAG, MAPQ=$MAPQ)"

  # Maximize alignments while staying clean:
  # - very-sensitive-local: soft-clip rescue for partials
  # - allow mixed/discordant: keep single-end & odd pairs (we'll count as single later)
  # - -X 2000: permissive insert size
  bowtie2 --very-sensitive-local -p "$THREADS" -k 1 -x "$IDX" -1 "$R1" -2 "$R2" -X 2000 \
    2> "$OUT/logs/${SAMPLE}_bowtie2_${RUN_TAG}.log" \
  | tee >(samtools view -h -b -F 0x904 \
           | samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.primary0.bam" \
           && samtools index "$OUT/bam/${SAMPLE}.primary0.bam") \
  | samtools view -h -b -q "$MAPQ" -F 0x904 \
  | samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam"

  samtools index "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam"

  # Add MD/NM tags for downstream identity/coverage calcs
  samtools calmd -bAr "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam" "$REF" > "$OUT/bam/${SAMPLE}.tmp.bam"
  mv -f "$OUT/bam/${SAMPLE}.tmp.bam" "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam"
  samtools index "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam"

  samtools idxstats "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam" > "$OUT/counts/${SAMPLE}_idxstats_${RUN_TAG}.tsv"
done

```

## 3. Map Cleaned Reads to MAGs - Q20

Use Bowtie2 and Samtools to map paired-end reads and sort/filter mapped reads.
```
#!/bin/bash
set -euo pipefail
# conda activate bowtie2_env2

BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS_DIR="$BASE/trimmed_reads"
OUT_BASE="$BASE/magmap_out_all"
RUN_TAG=q20
OUT="$OUT_BASE/$RUN_TAG"
THREADS=12

INDEX_DIR="$BASE/PRT_MAGs/mags_all_completeness_contam/renamed_all_mags"
IDX="$INDEX_DIR/all_MAGs_index"            # bowtie2 prefix (no extension)
REF="$INDEX_DIR/all_MAGs_unique.fa"        # concatenated MAGs FASTA

mkdir -p "$OUT"/{logs,bam,counts}

# sanity checks
ls -lh "${IDX}."*.bt2* >/dev/null
[ -f "$REF" ] || { echo "REF not found: $REF"; exit 1; }
[ -f "${REF}.fai" ] || samtools faidx "$REF"

for R1 in "$READS_DIR"/*_paired_R1.fastq.gz; do
  [ -e "$R1" ] || { echo "No R1 files in $READS_DIR"; break; }
  R2="${R1/_paired_R1/_paired_R2}"
  [ -f "$R2" ] || { echo "Missing mate for $R1"; exit 1; }
  SAMPLE="$(basename "$R1" | sed 's/_paired_R1\.fastq\.gz$//')"
  echo ">>> Mapping $SAMPLE  (RUN_TAG=$RUN_TAG)"

  # Maximize mapping:
  #  - --very-sensitive-local : allow soft clipping to rescue partial matches
  #  - allow mixed/discordant : keep single-end and odd pairs (we'll count as single later)
  #  - -X 2000               : permit wide insert sizes (metagenomes can be messy)
  bowtie2 --very-sensitive-local -p "$THREADS" -k 1 -x "$IDX" -1 "$R1" -2 "$R2" -X 2000 \
    2> "$OUT/logs/${SAMPLE}_bowtie2_${RUN_TAG}.log" \
  | tee >(samtools view -h -b -F 0x904 | samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.primary0.bam" \
           && samtools index "$OUT/bam/${SAMPLE}.primary0.bam") \
  | samtools view -h -b -q 20 -F 0x904 \
  | samtools sort -@ "$THREADS" -o "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam"

  samtools index "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam"

  # Add MD/NM tags so downstream tools can compute identity/coverage
  samtools calmd -bAr "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam" "$REF" > "$OUT/bam/${SAMPLE}.tmp.bam"
  mv -f "$OUT/bam/${SAMPLE}.tmp.bam" "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam"
  samtools index "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam"

  samtools idxstats "$OUT/bam/${SAMPLE}.${RUN_TAG}.primary.bam" > "$OUT/counts/${SAMPLE}_idxstats_${RUN_TAG}.tsv"
done
```

Summary of mapping rate (per sample): Q20 and Q30

```bash
OUT=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out_all/q20
DEST="$OUT/mapping_from_bowtie2_20.tsv"

echo -e "sample\ttotal_reads\tmapped_reads\toverall_alignment_rate" > "$DEST"

for f in "$OUT"/logs/*_bowtie2_q20.log; do
  sample=$(basename "$f" _bowtie2.log)

  # total reads
  total=$(awk '/ reads; of these:/{print $1; exit}' "$f")

  # mapped reads (concordant pairs only, matches --no-mixed --no-discordant)
  exact1=$(awk '/aligned concordantly exactly 1 time/{print $1; exit}' "$f")
  gt1=$(awk '/aligned concordantly >1 times/{print $1; exit}' "$f")
  mapped=$(( (${exact1:-0}) + (${gt1:-0}) ))

  # overall alignment rate string (e.g. "0.35% overall alignment rate")
  rate=$(awk '/overall alignment rate/{print $1; exit}' "$f")   # keeps the %

  printf "%s\t%s\t%s\t%s\n" "$sample" "${total:-0}" "${mapped:-0}" "${rate:-0%}" >> "$DEST"
done

# pretty view
column -t -s$'\t' "$DEST" | less -S
```

## 4. Verify MD/NM Tags for CoverM

CoverM requires NM (edit distance) and MD (mismatch string) tags in BAMs.

```bash
# Test one BAM
OUT=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out_all
samtools view -h $OUT/bam/260ATALL_S4_L003.q30.primary.bam | grep -m1 -E "NM:i|MD:Z" || echo "No tags found"
```
---

## 5. Quantify MAG Activity with CoverM

Estimate MAG coverage, RPKM, and relative abundance from mapped reads.
```
conda activate coverm_env
```
Q20 Strict setting 95 50
```bash
coverm genome \
  --bam-files /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out_all/q20/bam/*.bam \
  --genome-fasta-directory /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/PRT_MAGs/mags_all_completeness_contam/renamed_all_mags \
  --genome-fasta-extension fa \
  --methods covered_bases rpkm relative_abundance \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 50 \
  --output-file /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out_all/q20/mag_coverage_summary_q20strict.tsv \
  --threads 12
```
Q20 Strict avg 90 75
```bash
coverm genome \
  --bam-files /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out_all/q20/bam/*.bam \
  --genome-fasta-directory /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/PRT_MAGs/mags_all_completeness_contam/renamed_all_mags \
  --genome-fasta-extension fa \
  --methods covered_bases rpkm relative_abundance \
  --min-read-percent-identity 90 \
  --min-read-aligned-percent 75 \
  --output-file /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out_all/q20/mag_coverage_summary_q20_90_75.tsv \
  --threads 12
```

Relaxed settings 85 50
```bash
coverm genome \
  --bam-files /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out/bam/*.bam \
  --genome-fasta-directory /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/PRT_MAGs/renamed_derep_mags \
  --genome-fasta-extension fa \
  --min-read-percent-identity 85 \
  --min-read-aligned-percent 50 \
  --methods covered_bases rpkm relative_abundance \
  --threads 12 \
  --output-file /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out/mag_coverage_summary_85-50.tsv
```

Relaxed settings 80 50
```bash
coverm genome \
  --bam-files /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out/bam/*.bam \
  --genome-fasta-directory /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/PRT_MAGs/renamed_derep_mags \
  --genome-fasta-extension fa \
  --min-read-percent-identity 80 \
  --min-read-aligned-percent 50 \
  --methods covered_bases rpkm relative_abundance \
  --threads 12 \
  --output-file /scratch/mdesmarais/PRT_BONCAT-FACS-SEQ/magmap_out/mag_coverage_summary_80-50.tsv
```











