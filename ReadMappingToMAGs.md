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










## 5. Assign taxonomy to trimmed reads to check how much contamination we have.
### Got low mapping rate, contaminants?

```
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# -------- config --------
BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS_DIR="$BASE/trimmed_reads"
OUTDIR="$BASE/kraken"                 # all outputs live here
DB=/data_store/kraken_database        # adjust to your DB path
THREADS=16
CONF=0.10                             # Kraken2 confidence
READLEN=150                           # for Bracken k-mer model
MINREADS=10                           # Bracken filter
MINFRAC=0.0001                        # 0.01% fraction filter

mkdir -p "$OUTDIR"/{logs,bracken/genus,bracken/phylum}

# (activate env that has kraken2)
conda activate kraken_env

cd "$READS_DIR"

echo ">> Discovering samples in $(pwd)"
mapfile -t SAMPLES < <(
  printf '%s\n' *_paired_R1.fastq.gz \
  | sed -E 's/_L[0-9]{3}_paired_R1\.fastq\.gz$//; s/_paired_R1\.fastq\.gz$//' \
  | sort -u
)
echo "   Found ${#SAMPLES[@]} sample keys"

# ------------------------
# 1) Kraken2 per sample
# ------------------------
for s in "${SAMPLES[@]}"; do
  echo "== Kraken2: $s"
  # collect lane-chunked + singleton files
  r1=( "${s}_L"*"_paired_R1.fastq.gz" )
  r2=( "${s}_L"*"_paired_R2.fastq.gz" )
  [[ -e ${s}_paired_R1.fastq.gz ]] && r1+=( "${s}_paired_R1.fastq.gz" )
  [[ -e ${s}_paired_R2.fastq.gz ]] && r2+=( "${s}_paired_R2.fastq.gz" )

  if ((${#r1[@]}==0)) || ((${#r2[@]}==0)); then
    echo "   !! No paired inputs for $s — skipping"
    continue
  fi

  base="$(basename "$s")"
  REP="${OUTDIR}/${base}.kraken.report"
  OUT="${OUTDIR}/${base}.kraken.out"
  LOG="${OUTDIR}/logs/${base}.kraken2.log"

  (
    set -x
    kraken2 --db "$DB" --threads "$THREADS" --paired --use-names \
            --confidence "$CONF" \
            --report "$REP" --output "$OUT" \
            /dev/fd/63 /dev/fd/62 \
            63< <(zcat -f "${r1[@]}") \
            62< <(zcat -f "${r2[@]}")
  ) &> "$LOG"

  echo "   -> wrote: $REP  and  $OUT"
done

# ------------------------
# 2) Bracken (genus + phylum) + per-sample filters
# ------------------------
conda activate bracken_env

echo "== Bracken (genus & phylum)"
for rep in "${OUTDIR}"/*.kraken.report; do
  [[ -e "$rep" ]] || { echo "   (no reports found)"; break; }
  base=$(basename "$rep" .kraken.report)

  # Genus
  bracken -d "$DB" -i "$rep" -o "$OUTDIR/bracken/genus/${base}.bracken.G" \
          -l G -r $READLEN -t $MINREADS \
          -w "$OUTDIR/bracken/genus/${base}.bracken.G.report"

  awk -F'\t' -v OFS='\t' -v mr=$MINREADS -v mf=$MINFRAC -v S="$base" '
    NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
    ($h["new_est_reads"]+0 >= mr) && ($h["fraction_total_reads"]+0 >= mf) {
      print S, $h["name"], $h["new_est_reads"], $h["fraction_total_reads"]
    }' "$OUTDIR/bracken/genus/${base}.bracken.G" \
    > "$OUTDIR/bracken/genus/${base}.genus_filtered.tsv"

  # Phylum
  bracken -d "$DB" -i "$rep" -o "$OUTDIR/bracken/phylum/${base}.bracken.P" \
          -l P -r $READLEN -t $MINREADS \
          -w "$OUTDIR/bracken/phylum/${base}.bracken.P.report"

  awk -F'\t' -v OFS='\t' -v mr=$MINREADS -v mf=$MINFRAC -v S="$base" '
    NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
    ($h["new_est_reads"]+0 >= mr) && ($h["fraction_total_reads"]+0 >= mf) {
      print S, $h["name"], $h["new_est_reads"], $h["fraction_total_reads"]
    }' "$OUTDIR/bracken/phylum/${base}.bracken.P" \
    > "$OUTDIR/bracken/phylum/${base}.phylum_filtered.tsv"
done

# ------------------------
# 3) Simple aggregations (one file per level)
# ------------------------
{
  echo -e "sample\ttaxon\test_reads\tfraction"
  cat "$OUTDIR"/bracken/genus/*.genus_filtered.tsv 2>/dev/null || true
} > "$OUTDIR/bracken/genus_summary.tsv"

{
  echo -e "sample\ttaxon\test_reads\tfraction"
  cat "$OUTDIR"/bracken/phylum/*.phylum_filtered.tsv 2>/dev/null || true
} > "$OUTDIR/bracken/phylum_summary.tsv"

echo "All done."
echo "Kraken reports:         $OUTDIR/*.kraken.report"
echo "Bracken genus summary:  $OUTDIR/bracken/genus_summary.tsv"
echo "Bracken phylum summary: $OUTDIR/bracken/phylum_summary.tsv"

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
