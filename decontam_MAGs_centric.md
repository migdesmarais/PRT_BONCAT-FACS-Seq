# Conda install
# Concatenate reads (Samples and NC, separately)
```
# where your reads live
BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS=$BASE/trimmed_reads
OUT=$BASE/coassembly_reads
mkdir -p "$OUT"
cd "$READS"

# ---- build the two merged libraries ----
# ALL (everything paired, excluding NC_)
ls *_paired_R1.fastq.gz | grep -v '^NC_' > "$OUT/all_R1.list"
sed 's/_R1/_R2/' "$OUT/all_R1.list" > "$OUT/all_R2.list"
cat $(cat "$OUT/all_R1.list") > "$OUT/ALL_R1.fastq.gz"
cat $(cat "$OUT/all_R2.list") > "$OUT/ALL_R2.fastq.gz"

# NC controls
ls NC_*_paired_R1.fastq.gz > "$OUT/nc_R1.list"
sed 's/_R1/_R2/' "$OUT/nc_R1.list" > "$OUT/nc_R2.list"
cat $(cat "$OUT/nc_R1.list") > "$OUT/NC_R1.fastq.gz"
cat $(cat "$OUT/nc_R2.list") > "$OUT/NC_R2.fastq.gz"

echo "[done] Merged to:"
ls -lh "$OUT"/ALL_* "$OUT"/NC_*
```

# Coassemble with megahit
```
conda activate megahit
BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
ASM=$BASE/assemblies
QC=$BASE/qc
READS=$BASE/trimmed_reads
mkdir -p $QC
mkdir -p "$ASM"

megahit -1 "$OUT/ALL_R1.fastq.gz" -2 "$OUT/ALL_R2.fastq.gz" -o "$ASM/ALL" -t 24 --min-contig-len 1000
megahit -1 "$OUT/NC_R1.fastq.gz"  -2 "$OUT/NC_R2.fastq.gz"  -o "$ASM/NC"  -t 24 --min-contig-len 1000

conda activate quast_env

# NC
quast.py -t 16 -m 1000 -o $QC/quast_NC  $ASM/NC/final.contigs.fa

grep -E 'contigs \(>= 1000 bp\)|Total length|Largest contig|N50|L50|GC %' \
  $QC/quast_NC/report.txt

# contigs (>= 1000 bp)      70           
Total length (>= 0 bp)      490538       
Total length (>= 1000 bp)   490538       
Total length (>= 5000 bp)   411806       
Total length (>= 10000 bp)  321225       
Total length (>= 25000 bp)  77085        
Total length (>= 50000 bp)  0            
Largest contig              49808        
Total length                490538       
N50                         12586        
L50                         13

# ALL
quast.py -t 16 -m 1000 -o $QC/quast_ALL  $ASM/ALL/final.contigs.fa

grep -E 'contigs \(>= 1000 bp\)|Total length|Largest contig|N50|L50|GC %' \
  $QC/quast_ALL/report.txt

# contigs (>= 1000 bp)      35566        
Total length (>= 0 bp)      195460572    
Total length (>= 1000 bp)   195460572    
Total length (>= 5000 bp)   140668365    
Total length (>= 10000 bp)  99697914     
Total length (>= 25000 bp)  42884908     
Total length (>= 50000 bp)  16095322     
Largest contig              996967       
Total length                195460572    
N50                         10334        
L50                         4533 

```

# Map + bin samples + NC
# NC
```
conda activate mags

# 0) Clean any partial BAMs
BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS=$BASE/trimmed_reads
ASM=$BASE/assemblies
RUN=$BASE/mag_results/NC
REF=$ASM/NC/final.contigs.fa
mkdir -p $RUN/map
rm -f $RUN/map/*.bam $RUN/map/*.bai

# 1) Rebuild file lists with FULL paths
ls -1 $READS/NC_*_paired_R1.fastq.gz > $RUN/NC_R1.list
sed 's/_R1/_R2/' $RUN/NC_R1.list > $RUN/NC_R2.list
head $RUN/NC_R1.list   # should show /scratch/... full paths

# 2) Index (message “read 0 ALT contigs” is normal)
bwa index "$REF"

# 3) Map each NC library
paste $RUN/NC_R1.list $RUN/NC_R2.list | while read -r R1 R2; do
  S=$(basename "$R1" | sed 's/_L003.*//')
  echo "[*] Mapping $S"
  bwa mem -t 16 "$REF" "$R1" "$R2" | samtools sort -@8 -o $RUN/map/${S}.bam -
  samtools index $RUN/map/${S}.bam
done

# Sanity check
ls -lh $RUN/map/*.bam
samtools idxstats $RUN/map/*.bam | head

# mapping yield (per-BAM), how many reads mapped?
for B in $BASE/mag_results/NC/map/*.bam; do
  echo "== $(basename "$B") =="
  samtools flagstat "$B" | sed -n '1,8p'
done

# depths
jgi_summarize_bam_contig_depths --outputDepth $RUN/map/depth.txt $RUN/map/*.bam
head -n 3 $RUN/map/depth.txt

# Average contig depth
DEP=$BASE/mag_results/NC/map/depth.txt
awk 'NR>1 {n++; d+=$3} END{print "mean_totalAvgDepth =", (n? d/n : 0)}' "$DEP"

# THE ASSEMBLY IS VERY POOR, SO WHAT WE WILL DO NEXT IS JUST CLASSIFY READS WITH KRAKEN2 ONCE THE DB IS DOWNLOADED (GTDB-TK) AND THEN USE THE NC TAXONOMY TO REMOVE MAGS THAT MATCH IN THE SAMPLES "ALL".
```

# Map + bin samples + NC
# ALL
```
BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS=$BASE/trimmed_reads
ASM=$BASE/assemblies/ALL/final.contigs.fa
RUN=$BASE/mag_results/ALL
mkdir -p $RUN/{map,bins/metabat2,bins/maxbin2,concoct,checkm,gtdbtk}


# 1) map each non-NC library separately
conda activate samtools_env

#!/usr/bin/env bash
set -euo pipefail

# ---- config ----
BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
READS=$BASE/trimmed_reads
ASM=$BASE/assemblies/ALL/final.contigs.fa
RUN=$BASE/mag_results/ALL
MAPDIR=$RUN/map
LOGDIR=$RUN/logs
MASTER=$LOGDIR/mapping_master.log
mkdir -p "$MAPDIR" "$LOGDIR"

log(){ echo "[$(date '+%F %T')] $*" | tee -a "$MASTER" ; }

# Build list of ALL (non-NC) R1 files with full paths (robust)
find "$READS" -maxdepth 1 -name "*_paired_R1.fastq.gz" -printf "%f\n" \
  | grep -v '^NC_' \
  | while read -r f; do echo "$READS/$f"; done > "$RUN/ALL_R1.list"

N=$(wc -l < "$RUN/ALL_R1.list")
log "Found $N libraries for ALL"

# Index once
log "Indexing assembly: $ASM"
bwa index "$ASM" 2> >(tee -a "$MASTER" >&2)

# Map each library with live stderr and per-sample logs
i=0
paste "$RUN/ALL_R1.list" <(sed 's/_R1/_R2/' "$RUN/ALL_R1.list") | while read -r R1 R2; do
  i=$((i+1))
  S=$(basename "$R1" | sed 's/_L003.*//')
  BAM="$MAPDIR/${S}.bam"
  BWA_LOG="$LOGDIR/${S}.bwa.log"
  SORT_LOG="$LOGDIR/${S}.sort.log"
  FLAGSTAT="$LOGDIR/${S}.flagstat.txt"

  if [[ -s "$BAM" ]]; then
    log "($i/$N) [SKIP] $S (BAM exists)"
    continue
  fi

  log "($i/$N) Mapping $S"
  # bwa mem + samtools sort; stream stderr to both console and log files
  bwa mem -t 16 "$ASM" "$R1" "$R2" \
    2> >(tee -a "$BWA_LOG" >> "$MASTER" >&2) \
  | samtools sort -@8 -o "$BAM" - \
    2> >(tee -a "$SORT_LOG" >> "$MASTER" >&2)

  samtools index "$BAM"
  samtools flagstat "$BAM" > "$FLAGSTAT"

  # One-line summary to console+master
  MAPLINE=$(awk '/ mapped \(/ && $0 !~ /primary/{print; exit}' "$FLAGSTAT")
  log "($i/$N) Done $S :: $MAPLINE"
done

log "All mapping finished. BAMs in: $MAPDIR"
log "Tail the live log with: tail -f $MASTER"

# 2) contig depths for binners
conda activate binning_tools_env
jgi_summarize_bam_contig_depths --outputDepth $RUN/map/depth.txt $RUN/map/*.bam

# 3) binners
metabat2 -i "$ASM" -a $RUN/map/depth.txt -o $RUN/bins/metabat2/bin -m 1500 -t 24

run_MaxBin.pl -contig "$ASM" -abund $RUN/map/depth.txt -out $RUN/bins/maxbin2/bin -thread 24 -min_contig_length 1500

cut_up_fasta.py "$ASM" -c 10000 -o 0 -m -b $RUN/concoct/contigs_10K.bed > $RUN/concoct/contigs_10K.fa
concoct_coverage_table.py $RUN/concoct/contigs_10K.bed $RUN/map/*.bam > $RUN/concoct/coverage_table.tsv
concoct --composition_file $RUN/concoct/contigs_10K.fa --coverage_file $RUN/concoct/coverage_table.tsv -b $RUN/concoct/
merge_cutup_clustering.py $RUN/concoct/clustering_gt1000.csv > $RUN/concoct/clustering_merged.csv
mkdir -p $RUN/bins/concoct
extract_fasta_bins.py "$ASM" $RUN/concoct/clustering_merged.csv --output_path $RUN/bins/concoct

# 4) refine with DAS_Tool

# on your system:
source ~/miniforge3/bin/activate dastool   # or: conda activate dastool
which Rscript                              # should point to .../envs/dastool/bin/Rscript
DAS_Tool -h
# 0) use the env you just made
conda activate dastool

# use the env that has DAS_Tool
conda activate dastool

# paths
BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
RUN=$BASE/mag_results/ALL
ASM=$BASE/assemblies/ALL/final.contigs.fa
mkdir -p "$RUN/das_tool"

# make the contig↔bin map (MetaBAT2; change extension if needed)
BIN_DIR=$RUN/bins/metabat2
S2B=$RUN/das_tool/metabat2.s2b.tsv
: > "$S2B"
for f in "$BIN_DIR"/*.fa "$BIN_DIR"/*.fasta "$BIN_DIR"/*.fna; do
  [ -e "$f" ] || continue
  bin=$(basename "$f"); bin=${bin%.*}
  awk -v b="$bin" '/^>/{h=$1; sub(/^>/,"",h); print h "\t" b}' "$f" >> "$S2B"
done
wc -l "$S2B"; head -3 "$S2B"

# run DAS_Tool (explicit output prefix!)
OUTBASE=$RUN/das_tool/ALL_DASTool
DAS_Tool \
  -i "$S2B" -l metabat2 \
  -c "$ASM" \
  -o "$OUTBASE" \
  -t 24 --search_engine diamond --write_bins

# check outputs
ls -l ${OUTBASE}_DASTool_bins | head









conda activate dastool

BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
RUN=$BASE/mag_results/ALL
ASM=$BASE/assemblies/ALL/final.contigs.fa
S2B=$RUN/das_tool/metabat2.s2b.tsv     # the TSV you already built
OUTBASE=$RUN/das_tool/ALL_DASTool_lowt

DAS_Tool \
  -i "$S2B" -l metabat2 \
  -c "$ASM" \
  -o "$OUTBASE" \
  -t 24 --search_engine diamond \
  --score_threshold 0 \
  --duplicate_penalty 0.0 \
  --write_bins --write_unbinned

# bins should now be here:
ls -l ${OUTBASE}_DASTool_bins | head

# 5) CheckM on dastool and metabat2 bins
conda activate checkm_env   # or whatever env you installed checkm-genome in
THREADS=24

BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
RUN=$BASE/mag_results/ALL

# (A) DAS_Tool output you showed in the screenshot
DASTOOL_BINS=$RUN/das_tool/ALL_DASTool_lowt_DASTool_bins

# (B) original MetaBAT2 bins
MB2_BINS=$RUN/bins/metabat2

run_checkm () {
  local BIN_DIR="$1" LABEL="$2"
  mkdir -p "$RUN/checkm_${LABEL}/genomes"

  # detect extension: fa > fasta > fna
  local EXT=fa
  ls "$BIN_DIR"/*.fa >/dev/null 2>&1 || EXT=fasta
  ls "$BIN_DIR"/*.${EXT} >/dev/null 2>&1 || EXT=fna

  # symlink genomes except 'unbinned.*'
  find "$BIN_DIR" -maxdepth 1 -type f -name "*.${EXT}" ! -name "unbinned.*" -exec ln -sf {} "$RUN/checkm_${LABEL}/genomes/" \;

  echo "[*] Running CheckM on ${LABEL} (${EXT})"
  checkm lineage_wf -x "${EXT}" "$RUN/checkm_${LABEL}/genomes" "$RUN/checkm_${LABEL}" -t "$THREADS"
  checkm qa -o 2 -f "$RUN/checkm_${LABEL}/qa.tsv" "$RUN/checkm_${LABEL}/lineage.ms" "$RUN/checkm_${LABEL}"
  echo "[done] $RUN/checkm_${LABEL}/qa.tsv"
}

# (A) DAS_Tool bins
run_checkm "$DASTOOL_BINS" "dastool_lowt"

# (B) MetaBAT2 bins
run_checkm "$MB2_BINS" "metabat2"

awk 'BEGIN{FS=OFS="\t"} NR==1{h=$0; next} {print $0,"\tdastool_lowt"}' \
  "$RUN/checkm_dastool_lowt/qa.tsv" > "$RUN/checkm_all.tsv"

awk 'BEGIN{FS=OFS="\t"} FNR==1{next} {print $0,"\tmetabat2"}' \
  "$RUN/checkm_metabat2/qa.tsv" >> "$RUN/checkm_all.tsv"

echo -e "$(head -1 $RUN/checkm_dastool_lowt/qa.tsv)\tsource" | cat - "$RUN/checkm_all.tsv" > "$RUN/checkm_all.tmp" && mv "$RUN/checkm_all.tmp" "$RUN/checkm_all.tsv"
echo "[combined] $RUN/checkm_all.tsv"








# 5) CheckM
conda activate checkm_env

checkm lineage_wf -x fa $RUN/dastool_DASTool_bins $RUN/checkm -t 24
checkm qa -o 2 -f $RUN/checkm/qa.tsv $RUN/checkm/lineage.ms $RUN/checkm

# 6) GTDB-Tk (taxonomy for bins)
conda activate gtdbtk_env
export GTDBTK_DATA_PATH=/data_store/gtdbtk_files   # your DB

BASE=/scratch/mdesmarais/PRT_BONCAT-FACS-SEQ
RUN=$BASE/mag_results/ALL
MB2_BINS=$RUN/bins/metabat2
EXT=fa                          # change if your bins are .fasta or .fna
OUT=$RUN/gtdbtk_metabat2
mkdir -p "$OUT"

gtdbtk classify_wf \
  --genome_dir "$MB2_BINS" \
  -x "$EXT" \
  --out_dir "$OUT" \
  --cpus 24 \
  --skip_ani_screen








