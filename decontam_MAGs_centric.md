# Conda install
conda create -n mags -c conda-forge -c bioconda \
  spades=3.15.5 megahit=1.2.9 \
  bwa-mem2=2.2.1 samtools=1.20 bbmap=39.01 \
  metabat2=2.17 maxbin2=2.2.7 concoct=1.1.0 \
  dastool=1.1.6 checkm-genome=1.2.2 \
  gtdbtk=2.4.0 prodigal=2.6.3 hmmer=3.3.2 \
  fastani=1.33 parallel=20240722 seqkit=2.8.2
conda activate mags

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
<img width="252" height="154" alt="Capture d’écran, le 2025-09-23 à 15 57 09" src="https://github.com/user-attachments/assets/0118dfda-5642-4abd-8eb2-cc7631bfa892" />

# (when ready) ALL
quast.py -t 16 -m 1000 -o $QC/quast_ALL $ASM/ALL/final.contigs.fa

grep -E 'contigs \(>= 1000 bp\)|Total length|Largest contig|N50|L50|GC %' \
  $QC/quast_ALL/report.txt






# Diagnostic
# 1) How many reads in the NC trimmed files?
seqkit stats -T -a $READS/NC_*_paired_R1.fastq.gz $READS/NC_*_paired_R2.fastq.gz \
  > $OUT/nc_trimmed_seqstats.tsv
column -t -s $'\t' $OUT/nc_trimmed_seqstats.tsv | sed -n '1,5p'


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






