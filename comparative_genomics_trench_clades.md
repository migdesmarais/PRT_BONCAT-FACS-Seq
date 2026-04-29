1. Run prodigal on all MAGs/SAGs
```
conda activate prodigal

which prodigal
prodigal -v

for fna in genomes_fna_all/*.fna
do
    base=$(basename "$fna" .fna)

    echo "Running Prodigal on $base"

    prodigal \
        -i "$fna" \
        -a "prodigal_out/faa/${base}.faa" \
        -d "prodigal_out/genes_fna/${base}.genes.fna" \
        -o "prodigal_out/gff/${base}.gff" \
        -p meta \
        -f gff
done

ls prodigal_out/faa | wc -l
ls prodigal_out/faa | head

grep -c "^>" prodigal_out/faa/*.faa | head
grep -c "^>" prodigal_out/faa/*.faa | sort -t: -k2,2n | head

grep -c "^>" prodigal_out/faa/*.faa \
  | sed 's#prodigal_out/faa/##; s#.faa:#,#' \
  > results/protein_counts_by_genome.csv
```

# Then run R code for physicochemical analysis.

# 0 Project structure
```
# Main project directory

BASE=/scratch/mdesmarais/JS1_comparative_genomics
cd $BASE

# Expected input folders
# clades/clade_deep/     = trench-associated JS1 genomes
# clades/clade_shallow/  = other/non-trench JS1 genomes
# mags_metadata.xlsx     = metadata table with genome IDs and clade labels

# Create working folders
mkdir -p genomes_fna_all
mkdir -p prodigal_out/faa
mkdir -p prodigal_out/genes_fna
mkdir -p prodigal_out/gff
mkdir -p results
mkdir -p figures
mkdir -p logs
```

# 1 run eggnog on prodigal ORF
```
conda activate eggnog

BASE_DIR="/scratch/mdesmarais/JS1_comparative_genomics"
FAA_DIR="${BASE_DIR}/prodigal_out/faa"
OUT_DIR="${BASE_DIR}/eggnog_out_per_genome"
DATA_DIR="/home/mdesmarais/.conda/envs/eggnog/lib/python2.7/site-packages/data"

mkdir -p "${OUT_DIR}"

rm -f "${OUT_DIR}"/*.eggnog.stdout.txt
rm -f "${OUT_DIR}"/*.eggnog.stderr.txt
rm -f "${OUT_DIR}"/*.emapper.*

for faa in "${FAA_DIR}"/*.faa; do
  genome=$(basename "${faa}" .faa)

  echo "============================================================"
  echo "Running eggNOG-mapper: ${genome}"
  echo "Input: ${faa}"
  echo "============================================================"

  emapper.py \
    -i "${faa}" \
    -o "${OUT_DIR}/${genome}" \
    --cpu 8 \
    -m diamond \
    --data_dir "${DATA_DIR}" \
    > "${OUT_DIR}/${genome}.eggnog.stdout.txt" \
    2> "${OUT_DIR}/${genome}.eggnog.stderr.txt"
done

echo "Done."
echo "Annotation files:"
find "${OUT_DIR}" -name "*.emapper.annotations" | wc -l
find "${OUT_DIR}" -name "*.emapper.annotations" | head
```

# Run eggnog
```
conda activate eggnog

emapper.py \
  -i eggnog_out/JS1_all_predicted_proteins.prefixed.faa \
  --qtype seq \
  -m diamond \
  -o JS1_eggnog \
  --output_dir eggnog_out \
  --cpu 16 \
  2>&1 | tee logs/JS1_eggnog_mapper.log
```

# Run kofamscan: KO-focused annotation using KEGG HMM profiles.
```
conda activate kofamscan

KOFAM_DB_DIR="/scratch/mdesmarais/databases/kofam"
mkdir -p "${KOFAM_DB_DIR}"
cd "${KOFAM_DB_DIR}"

# Download database files
wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz

conda activate kofamscan

BASE_DIR="/scratch/mdesmarais/JS1_comparative_genomics"
FAA_DIR="${BASE_DIR}/prodigal_out/faa"
KOFAM_DB_DIR="/scratch/mdesmarais/databases/kofam"
OUT_DIR="${BASE_DIR}/kofamscan_out"

mkdir -p "${OUT_DIR}"

THREADS=8

# Optional: remove failed old outputs
rm -f "${OUT_DIR}"/*.kofamscan.stderr.txt
rm -f "${OUT_DIR}"/*.kofamscan.stdout.txt
rm -f "${OUT_DIR}"/*.kofamscan.detail.tsv

for faa in "${FAA_DIR}"/*.faa; do
  genome=$(basename "${faa}" .faa)

  echo "============================================================"
  echo "Running KOfamScan: ${genome}"
  echo "============================================================"

  exec_annotation \
    -p "${KOFAM_DB_DIR}/profiles" \
    -k "${KOFAM_DB_DIR}/ko_list" \
    --cpu "${THREADS}" \
    -f detail-tsv \
    -o "${OUT_DIR}/${genome}.kofamscan.detail.tsv" \
    "${faa}" \
    > "${OUT_DIR}/${genome}.kofamscan.stdout.txt" \
    2> "${OUT_DIR}/${genome}.kofamscan.stderr.txt"
done

echo "Done."
ls -lh "${OUT_DIR}" | head
```

# dbCAN: CAZyme/carbohydrate-active enzyme annotation.
```
conda create -n dbcan_new -c conda-forge -c bioconda python=3.11 diamond hmmer prodigal -y
conda activate dbcan_new

pip install dbcan
which run_dbcan
run_dbcan --help | head -n 30

conda activate dbcan_new

DBCAN_DB_DIR="/scratch/mdesmarais/databases/dbcan"
mkdir -p "$DBCAN_DB_DIR"

run_dbcan database --db_dir "$DBCAN_DB_DIR" --no-cgc --aws_s3

# ============================================================
# Run dbCAN on all JS1 protein FASTA files — correct syntax
# ============================================================

conda activate dbcan_new

BASE_DIR="/scratch/mdesmarais/JS1_comparative_genomics"
FAA_DIR="${BASE_DIR}/prodigal_out/faa"
DBCAN_DB_DIR="/scratch/mdesmarais/databases/dbcan"
OUT_DIR="${BASE_DIR}/dbcan_out"

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

for faa in "${FAA_DIR}"/*.faa; do

  genome=$(basename "${faa}" .faa)

  echo "============================================================"
  echo "Running dbCAN: ${genome}"
  echo "============================================================"

  run_dbcan CAZyme_annotation \
    --input_raw_data "${faa}" \
    --mode protein \
    --output_dir "${OUT_DIR}/${genome}" \
    --db_dir "${DBCAN_DB_DIR}" \
    --methods diamond,hmm,dbCANsub \
    --threads 8 \
    > "${OUT_DIR}/${genome}.dbcan.stdout.txt" \
    2> "${OUT_DIR}/${genome}.dbcan.stderr.txt"

done

echo "Done."

echo
echo "Genome output directories:"
find "${OUT_DIR}" -mindepth 1 -maxdepth 1 -type d | wc -l

echo
echo "Overview files:"
find "${OUT_DIR}" -name "overview.tsv" -o -name "overview.txt" | head

echo
echo "Errors:"
grep -i "error\|failed\|cannot\|traceback\|no such option" "${OUT_DIR}"/*.stderr.txt | head -n 50
```

# Merops: Peptidase/protease family annotation.
```
conda create -n merops -c conda-forge -c bioconda diamond blast -y
conda activate merops

MEROPS_DB_DIR="/scratch/mdesmarais/databases/merops"
mkdir -p "$MEROPS_DB_DIR"
cd "$MEROPS_DB_DIR"

# Download MEROPS peptidase-unit FASTA from current release
wget -c ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib

# Check it
ls -lh pepunit.lib
head -n 5 pepunit.lib

diamond makedb \
  --in "${MEROPS_DB_DIR}/pepunit.lib" \
  -d "${MEROPS_DB_DIR}/merops_pepunit"

ls -lh "${MEROPS_DB_DIR}"

# ============================================================
# Clean MEROPS pepunit.lib for DIAMOND
# ============================================================

conda activate merops

MEROPS_DB_DIR="/scratch/mdesmarais/databases/merops"
cd "${MEROPS_DB_DIR}"

awk '
  /^>/ {
    print
    next
  }
  {
    gsub(/[^A-Za-z*.-]/, "", $0)
    gsub(/\*/, "", $0)
    gsub(/\./, "", $0)
    gsub(/-/, "", $0)
    print toupper($0)
  }
' pepunit.lib > pepunit.clean.faa

echo "Check cleaned file:"
head -n 8 pepunit.clean.faa

echo "Any spaces left?"
grep -n " " pepunit.clean.faa | head

# BUILD DATABASE
conda activate merops

MEROPS_DB_DIR="/scratch/mdesmarais/databases/merops"

diamond makedb \
  --in "${MEROPS_DB_DIR}/pepunit.clean.faa" \
  -d "${MEROPS_DB_DIR}/merops_pepunit"

ls -lh "${MEROPS_DB_DIR}"

# RUN ALL GENOMES
conda activate merops

BASE_DIR="/scratch/mdesmarais/JS1_comparative_genomics"
FAA_DIR="${BASE_DIR}/prodigal_out/faa"
MEROPS_DB_DIR="/scratch/mdesmarais/databases/merops"
OUT_DIR="${BASE_DIR}/merops_out"

mkdir -p "${OUT_DIR}"

for faa in "${FAA_DIR}"/*.faa; do
  genome=$(basename "${faa}" .faa)

  echo "============================================================"
  echo "Running MEROPS: ${genome}"
  echo "============================================================"

  diamond blastp -q "${faa}" -d "${MEROPS_DB_DIR}/merops_pepunit" -o "${OUT_DIR}/${genome}.merops.diamond.tsv" --threads 8 --evalue 1e-10 --query-cover 50 --subject-cover 50 --max-target-seqs 5 --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore stitle > "${OUT_DIR}/${genome}.merops.stdout.txt" 2> "${OUT_DIR}/${genome}.merops.stderr.txt"
done

echo "Done."
ls -lh "${OUT_DIR}" | head
```

# Inteprorscan: Protein family/domain annotation.
```
# ============================================================
# Run InterProScan on all JS1 protein FASTA files
# ============================================================

# ============================================================
# Run InterProScan on all cleaned JS1 protein FASTA files overnight
# ============================================================

conda activate java11-env

BASE_DIR="/scratch/mdesmarais/JS1_comparative_genomics"
CLEAN_FAA_DIR="${BASE_DIR}/prodigal_out/faa_interproscan_clean"
IPRSCAN_DIR="/scratch/mdesmarais/PRT_DGE/MAGs/prodigal/my_interproscan/interproscan-5.74-105.0"
OUT_DIR="${BASE_DIR}/interproscan_out"

mkdir -p "${OUT_DIR}"

THREADS=8

echo "Starting InterProScan all-genome run:"
date
echo "Input directory: ${CLEAN_FAA_DIR}"
echo "Output directory: ${OUT_DIR}"

for faa in "${CLEAN_FAA_DIR}"/*.clean.faa; do

  genome=$(basename "${faa}" .clean.faa)

  echo "============================================================"
  echo "Running InterProScan: ${genome}"
  echo "Started: $(date)"
  echo "============================================================"

  if [ -s "${OUT_DIR}/${genome}.tsv" ]; then
    echo "Already exists, skipping: ${OUT_DIR}/${genome}.tsv"
    continue
  fi

  "${IPRSCAN_DIR}/interproscan.sh" \
    -i "${faa}" \
    -b "${OUT_DIR}/${genome}" \
    -f TSV,GFF3 \
    -goterms \
    -pa \
    -exclappl MobiDBLite \
    -cpu "${THREADS}" \
    > "${OUT_DIR}/${genome}.interproscan.log" \
    2>&1

  echo "Finished: $(date)"
  echo

done

echo "============================================================"
echo "InterProScan all-genome run complete"
date
echo "============================================================"

echo "Number of TSV outputs:"
ls "${OUT_DIR}"/*.tsv 2>/dev/null | wc -l

echo "Check for errors:"
grep -i "error\|failed\|exception\|killed\|asterisk" "${OUT_DIR}"/*.interproscan.log | head -n 100
```

# TCDB / transporter-specific annotation
```
# ============================================================
# Transporter annotation environment
# ============================================================

conda create -n tcdb -c conda-forge -c bioconda diamond seqkit -y
conda activate tcdb

TCDB_DB_DIR="/scratch/mdesmarais/databases/tcdb"
mkdir -p "${TCDB_DB_DIR}"
cd "${TCDB_DB_DIR}"

# Download all TCDB protein sequences
wget -c http://www.tcdb.org/public/tcdb -O tcdb_proteins.raw.faa

# Check it is FASTA
ls -lh tcdb_proteins.raw.faa
head -n 10 tcdb_proteins.raw.faa

# Build database
TCDB_DB_DIR="/scratch/mdesmarais/databases/tcdb"
cd "${TCDB_DB_DIR}"

awk '
  /^>/ {
    print
    next
  }
  {
    gsub(/\*/, "", $0)
    gsub(/[^A-Za-z]/, "", $0)
    print toupper($0)
  }
' tcdb_proteins.raw.faa > tcdb_proteins.clean.faa

echo "Check cleaned FASTA:"
head -n 10 tcdb_proteins.clean.faa

diamond makedb \
  --in "${TCDB_DB_DIR}/tcdb_proteins.clean.faa" \
  -d "${TCDB_DB_DIR}/tcdb_proteins"

ls -lh "${TCDB_DB_DIR}"

# Run
conda activate tcdb

BASE_DIR="/scratch/mdesmarais/JS1_comparative_genomics"
FAA_DIR="${BASE_DIR}/prodigal_out/faa"
TCDB_DB_DIR="/scratch/mdesmarais/databases/tcdb"
OUT_DIR="${BASE_DIR}/tcdb_out"

mkdir -p "${OUT_DIR}"

for faa in "${FAA_DIR}"/*.faa; do

  genome=$(basename "${faa}" .faa)

  echo "============================================================"
  echo "Running TCDB transporter search: ${genome}"
  echo "============================================================"

  diamond blastp -q "${faa}" -d "${TCDB_DB_DIR}/tcdb_proteins" -o "${OUT_DIR}/${genome}.tcdb.diamond.tsv" --threads 8 --evalue 1e-20 --query-cover 50 --subject-cover 40 --max-target-seqs 5 --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore stitle > "${OUT_DIR}/${genome}.tcdb.stdout.txt" 2> "${OUT_DIR}/${genome}.tcdb.stderr.txt"

done

echo "Done."
ls -lh "${OUT_DIR}" | head
```





