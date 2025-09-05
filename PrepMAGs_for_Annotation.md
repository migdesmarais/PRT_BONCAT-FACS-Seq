# MAG Preparation and Annotation Workflow

This workflow outlines essential steps for preparing Metagenome-Assembled Genomes (MAGs) for annotation and downstream analysis. It includes ORF prediction, basic stats, and common annotation tools.

---

## 1. Predict Open Reading Frames (ORFs) with Prodigal

Set up environment and run Prodigal for all MAGs:

```bash
conda create --name prodigal
conda activate prodigal
conda install -c bioconda -c conda-forge prodigal seqkit
mkdir prodigal

for MAG in *.fasta; do
  BASENAME=$(basename "$MAG" .fasta)
  prodigal -i "$MAG" \
    -a "prodigal/${BASENAME}_proteins.faa" \
    -d "prodigal/${BASENAME}_nucleotides.fna" \
    -f gff \
    -o "prodigal/${BASENAME}_genes.gff" \
    -p meta
  seqkit stats "prodigal/${BASENAME}_nucleotides.fna"
  seqkit stats "prodigal/${BASENAME}_proteins.faa"
done
```

---

## 2. DRAM Annotation

> DRAM is easiest to run on KBase, Docker, or Singularity. Manual install can be difficult.

```bash
git clone https://github.com/WrightonLabCSU/DRAM.git
cd DRAM
conda env create -f environment.yaml -n DRAM
conda activate DRAM
pip install .

DRAM-setup.py prepare_databases --output_dir DRAM_data --threads 8 &> dram_setup.log

# Optionally filter proteins by minimum length:
seqkit seq -m 100 prodigal/flavo_proteins.faa > prodigal/flavo_proteins_min100.faa

DRAM.py annotate \
  -i prodigal/flavo_proteins_min100.faa \
  -o DRAM/flavo_dram_output \
  --threads 8
```

---

## 3. Eggnog-mapper Annotation

```bash
conda create --name eggnog
conda activate eggnog
conda install -c bioconda -c conda-forge eggnog-mapper
download_eggnog_data.py

emapper.py \
  -i prodigal/flavo_proteins.faa \
  -o eggnog/flavo_emapper \
  --cpu 8 \
  -m diamond
# Repeat for other MAGs as needed
```

---

## 4. Prokka Annotation

```bash
conda create --name prokka
conda activate prokka
conda install -c conda-forge -c bioconda prokka

prokka --outdir prokka/flavo --prefix flavo flavo.fasta
# Repeat for other MAGs as needed
```

---

## 5. InterProScan Annotation

```bash
conda create --name interproscan
conda activate interproscan
conda install -c bioconda -c conda-forge interproscan openjdk=11

mkdir interproscan
cd interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.74-105.0/interproscan-5.74-105.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.74-105.0/interproscan-5.74-105.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.74-105.0-64-bit.tar.gz.md5
tar -pxvzf interproscan-5.74-105.0-*-bit.tar.gz

# Clean protein FASTA headers and remove stop codons (optional):
awk '/^>/ {match($0, /ID=([^;]+)/, a); print ">" a[1]; next} {gsub(/\*$/, "", $0); print}' prodigal/flavo_proteins.faa > prodigal/flavo_clean.faa

./interproscan.sh -i prodigal/flavo_clean.faa -exclappl MobiDBLite -goterms &> flavo.log
# Repeat for other MAGs as needed
```

---

## Notes

- Adjust file paths and sample names as needed.
- Some annotation tools require significant compute/memory resources.
- Always check logs and outputs for errors.
- Modify and extend this workflow for your specific needs or additional tools.

---

## References

- [Prodigal](https://github.com/hyattpd/Prodigal)
- [DRAM](https://github.com/WrightonLabCSU/DRAM)
- [Eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)
- [Prokka](https://github.com/tseemann/prokka)
- [InterProScan](https://github.com/ebi-pf-team/interproscan)