# ------------------------------
# Run CheckM on your cleaned MAGs
# ------------------------------
# Folder expected:
#   /scratch/mdesmarais/Atribacterota_tree/genomes_clean/*.fna
#
# Output:
#   /scratch/mdesmarais/Atribacterota_tree/qc/checkm1/
#   /scratch/mdesmarais/Atribacterota_tree/qc/checkm1/qa.tsv  (summary table)
#
# Notes:
# - CheckM can take hours depending on number of genomes + threads.
# - For phylogenomics, a common “paper” filter is completeness >= 70% and contamination <= 10%
#   (you can be looser for a supplement, e.g., >=50% and <=15%).
# - One genome per file; headers can be simple (>contig_1 ...) and no spaces is safest.

# 1) activate env
```
conda activate checkm_env
which checkm
```
# 2) run CheckM
```
cd /scratch/mdesmarais/Atribacterota_tree
mkdir -p qc/checkm1

checkm lineage_wf \
  -x fna \
  -t 24 \
  genomes_clean \
  qc/checkm1
```

# 3) write a simple QC table (tab-delimited)
```
# if you're not already there:
cd /scratch/mdesmarais/Atribacterota_tree/qc/checkm1

# write a tab-delimited summary table
checkm qa lineage.ms . -o 2 > qa.tsv

# quick look
head -n 25 qa.tsv
```

###############################################################################
# Assign taxonomy with GTDB-Tk after CheckM
# Project: Atribacterota_tree
#
# Assumptions:
# - Your genome FASTAs are one-per-file in:
#     /scratch/mdesmarais/Atribacterota_tree/genomes_clean/*.fna
# - You already ran CheckM (optional, but recommended).
#
# Notes:
# - GTDB-Tk needs a GTDB release database. Many clusters keep this in a shared
#   path and you point to it with GTDBTK_DATA_PATH.
# - Output taxonomy summary files will be in:
#     qc/gtdbtk/
###############################################################################

# 1) activate the conda env that has gtdbtk (replace with your actual env name)
#    (Do NOT use checkm_env unless gtdbtk is installed there.)
```
conda activate gtdbtk_env
which gtdbtk
gtdbtk --version
```

# 2) (only if needed) point GTDB-Tk to the database
#    If your system already has GTDBTK_DATA_PATH set, you can skip this.
```
export GTDBTK_DATA_PATH=/data_store/gtdbtk_db/release226
echo "GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"
```

# 3) run GTDB-Tk classification
```
mkdir -p qc/gtdbtk

gtdbtk classify_wf \
  --genome_dir genomes_clean \
  --out_dir qc/gtdbtk \
  --extension fna \
  --cpus 24
```

# 4) write a SMALL tab-delimited summary (easy to open in Excel)
#    Keeps only the most useful columns for manual filtering.
```
python3 - <<'PY'
import csv

in_tsv  = "qc/gtdbtk/gtdbtk.bac120.summary.tsv"
out_tsv = "qc/gtdbtk/gtdbtk_taxonomy_summary.tsv"

with open(in_tsv, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    # columns we want (only keep those that actually exist)
    wanted = ["user_genome", "classification", "warnings", "note"]
    wanted = [c for c in wanted if c in reader.fieldnames]

    with open(out_tsv, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=["genome_id"] + [c for c in wanted if c != "user_genome"],
                           delimiter="\t")
        w.writeheader()
        for row in reader:
            genome_id = row.get("user_genome", "")
            out_row = {"genome_id": genome_id}
            for c in wanted:
                if c == "user_genome": 
                    continue
                out_row[c] = row.get(c, "")
            w.writerow(out_row)

print(f"Wrote: {out_tsv}")
PY
```

# =============================================================================
# Step 3 — Build a phylogenomic tree from the filtered genomes (GToTree → IQ-TREE)
# =============================================================================
# Assumptions:
# - You have already created a folder of genomes you want in the tree:
#     /scratch/mdesmarais/Atribacterota_tree/genomes_filtered/
#   (i.e., after CheckM QC + GTDB-Tk phylum filtering / manual curation)
#
# Outputs:
#   phylo/gtotree_out/                 (marker extraction + alignments)
#   phylo/atrib_iqtree.treefile        (final ML tree)
#
# Notes:
# - GToTree builds a concatenated marker alignment (protein) using MAFFT + trimAl.
# - IQ-TREE then infers a maximum-likelihood phylogeny with bootstrap support.
# - After this, upload the .treefile to iTOL and annotate with metadata rings.
# =============================================================================

# 1) activate tree-building env
```
conda activate phylo_env
which iqtree
which mafft
which trimal
which GToTree
```

# 2) Run GToTree to extract/align/trim/concatenate conserved markers in best-hit mode.
#    Marker set: "Bacteria" (GToTree built-in HMM set).
#    If your install uses different names, run: GToTree -h
```
# Make an explicit list of genome FASTAs (absolute paths are safest), fix naming _ :, and remove MAGs that are only gaps.
# (If you have mixed extensions, see the commented "mixed ext" option below.)
cd /scratch/mdesmarais/Atribacterota_tree/genomes_filtered_sediment

# replace ":" with "_" in filenames (safe; only affects names, not sequences)
for f in *:*; do
  mv "$f" "${f//:/_}"
done

cd /scratch/mdesmarais/Atribacterota_tree
conda activate phylo_env

ls -1 "$PWD"/Sediments_all/*.fna > Sediments_all/sediments_all.list

grep -v -E "Amon_Mud_Volcano_CSMAG-886_1250mbsl_N_Acmbsf_GCA_030666935\.1|Laptev_Sea_Cold_Seep_CSMAG-808_N_Ambsl_14-18cmbsf_GCA_030668465\.1|Scotian_Basin_Oil_Gas_Seep_CSMAG-2406_2306mbsl_60cmbsf_GCA_030612975\.1" \
  phylo_sed_only/genomes_filtered_sediment.list > phylo_sed_only/genomes_filtered_sediment.no_allgap.list

wc -l phylo_sed_only/genomes_filtered_sediment.no_allgap.list

# Run GToTree in best-hit mode, FORCE overwriting output directory if it exists
rm -rf phylo/gtotree_out_besthit_no_allgap_sediment

GToTree \
  -f Sediments_all/sediments_all.list \
  -H Bacteria \
  -B \
  -o Sediments_all/gtotree_out_besthit_sediments_all \
  -t 24
```


GToTree \
  -f genomes_filtered_sediment.list \
  -H Bacteria \
  -B \
  -o gtotree_out_besthit \
  -t 24
  
GToTree \
  -f Sediments_all/sediments_all.list \
  -H Bacteria \
  -o Sediments_all/gtotree_out_sediments_all_noBestHit \
  -t 24

# 4) Infer the final maximum-likelihood tree with IQ-TREE
#    -m MFP  : ModelFinder picks best-fit amino acid model
#    -B 1000 : ultrafast bootstrap support
#    -T 24   : threads
#    Output: phylo/atrib_iqtree.treefile
```
cd /scratch/mdesmarais/Atribacterota_tree

iqtree \
  -s Sediments_all/gtotree_out_sediments_all_noBestHit/Aligned_SCGs.faa \
  -m MFP \
  -B 1000 \
  -T 24 \
  --prefix Sediments_all/atrib_iqtree
```



iqtree \
  -s Aligned_SCGs.faa \
  -m MFP \
  -B 1000 \
  -T 24 \
  --prefix KB1_iqtree

