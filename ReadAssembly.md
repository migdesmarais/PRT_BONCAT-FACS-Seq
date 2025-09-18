## Assemble BONCAT-FACS-SEQ reads to hopefully recover some partial MAGs since we're getting very low mapping of reads to in stiu MAGs from Corinna.

## Make script
```
nano run_megahit_all.sh
```

## Script:
```
for r1 in *paired_R1.fastq.gz; do
  base=$(echo "$r1" | sed 's/paired_R1.fastq.gz//')
  r2="${base}_paired_R2.fastq.gz"
  outdir="megahit/${base}"

  echo "Assembling $base..."
  megahit \
    -1 "$r1" \
    -2 "$r2" \
    -o "$outdir" \
    --min-contig-len 1000 \
    --presets meta-sensitive \
    --num-cpu-threads 20
done
```
## Make it executable
```
chmod +x run_megahit_all.sh
```
## Run it
```
conda activate megahit
mkdir -p megahit
./run_megahit_all.sh
```

## Check assembly quality with quast
```
conda create --name quast_env
conda activate quast_env
conda install -c bioconda quast
quast --version
```

```
quast -o QC *_contigs.fa
scp -r migdesmarais@fram.ucsd.edu:/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/megahit_assemblies/contigs/QC /Users/migueldesmarais/Downloads
```

