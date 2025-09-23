# REMOVE CONTAMINANTS FROM BRACKEN RESULTS
```
screen
export DB=/scratch/mdesmarais/kraken_gtdb_r226_128g
echo "$DB"    # sanity check
URL=https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz
wget -c -O "$DB/gtdb_genomes_reps.tar.gz" "$URL" --progress=dot:giga

DB=/scratch/mdesmarais/kraken_gtdb_r226_128g
k2 build --db "$DB" --special gtdb \
  --gtdb-server data.gtdb.ecogenomic.org \
  --gtdb-files 'https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz' \
  --max-db-size 128GiB \
  --threads 24 \
  --log "$DB/build.log"

