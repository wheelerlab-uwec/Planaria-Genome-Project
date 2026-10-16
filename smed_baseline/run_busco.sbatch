#!/bin/bash
#SBATCH --job-name=smed_busco
#SBATCH --partition=week
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=smed_busco-%j.out
# Baseline for the internal-stop-codon rate: run the same BUSCO setup as the
# pipeline (BUSCO 6, metazoa_odb10, miniprot) on S. mediterranea references.
# If they show a similar ~33% rate, the stop codons are a method artifact on
# planarian genes rather than errors in our assembly.
set -eo pipefail
module load python-libs
source activate planaria-env

OUT=/data/users/wheelenj/GitHub/Planaria-Genome-Project/smed_baseline
# run from the snakemake dir to reuse busco_downloads
cd /data/users/wheelenj/GitHub/Planaria-Genome-Project/snakemake
for acc in GCA_045838265.1 GCA_002600895.1; do
    busco -i $OUT/$acc.fna -l metazoa_odb10 -m genome --offline \
          -o busco_$acc --out_path $OUT -f --cpu 32 > $OUT/busco_$acc.log 2>&1
done
grep -H -A2 "C:" $OUT/busco_*/short_summary*.txt
