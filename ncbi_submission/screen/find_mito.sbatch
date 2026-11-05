#!/bin/bash
#SBATCH --job-name=find_mito
#SBATCH --partition=week
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=find_mito-%j.out
# Locate mitochondrial contigs by aligning published Girardia mitogenomes
# (G. tigrina MW972220; Girardia sp. KP090061) to the primary assembly and
# the purged haplotigs. Read depth (from purge_dups' PB.base.cov) is
# reported for each hit, since mtDNA should be far above the ~100x nuclear peak.
set -eo pipefail
module load python-libs
source activate planaria-env
FIN=/data/users/wheelenj/GitHub/Planaria-Genome-Project/snakemake/results/final
OUT=/data/users/wheelenj/GitHub/Planaria-Genome-Project/ncbi_submission/screen/results/mito
mkdir -p $OUT && cd $OUT

curl -sSf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MW972220,KP090061&rettype=fasta&retmode=text" > girardia_mito.fa
grep ">" girardia_mito.fa

for asm in purged_assembly haplotigs; do
    minimap2 -x asm20 -t 16 $FIN/$asm.fasta girardia_mito.fa > mito_vs_$asm.paf 2> /dev/null
done
# query, qlen, qstart, qend, target, tlen, tstart, tend, matches, alnlen
cat mito_vs_*.paf | awk -v OFS="\t" '$11>=1000 {print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11}' \
    | sort -k5,5 -k7,7n > mito_hits.tsv
cut -f5,6 mito_hits.tsv | sort -u > mito_contigs.txt
cat mito_contigs.txt
# mean read depth of each hit contig
awk 'NR==FNR{k[$1];next} /^>/{c=substr($1,2); next} (c in k){s[c]+=$1; n[c]++} END{for(x in s) printf "%s\t%.0f\n", x, s[x]/n[x]}' \
    <(cut -f1 mito_contigs.txt) $FIN/purge_dups/PB.base.cov | sort -k2,2nr | tee mito_contig_depth.tsv
