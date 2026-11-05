#!/bin/bash
#SBATCH --job-name=ploidy_af
#SBATCH --partition=week
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=ploidy_af-%j.out
# One-off ploidy / pool-mixture check: allele frequencies at SNP sites in
# collapsed (~100x) vs haplotype-level (~25x) contigs, using the reads that
# Medaka already aligned to the racon assembly.
set -eo pipefail
module load python-libs
source activate planaria-env

RES=/data/users/wheelenj/GitHub/Planaria-Genome-Project/snakemake/results/polished
BAM=$RES/medaka/calls_to_draft.bam
REF=$RES/racon/final_racon_polished.fasta
MEAN_READ_LEN=13173  # filtered reads, NanoStats

mkdir -p sites
# Rough per-contig depth from mapped read counts; contigs 100-500 kb only
samtools idxstats $BAM | awk -v rl=$MEAN_READ_LEN '$2>=100000 && $2<=500000 {printf "%s\t%d\t%.1f\n", $1, $2, $3*rl/$2}' > contig_depth.tsv
awk '$3>=85 && $3<=125' contig_depth.tsv | shuf -n 200 --random-source=<(yes) | sed 's/$/\tcollapsed/' > selected.tsv
awk '$3>=18 && $3<=35'  contig_depth.tsv | shuf -n 200 --random-source=<(yes) | sed 's/$/\thaplotype/' >> selected.tsv
cut -f4 selected.tsv | sort | uniq -c

export BAM REF
cut -f1 selected.tsv | xargs -P $SLURM_CPUS_PER_TASK -I{} bash -c '
  bcftools mpileup -X ont -I -Q 20 -q 20 -d 1000 -a FORMAT/AD -f "$REF" -r "{}" "$BAM" 2>/dev/null |
  bcftools query -f "%CHROM\t%POS\t%ALT\t[%AD]\n" |
  awk -F"\t" "\$3 != \"<*>\"" | cut -f1,2,4 > "sites/{}.tsv"'

python plot_af.py selected.tsv sites af_hist.png > af_summary.txt
cat af_summary.txt
