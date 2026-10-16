#!/bin/bash
#SBATCH --job-name=polish_test
#SBATCH --partition=GPU
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=1-12:00:00
#SBATCH --output=polish_test-%j.out
# One-off comparison of Dorado polish vs Medaka on the same input:
# the 3x Racon assembly and the reads Medaka already aligned to it.
# Only contigs carrying a BUSCO hit in the Medaka assembly are compared
# (contigs.txt), so BUSCO on each subset is directly comparable.
set -eo pipefail
module load python-libs
source activate planaria-env

PROJ=/data/users/wheelenj/GitHub/Planaria-Genome-Project
SMK=$PROJ/snakemake
OUT=$PROJ/polish_test
DORADO=/data/users/wheelenj/software/dorado-2.1.2-linux-x64/bin/dorado
MODEL=/data/users/wheelenj/software/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0_polish
BAM=$SMK/results/polished/medaka/calls_to_draft.bam
RACON=$SMK/results/polished/racon/final_racon_polished.fasta
MEDAKA=$SMK/results/polished/medaka/consensus.fasta

cd $OUT
nvidia-smi
awk 'NR==FNR{k[$1];next} ($1 in k){print $1"\t0\t"$2}' contigs.txt $RACON.fai > contigs.bed

# Dorado polish only accepts BAMs from dorado aligner, so pull the reads
# Medaka aligned to the test contigs and realign them to the whole Racon
# assembly (reads can still land on their best contig anywhere).
samtools view -u -@ 8 -L contigs.bed $BAM | samtools fastq -@ 8 - > test_reads.fastq 2> fastq.log
/usr/bin/time -v $DORADO aligner -t 32 $RACON test_reads.fastq 2> aligner.log \
    | samtools sort -@ 16 -m 2G -o dorado_aligned.bam -
samtools index dorado_aligned.bam

# Dorado polish (reads are basecalled with sup@v4.3.0; the BAM has no
# read groups, so the matching model is given explicitly)
/usr/bin/time -v $DORADO polish dorado_aligned.bam $RACON \
    --regions contigs.bed \
    --model-override $MODEL \
    --ignore-read-groups \
    -x cuda:0 -t 32 \
    -o dorado 2> dorado.log

# Same contigs from each assembly
seqkit grep -f contigs.txt $RACON > racon.subset.fasta
seqkit grep -f contigs.txt $MEDAKA > medaka.subset.fasta
seqkit grep -r -f <(sed 's/^/^/; s/$/([:_ ]|$)/' contigs.txt) dorado/consensus.fasta > dorado.subset.fasta
seqkit stats -a -T racon.subset.fasta medaka.subset.fasta dorado.subset.fasta > subset_stats.tsv

# BUSCO on each subset (run from the snakemake dir to reuse busco_downloads)
cd $SMK
for s in racon medaka dorado; do
    busco -i $OUT/$s.subset.fasta -l metazoa_odb10 -m genome --offline \
          -o busco_$s --out_path $OUT -f --cpu 32 > $OUT/busco_$s.log 2>&1
done
grep -H "C:" $OUT/busco_*/short_summary*.txt
