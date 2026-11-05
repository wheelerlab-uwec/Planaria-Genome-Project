#!/bin/bash
#SBATCH --job-name=sra_prep
#SBATCH --partition=week
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=sra_prep-%j.out
# Build one FASTQ per flow cell for SRA from the MinKNOW fastq_pass chunks.
# 1) verify every chunk against the core's md5sum.txt
# 2) concatenate chunks (gzip members concatenate into a valid .fastq.gz)
# 3) record md5, size and read stats of the upload files
set -eo pipefail
module load python-libs
source activate planaria-env

RAW=/data/groups/wheelenj/sequencing/20250916_M009242
OUT=/data/users/wheelenj/ncbi_upload/sra
mkdir -p $OUT

for run in 4650P/Planaria1/20250909_1219_3B_PBE82635_fbdcc6f4 \
           4651P/Planaria2/20250909_1354_3C_PBE69031_20f7cd5b; do
    fc=$(basename $run | cut -d_ -f4)
    echo "== $fc"
    (cd $RAW/$run && grep "fastq_pass/" md5sum.txt | md5sum -c --quiet) \
        && echo "$fc: all fastq_pass chunks match md5sum.txt"
    ls $RAW/$run/fastq_pass/*.fastq.gz | sort -V | xargs cat > $OUT/$fc.fastq.gz
done

cd $OUT
md5sum PBE82635.fastq.gz PBE69031.fastq.gz | tee md5sums.txt
seqkit stats -a -T -j 8 PBE82635.fastq.gz PBE69031.fastq.gz | tee read_stats.tsv
