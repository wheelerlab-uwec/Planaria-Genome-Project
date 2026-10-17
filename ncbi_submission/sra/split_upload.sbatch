#!/bin/bash
#SBATCH --job-name=sra_split
#SBATCH --partition=week
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=sra_split-%j.out
# SRA rejects files over 100 GiB, and PBE69031.fastq.gz is 100.3 GiB.
# Split it into two halves by read count, upload both to the preload folder,
# then move the oversized original out of that folder (an FTP rename; the
# original stays on BOSE).
# Usage: sbatch split_upload.sbatch <account folder> <subfolder>
# Credentials come from ~/.netrc, as for upload_reads.sbatch.
set -eo pipefail
module load python-libs
source activate planaria-env
FOLDER=${1:?account folder (uploads/...) required}
SUB=${2:?subfolder name required}
SRC=/data/users/wheelenj/ncbi_upload/sra
FTP=ftp://ftp-private.ncbi.nlm.nih.gov/$FOLDER

[ "$(stat -c %a ~/.netrc)" = 600 ] || { echo "~/.netrc must be chmod 600"; exit 1; }
cd $SRC
echo "$(date) splitting PBE69031.fastq.gz"
seqkit split2 -p 2 -j $SLURM_CPUS_PER_TASK -O split PBE69031.fastq.gz
mv split/PBE69031.part_001.fastq.gz PBE69031_part1.fastq.gz
mv split/PBE69031.part_002.fastq.gz PBE69031_part2.fastq.gz
rmdir split
seqkit stats -T -j 2 PBE69031.fastq.gz PBE69031_part1.fastq.gz PBE69031_part2.fastq.gz | tee split_stats.tsv
# the two parts must hold every read of the original
awk -F'\t' 'NR==2{n=$4; b=$5} NR>2{m+=$4; c+=$5} END{if(n!=m||b!=c){print "read/base count mismatch"; exit 1} else print "counts match:", n, "reads,", b, "bases"}' split_stats.tsv
md5sum PBE69031_part1.fastq.gz PBE69031_part2.fastq.gz | tee split_md5.txt
ls -l PBE69031_part*.fastq.gz

for f in PBE69031_part1.fastq.gz PBE69031_part2.fastq.gz; do
    echo "$(date) uploading $f"
    curl --netrc --ftp-create-dirs --retry 5 --retry-delay 60 -C - -T $f "$FTP/$SUB/$f"
done

echo "$(date) moving the oversized original out of $SUB"
curl --netrc -s -Q "MKD ${SUB}_oversize" "$FTP/" -l > /dev/null || true
curl --netrc -s -Q "RNFR $SUB/PBE69031.fastq.gz" -Q "RNTO ${SUB}_oversize/PBE69031.fastq.gz" "$FTP/" -l > /dev/null
echo "$(date) listing $SUB:"
curl --netrc -s -l "$FTP/$SUB/"
