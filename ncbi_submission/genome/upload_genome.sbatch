#!/bin/bash
#SBATCH --job-name=genome_upload
#SBATCH --partition=week
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=genome_upload-%j.out
# Compress the cleaned primary assembly (screen/clean_fcs.sbatch) and upload it
# to the NCBI preload area for the Genome submission portal.
# Usage: sbatch upload_genome.sbatch <account folder> <subfolder>
#   <account folder>: the "uploads/..." folder shown in the submission wizard
#   <subfolder>:      any new name, e.g. gdor_ma-c2_genome
# Credentials come from ~/.netrc, as for sra/upload_reads.sbatch.
set -eo pipefail
FOLDER=${1:?account folder (uploads/...) required}
SUB=${2:?subfolder name required}
DIR=/data/users/wheelenj/GitHub/Planaria-Genome-Project/ncbi_submission/genome
F=Gdor_MA-C2_primary.fsa

[ "$(stat -c %a ~/.netrc)" = 600 ] || { echo "~/.netrc must be chmod 600"; exit 1; }
cd $DIR
md5sum -c $F.md5
gzip -c $F > $F.gz
md5sum $F.gz | tee $F.gz.md5
echo "$(date) uploading $F.gz"
curl --netrc --ftp-create-dirs --retry 5 --retry-delay 60 -C - \
     -T $F.gz "ftp://ftp-private.ncbi.nlm.nih.gov/$FOLDER/$SUB/$F.gz"
echo "$(date) listing:"
curl --netrc -s --list-only "ftp://ftp-private.ncbi.nlm.nih.gov/$FOLDER/$SUB/"
