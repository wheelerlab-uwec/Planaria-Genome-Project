#!/bin/bash
#SBATCH --job-name=annot_upload
#SBATCH --partition=week
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=annot_upload-%j.out
# Upload the annotated genome (.sqn from snakemake rule annot_table2asn) to the
# NCBI preload area, for the portal FIX of WGS submission SUB16566354.
# Usage: sbatch upload_annotation.sbatch <account folder> <subfolder>
#   <account folder>: the "uploads/..." folder shown in the submission wizard
#   <subfolder>:      any new name, e.g. gdor_ma-c2_annotation
# Credentials come from ~/.netrc, as for upload_genome.sbatch.
set -eo pipefail
FOLDER=${1:?account folder (uploads/...) required}
SUB=${2:?subfolder name required}
SQN=/data/users/wheelenj/GitHub/Planaria-Genome-Project/snakemake/results/annotation/ncbi/Gdor.sqn
DIR=/data/users/wheelenj/GitHub/Planaria-Genome-Project/ncbi_submission/genome
F=Gdor_MA-C2_annotated.sqn

[ "$(stat -c %a ~/.netrc)" = 600 ] || { echo "~/.netrc must be chmod 600"; exit 1; }
cd $DIR
gzip -c $SQN > $F.gz
md5sum $F.gz | tee $F.gz.md5
echo "$(date) uploading $F.gz"
curl --netrc --ftp-create-dirs --retry 5 --retry-delay 60 -C - \
     -T $F.gz "ftp://ftp-private.ncbi.nlm.nih.gov/$FOLDER/$SUB/$F.gz"
echo "$(date) listing:"
curl --netrc -s --list-only "ftp://ftp-private.ncbi.nlm.nih.gov/$FOLDER/$SUB/"
