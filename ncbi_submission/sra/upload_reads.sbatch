#!/bin/bash
#SBATCH --job-name=sra_upload
#SBATCH --partition=week
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH --output=sra_upload-%j.out
# Upload the per-flow-cell FASTQs to the NCBI SRA preload area.
# Usage: sbatch upload_reads.sbatch <account folder> <subfolder>
#   <account folder>: the "uploads/..." folder shown in the SRA submission
#                     wizard (FTP preload option)
#   <subfolder>:      any new name, e.g. gdor_ma-c2_ont
# Credentials: put the FTP username/password from the wizard in ~/.netrc
# (chmod 600), as
#   machine ftp-private.ncbi.nlm.nih.gov login subftp password <password>
# so they never appear in this script, the command line or the log.
set -eo pipefail
FOLDER=${1:?account folder (uploads/...) required}
SUB=${2:?subfolder name required}
SRC=/data/users/wheelenj/ncbi_upload/sra

[ "$(stat -c %a ~/.netrc)" = 600 ] || { echo "~/.netrc must be chmod 600"; exit 1; }
for f in PBE82635.fastq.gz PBE69031.fastq.gz; do
    echo "$(date) uploading $f"
    curl --netrc --ftp-create-dirs --retry 5 --retry-delay 60 -C - \
         -T $SRC/$f "ftp://ftp-private.ncbi.nlm.nih.gov/$FOLDER/$SUB/$f"
done
echo "$(date) listing:"
curl --netrc -s --list-only "ftp://ftp-private.ncbi.nlm.nih.gov/$FOLDER/$SUB/"
