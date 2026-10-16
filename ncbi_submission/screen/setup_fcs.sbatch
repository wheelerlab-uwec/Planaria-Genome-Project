#!/bin/bash
#SBATCH --job-name=fcs_setup
#SBATCH --partition=week
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=fcs_setup-%j.out
# Fetch NCBI FCS (v0.5.5) containers and the FCS-GX database (~470 GiB).
set -eo pipefail
module load apptainer python-libs
source activate planaria-env  # fcs.py needs Python >= 3.7
FCS=/data/users/wheelenj/software/fcs
mkdir -p $FCS && cd $FCS

curl -sSfLO https://github.com/ncbi/fcs/raw/v0.5.5/dist/fcs.py
curl -sSfLO https://github.com/ncbi/fcs/raw/v0.5.5/dist/run_fcsadaptor.sh
chmod +x run_fcsadaptor.sh
curl -sSfL -o fcs-gx.sif https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.5.5/fcs-gx.sif
curl -sSfL -o fcs-adaptor.sif https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.5.5/fcs-adaptor.sif

export FCS_DEFAULT_IMAGE=$FCS/fcs-gx.sif
python3 fcs.py db get --mft https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest --dir $FCS/gxdb
python3 fcs.py db check --mft https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest --dir $FCS/gxdb
du -sh $FCS/gxdb
