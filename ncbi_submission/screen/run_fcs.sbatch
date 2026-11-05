#!/bin/bash
#SBATCH --job-name=fcs_screen
#SBATCH --partition=highmemory
#SBATCH --cpus-per-task=48
#SBATCH --mem=600G
#SBATCH --time=1-00:00:00
#SBATCH --output=fcs_screen-%j.out
# Screen the purged primary assembly as NCBI will: FCS-adaptor (vector/adaptor)
# and FCS-GX (foreign organisms), with taxid 135777 (G. dorotocephala).
# Reports only; nothing is removed here.
set -eo pipefail
module load apptainer python-libs
source activate planaria-env  # fcs.py needs Python >= 3.7
FCS=/data/users/wheelenj/software/fcs
ASM=/data/users/wheelenj/GitHub/Planaria-Genome-Project/snakemake/results/final/purged_assembly.fasta
OUT=/data/users/wheelenj/GitHub/Planaria-Genome-Project/ncbi_submission/screen/results
mkdir -p $OUT/adaptor $OUT/gx

$FCS/run_fcsadaptor.sh --fasta-input $ASM --output-dir $OUT/adaptor --euk \
    --container-engine singularity --image $FCS/fcs-adaptor.sif

export FCS_DEFAULT_IMAGE=$FCS/fcs-gx.sif
export GX_NUM_CORES=48
python3 $FCS/fcs.py screen genome --fasta $ASM --out-dir $OUT/gx \
    --gx-db $FCS/gxdb/all --tax-id 135777

ls -la $OUT/adaptor $OUT/gx
