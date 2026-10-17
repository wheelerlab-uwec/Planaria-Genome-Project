#!/bin/bash
#SBATCH --job-name=fcs_clean
#SBATCH --partition=highmemory
#SBATCH --cpus-per-task=48
#SBATCH --mem=600G
#SBATCH --time=1-00:00:00
#SBATCH --output=fcs_clean-%j.out
# Apply the FCS results (run_fcs.sbatch) to the purged primary assembly (dropping
# sequences under 1 kb, the INSDC minimum from 2027) and
# rescreen the cleaned file with FCS-adaptor and FCS-GX to confirm it is clean.
set -eo pipefail
module load apptainer python-libs
source activate planaria-env
FCS=/data/users/wheelenj/software/fcs
REPO=/data/users/wheelenj/GitHub/Planaria-Genome-Project
ASM=$REPO/snakemake/results/final/purged_assembly.fasta
SCR=$REPO/ncbi_submission/screen
RES=$SCR/results
OUT=$REPO/ncbi_submission/genome
mkdir -p $OUT $RES/recheck/adaptor $RES/recheck/gx

python3 $SCR/clean_assembly.py --fasta $ASM \
    --gx-report $RES/gx/purged_assembly.135777.fcs_gx_report.txt \
    --adaptor-report $RES/adaptor/fcs_adaptor_report.txt \
    --skipped-trims $RES/adaptor/skipped_trims.jsonl --min-len 1000 \
    --out $OUT/Gdor_MA-C2_primary.fsa --contam-out $OUT/removed_contaminants.fasta \
    --log $OUT/cleaning_log.tsv
head -3 $OUT/cleaning_log.tsv
cut -f1 $OUT/cleaning_log.tsv | grep -v "^#" | sort | uniq -c

seqkit stats -a -T $ASM $OUT/Gdor_MA-C2_primary.fsa > $OUT/assembly_stats.tsv
cat $OUT/assembly_stats.tsv

# rescreen
$FCS/run_fcsadaptor.sh --fasta-input $OUT/Gdor_MA-C2_primary.fsa --output-dir $RES/recheck/adaptor --euk \
    --container-engine singularity --image $FCS/fcs-adaptor.sif
export FCS_DEFAULT_IMAGE=$FCS/fcs-gx.sif
export GX_NUM_CORES=48
python3 $FCS/fcs.py screen genome --fasta $OUT/Gdor_MA-C2_primary.fsa --out-dir $RES/recheck/gx \
    --gx-db $FCS/gxdb/all --tax-id 135777

echo "== recheck: adaptor calls"; grep -vc "^#" $RES/recheck/adaptor/fcs_adaptor_report.txt || true
echo "== recheck: skipped trims"; wc -l < $RES/recheck/adaptor/skipped_trims.jsonl
echo "== recheck: GX calls"; grep -v "^#" $RES/recheck/gx/*.fcs_gx_report.txt || true
md5sum $OUT/Gdor_MA-C2_primary.fsa > $OUT/Gdor_MA-C2_primary.fsa.md5
