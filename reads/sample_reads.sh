#!/bin/bash
#SBATCH --job-name=subsample_planaria
#SBATCH --output=subsample_%j.log
#SBATCH --error=subsample_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=week

# Load modules if needed
module load python-libs
source $(conda info --base)/etc/profile.d/conda.sh
conda activate biol343

# Set variables
FC1="/data/groups/wheelenj/sequencing/20250916_M009242/4650P/Planaria1/20250909_1219_3B_PBE82635_fbdcc6f4/fastq_pass"
FC2="/data/groups/wheelenj/sequencing/20250916_M009242/4651P/Planaria2/20250909_1354_3C_PBE69031_20f7cd5b/fastq_pass"
OUTFILE="planaria_test_subset.fastq.gz"
NREADS=1000000
SEED=100

# Run sampling
echo "Starting subsampling at $(date)"
cat ${FC1}/*.fastq.gz ${FC2}/*.fastq.gz | \
    pigz -dc -p ${SLURM_CPUS_PER_TASK} | \
    seqtk sample -s ${SEED} - ${NREADS} | \
    pigz -p ${SLURM_CPUS_PER_TASK} > ${OUTFILE}

echo "Finished subsampling at $(date)"
echo "Output file size:"
ls -lh ${OUTFILE}

# Quick QC check
echo "Number of reads in output:"
zcat ${OUTFILE} | wc -l | awk '{print $1/4}'
