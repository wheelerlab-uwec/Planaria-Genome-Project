#!/bin/bash
#SBATCH --job-name=combine_planaria
#SBATCH --output=combine_%j.log
#SBATCH --error=combine_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=week

# Load modules if needed
module load python-libs
source $(conda info --base)/etc/profile.d/conda.sh
conda activate biol343

# Set variables
FC1="/data/groups/wheelenj/sequencing/20250916_M009242/4650P/Planaria1/20250909_1219_3B_PBE82635_fbdcc6f4/fastq_pass/*.fastq.gz"
FC2="/data/groups/wheelenj/sequencing/20250916_M009242/4651P/Planaria2/20250909_1354_3C_PBE69031_20f7cd5b/fastq_pass/*.fastq.gz"
OUTFILE="planaria_all_reads.fastq.gz"

# Run concatenation
echo "Starting concatenation at $(date)"
echo "Combining all reads from both flow cells..."

cat ${FC1}/*.fastq.gz ${FC2}/*.fastq.gz > ${OUTFILE}

echo "Finished concatenation at $(date)"
echo "Output file size:"
ls -lh ${OUTFILE}

# Quick QC check
echo "Number of reads in output:"
zcat ${OUTFILE} | wc -l | awk '{print $1/4}'