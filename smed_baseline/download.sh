#!/bin/bash
# Download S. mediterranea reference assemblies from NCBI (run on the login
# node; compute nodes may not have internet access).
#   GCA_045838265.1  schMedS3h1  (2024; PacBio HiFi + Illumina, chromosome-level)
#   GCA_002600895.1  ASM260089v1 (Grohme et al. 2018; PacBio CLR)
set -eo pipefail
cd "$(dirname "$0")"
for acc in GCA_045838265.1 GCA_002600895.1; do
    [ -s $acc.fna ] && continue
    curl -sSfL -o $acc.zip \
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA"
    unzip -p $acc.zip "ncbi_dataset/data/$acc/*.fna" > $acc.fna
    rm $acc.zip
done
grep -c ">" *.fna
