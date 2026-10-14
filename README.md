# Planaria-Genome-Project
RNA-seq + Genome assembly of Girardia dorotocephala.

1 million subsampled reads can be found at `/data/groups/wheelenj/sequencing/20250916_M009242/planaria_test_subset.fastq.gz`

## Build the env
`mamba env create`

## Install racon
`mamba install gtest`
`git clone https://github.com/lbcb-sci/racon && cd racon && mkdir build && cd build`
`cmake -DCMAKE_BUILD_TYPE=Release .. && make`

Run with `./racon/build/bin/racon`

## Sample provenance
- Putatively clonal *G. dorotocephala* from **[TK: collaborator]**, mucus-stripped (NAC) and flash-frozen 2025-06-25 in 1/10/20/30/50 mg aliquots (~5 mg per animal).
- HMW DNA was extracted by the UWBC from **two 50 mg tubes (20 animals)** using Grohme et al. 2018 S1.3–S1.4 (GuSCN/sarkosyl lysis + CTAB cleanup). The first attempt, with the Nanobind kit on the 10–50 mg tubes, failed and was not sequenced.
- **One** SQK-LSK114 library was run on **two** PromethION R10.4.1 flow cells (`Planaria1` = 4650P, `Planaria2` = 4651P), which are technical replicates of the same DNA. Basecalling used SUP v4.3.0 (Dorado in MinKNOW 24.11).
- The UWBC may still hold unused aliquots (3× 10 mg, 3× 20 mg, 3× 30 mg, 1× 50 mg, 4× 1 mg).

Full methods: [METHODS.md](METHODS.md). One-off ploidy check: `ploidy_check/`.
