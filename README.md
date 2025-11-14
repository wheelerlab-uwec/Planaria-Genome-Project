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