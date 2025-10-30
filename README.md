# Planaria-Genome-Project
RNA-seq + Genome assembly of Girardia dorotocephala.

1 million subsampled reads can be found at `/data/groups/wheelenj/sequencing/20250916_M009242/planaria_test_subset.fastq.gz`

<<<<<<< HEAD
## Install the uv env
`uv sync`

## Install Filtlong
```
mkdir tools
cd tools
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j
```

## Install Flye
```
cd tools
git clone https://github.com/fenderglass/Flye
cd Flye
make
```
## Installing other tools
Use pypi.org to find tool
```
uv add 'tool'
```
=======
## Build the env
`mamba env create`
>>>>>>> 745063808ae49bedbbfe479716bae78ef3ed8e88
