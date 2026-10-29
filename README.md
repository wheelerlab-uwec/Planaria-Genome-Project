# Planaria-Genome-Project
Genome assembly and annotation of the planarian *Girardia dorotocephala* (clonal line MA-C2), from Oxford Nanopore reads, with RNA-seq evidence for gene prediction.

1 million subsampled reads can be found at `/data/groups/wheelenj/sequencing/20250916_M009242/planaria_test_subset.fastq.gz`

## Status
- Nuclear assembly: 1.47 Gb, 21,035 contigs, N50 131.7 kb. Submitted as WGS JCFHDI000000000.
- Annotation: 21,116 protein-coding genes and 2,245 tRNA genes. Added to the WGS submission on 2026-10-29.
- Mitochondrial genome: 22,455 bp, GenBank QB207629.
- Manuscript: drafting for a GigaByte Data Release ([METHODS.md](METHODS.md); checklist in its internal notes).

## Workflow
Everything runs from `snakemake/` with Snakemake on the BOSE Slurm cluster:

```
cd snakemake
snakemake --profile executor/            # everything in rule all
snakemake --profile executor/ results/annotation/ncbi/Gdor.sqn   # one target
```

| Rules | Does |
|---|---|
| `Snakefile` | read QC and filtering, Flye, Racon, Medaka, purge_dups, BUSCO, reports |
| `rules/mito.smk` | mitochondrial genome assembly, annotation and `.sqn` |
| `rules/annotation.smk` | repeats (RepeatModeler, RepeatMasker), RNA-seq download and STAR alignment |
| `rules/annotation_genes.smk` | BRAKER3, BUSCO and OMArk on the proteins, eggNOG-mapper, InterProScan |
| `rules/annotation_ncbi.smk` | product names, tRNAscan-SE, locus tags, `table2asn` |

Databases (OrthoDB, OMA, eggNOG, InterProScan, Swiss-Prot names, `table2asn`) are downloaded by rules into `snakemake/resources/`, at the versions pinned in `snakemake/config.yaml`.

## Environments
- `planaria-env` (`environment.yml`): `mamba env create`. Most rules use it.
- `braker-env` (`snakemake/envs/braker.yaml`): BRAKER3. GeneMark-ETP is not on bioconda; the file explains how it is cloned into the env.
- `omark-env` (`snakemake/envs/omark.yaml`): OMArk 0.5.0, which does not solve into `planaria-env`.

### Install racon
`mamba install gtest`
`git clone https://github.com/lbcb-sci/racon && cd racon && mkdir build && cd build`
`cmake -DCMAKE_BUILD_TYPE=Release .. && make`

The workflow runs it from `./racon/build/bin/racon`.

## Sample provenance
- *G. dorotocephala* asexual clonal line MA-C2 (2n = 24; Almazan et al. 2018) from Labib Rouhana (UMass Boston), mucus-stripped (NAC) and flash-frozen 2025-06-25 in 1/10/20/30/50 mg aliquots (~5 mg per animal).
- HMW DNA was extracted by the UWBC from **two 50 mg tubes (20 animals)** using Grohme et al. 2018 S1.3–S1.4 (GuSCN/sarkosyl lysis + CTAB cleanup). The first attempt, with the Nanobind kit on the 10–50 mg tubes, failed and was not sequenced.
- **One** SQK-LSK114 library was run on **two** PromethION R10.4.1 flow cells (`Planaria1` = 4650P, `Planaria2` = 4651P), which are technical replicates of the same DNA. Basecalling used SUP v4.3.0 (Dorado in MinKNOW 24.11).
- The UWBC may still hold unused aliquots (3× 10 mg, 3× 20 mg, 3× 30 mg, 1× 50 mg, 4× 1 mg).
- RNA-seq for annotation: four MA-C2 runs from BioProject PRJNA317859 (Almazan et al. 2018).

## NCBI
BioProject **PRJNA1545025**, BioSample **SAMN64250158**, locus_tag prefix **ADAYWY**. Genome assembly: WGS **JCFHDI000000000** (this version JCFHDI010000000), with annotation. Mitochondrial genome: GenBank **QB207629**. SRA accessions are pending. All records are private until 2027-10-16 or publication. Submission steps: [ncbi_submission/README.md](ncbi_submission/README.md).

Full methods: [METHODS.md](METHODS.md). One-off ploidy check: `ploidy_check/`.
