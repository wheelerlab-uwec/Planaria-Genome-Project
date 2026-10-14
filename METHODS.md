# Methods: *Girardia dorotocephala* genome assembly

Draft methods covering work completed through October 2026. Items marked **[TK]** still need a source. Steps marked *(in progress)* were running when this was written. The "Notes for revision" section at the end is for internal use, not for publication.

## Animals and tissue preparation

*Girardia dorotocephala* from a putatively clonal (asexual) laboratory line were provided by **[TK: collaborator, institution, line designation; the core's sample sheet calls the samples "MA-C2 Planaria Clonal Samples"]**. Animals arrived on approximately 5 June 2025. They were held in the shipping water for about one week to acclimate, and any unhealthy animals were removed. They were then kept in artificial pond water (APW) at **[TK: temperature]**. Animals were last fed on 20 June 2025 and starved for 5 days before tissue collection.

On 25 June 2025, mucus was removed with N-acetyl-L-cysteine (NAC), following Grohme et al. (2018). The stripping solution was 0.5% (w/v) NAC in 20 mM HEPES-NaOH (pH 7.25–7.4), with phenol red as a pH indicator, adjusted to pH > 7 with 1 M NaOH. Animals were inverted in the solution every 10–15 s for 10–15 min at room temperature and then rinsed in APW. They were divided into aliquots by approximate mass (~5 mg per animal):

| Aliquot | Animals per tube | Tubes |
|---|---|---|
| 10 mg | 2 | 4 |
| 20 mg | 4 | 4 |
| 30 mg | 6 | 4 |
| 50 mg | 10 | 4 |
| 1 mg | ¼ of one animal | 4 |

Residual water was removed. Tubes were flash-frozen in liquid nitrogen, stored at −80 °C, and shipped on dry ice to the University of Wisconsin Biotechnology Center (UWBC) DNA Sequencing Facility (RRID:SCR_017759).

## High-molecular-weight DNA extraction

All extractions were done by the UWBC. The first attempt (17 July 2025) used the Nanobind Tissue kit (PacBio) on one tube each of the 10, 20, 30 and 50 mg aliquots, with the eluates pooled. Lysates were not viscous and the yield was low (~6.6 ng/µl, ~660 ng total), so this extract was not used.

The second attempt followed the planarian-specific protocol of Grohme et al. (2018; Supplementary Methods S1.3–S1.4). It used two 50 mg aliquots combined (**20 animals**). S1.3 is a guanidinium thiocyanate / N-lauroylsarcosine lysis. S1.4 is a CTAB polysaccharide cleanup with RNase treatment. Lysis began 25 August 2025.

The final extract measured 187 ng/µl by NanoDrop and 140 ng/µl by Qubit. Femto Pulse showed a modal fragment length of about 20 kb, with some fragments over 165 kb. **All sequence data come from this single 20-animal extract.**

## Library preparation and sequencing

One ligation library (SQK-LSK114-XL, Oxford Nanopore Technologies) was prepared from the 20-animal extract. It was sequenced on two PromethION R10.4.1 flow cells (FLO-PRO114M; PBE82635 and PBE69031), loaded on 9 September 2025. Each flow cell ran for 72 h with a nuclease flush and library top-off. The core labeled the two runs `Planaria1` (4650P) and `Planaria2` (4651P), but both used the same library.

Data were collected at 5 kHz in MinKNOW 24.11.11 (core 6.2.8). Basecalling was done during the run with the Dorado basecaller (7.6.x) and the super-accuracy model `dna_r10.4.1_e8.2_400bps_sup@v4.3.0`. Reads with mean Q ≥ 10 were kept as "pass" reads.

## Read processing

All computation ran on the UWEC BOSE HPC cluster. Steps were orchestrated with Snakemake 9.14.0 and submitted to Slurm through a cluster-generic executor profile. Software was installed with conda/bioconda (Python 3.12).

Pass reads from both flow cells were concatenated: 31.5 M reads, 192.6 Gb, read N50 13.0 kb, mean read Q 16.7. Reads were filtered with Filtlong 0.3.1 (`--min_length 1000 --keep_percent 80 --min_mean_q 8`). This kept 11.7 M reads and 154.1 Gb (read N50 16.4 kb, mean Q 18.4). Read quality before and after filtering was summarized with NanoPlot 1.46.1.

## Assembly

Filtered reads were assembled with Flye 2.9.6 using these options:
- `--nano-raw`
- `--genome-size 1g`
- `--asm-coverage 80` (the longest reads, up to 80× of the 1 Gb estimate, were used for the initial disjointig assembly)
- `--iterations 2`
- `--scaffold`

Flye chose a minimum overlap of 6 kb automatically. The final run used 48 threads and 500 GB of memory and took about 28 days of wall time.

## Polishing

The Flye assembly was polished with three rounds of Racon 1.5.0 (built from source), using default parameters. In each round, all filtered reads were aligned to the current assembly with minimap2 2.30 (`-ax map-ont`), and the alignments were passed to Racon.

The Racon-polished assembly was then polished once with Medaka 2.1.1 (`medaka_consensus`) using all filtered reads, model `r1041_e82_400bps_sup_g615`, and batch size 100.

## Haplotig and overlap purging *(in progress)*

Duplicate haplotypes and overlapping contigs were removed from the Medaka assembly with purge_dups 1.2.6:
1. Filtered reads were mapped to the unsplit assembly with minimap2 (`-x map-ont`).
2. Per-base read depth and a depth histogram were computed with `pbcstat`.
3. Depth cutoffs were set manually with `calcuts -l 5 -m 86 -u 306`. The values come from the read-depth histogram, which has peaks at about 26× and 101× and a minimum between them at about 86×.
4. The assembly was split at gaps with `split_fa` and aligned to itself with minimap2 (`-x asm5 -DP`).
5. Duplicates were identified with `purge_dups -2` using the self-alignment, depth file and cutoffs.
6. Haplotigs and overlaps at contig ends were removed with `get_seqs -e`, giving the primary (purged) assembly and a separate haplotig set.

## Assembly evaluation

The Flye, Racon, Medaka and purged assemblies were each evaluated with:
- BUSCO 6.0.0 in genome mode against the `metazoa_odb10` lineage (n = 954), with miniprot 0.18 and HMMER 3.4
- QUAST 5.3.0
- SeqKit 2.11.0 (`seqkit stats -a`)

Reports were compiled with MultiQC 1.31.

## Exploratory analysis of allele frequencies

To investigate the two read-depth peaks, we computed allele frequencies with bcftools 1.22 / samtools 1.22.1. Reads came from Medaka's read alignment to the Racon-polished assembly (minimap2, primary alignments only).

Contigs of 100–500 kb were sorted into two classes by approximate depth (mapped reads × mean read length ÷ contig length):
- "collapsed" (85–125×)
- "haplotype-level" (18–35×)

Two hundred contigs were randomly sampled from each class. Pileups were generated with `bcftools mpileup -X ont -I -Q 20 -q 20 -d 1000 -a FORMAT/AD`.

For each SNP site, the alternate-allele frequency was computed from the most common non-reference allele. Sites were kept if their depth was within 0.6–1.4× of the class median and they had ≥ 3 alternate reads. Scripts are in `ploidy_check/`.

## Code availability

The workflow (Snakefile, configuration and executor profile) is at https://github.com/wheelerlab-uwec/Planaria-Genome-Project (branch `nic-pipeline`).

## References

- Grohme MA, Schloissnig S, Rozanski A, et al. (2018). The genome of *Schmidtea mediterranea* and the evolution of core cellular mechanisms. *Nature* 554:56–61.
- Kolmogorov M, Yuan J, Lin Y, Pevzner PA (2019). Assembly of long, error-prone reads using repeat graphs. *Nat Biotechnol* 37:540–546. (Flye)
- Vaser R, Sović I, Nagarajan N, Šikić M (2017). Fast and accurate de novo genome assembly from long uncorrected reads. *Genome Res* 27:737–746. (Racon)
- Li H (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34:3094–3100.
- Guan D, McCarthy SA, Wood J, et al. (2020). Identifying and removing haplotypic duplication in primary genome assemblies. *Bioinformatics* 36:2896–2898. (purge_dups)
- Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM (2021). BUSCO update. *Mol Biol Evol* 38:4647–4654.
- Mölder F, Jablonski KP, Letcher B, et al. (2021). Sustainable data analysis with Snakemake. *F1000Research* 10:33.
- Medaka, Oxford Nanopore Technologies, https://github.com/nanoporetech/medaka

---

## Notes for revision (internal; not for publication)

- **Medaka model mismatch.** Reads were basecalled with `sup@v4.3.0`, but Medaka used `r1041_e82_400bps_sup_g615`. If Medaka is rerun, use the model that matches v4.3.0. Medaka also ran on the CPU despite being on a GPU node.
- **Flye read mode.** `--nano-raw` was used. `--nano-hq` is the recommended mode for Q18 SUP reads. Rerunning Flye is not planned because of the ~28-day runtime.
- **Racon** had no measurable effect on BUSCO and can probably be dropped or described briefly.
- **Ploidy and genotype mixture are unresolved.** The two read-depth peaks are at a ~1:4 ratio rather than the 1:2 expected for a diploid. The allele-frequency analysis could not separate real variants from ONT systematic error. A single-animal Illumina k-mer analysis (GenomeScope/Smudgeplot) would resolve this, if unused aliquots (10, 20, 30, 1 mg) are still at the UWBC.
- **No genome-size estimate exists for *G. dorotocephala*.** Related species range from 0.78 to 1.56 Gb in assemblies, and *G. tigrina* has C-values of 1.40–1.88 pg.
