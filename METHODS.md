# Methods: *Girardia dorotocephala* genome assembly

Draft methods covering work completed through October 2026. Items marked **[TK]** still need a source. The "Notes for revision" section at the end is for internal use, not for publication.

## Animals and tissue preparation

*Girardia dorotocephala* of the asexual clonal line MA-C2 were provided by Labib Rouhana (University of Massachusetts Boston). Almazan et al. (2018) established the line from a single animal obtained from Carolina Biological Supply Company (item #132970) and reported a karyotype of 2n = 24 (twelve small chromosome pairs). The same line was used for the published *G. dorotocephala* transcriptome (BioProject PRJNA317859). Animals arrived on approximately 5 June 2025. They were held in the shipping water for about one week to acclimate, and any unhealthy animals were removed. They were then kept in artificial pond water (APW) at **[TK: temperature]**. Animals were last fed on 20 June 2025 and starved for 5 days before tissue collection.

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

## Haplotig and overlap purging

Duplicate haplotypes and overlapping contigs were removed from the Medaka assembly with purge_dups 1.2.6:
1. Filtered reads were mapped to the unsplit assembly with minimap2 (`-x map-ont`).
2. Per-base read depth and a depth histogram were computed with `pbcstat`.
3. Depth cutoffs were set manually with `calcuts -l 5 -m 86 -u 306`. The values come from the read-depth histogram, which has peaks at about 26× and 101× and a minimum between them at about 86×.
4. The assembly was split at gaps with `split_fa` and aligned to itself with minimap2 (`-x asm5 -DP`).
5. Duplicates were identified with `purge_dups -2` using the self-alignment, depth file and cutoffs.
6. Haplotigs and overlaps at contig ends were removed with `get_seqs -e`, giving the primary (purged) assembly and a separate haplotig set.

purge_dups flagged 1.17 Gb of sequence, as follows:

| Class | Regions | Sequence |
|---|---|---|
| Haplotigs | 26,209 | 1,024 Mb |
| Overlapping contig ends | 3,172 | 122 Mb |
| Repeats | 660 | 23 Mb |
| Low-coverage junk | 417 | 3.4 Mb |
| High-coverage | 54 | 0.7 Mb |

The removed sequence (1.16 Gb) was kept as a separate haplotig set. The high-coverage class included the mitochondrial genome, which Flye had assembled as a 47.5 kb contig at 558× depth; it was reassembled separately (see *Mitochondrial genome*).

## Assembly results

| Stage | Size (Gb) | Contigs | N50 (kb) | BUSCO C | S | D | F | M |
|---|---|---|---|---|---|---|---|---|
| Flye | 2.60 | 52,321 | 91.1 | 68.2% | 27.4% | 40.9% | 3.1% | 28.6% |
| Racon (3 rounds) | 2.62 | 48,551 | 93.4 | 68.2% | 27.4% | 40.9% | 3.1% | 28.6% |
| Medaka | 2.63 | 48,551 | 93.5 | 68.4% | 27.1% | 41.3% | 3.1% | 28.4% |
| Purged (primary) | 1.47 | 21,045 | 131.7 | 66.1% | 54.1% | 12.1% | 3.9% | 30.0% |

BUSCO columns: C = complete, S = complete and single-copy, D = complete and duplicated, F = fragmented, M = missing; metazoa_odb10, n = 954.

The primary assembly is 1,465,713,240 bp in 21,045 contigs:
- N50 131.7 kb (L50 = 3,399); N90 32.8 kb
- largest contig 1.18 Mb
- GC content 30.5%
- 0.01 N per 100 kb

Purging cut duplicated BUSCOs from 41.3% to 12.1%, while complete BUSCOs fell by only 2.3 percentage points. Across all stages, about one-third of BUSCO matches (32–34%) contained internal stop codons. This is consistent with residual indel errors in the consensus.

## Assembly evaluation

The Flye, Racon, Medaka and purged assemblies were each evaluated with:
- BUSCO 6.0.0 in genome mode against the `metazoa_odb10` lineage (n = 954), with miniprot 0.18 and HMMER 3.4
- QUAST 5.3.0
- SeqKit 2.11.0 (`seqkit stats -a`)

Reports were compiled with MultiQC 1.31.

## Contamination screening

The purged primary assembly was screened with NCBI's Foreign Contamination Screen tools (`ncbi_submission/screen/`), as NCBI does for genome submissions:
- **FCS-adaptor** (`--euk`) screened for adaptor and vector sequence.
- **FCS-GX 0.5.5** screened for foreign organisms (database build 2023-01-24, taxid 135777).

The results were applied with `clean_assembly.py`:
- **FCS-GX** flagged 7 contigs (79.3 kb) for removal. One was a 37.8 kb contig with bacteriophage hits, one an 18.3 kb contig with ciliate (*Paramecium*) hits, and five were 1.8–7.5 kb contigs with arthropod hits.
- **FCS-adaptor** found Oxford Nanopore adaptor sequence at 37 contig ends (25–45 bp each), which was trimmed. It also found 8 internal adaptor spans in 6 contigs. Each contig was split at the adaptor, which was removed, because an internal adaptor marks a chimeric join; the resulting pieces were 1.7–292 kb.
- 100 leading Ns on one contig, left by purging at a scaffold gap, were trimmed.
- 11 sequences shorter than 1 kb (517–999 bp; 8.2 kb in total) were removed, following the INSDC minimum sequence length that applies to genome submissions from 2027.

The cleaned assembly has 21,035 sequences totalling 1,465,624,171 bp (N50 131.7 kb; shortest 1,037 bp). Rescreening it with both tools found no adaptor or foreign sequence. (The rescreen was run before the 1 kb filter; removing sequences cannot add contamination.)

## Exploratory analysis of allele frequencies

To investigate the two read-depth peaks, we computed allele frequencies with bcftools 1.22 / samtools 1.22.1. Reads came from Medaka's read alignment to the Racon-polished assembly (minimap2, primary alignments only).

Contigs of 100–500 kb were sorted into two classes by approximate depth (mapped reads × mean read length ÷ contig length):
- "collapsed" (85–125×)
- "haplotype-level" (18–35×)

Two hundred contigs were randomly sampled from each class. Pileups were generated with `bcftools mpileup -X ont -I -Q 20 -q 20 -d 1000 -a FORMAT/AD`.

For each SNP site, the alternate-allele frequency was computed from the most common non-reference allele. Sites were kept if their depth was within 0.6–1.4× of the class median and they had ≥ 3 alternate reads. Scripts are in `ploidy_check/`.

## Mitochondrial genome

The mitochondrial genome was assembled separately from the nuclear assembly (`snakemake/rules/mito.smk`).

**Seed contig.** The two published *Girardia* mitogenomes (MW972220.1, *G. tigrina*; KP090061.1, *Girardia* sp.) were compared with BLAST+ 2.16.0 (`blastn -task dc-megablast`, E ≤ 1e-10) against the 61 nuclear Flye contigs with depth ≥ 300×. The contig with the highest summed bitscore (contig_110383; 47.5 kb, 558×) was used as the seed.

**Read selection and assembly.**
1. Filtered reads were mapped to the seed with minimap2 2.30 (`-x map-ont`).
2. Reads were kept if they were ≥ 3 kb and ≥ 80% of their length aligned to the seed. This excludes reads from nuclear mitochondrial insertions (NUMTs) that carry nuclear flanking sequence. The resulting set had 19,782 reads (175.6 Mb).
3. Filtlong 0.3.1 selected the longest, highest-quality reads up to 4 Mb (152 reads, mean 26.3 kb).
4. These were assembled with Flye 2.9.6 (`--nano-hq --genome-size 16k --iterations 3`).

Flye produced two circular contigs:
- a 21.8 kb contig
- a 44.9 kb contig that is a head-to-tail dimer of the first (~99% identity over both copies)

The monomer was kept.

**Polishing and orientation.**
1. The monomer was polished with Medaka 2.1.1 (model `r1041_e82_400bps_sup_v4.3.0`) using the 4 Mb read subset.
2. It was rotated to start at *cox1* on the forward strand. The *cox1* position was found by BLAST of the reference *cox1* genes against the sequence doubled end-to-end (`scripts/mito_rotate.py`).
3. It was polished again with Medaka so that the original circle junction was polished in the interior of the sequence.

**Annotation.** The genome was annotated with MITOS2 2.1.9 (RefSeq 89 Metazoa reference data, genetic code 9). The annotation was then curated (rules `mito_rrna`, `mito_trna`, `mito_feature_table`; `scripts/mito_feature_table.py`):
- **tRNAs.** tRNAs were also predicted with ARAGORN 1.2.41 (`-mt -gc9 -c`) and tRNAscan-SE 2.0 (`-O`). A tRNA that MITOS2 missed was added from these predictions if it overlapped no protein-coding or rRNA gene, or other tRNA, by more than 10 bp. MITOS2 tRNAs longer than 90 bp were replaced by the overlapping ARAGORN prediction for the same tRNA. A tRNA still missing after this was taken from a relaxed ARAGORN run (`-ps80`) only if it had the expected anticodon, lay on the same strand as the other genes and fell in a gap between genes; the highest-scoring such candidate was used.
- **rRNAs.** The rRNA genes of three annotated Dugesiidae mitogenomes (*S. mediterranea* NC_022448.1, *Dugesia japonica* AB618487.1, *D. ryukyuensis* AB618488.1) were aligned to the sequence (`blastn -task dc-megablast`). Each hit was extended to the full length of the reference gene. Where at least two references hit, the boundary is the median of the projections; otherwise it is the union of the MITOS2 call and the projection. rRNAs were then trimmed so they do not overlap neighbouring genes. The two *Girardia* mitogenomes were not used: their 12S and 16S labels are swapped relative to the other planarian records.
- **Start codons.** Protein-coding genes whose start codon is not ATG or GTG (the start codons of genetic code 9) were annotated as 5′ partial with the note "start codon not determined", as in the *S. mediterranea* and *Dugesia* records.
- The GenBank submission file was built with table2asn 1.28 and checked with its validator and discrepancy report (rule `mito_table2asn`).

**Validation** (rule `mito_qc`):
- All 19,782 mitochondrial reads were mapped to the final sequence with minimap2, and depth was computed with samtools 1.22.1.
- Very long reads (≥ 50 kb) were searched for repeated *cox1* copies; the distance between consecutive copies estimates the genome length directly from single reads.
- A third Medaka round was run to test whether the consensus had converged.
- For species confirmation (rule `mito_barcode`), the *cox1* gene was compared by BLAST with all *Girardia* COI records in GenBank (318 records; retrieved 2026-10-16).

### Mitochondrial genome results

The mitochondrial genome is a circular molecule of 22,455 bp with 31.5% GC.

**Read support for the length.**
- Mean read depth was 7,312× (median of 500 bp windows; range 5,937–8,345×), with no gaps or drops in depth.
- 174 reads aligned contiguously across ≥ 20 kb of the sequence.
- In reads ≥ 50 kb, the most common spacing between consecutive *cox1* copies was 22.0–23.0 kb (6 of 13 spacings).

The genome is therefore substantially longer than the two published *Girardia* mitogenomes (15.9 kb). The extra length is in non-coding regions.

**Annotation.** The curated annotation contains:
- all 12 protein-coding genes expected for flatworms (*atp8* is absent in Platyhelminthes). All 12 are open reading frames with no internal stop codons. *cox1* (ATT) and *nad5* (TTG) have no ATG or GTG start nearby and are annotated as 5′ partial.
- *rrnS* (14,257–14,976; 720 bp) and *rrnL* (15,643–16,601; 959 bp). MITOS2 had split *rrnS* into two fragments and placed the start of *rrnL* about 280 bp too far downstream. All three Dugesiidae references agree on the *rrnL* ends; *rrnS* is supported by only one reference and its ends are less certain.
- all 22 tRNAs:
  - 19 from MITOS2
  - *trnD* (4,948–5,003) and *trnS1* (15,206–15,279), added from ARAGORN; both lie in gaps between other genes
  - *trnY* (5,150–5,224; anticodon GTA), from the relaxed ARAGORN run. It lies in the gap between *trnD* and *trnL2* and folds into a complete cloverleaf, but was missed at the default threshold, by MITOS2 and by tRNAscan-SE, so it is the least certain tRNA.
  - *trnR*, shortened from MITOS2's 103 bp to ARAGORN's 65 bp

All genes are on one strand. table2asn reported no errors. Its only warnings concern the two 5′-partial genes.

**Consensus stability.** A third Medaka round changed 12 substitutions and 12 indels (+14 bp). These changes, and most positions where < 80% of reads agree with the consensus, lie in two non-coding regions (~2.0–3.5 kb and ~12.0–12.5 kb). These regions share a 570 bp repeat (85% identity). Their sequence is less certain than the rest of the genome and may include length heteroplasmy.

**Species identity.** The *cox1* gene is 100% identical to 29 *G. dorotocephala* COI records (Benítez-Álvarez et al. 2023), and ≤ 92% identical to all other *Girardia* species (*G. sinensis*, *G. tigrina* and *Girardia* spp.: 91.4–91.6%). This confirms the species identity of the sequenced line. Identity to the two published *Girardia* mitogenomes is low (66–72% over local alignments of 0.6–1.4 kb). Both belong to the *G. tigrina* lineage, so this is, to our knowledge, the first mitochondrial genome for *G. dorotocephala*.

## Data availability

The project is registered at NCBI as BioProject PRJNA1545025, and the sequenced sample as BioSample SAMN64250158. The locus_tag prefix is ADAYWY. Accessions for the raw reads (SRA), the genome assembly and the mitochondrial genome are pending; all records stay private until publication.

## Code availability

The workflow (Snakefile, configuration and executor profile) is at https://github.com/wheelerlab-uwec/Planaria-Genome-Project (branch `nic-pipeline`).

## References

- Almazan EMH, Lesko SL, Markey MP, Rouhana L (2018). *Girardia dorotocephala* transcriptome sequence, assembly, and validation through characterization of piwi homologs and stem cell progeny markers. *Dev Biol* 433:433–447.
- Grohme MA, Schloissnig S, Rozanski A, et al. (2018). The genome of *Schmidtea mediterranea* and the evolution of core cellular mechanisms. *Nature* 554:56–61.
- Ivanković M, Brand JN, Pandolfini L, et al. (2024). A comparative analysis of planarian genomes reveals regulatory conservation in the face of rapid structural divergence. *Nat Commun* 15:8215.
- Kolmogorov M, Yuan J, Lin Y, Pevzner PA (2019). Assembly of long, error-prone reads using repeat graphs. *Nat Biotechnol* 37:540–546. (Flye)
- Vaser R, Sović I, Nagarajan N, Šikić M (2017). Fast and accurate de novo genome assembly from long uncorrected reads. *Genome Res* 27:737–746. (Racon)
- Li H (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34:3094–3100.
- Guan D, McCarthy SA, Wood J, et al. (2020). Identifying and removing haplotypic duplication in primary genome assemblies. *Bioinformatics* 36:2896–2898. (purge_dups)
- Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM (2021). BUSCO update. *Mol Biol Evol* 38:4647–4654.
- Mölder F, Jablonski KP, Letcher B, et al. (2021). Sustainable data analysis with Snakemake. *F1000Research* 10:33.
- Medaka, Oxford Nanopore Technologies, https://github.com/nanoporetech/medaka
- Benítez-Álvarez L, Sluys R, Leal-Zanchet AM, Leria L, Riutort M (2023). First molecular phylogeny of the freshwater planarian genus *Girardia* (Platyhelminthes: Tricladida) unveils hidden taxonomic diversity and initiates resolution of its historical biogeography. *Zool J Linn Soc*. doi:10.1093/zoolinnean/zlac065
- Donath A, Jühling F, Al-Arab M, et al. (2019). Improved annotation of protein-coding genes boundaries in metazoan mitochondrial genomes. *Nucleic Acids Res* 47:10543–10552. (MITOS2)
- Camacho C, Coulouris G, Avagyan V, et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics* 10:421.
- Filtlong, R. Wick, https://github.com/rrwick/Filtlong
- Laslett D, Canbäck B (2004). ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. *Nucleic Acids Res* 32:11–16.
- Chan PP, Lin BY, Mak AJ, Lowe TM (2021). tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. *Nucleic Acids Res* 49:9077–9096.

---

## Notes for revision (internal; not for publication)

- **Medaka model mismatch.** Reads were basecalled with `sup@v4.3.0`, but Medaka used `r1041_e82_400bps_sup_g615`. If Medaka is rerun, use the model that matches v4.3.0. Medaka also ran on the CPU: Slurm job 349398 held a GPU for 2 d 13 h but did not use it.
- **BUSCO internal stop codons are a method artifact, not assembly error.** About 33% of complete BUSCOs contain internal stop codons at every stage. Re-polishing a 164 Mb BUSCO-bearing subset with Dorado polish (sup@v4.3.0 model; `polish_test/`) gave the same BUSCO result as Medaka. The *S. mediterranea* references run through the same BUSCO setup (`smed_baseline/`) show an even higher rate:

  | Reference | BUSCO C | S | D | M | Complete with internal stops |
  |---|---|---|---|---|---|
  | schMedS3h1 (HiFi + Illumina, chromosome-level) | 70.4% | 65.2% | 5.2% | 27.6% | 41.8% |
  | ASM260089v1 (PacBio CLR, 2018) | 69.7% | 51.4% | 18.3% | 28.2% | 40.8% |

  This is consistent with the stop-codon flag reflecting miniprot alignment on planarian genes, so it should not be read as a measure of base accuracy, and planarian BUSCO completeness tops out near 70%. A true base-level QV needs short reads (Merqury).

  These reference numbers are our own runs; neither *S. mediterranea* genome paper reports a stop-codon rate (BUSCO v5.0, which they used, predates the miniprot-based count). Ivanković et al. (2024) do report that "transcriptome or gene model-based BUSCO scores were consistently better than genome-mode BUSCO assessments … indicating that the detection method used by BUSCO is sub-optimal for planarian genomes," and that the "comparatively high number of 'missing' BUSCOs reflects a high proportion of genuine gene losses in planarians." Their genome-mode scores for *S. polychroa*, *S. nova* and *S. lugubris* were 602, 592 and 597 complete of 954 (~62–63%; BUSCO v5.0.0, metazoa_odb10; their Fig. 3c). Grohme et al. (2018) did not use BUSCO.
- **Residual duplication.** 12.1% of BUSCOs are still duplicated after purging. These may be divergent haplotigs or real gene copies; this has not been checked.
- **Flye read mode.** `--nano-raw` was used. `--nano-hq` is the recommended mode for Q18 SUP reads. Rerunning Flye is not planned because of the ~28-day runtime.
- **Racon** had no measurable effect on BUSCO and can probably be dropped or described briefly.
- **Ploidy and genotype mixture are unresolved.** The two read-depth peaks are at a ~1:4 ratio rather than the 1:2 expected for a diploid. The allele-frequency analysis could not separate real variants from ONT systematic error. A single-animal Illumina k-mer analysis (GenomeScope/Smudgeplot) would resolve this, if unused aliquots (10, 20, 30, 1 mg) are still at the UWBC.
- **No genome-size estimate exists for *G. dorotocephala*.** Related species range from 0.78 to 1.56 Gb in assemblies, and *G. tigrina* has C-values of 1.40–1.88 pg.
- **Mitochondrial annotation.** *trnY* comes from a relaxed ARAGORN search (added after GenBank flagged it as missing), and the *rrnS* ends rest on one reference. The submission template (`ncbi_submission/template.sbt`) lists N.J. Wheeler as the only author; update it before submission. Check volume and pages for Benítez-Álvarez et al. (published online 2022).
