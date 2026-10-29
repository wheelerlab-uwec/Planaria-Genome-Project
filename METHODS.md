# Methods: *Girardia dorotocephala* genome assembly

Draft methods, organised as the Methods and Data Validation sections of a GigaByte Data Release. The "Notes for revision" section at the end is internal.

## Methods

### Sampling strategy

*Girardia dorotocephala* of the asexual clonal line MA-C2 were provided by Labib Rouhana (University of Massachusetts Boston). Almazan et al. (2018) established the line from a single animal obtained from Carolina Biological Supply Company (item #132970), reported a karyotype of 2n = 24, and used it for the published *G. dorotocephala* transcriptome (BioProject PRJNA317859).

Animals were kept in artificial pond water at room temperature and starved for 5 days. Mucus was then removed with 0.5% (w/v) N-acetyl-L-cysteine in 20 mM HEPES-NaOH (pH > 7) for 10–15 min, following Grohme et al. (2018). Animals were rinsed, flash-frozen in liquid nitrogen in aliquots of about 5 mg per animal, and stored at −80 °C. All sequence data come from a single DNA extract of 20 pooled animals.

### DNA extraction and sequencing

High-molecular-weight DNA was extracted by the University of Wisconsin Biotechnology Center (UWBC; RRID:SCR_017759) using the planarian protocol of Grohme et al. (2018; Supplementary Methods S1.3–S1.4). That protocol uses guanidinium thiocyanate/N-lauroylsarcosine lysis, followed by CTAB polysaccharide cleanup with RNase treatment. The extract measured 140 ng/µl by Qubit. On Femto Pulse it had a modal fragment length of about 20 kb, with fragments over 165 kb.

One ligation library (SQK-LSK114-XL, Oxford Nanopore Technologies) was sequenced on two PromethION R10.4.1 flow cells (FLO-PRO114M) for 72 h each, with a nuclease flush and library top-off. Data were collected at 5 kHz in MinKNOW 24.11.11. Reads were basecalled during the run with Dorado 7.6 and the super-accuracy model `dna_r10.4.1_e8.2_400bps_sup@v4.3.0`, keeping reads with mean Q ≥ 10.

### Read processing and assembly

All analyses were run with Snakemake 9.14.0 on a Slurm cluster, with software from bioconda.

1. **Reads.** Reads from both flow cells (31.5 M reads, 192.6 Gb, read N50 13.0 kb, mean Q 16.7) were filtered with Filtlong 0.3.1 (`--min_length 1000 --keep_percent 80 --min_mean_q 8`). This left 11.7 M reads and 154.1 Gb (read N50 16.4 kb, mean Q 18.4). Read quality was summarised with NanoPlot 1.46.1.
2. **Assembly.** Filtered reads were assembled with Flye 2.9.6 (`--nano-raw --genome-size 1g --asm-coverage 80 --iterations 2 --scaffold`).
3. **Polishing.** The assembly was polished with three rounds of Racon 1.5.0, each using minimap2 2.30 alignments of all filtered reads (`-ax map-ont`). It was then polished once with Medaka 2.1.1 (model `r1041_e82_400bps_sup_g615`).
4. **Purging.** Haplotigs and overlapping contig ends were removed with purge_dups 1.2.6. The read-depth histogram had peaks at about 26× and 101× and a minimum at about 86×, so depth cutoffs were set manually (`calcuts -l 5 -m 86 -u 306`). Duplicates were identified with `purge_dups -2` from a self-alignment of the split assembly (minimap2 `-x asm5 -DP`) and removed with `get_seqs -e`. purge_dups removed 1.16 Gb, mostly haplotigs (1.02 Gb) and overlapping ends (0.12 Gb). The removed sequence was kept as a separate haplotig set.
5. **Contamination screening.** The purged assembly was screened with NCBI's FCS-adaptor (`--euk`) and FCS-GX 0.5.5 (database 2023-01-24, taxid 135777). Then:
   - the 7 contigs flagged by FCS-GX were removed
   - terminal Oxford Nanopore adaptors were trimmed
   - contigs with internal adaptors were split at the adaptor
   - terminal Ns were trimmed
   - sequences shorter than 1 kb were removed (the INSDC minimum for genome submissions from 2027)

### Mitochondrial genome

The mitochondrial genome is absent from the primary assembly: purge_dups removed it as high-coverage sequence, a 47.5 kb contig at 558× depth. It was reassembled from reads in the following steps.

1. **Seed.** That contig was identified by BLAST+ 2.16.0 (`dc-megablast`) of two published *Girardia* mitogenomes (MW972220.1, KP090061.1) against high-depth Flye contigs.
2. **Read selection.** Filtered reads ≥ 3 kb with ≥ 80% of their length aligned to the seed (minimap2) were kept, which excludes nuclear mitochondrial insertions with nuclear flanks. This gave 19,782 reads. Filtlong selected the longest of them, up to 4 Mb.
3. **Assembly.** The selected reads were assembled with Flye 2.9.6 (`--nano-hq`). Flye returned a 21.8 kb circular contig, and a 44.9 kb head-to-tail dimer of it, which was discarded.
4. **Polishing and rotation.** The contig was polished with Medaka 2.1.1 (`r1041_e82_400bps_sup_v4.3.0`), rotated to start at *cox1*, and polished again.
5. **Annotation.** The genome was annotated with MITOS2 2.1.9 (RefSeq 89 Metazoa, genetic code 9) and then curated:
   - tRNAs missed by MITOS2 were added from ARAGORN 1.2.41 (`-mt -gc9`) or, as a last resort, a relaxed ARAGORN search (`-ps80`). Only candidates with the expected anticodon, lying in gaps between genes, were accepted.
   - rRNA ends were projected from the rRNA genes of three Dugesiidae mitogenomes (*S. mediterranea* NC_022448.1, *Dugesia japonica* AB618487.1, *D. ryukyuensis* AB618488.1).
   - Protein-coding genes without an ATG or GTG start were marked 5′ partial.
   - The GenBank file was built and validated with table2asn 1.28.

### Nuclear genome annotation

The submitted primary assembly was annotated with the workflow in `snakemake/rules/annotation*.smk`.

1. **Repeats.** A species-specific repeat library was built with RepeatModeler 2.0.8 (`-LTRStruct`). The genome was soft-masked with RepeatMasker 4.2.3 (`-xsmall`) using that library.
2. **RNA-seq.** Four paired-end Illumina RNA-seq runs from the same clonal line (MA-C2; BioProject PRJNA317859, Almazan et al. 2018) were used: intact animals (SRR3479052), regenerates at 1 and 4 days post-amputation (SRR3479045, SRR3479046) and auricles (SRR3479048). Reads were downloaded with SRA Toolkit 3.4.1 and trimmed with fastp 1.4.0 (`--detect_adapter_for_pe`). They were aligned with STAR 2.7.11b in two-pass mode (`--twopassMode Basic`, introns and mate gaps ≤ 100 kb, non-canonical junctions removed), and the four alignments were merged.
3. **Gene prediction.** Protein-coding genes were predicted with BRAKER 3.0.8 (GeneMark-ETP, AUGUSTUS 3.5.0) on the soft-masked genome, using the merged RNA-seq alignments and the Metazoa partition of OrthoDB 12 as protein evidence.
4. **Functional annotation.** Proteins were annotated with eggNOG-mapper 2.1.15 (eggNOG 5.0.2, DIAMOND 2.1.13) and InterProScan 5.78-109.0 (InterPro 109.0, with GO terms). Product names follow the International Protein Nomenclature Guidelines used by NCBI, assigned in this order:
   1. the Swiss-Prot recommended name (UniProtKB 2026_03, reviewed Metazoa entries) for the eggNOG gene symbol, when the eggNOG ortholog group was at the Metazoa level or below and no more than three genes shared the symbol;
   2. the eggNOG description, if it was a protein name;
   3. the name of the InterPro family hit;
   4. the name of the InterPro domain or repeat hit with the greatest coverage (as "… domain-containing protein");
   5. otherwise, "hypothetical protein".

   Names were then cleaned to NCBI style (EC numbers, species, "homolog" and similar words removed; American spelling; lower case except for acronyms).
5. **tRNAs.** tRNA genes were predicted with tRNAscan-SE 2.0.13 in eukaryotic mode (`-E -I -H`, Infernal 1.1.5) on the unmasked genome. Pseudogenes and tRNA-derived repeats were removed with EukHighConfidenceFilter, keeping only the high-confidence set. Undetermined and suppressor tRNAs, and any tRNA overlapping a coding sequence, were excluded.
6. **GenBank submission.** Models were given ADAYWY locus tags. Transcripts were removed when the protein was shorter than 30 aa, when a missing start or stop codon was not at a contig end, or when they repeated another gene's coding structure. Genes running off a contig end were marked partial. The GenBank file was built and validated with NCBI table2asn 1.29.324 (release 2025-09-16).
7. **Completeness.** The longest protein per gene was assessed with BUSCO 6.0.0 (protein mode, metazoa_odb10) and OMArk 0.5.0 (OMAmer 2.1.2, OMA All.May2026 database).

## Data validation and quality control

### Nuclear assembly

| Stage | Size (Gb) | Contigs | N50 (kb) | BUSCO C | S | D | F | M |
|---|---|---|---|---|---|---|---|---|
| Flye | 2.60 | 52,321 | 91.1 | 68.2% | 27.4% | 40.9% | 3.1% | 28.6% |
| Racon (3 rounds) | 2.62 | 48,551 | 93.4 | 68.2% | 27.4% | 40.9% | 3.1% | 28.6% |
| Medaka | 2.63 | 48,551 | 93.5 | 68.4% | 27.1% | 41.3% | 3.1% | 28.4% |
| Purged (primary) | 1.47 | 21,045 | 131.7 | 66.1% | 54.1% | 12.1% | 3.9% | 30.0% |

BUSCO 6.0.0, genome mode, metazoa_odb10 (n = 954; miniprot 0.18, HMMER 3.4). C = complete, S = single-copy, D = duplicated, F = fragmented, M = missing. Assembly statistics are from QUAST 5.3.0 and SeqKit 2.11.0.

**Primary assembly as submitted** (after screening): 1,465,624,171 bp in 21,035 contigs; N50 131.7 kb (L50 3,399); N90 32.8 kb; longest contig 1.18 Mb; GC 30.5%; no gaps.

- **Purging** cut duplicated BUSCOs from 41.3% to 12.1%, while complete BUSCOs fell by only 2.3 percentage points.
- **BUSCO completeness and stop codons.** About one-third of complete BUSCOs contained internal stop codons at every stage. This does not indicate consensus errors:
  - Re-polishing a 164 Mb BUSCO-bearing subset with Dorado polish (sup@v4.3.0) gave the same result as Medaka.
  - Two *S. mediterranea* assemblies run through the same BUSCO setup scored 70.4% and 69.7% complete, with 41.8% and 40.8% of complete BUSCOs containing internal stops (schMedS3h1 and ASM260089v1).

  Genome-mode BUSCO is known to underperform on planarian genomes (Ivanković et al. 2024), so gene-level completeness will be assessed on the annotation.
- **Contamination screening** removed 79.3 kb in 7 contigs (bacteriophage, ciliate and arthropod hits). It also trimmed adaptors at 37 contig ends and split 6 contigs at 8 internal adaptors. Rescreening the cleaned assembly with both tools found no adaptor or foreign sequence.

### Mitochondrial genome

The mitochondrial genome is a circular molecule of 22,455 bp (31.5% GC). It is longer than the published *Girardia* mitogenomes (15.9 kb); the difference is in non-coding sequence.

- **Read support for the length.**
  - Mean read depth was 7,312× (range of 500 bp windows 5,937–8,345×), with no gaps.
  - 174 reads aligned contiguously across ≥ 20 kb of the sequence.
  - In reads ≥ 50 kb, consecutive *cox1* copies were most often 22–23 kb apart.
- **Annotation.**
  - All 12 protein-coding genes expected for flatworms are present as open reading frames; *cox1* and *nad5* are 5′ partial.
  - Both rRNAs and all 22 tRNAs are annotated, all on one strand. *trnY* comes from the relaxed search and is the least certain tRNA, and the *rrnS* ends rest on a single reference.
- **Consensus stability.** A third Medaka round changed 24 positions. Most fall in two non-coding regions that share a 570 bp repeat, so the sequence there is less certain.
- **Species identity.** *cox1* is 100% identical to 29 *G. dorotocephala* COI records (Benítez-Álvarez et al. 2023) and ≤ 92% identical to other *Girardia* species. To our knowledge this is the first *G. dorotocephala* mitogenome.

### Nuclear annotation

**Repeats.** 69.8% of the assembly (1.02 Gb) was masked.

| Class | % of genome |
|---|---|
| LTR retrotransposons | 24.5 (Gypsy 16.6, BEL/Pao 1.2) |
| DNA transposons | 21.9 (Tc1/mariner 10.5, hAT 3.8, MULE 1.3) |
| LINEs | 5.9 |
| Unclassified | 15.2 |

**RNA-seq alignment.** 69.2–75.8% of read pairs from each run mapped uniquely and 8.7–9.0% mapped to multiple loci.

**Gene models.**
- BRAKER3 predicted 21,181 genes (25,016 transcripts; 1.18 per gene). Transcripts average 5.55 exons; 5,470 are single-exon. The median coding sequence is 909 bp and the median intron 123 bp.
- After filtering, the GenBank annotation has 21,116 protein-coding genes (24,939 mRNAs, 95 partial at contig ends). 62 transcripts were removed as < 30 aa and 15 as incomplete away from a contig end.
- 15,656 genes (74.1%) received a functional name; the rest are "hypothetical protein".
- table2asn reported no errors. Its 211 warnings are 207 short internal coding exons and 4 pairs of genes with identical spans.

**Completeness of the gene set.**

| Assessment | Result |
|---|---|
| BUSCO, metazoa_odb10 (n = 954), proteins | C 75.5% [S 63.7%, D 11.7%], F 4.3%, M 20.2% |
| BUSCO, same lineage, genome mode (above) | C 66.1% |
| OMArk, Platyhelminthes conserved HOGs (n = 2,828) | S 60.5%, D 14.8% [unexpected 12.8%, expected 2.0%], M 24.7% |
| OMArk proteome consistency (21,174 proteins) | consistent 46.4%, inconsistent 11.7%, contaminant 1.5%, unknown 40.4% |

OMArk placed the proteome closest to *Schmidtea mediterranea*. The 40.4% "unknown" proteins have no OMA family hit, as is common for planarian-specific genes.

**tRNAs.** tRNAscan-SE predicted 2,584 tRNAs, including 68 pseudogenes. The high-confidence set has 2,245, all added to the GenBank annotation (none overlapped a coding sequence). Isotype counts are very uneven (Arg 376, Ser 341, Leu 307, Asn 268 [all GTT], Pro 255, Trp 153, but Ile 9, Phe 5, Met/iMet 4, Thr 2, Cys 1), and the raw predictions show the same pattern, so the large families probably include tRNA-derived repeat copies that pass the filter.

## Data availability

The project is registered at NCBI as BioProject PRJNA1545025, and the sequenced sample as BioSample SAMN64250158 (locus_tag prefix ADAYWY). This Whole Genome Shotgun project has been deposited at DDBJ/ENA/GenBank under the accession JCFHDI000000000; the version described here is JCFHDI010000000. The gene and tRNA annotation was submitted with the WGS record. The mitochondrial genome is GenBank QB207629. The raw reads are in the Sequence Read Archive as run SRR41252593. All records stay private until publication.

## Code availability

The workflow (Snakemake rules, scripts, configuration and environment files) is at https://github.com/wheelerlab-uwec/Planaria-Genome-Project.

## References

- Almazan EMH, Lesko SL, Markey MP, Rouhana L (2018). *Girardia dorotocephala* transcriptome sequence, assembly, and validation through characterization of piwi homologs and stem cell progeny markers. *Dev Biol* 433:433–447.
- Benítez-Álvarez L, Sluys R, Leal-Zanchet AM, Leria L, Riutort M (2023). First molecular phylogeny of the freshwater planarian genus *Girardia* (Platyhelminthes: Tricladida) unveils hidden taxonomic diversity and initiates resolution of its historical biogeography. *Zool J Linn Soc*. doi:10.1093/zoolinnean/zlac065
- Camacho C, Coulouris G, Avagyan V, et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics* 10:421.
- Donath A, Jühling F, Al-Arab M, et al. (2019). Improved annotation of protein-coding genes boundaries in metazoan mitochondrial genomes. *Nucleic Acids Res* 47:10543–10552. (MITOS2)
- Grohme MA, Schloissnig S, Rozanski A, et al. (2018). The genome of *Schmidtea mediterranea* and the evolution of core cellular mechanisms. *Nature* 554:56–61.
- Guan D, McCarthy SA, Wood J, et al. (2020). Identifying and removing haplotypic duplication in primary genome assemblies. *Bioinformatics* 36:2896–2898. (purge_dups)
- Ivanković M, Brand JN, Pandolfini L, et al. (2024). A comparative analysis of planarian genomes reveals regulatory conservation in the face of rapid structural divergence. *Nat Commun* 15:8215.
- Kolmogorov M, Yuan J, Lin Y, Pevzner PA (2019). Assembly of long, error-prone reads using repeat graphs. *Nat Biotechnol* 37:540–546. (Flye)
- Laslett D, Canbäck B (2004). ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. *Nucleic Acids Res* 32:11–16.
- Li H (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34:3094–3100.
- Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM (2021). BUSCO update. *Mol Biol Evol* 38:4647–4654.
- Mölder F, Jablonski KP, Letcher B, et al. (2021). Sustainable data analysis with Snakemake. *F1000Research* 10:33.
- Vaser R, Sović I, Nagarajan N, Šikić M (2017). Fast and accurate de novo genome assembly from long uncorrected reads. *Genome Res* 27:737–746. (Racon)
- Filtlong, R. Wick, https://github.com/rrwick/Filtlong
- Medaka, Oxford Nanopore Technologies, https://github.com/nanoporetech/medaka
- Still to cite: FCS (Astashyn et al. 2024), NanoPlot, QUAST, SeqKit, Dorado; for the annotation, RepeatModeler2 (Flynn et al. 2020), RepeatMasker, SRA Toolkit, fastp (Chen et al. 2018), STAR (Dobin et al. 2013), BRAKER3 (Gabriel et al. 2024), GeneMark-ETP (Brůna et al. 2024), AUGUSTUS (Stanke et al. 2008), OrthoDB (Tegenfeldt et al. 2025), eggNOG-mapper 2 (Cantalapiedra et al. 2021), eggNOG 5 (Huerta-Cepas et al. 2019), DIAMOND (Buchfink et al. 2021), InterProScan (Jones et al. 2014), InterPro (Blum et al. 2025), UniProt (UniProt Consortium 2025), tRNAscan-SE 2 (Chan et al. 2021), Infernal (Nawrocki & Eddy 2013), OMArk (Nevers et al. 2025), OMAmer (Rossier et al. 2021). Check years and add full entries.

---

## Notes for revision (internal; not for publication)

- **Medaka model mismatch.** Reads were basecalled with `sup@v4.3.0`, but Medaka used `r1041_e82_400bps_sup_g615`. Rerun with the matching model, or justify the choice with the Dorado polish comparison (identical BUSCO).
- **Flye read mode.** `--nano-raw` was used; `--nano-hq` is recommended for Q18 SUP reads. Rerunning Flye (~28 days) is not planned; be ready to justify this.
- **Racon** had no measurable effect on BUSCO; consider describing it in one line or dropping it.
- **Residual duplication.** 12.1% of BUSCOs are still duplicated after purging, either divergent haplotigs or real gene copies. Not yet checked.
- **Ploidy and genotype mixture are unresolved.** The two read-depth peaks are at ~1:4 rather than the 1:2 expected for a diploid. An exploratory allele-frequency analysis (bcftools on contigs from both depth classes; `ploidy_check/`) could not separate real variants from ONT systematic error. A single-animal Illumina k-mer analysis (GenomeScope/Smudgeplot) would resolve this, and would also give a QV (Merqury) and a genome-size estimate, if unused aliquots (10, 20, 30, 1 mg) are still at the UWBC.
- **No genome-size estimate exists for *G. dorotocephala*.** Related species range from 0.78 to 1.56 Gb in assemblies; *G. tigrina* C-values are 1.40–1.88 pg.
- **BUSCO baseline detail** (for a supplementary table): schMedS3h1 C 70.4 / S 65.2 / D 5.2 / M 27.6%; ASM260089v1 C 69.7 / S 51.4 / D 18.3 / M 28.2%. These are our own runs; neither *S. mediterranea* paper reports a stop-codon rate. Ivanković et al. (2024) report genome-mode scores of ~62–63% for *S. polychroa*, *S. nova* and *S. lugubris* (BUSCO 5.0.0, their Fig. 3c).
- **Mitochondrial annotation.** *trnY* was added after GenBank flagged it as missing.
- **Removed from the methods text** (kept in the lab records and repo): sample-collection dates, aliquot table, the failed Nanobind extraction, flow-cell/run labels, the Flye runtime, the purge_dups class table, the allele-frequency method, and the step-by-step mito curation rules.
- Check volume and pages for Benítez-Álvarez et al. (published online 2022).
- The submission template (`ncbi_submission/template.sbt`) lists N.J. Wheeler as the only author.
- **Annotation naming choices.** A Swiss-Prot name is used only when ≤ 3 genes share the eggNOG symbol, so large paralog families get family-level names instead of one gene's name (e.g. 58 genes were all "TRAF6"). Be ready to explain this.
- **OMArk contamination line.** OMArk flags 342 proteins (1.6%) as Formicidae (ants); the FCS-GX rescreen of the cleaned assembly found no foreign sequence. Not investigated.
- **tRNA isotype skew** (see validation): mention if a tRNA count is quoted.
- **GigaByte checklist (2026-10-28).** Still needed: GigaDB dataset (readme.txt + md5) for reviewer access; LICENSE, merge `nic-pipeline` to `main`, tagged release with a DOI; abstract, introduction, Context, Re-use potential, abbreviations, declarations (funding, CRediT, acknowledgements), author list, cover letter; numbered references; specimen voucher. APC $535.
