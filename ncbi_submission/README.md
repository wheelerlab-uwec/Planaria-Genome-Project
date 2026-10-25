# NCBI submissions

Log of the *Girardia dorotocephala* MA-C2 submissions to NCBI: what was sent, how
it was made, and the steps done by hand in the NCBI Submission Portal. All
records are held until 2027-10-16 or publication, whichever is first.

## Accessions

| Record | Accession | Submission | Date |
|---|---|---|---|
| BioProject | PRJNA1545025 | | 2026-10-16 |
| BioSample | SAMN64250158 (Invertebrate.1.0, isolate MA-C2) | | 2026-10-16 |
| Locus tag prefix | ADAYWY | | 2026-10-16 |
| SRA (ONT reads) | pending | SUB16565902 | 2026-10-19 |
| Nuclear genome (WGS) | JCFHDI000000000 (cite JCFHDI010000000) | SUB16566354 | 2026-10-19 |
| Mitochondrial genome | QB207629 | SUB16566882 | 2026-10-16, revised 2026-10-19 |

## SRA (`sra/`)

- `prepare_reads.sbatch`: one FASTQ per flow cell from the MinKNOW fastq_pass
  chunks, checked against the sequencing core's md5 sums.
- `split_upload.sbatch`: splits PBE69031 into `_part1`/`_part2`, because SRA
  caps files at 100 GiB (it was 100.3 GiB).
- `upload_reads.sbatch`: FTP upload to the account's preload folder.
- `sra_metadata.tsv`: the metadata table uploaded to the portal (one run,
  three files).
- Portal: SRA wizard, files picked from the preload folder, metadata TSV
  uploaded, submitted 2026-10-19.

## Nuclear genome (`screen/`, `genome/`)

- `screen/run_fcs.sbatch`: NCBI FCS-GX and FCS-adaptor on the purged primary
  assembly (`setup_fcs.sbatch` fetches the databases).
- `screen/clean_assembly.py` via `clean_fcs.sbatch` (`--min-len 1000`): drops
  the 7 FCS-GX contaminant contigs, trims 37 terminal adaptors, splits 6 contigs
  at 8 internal adaptors (`_s1`/`_s2` names), trims terminal Ns and drops
  sequences under 1 kb. Output `genome/Gdor_MA-C2_primary.fsa`: 21,035
  sequences, 1,465,624,171 bp (`cleaning_log.tsv`, `assembly_stats.tsv`).
- `genome/upload_genome.sbatch`: FTP upload to preload folder
  `gdor_ma-c2_genome_v2`.
- Portal: Genome wizard, WGS, primary assembly only, file from the preload
  folder, submitted 2026-10-19. The portal looped without an error message
  while the file still had contigs under 1 kb; the 1 kb cutoff fixed it.

## Mitochondrial genome

- Built by `snakemake/rules/mito.smk`; `mito_table2asn` writes the `.sqn`
  under `snakemake/results/mito/submission/`.
- Portal: GenBank (SP-GenBank) with FASTA plus feature table, submitted
  2026-10-16. GenBank flagged a missing *trnY*. It was found with a relaxed
  ARAGORN search (`-ps80`, anticodon GTA, positions 5150–5224; commit d8d3713)
  and the corrected files were sent through the portal's Fix option on
  2026-10-19.

## Nuclear annotation (not yet submitted)

- Built by `snakemake/rules/annotation_ncbi.smk`: eggNOG-derived product
  names, ADAYWY locus tags and `table2asn` (`.sqn`, `.val`, `.dr`) under
  `snakemake/results/annotation/ncbi/`.
- To be sent as annotation for JCFHDI after checking the submission status with
  NCBI (genomes@ncbi.nlm.nih.gov).
