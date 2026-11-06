#!/usr/bin/env python3
"""Tabulate our assembly against published triclad genomes.

Size and contiguity of the published assemblies come from their NCBI dataset
reports; ours are computed from the .fai (no gaps, so contigs = sequences).
BUSCO columns come from the short summaries of the same BUSCO run on every
assembly.
Repeat % for published assemblies is the figure from their papers (methods
differ between studies); ours is the masked share from RepeatMasker.

Usage: cmp_table.py META_JSON genome.fa.fai genome.fa.tbl REPORT_DIR BUSCO_DIR
  META_JSON: {"GCA_...": ["species", "family", repeat % or null, "source"], ...}
"""
import json
import re
import sys

BUSCO_RE = re.compile(
    r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%"
)
COLS = ["species", "family", "assembly", "accession", "released", "technology",
        "level", "size_Mb", "contigs", "contig_N50_kb", "scaffolds",
        "scaffold_N50_kb", "BUSCO_C", "BUSCO_S", "BUSCO_D", "BUSCO_F",
        "BUSCO_M", "repeat_pct", "repeat_source"]


def busco(path):
    with open(path) as fh:
        m = BUSCO_RE.search(fh.read())
    return list(m.groups())


def n50(lengths):
    half, cum = sum(lengths) / 2, 0
    for n in sorted(lengths, reverse=True):
        cum += n
        if cum >= half:
            return n


def masked_pct(tbl):
    with open(tbl) as fh:
        return re.search(r"bases masked:.*\(\s*([\d.]+) %\)", fh.read()).group(1)


def ours(fai):
    with open(fai) as fh:
        lengths = [int(line.split("\t")[1]) for line in fh]
    k = n50(lengths)
    return ["Girardia dorotocephala", "Dugesiidae", "Gdor MA-C2 (this study)",
            "JCFHDI000000000", "2026", "Oxford Nanopore PromethION", "Contig",
            f"{sum(lengths) / 1e6:.1f}", len(lengths), f"{k / 1e3:.1f}",
            len(lengths), f"{k / 1e3:.1f}"]


def published(acc, species, family, report_dir):
    with open(f"{report_dir}/{acc}.json") as fh:
        r = json.load(fh)["reports"][0]
    a, s = r["assembly_info"], r["assembly_stats"]
    contigs = int(s["number_of_contigs"])
    scaffolds = int(s.get("number_of_scaffolds", contigs))
    contig_n50 = int(s["contig_n50"])
    scaffold_n50 = int(s.get("scaffold_n50", contig_n50))
    return [species, family, a["assembly_name"], acc, a.get("release_date", "")[:4],
            a.get("sequencing_tech", ""), a["assembly_level"],
            f"{int(s['total_sequence_length']) / 1e6:.1f}", contigs,
            f"{contig_n50 / 1e3:.1f}", scaffolds, f"{scaffold_n50 / 1e3:.1f}"]


def main(meta_json, fai, tbl, report_dir, busco_dir):
    meta = json.loads(meta_json)
    rows = [ours(fai) + busco(f"{busco_dir}/Gdor/short_summary.txt")
            + [masked_pct(tbl), "this study (RepeatMasker)"]]
    for acc, (species, family, repeat_pct, source) in meta.items():
        rows.append(published(acc, species, family, report_dir)
                    + busco(f"{busco_dir}/{acc}/short_summary.txt")
                    + ["" if repeat_pct is None else repeat_pct, source])
    print("\t".join(COLS))
    for row in rows:
        print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    main(*sys.argv[1:6])
