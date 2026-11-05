#!/usr/bin/env python3
"""Summary statistics for a BRAKER GTF (genes, transcripts, exons, introns, lengths).
Usage: annot_gene_stats.py braker.gtf > gene_stats.tsv
"""
import re
import statistics
import sys
from collections import defaultdict


def main():
    cds = defaultdict(list)      # transcript -> [(start, end)]
    exons = defaultdict(list)
    tx_gene = {}
    with open(sys.argv[1]) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            g = re.search(r'gene_id "([^"]+)"', f[8])
            if not m:
                continue
            tx = m.group(1)
            if g:
                tx_gene[tx] = g.group(1)
            s, e = int(f[3]), int(f[4])
            if f[2] == "CDS":
                cds[tx].append((s, e))
            elif f[2] == "exon":
                exons[tx].append((s, e))

    feats = exons if exons else cds
    genes = set(tx_gene.values()) or {t.rsplit(".t", 1)[0] for t in feats}
    tx_len, cds_len, n_exon, intron_len = [], [], [], []
    single = 0
    for tx, iv in feats.items():
        iv.sort()
        n_exon.append(len(iv))
        tx_len.append(iv[-1][1] - iv[0][0] + 1)
        intron_len += [b[0] - a[1] - 1 for a, b in zip(iv, iv[1:])]
        if len(iv) == 1:
            single += 1
    for tx, iv in cds.items():
        cds_len.append(sum(e - s + 1 for s, e in iv))

    rows = [
        ("genes", len(genes)),
        ("transcripts", len(feats)),
        ("transcripts_per_gene", round(len(feats) / max(1, len(genes)), 2)),
        ("mean_exons_per_transcript", round(statistics.mean(n_exon), 2) if n_exon else 0),
        ("single_exon_transcripts", single),
        ("median_transcript_span_bp", int(statistics.median(tx_len)) if tx_len else 0),
        ("median_cds_length_bp", int(statistics.median(cds_len)) if cds_len else 0),
        ("mean_cds_length_bp", int(statistics.mean(cds_len)) if cds_len else 0),
        ("introns", len(intron_len)),
        ("median_intron_length_bp", int(statistics.median(intron_len)) if intron_len else 0),
        ("mean_intron_length_bp", int(statistics.mean(intron_len)) if intron_len else 0),
    ]
    print("metric\tvalue")
    for k, v in rows:
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
