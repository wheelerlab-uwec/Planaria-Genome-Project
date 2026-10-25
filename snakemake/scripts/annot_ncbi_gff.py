#!/usr/bin/env python3
"""Rewrite BRAKER's GFF3 into the form table2asn expects for a GenBank genome.

- genes get locus tags <prefix>_000010, _000020, ... in file order (steps of 10
  leave room for genes added later)
- mRNAs and CDSs get products (from annot_ncbi_products.py) and NCBI ids:
  transcript_id gnl|<db>|mrna.<locus_tag>-T<n>, protein_id gnl|<db>|<locus_tag>-T<n>
- all segments of a CDS share one ID, as table2asn requires
- models without a start or stop codon (genes running off a contig end) are
  marked partial with start_range/end_range
- start_codon, stop_codon and intron rows are dropped
- transcripts are dropped (and genes with none left) when the protein is
  shorter than --min-aa, when a missing start/stop is not within 3 bp of a
  contig end (an incomplete model, not a gene cut off by the assembly), or
  when a gene repeats another gene's exact coding structure

Usage: annot_ncbi_gff.py braker.gff3 products.tsv PREFIX DB contig_lengths.tsv
                         out.gff3 locus_tags.tsv [min_aa]
"""
import re
import sys
from collections import defaultdict


def attrs(col):
    return dict(kv.split("=", 1) for kv in col.strip().strip(";").split(";") if kv)


def esc(v):
    # GFF3 reserves , ; = in attribute values (a comma splits multiple values)
    return v.replace("%", "%25").replace(",", "%2C").replace(";", "%3B").replace("=", "%3D")


def fmt(a):
    return ";".join(f"{k}={esc(v) if k == 'product' else v}" for k, v in a.items())


def main():
    gff, products, prefix, db, lengths, out, mapping = sys.argv[1:8]
    min_aa = int(sys.argv[8]) if len(sys.argv) > 8 else 30
    contig_len = {}
    with open(lengths) as fh:
        for line in fh:
            name, n = line.split("\t")[:2]
            contig_len[name.split()[0]] = int(n)
    product = {}
    with open(products) as fh:
        next(fh)
        for line in fh:
            g, p = line.split("\t")[:2]
            product[g] = p

    genes, mrnas = {}, {}
    order = []
    children = defaultdict(list)         # mRNA -> exon/CDS rows
    codons = defaultdict(set)            # mRNA -> {"start_codon", "stop_codon"}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            a = attrs(f[8])
            if f[2] == "gene":
                genes[a["ID"]] = f
                order.append(a["ID"])
            elif f[2] == "mRNA":
                mrnas[a["ID"]] = (a["Parent"], f)
            elif f[2] in ("exon", "CDS"):
                children[a["Parent"]].append(f)
            elif f[2] in ("start_codon", "stop_codon"):
                codons[a["Parent"]].add(f[2])

    # Filter transcripts, then genes
    dropped = defaultdict(int)
    by_gene = defaultdict(list)
    seen_cds = set()
    for mid, (gid, f) in mrnas.items():
        cds = sorted((int(r[3]), int(r[4])) for r in children[mid] if r[2] == "CDS")
        aa = sum(e - s + 1 for s, e in cds) // 3 - ("stop_codon" in codons[mid])
        if aa < min_aa:
            dropped["short protein"] += 1
            continue
        strand, start, end = f[6], int(f[3]), int(f[4])
        no_start = "start_codon" not in codons[mid]
        no_stop = "stop_codon" not in codons[mid]
        low_open = no_start if strand == "+" else no_stop
        high_open = no_stop if strand == "+" else no_start
        if (low_open and start > 3) or (high_open and end < contig_len[f[0]] - 3):
            dropped["incomplete, not at a contig end"] += 1
            continue
        key = (f[0], strand, tuple(cds))
        if key in seen_cds:
            dropped["duplicate coding structure"] += 1
            continue
        seen_cds.add(key)
        by_gene[gid].append(mid)
    n_genes = len(order)
    order = [g for g in order if by_gene[g]]
    # Gene spans follow the transcripts that are kept
    for gid in order:
        kept = [mrnas[m][1] for m in by_gene[gid]]
        genes[gid][3] = str(min(int(k[3]) for k in kept))
        genes[gid][4] = str(max(int(k[4]) for k in kept))

    n_partial = 0
    with open(out, "w") as o, open(mapping, "w") as m:
        o.write("##gff-version 3\n")
        m.write("braker_gene\tlocus_tag\tproduct\n")
        for i, gid in enumerate(order, 1):
            locus = f"{prefix}_{i * 10:06d}"
            prod = product.get(gid, "hypothetical protein")
            m.write(f"{gid}\t{locus}\t{prod}\n")
            g = genes[gid]
            gstart, gend = int(g[3]), int(g[4])
            lines = []
            g_ranges = {}
            for mid in sorted(by_gene[gid], key=lambda x: int(x.rsplit(".t", 1)[1])):
                n = mid.rsplit(".t", 1)[1]
                f = mrnas[mid][1]
                strand, start, end = f[6], int(f[3]), int(f[4])
                no_start = "start_codon" not in codons[mid]
                no_stop = "stop_codon" not in codons[mid]
                # start_range/end_range refer to the low/high coordinate
                low_open = no_start if strand == "+" else no_stop
                high_open = no_stop if strand == "+" else no_start
                ranges = {}
                if low_open:
                    ranges["start_range"] = f".,{start}"
                    if start == gstart:
                        g_ranges["start_range"] = f".,{gstart}"
                if high_open:
                    ranges["end_range"] = f"{end},."
                    if end == gend:
                        g_ranges["end_range"] = f"{gend},."
                n_partial += bool(ranges)
                tid = f"gnl|{db}|mrna.{locus}-T{n}"
                pid = f"gnl|{db}|{locus}-T{n}"
                ma = {"ID": mid, "Parent": gid, "product": prod, "transcript_id": tid,
                      "protein_id": pid, **ranges}
                lines.append(f[:8] + [fmt(ma)])
                rows = sorted(children[mid], key=lambda r: (int(r[3]), r[2]))
                for k, r in enumerate([r for r in rows if r[2] == "exon"], 1):
                    lines.append(r[:7] + ["."] + [fmt({"ID": f"{mid}.exon{k}", "Parent": mid})])
                cds = [r for r in rows if r[2] == "CDS"]
                for r in cds:
                    ca = {"ID": f"cds.{locus}-T{n}", "Parent": mid, "product": prod,
                          "protein_id": pid}
                    if low_open and r is cds[0]:
                        ca["start_range"] = f".,{r[3]}"
                    if high_open and r is cds[-1]:
                        ca["end_range"] = f"{r[4]},."
                    lines.append(r[:8] + [fmt(ca)])
            ga = {"ID": gid, "locus_tag": locus, **g_ranges}
            o.write("\t".join(g[:8] + [fmt(ga)]) + "\n")
            for r in lines:
                o.write("\t".join(r) + "\n")
    print(f"{len(order)} genes kept of {n_genes}; "
          f"{sum(len(v) for v in by_gene.values())} mRNAs kept of {len(mrnas)}; "
          f"{n_partial} partial models", file=sys.stderr)
    for why, n in dropped.items():
        print(f"  dropped {n} transcripts: {why}", file=sys.stderr)


if __name__ == "__main__":
    main()
