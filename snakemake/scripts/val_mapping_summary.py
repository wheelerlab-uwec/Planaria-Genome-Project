#!/usr/bin/env python3
"""Summarise read mapping back to the assembly as a two-column TSV.

From `samtools stats` (SN lines): reads and read bases, and the fraction
mapped. The NM-based "error rate" is not reported: NM counts every base of
long gaps, so it is inflated far above read error.

From the per-base depth histogram (depth, count; depths >= 1000 pooled):
mean and median depth, and the share of the assembly at zero depth, below
10x, and within half to twice the median.

From a read subset aligned to the primary assembly and to the purged
haplotigs (base-level PAF): per-read divergence from the minimap2 de tag
(gap-compressed, so each indel counts once), using each read's longest
primary alignment. For Q18-20 ONT reads this is mostly read error, so its
median is an upper bound on consensus error, not a QV. Also the share of
reads closer to a haplotig than to the primary assembly (de lower by more
than 0.5 percentage points), i.e. reads from the purged haplotype.

Usage: val_mapping_summary.py samtools_stats.txt depth_histogram.tsv
                              subset_vs_primary.paf.gz subset_vs_haplotigs.paf.gz
"""
import gzip
import sys


def read_sn(path):
    sn = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("SN\t"):
                f = line.rstrip("\n").split("\t")
                sn[f[1].rstrip(":")] = float(f[2])
    return sn


def read_hist(path):
    hist = []
    with open(path) as fh:
        for line in fh:
            d, n = line.split()
            hist.append((int(d), int(n)))
    return sorted(hist)


def read_divergence(path):
    """de of each read's longest primary alignment: {read: de}"""
    best = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if "tp:A:P" not in f[12:]:
                continue
            de = next(float(t[5:]) for t in f[12:] if t.startswith("de:f:"))
            aln = int(f[3]) - int(f[2])
            if f[0] not in best or aln > best[f[0]][1]:
                best[f[0]] = (de, aln)
    return {q: v[0] for q, v in best.items()}


def quantile(sorted_values, p):
    return sorted_values[int(p * (len(sorted_values) - 1))]


def main(stats_path, hist_path, primary_paf, haplotig_paf):
    sn = read_sn(stats_path)
    hist = read_hist(hist_path)

    reads = sn["raw total sequences"]
    bases = sn["total length"]
    rows = [
        ("reads", f"{reads:.0f}"),
        ("reads_mapped", f"{sn['reads mapped']:.0f}"),
        ("reads_mapped_pct", f"{100 * sn['reads mapped'] / reads:.2f}"),
        ("read_bases", f"{bases:.0f}"),
        ("bases_mapped_cigar", f"{sn['bases mapped (cigar)']:.0f}"),
        ("bases_mapped_pct", f"{100 * sn['bases mapped (cigar)'] / bases:.2f}"),
        ("supplementary_alignments", f"{sn['supplementary alignments']:.0f}"),
    ]

    total = sum(n for _, n in hist)
    mean = sum(d * n for d, n in hist) / total
    half, cum, median = total / 2, 0, None
    for d, n in hist:
        cum += n
        if cum >= half:
            median = d
            break
    zero = sum(n for d, n in hist if d == 0)
    low = sum(n for d, n in hist if d < 10)
    near = sum(n for d, n in hist if 0.5 * median <= d <= 2 * median)
    rows += [
        ("assembly_bp", f"{total}"),
        ("depth_mean", f"{mean:.1f}"),
        ("depth_median", f"{median}"),
        ("bp_depth0_pct", f"{100 * zero / total:.3f}"),
        ("bp_depth_lt10_pct", f"{100 * low / total:.3f}"),
        ("bp_half_to_2x_median_pct", f"{100 * near / total:.2f}"),
    ]

    prim = read_divergence(primary_paf)
    hap = read_divergence(haplotig_paf)
    de = sorted(prim.values())
    both = [q for q in prim if q in hap]
    closer = sum(1 for q in both if hap[q] + 0.005 < prim[q])
    rows += [
        ("subset_reads_aligned", f"{len(prim)}"),
        ("subset_divergence_median_pct", f"{100 * quantile(de, 0.5):.2f}"),
        ("subset_divergence_p10_pct", f"{100 * quantile(de, 0.1):.2f}"),
        ("subset_divergence_p90_pct", f"{100 * quantile(de, 0.9):.2f}"),
        ("subset_closer_to_haplotig_pct", f"{100 * closer / len(both):.1f}"),
    ]

    print("metric\tvalue")
    for k, v in rows:
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main(*sys.argv[1:5])
