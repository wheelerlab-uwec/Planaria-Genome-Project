"""Histogram of alt-allele frequency at SNP sites, split by contig depth class."""
import sys, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

selected, site_dir, out_png = sys.argv[1:4]
cls = {}
for line in open(selected):
    name, _, _, c = line.rstrip("\n").split("\t")
    cls[name] = c

freqs = {"collapsed": [], "haplotype": []}
depths = {"collapsed": [], "haplotype": []}
for f in glob.glob(os.path.join(site_dir, "*.tsv")):
    for line in open(f):
        chrom, pos, ad = line.rstrip("\n").split("\t")
        counts = [int(x) for x in ad.split(",") if x != "."]
        dp = sum(counts)
        if dp == 0:
            continue
        alt = max(counts[1:]) if len(counts) > 1 else 0
        depths[cls[chrom]].append(dp)
        freqs[cls[chrom]].append((alt, dp))

fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
for ax, c in zip(axes, ["collapsed", "haplotype"]):
    if not depths[c]:
        print(f"{c}: no sites")
        continue
    med = np.median(depths[c])
    af = np.array([a / d for a, d in freqs[c] if 0.6 * med <= d <= 1.4 * med and a >= 3])
    af = af[(af >= 0.08) & (af <= 0.98)]
    ax.hist(af, bins=np.arange(0.08, 1.0, 0.01), color="steelblue")
    for x in (0.25, 1/3, 0.5, 2/3, 0.75):
        ax.axvline(x, color="grey", ls=":", lw=0.8)
    ax.set_title(f"{c} contigs (median site depth {med:.0f}x, {len(af):,} SNP sites)")
    ax.set_ylabel("sites")
    hist, edges = np.histogram(af, bins=np.arange(0.05, 1.0, 0.05))
    print(f"== {c}: median depth {med:.0f}, sites {len(af)}")
    for h, e in zip(hist, edges):
        print(f"  {e:.2f}-{e+0.05:.2f}\t{h}")
axes[1].set_xlabel("alt allele frequency (dotted: 1/4, 1/3, 1/2, 2/3, 3/4)")
fig.tight_layout()
fig.savefig(out_png, dpi=150)
