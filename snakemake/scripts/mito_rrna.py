#!/usr/bin/env python3
"""Estimate rRNA gene boundaries in a mitogenome from related annotated mitogenomes.

rrnL and rrnS are taken from the reference GenBank records (rRNA features whose
product names the large or small subunit). Each reference gene is aligned to
the query (blastn, dc-megablast), and the best hit near the existing annotation
is extended to the reference gene's full length using its query coordinates.
The reported boundary is the median over references; the per-reference values
are written alongside it so the spread can be checked.
"""
import argparse
import os
import re
import statistics
import subprocess

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def rrna_from_genbank(path):
    """Yield (accession, organism, gene, sequence) for each rRNA feature."""
    with open(path) as fh:
        text = fh.read()
    for rec in text.split("\n//"):
        m = re.search(r"^ACCESSION\s+(\S+)", rec, re.M)
        if not m or "ORIGIN" not in rec:
            continue
        acc = m.group(1)
        org = re.search(r"^\s+ORGANISM\s+(.+)$", rec, re.M).group(1).strip()
        seq = re.sub(r"[^acgtnACGTN]", "", rec.split("ORIGIN", 1)[1]).upper()
        feats = rec.split("FEATURES", 1)[1].split("ORIGIN", 1)[0]
        for block in re.split(r"\n(?=     \S)", feats):
            if not block.strip().startswith("rRNA"):
                continue
            joined = re.sub(r"\n\s+", " ", block)
            prod = re.search(r'/product="([^"]+)"', joined)
            gene = re.search(r'/gene="([^"]+)"', joined)
            label = ((prod.group(1) if prod else "") + " " + (gene.group(1) if gene else "")).lower()
            if re.search(r"large|16s|rrnl|l-rrna|lrrna", label):
                name = "rrnL"
            elif re.search(r"small|12s|rrns|s-rrna|srrna", label):
                name = "rrnS"
            else:
                continue
            loc = joined.split()[1]
            nums = [int(x) for x in re.findall(r"\d+", loc)]
            s = seq[min(nums) - 1:max(nums)]
            if "complement" in loc:
                s = revcomp(s)
            yield acc, org, name, s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, help="mitogenome (single sequence)")
    ap.add_argument("--genbank", required=True, help="annotated reference mitogenomes")
    ap.add_argument("--bed", required=True, help="MITOS2 result.bed (for the expected location)")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True, help="TSV: gene, start, end, strand, n_refs, per-reference values")
    ap.add_argument("--window", type=int, default=1500, help="max distance of a hit from the MITOS2 annotation")
    a = ap.parse_args()

    os.makedirs(a.workdir, exist_ok=True)
    seq = read_fasta(a.fasta)
    n = len(seq)

    # MITOS2 location of each rRNA (pieces such as rrnS_0/rrnS_1 are merged)
    mitos = {}
    with open(a.bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            m = re.match(r"(rrn[LS])", f[3])
            if m:
                s, e = int(f[1]), int(f[2])
                lo, hi = mitos.get(m.group(1), (s, e))
                mitos[m.group(1)] = (min(lo, s), max(hi, e))

    refs = list(rrna_from_genbank(a.genbank))
    q = os.path.join(a.workdir, "ref_rrna.fasta")
    qlen = {}
    with open(q, "w") as fh:
        for acc, org, name, s in refs:
            qid = f"{name}|{acc}|{org.replace(' ', '_')}"
            qlen[qid] = len(s)
            fh.write(f">{qid}\n{s}\n")

    out = subprocess.run(
        ["blastn", "-task", "dc-megablast", "-query", q, "-subject", a.fasta, "-evalue", "1e-5",
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"],
        check=True, capture_output=True, text=True).stdout
    with open(os.path.join(a.workdir, "rrna_hits.tsv"), "w") as fh:
        fh.write(out)

    best = {}
    for line in out.strip().splitlines():
        f = line.split("\t")
        qid = f[0]
        gene = qid.split("|")[0]
        qs, qe, ss, se, bits = int(f[4]), int(f[5]), int(f[6]), int(f[7]), float(f[9])
        if ss > se:          # all genes are on the forward strand in these genomes
            continue
        lo, hi = mitos.get(gene, (0, n))
        if ss > hi + a.window or se < lo - a.window:
            continue
        # project the hit to the full reference gene length
        start = ss - (qs - 1)
        end = se + (qlen[qid] - qe)
        if qid not in best or bits > best[qid][2]:
            best[qid] = (start, end, bits, float(f[2]), qe - qs + 1)

    with open(a.out, "w") as fh:
        fh.write("gene\tstart\tend\tn_refs\tmitos_start\tmitos_end\tper_reference\n")
        for gene in ("rrnL", "rrnS"):
            vals = [(qid, v) for qid, v in best.items() if qid.startswith(gene + "|")]
            if not vals:
                continue
            start = int(statistics.median(v[0] for _, v in vals))
            end = int(statistics.median(v[1] for _, v in vals))
            lo, hi = mitos.get(gene, (None, None))
            per = "; ".join(f"{qid.split('|')[2]}:{v[0]}-{v[1]} ({v[3]:.0f}% over {v[4]} bp)" for qid, v in vals)
            fh.write(f"{gene}\t{start}\t{end}\t{len(vals)}\t{'' if lo is None else lo + 1}\t{'' if hi is None else hi}\t{per}\n")


if __name__ == "__main__":
    main()
