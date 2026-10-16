#!/usr/bin/env python3
"""Rotate a circular mitogenome to start at cox1 on the forward strand.

cox1 is taken from the reference GenBank records (CDS with gene cox1/COX1 or
product "cytochrome c oxidase subunit I"). The reference cox1 sequences are
BLASTed (dc-megablast) against the sequence doubled end-to-end, so a gene that
spans the old start/end is still found.
"""
import argparse
import os
import re
import subprocess
import sys

COMP = str.maketrans("ACGTacgtNnRYKMSWrykmswBDHVbdhv", "TGCAtgcaNnYRMKSWyrmkswVHDBvhdb")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    break  # single-sequence input
                name = line[1:]
            elif line:
                seq.append(line)
    return name, "".join(seq)


def cox1_from_genbank(path):
    """Yield (record id, cox1 nucleotide sequence) for each GenBank record."""
    with open(path) as fh:
        text = fh.read()
    for rec in text.split("\n//"):
        m = re.search(r"^ACCESSION\s+(\S+)", rec, re.M)
        if not m or "ORIGIN" not in rec:
            continue
        acc = m.group(1)
        origin = rec.split("ORIGIN", 1)[1]
        seq = re.sub(r"[^acgtnACGTN]", "", origin).upper()
        feats = rec.split("FEATURES", 1)[1].split("ORIGIN", 1)[0]
        # split features at their 5-space-indented keys
        for block in re.split(r"\n(?=     \S)", feats):
            if not block.strip().startswith("CDS"):
                continue
            joined = re.sub(r"\n\s+", " ", block)
            gene = re.search(r'/gene="([^"]+)"', joined)
            prod = re.search(r'/product="([^"]+)"', joined)
            g = (gene.group(1) if gene else "").lower()
            p = (prod.group(1) if prod else "").lower()
            if g not in ("cox1", "coi", "co1") and not re.fullmatch(r"cytochrome c oxidase subunit (i|1)", p):
                continue
            loc = joined.split()[1]
            nums = [int(x) for x in re.findall(r"\d+", loc)]
            start, end = min(nums), max(nums)
            s = seq[start - 1:end]
            if "complement" in loc:
                s = revcomp(s)
            yield acc, s
            break


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--genbank", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    os.makedirs(a.workdir, exist_ok=True)
    name, seq = read_fasta(a.fasta)
    n = len(seq)

    cox1 = list(cox1_from_genbank(a.genbank))
    if not cox1:
        sys.exit("No cox1 CDS found in the reference GenBank records")
    q = os.path.join(a.workdir, "ref_cox1.fasta")
    with open(q, "w") as fh:
        for acc, s in cox1:
            fh.write(f">{acc}_cox1\n{s}\n")
    s2 = os.path.join(a.workdir, "doubled.fasta")
    with open(s2, "w") as fh:
        fh.write(f">doubled\n{seq}{seq}\n")

    out = subprocess.run(
        ["blastn", "-task", "dc-megablast", "-query", q, "-subject", s2, "-evalue", "1e-10",
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"],
        check=True, capture_output=True, text=True).stdout
    with open(os.path.join(a.workdir, "cox1_hits.tsv"), "w") as fh:
        fh.write(out)
    hits = [l.split("\t") for l in out.strip().splitlines()]
    if not hits:
        sys.exit("No cox1 hit in the assembly")
    best = max(hits, key=lambda h: float(h[9]))
    qstart, sstart, send = int(best[4]), int(best[6]), int(best[7])
    print(f"best cox1 hit: {best}", file=sys.stderr)

    if sstart <= send:   # forward strand: extrapolate back to cox1 position 1
        start0 = (sstart - qstart) % n
        rotated = seq[start0:] + seq[:start0]
    else:                # reverse strand: flip, then do the same
        rc = revcomp(seq)
        pos_rc = 2 * n - sstart + 1          # hit start on the reverse complement of the doubled seq
        start0 = (pos_rc - qstart) % n
        rotated = rc[start0:] + rc[:start0]
    print(f"length {n}; strand {'+' if sstart <= send else '-'}; new start (0-based) {start0}",
          file=sys.stderr)

    with open(a.out, "w") as fh:
        fh.write(f">{name.split()[0]}\n")
        for i in range(0, n, 80):
            fh.write(rotated[i:i + 80] + "\n")


if __name__ == "__main__":
    main()
