#!/usr/bin/env python3
"""Keep the longest protein per gene from a BRAKER protein FASTA.

BRAKER names transcripts <gene>.t<n> (e.g. g123.t2). Stop-codon asterisks are
removed, since InterProScan rejects them.
Usage: annot_longest_isoform.py braker.aa out.faa
"""
import sys


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name:
        yield name, "".join(seq)


def main():
    src, dst = sys.argv[1:3]
    best = {}
    n = 0
    for name, seq in read_fasta(src):
        n += 1
        gene = name.rsplit(".t", 1)[0]
        seq = seq.replace("*", "")
        if gene not in best or len(seq) > len(best[gene][1]):
            best[gene] = (name, seq)
    with open(dst, "w") as out:
        for name, seq in best.values():
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
    print(f"{n} transcripts, {len(best)} genes", file=sys.stderr)


if __name__ == "__main__":
    main()
