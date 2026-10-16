#!/usr/bin/env python3
"""Apply the FCS results to the purged assembly to make the submission FASTA.

- FCS-GX: contigs with action EXCLUDE are removed (written to --contam-out).
- FCS-adaptor: adaptor spans at contig ends (fcs_adaptor_report.txt,
  ACTION_TRIM) are trimmed. Internal spans that FCS-adaptor skipped
  (skipped_trims.jsonl) are cut out and the contig is split there; an
  internal adaptor in a nanopore assembly marks a chimeric join.
- Leading and trailing Ns are trimmed (NCBI rejects sequences that start or
  end with N).
- Pieces shorter than --min-len are dropped.

Every change is listed in --log. Split pieces are named <contig>_s1, _s2, ...
"""
import argparse
import json
import re


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def parse_ranges(text):
    """'1..37,8000..8030' -> [(1, 37), (8000, 8030)] (1-based, inclusive)."""
    return [tuple(int(x) for x in r.split("..")) for r in text.split(",") if r]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gx-report", required=True)
    ap.add_argument("--adaptor-report", required=True)
    ap.add_argument("--skipped-trims", required=True)
    ap.add_argument("--min-len", type=int, default=200)
    ap.add_argument("--out", required=True)
    ap.add_argument("--contam-out", required=True)
    ap.add_argument("--log", required=True)
    a = ap.parse_args()

    exclude = {}
    for line in open(a.gx_report):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[4] == "EXCLUDE":
            exclude[f[0]] = f"{f[5]}; {f[7]}"
        elif f[4] != "REVIEW":
            raise SystemExit(f"Unhandled FCS-GX action {f[4]} for {f[0]}")

    cuts = {}   # contig -> list of (start, end) spans to remove
    for line in open(a.adaptor_report):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] == "ACTION_TRIM":
            cuts.setdefault(f[0], []).extend(parse_ranges(f[3]))
        elif f[2] == "ACTION_EXCLUDE":
            exclude[f[0]] = "adaptor: " + f[4]
        else:
            raise SystemExit(f"Unhandled FCS-adaptor action {f[2]} for {f[0]}")
    for line in open(a.skipped_trims):
        d = json.loads(line)
        for s in d["spans"]:
            cuts.setdefault(d["seqId"], []).append((int(s["range"]["start"]), int(s["range"]["stop"])))

    log = []
    n_in = n_out = bp_in = bp_out = 0
    with open(a.out, "w") as out, open(a.contam_out, "w") as contam:
        for name, seq in read_fasta(a.fasta):
            name = name.removeprefix("lcl|")
            n_in += 1
            bp_in += len(seq)
            if name in exclude:
                contam.write(f">{name}\n{seq}\n")
                log.append(f"EXCLUDE\t{name}\t{len(seq)}\t{exclude[name]}")
                continue
            pieces = [seq]
            if name in cuts:
                spans = sorted(set(cuts[name]))
                pieces, prev = [], 0
                for s, e in spans:
                    pieces.append(seq[prev:s - 1])
                    prev = e
                pieces.append(seq[prev:])
                log.append(f"ADAPTOR\t{name}\t{len(seq)}\tremoved {';'.join(f'{s}..{e}' for s, e in spans)}")
            kept = []
            for p in pieces:
                q = p.strip("Nn")
                if len(q) != len(p) and q:
                    log.append(f"NTRIM\t{name}\t{len(p)}\ttrimmed {len(p) - len(q)} terminal N")
                if len(q) >= a.min_len:
                    kept.append(q)
                elif q:
                    log.append(f"SHORT\t{name}\t{len(q)}\tpiece dropped (< {a.min_len} bp)")
            for i, q in enumerate(kept, 1):
                pid = name if len(kept) == 1 else f"{name}_s{i}"
                if len(kept) > 1:
                    log.append(f"SPLIT\t{pid}\t{len(q)}\tpiece {i} of {len(kept)}")
                out.write(f">{pid}\n")
                for j in range(0, len(q), 80):
                    out.write(q[j:j + 80] + "\n")
                n_out += 1
                bp_out += len(q)

    with open(a.log, "w") as fh:
        fh.write(f"# input: {n_in} sequences, {bp_in} bp\n")
        fh.write(f"# output: {n_out} sequences, {bp_out} bp\n")
        fh.write("# action\tsequence\tlength\tdetail\n")
        fh.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
