#!/usr/bin/env python3
"""Build a GenBank feature table (.tbl) for the mitogenome.

Sources:
- protein-coding genes and tRNAs: MITOS2 (result.bed)
- tRNAs that MITOS2 missed: ARAGORN (-mt) and tRNAscan-SE (-O), added only when
  the prediction does not overlap a protein-coding or rRNA gene, or another tRNA,
  by more than --max-overlap bp. Any tRNA still missing after that is taken
  from a relaxed ARAGORN run (--aragorn-relaxed; -ps80), only if it has the
  expected anticodon, lies on the same strand as the other genes and sits in
  a gap between genes
- rRNA boundaries: scripts/mito_rrna.py (reference-projected, median over
  references), merged with the MITOS2 calls and trimmed so they do not run into
  neighbouring genes

Protein-coding genes without a legal start codon for the genetic code are
marked 5' partial with the note "start codon not determined", as in the
Schmidtea and Dugesia mitogenome records. MITOS2 tRNAs longer than
--max-trna-len are replaced by an overlapping ARAGORN prediction of the same
tRNA. MITOS2 replication-origin (OH/OL) predictions are not carried over.
"""
import argparse
import re

PCG_PRODUCT = {
    "cox1": "cytochrome c oxidase subunit I",
    "cox2": "cytochrome c oxidase subunit II",
    "cox3": "cytochrome c oxidase subunit III",
    "cob": "cytochrome b",
    "atp6": "ATP synthase F0 subunit 6",
    "atp8": "ATP synthase F0 subunit 8",
    "nad1": "NADH dehydrogenase subunit 1",
    "nad2": "NADH dehydrogenase subunit 2",
    "nad3": "NADH dehydrogenase subunit 3",
    "nad4": "NADH dehydrogenase subunit 4",
    "nad4l": "NADH dehydrogenase subunit 4L",
    "nad5": "NADH dehydrogenase subunit 5",
    "nad6": "NADH dehydrogenase subunit 6",
}
GENE_NAME = {"nad4l": "nad4L"}
RRNA_PRODUCT = {"rrnS": "12S ribosomal RNA", "rrnL": "16S ribosomal RNA"}
AA3 = {"A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly", "H": "His",
       "I": "Ile", "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln",
       "R": "Arg", "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr"}
AA1 = {v.upper(): k for k, v in AA3.items()}
# the 22 metazoan mitochondrial tRNAs; L and S each have two, told apart by anticodon
EXPECTED = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L1", "L2", "M", "N", "P", "Q", "R",
            "S1", "S2", "T", "V", "W", "Y"]


def trna_label(aa, anticodon):
    anticodon = anticodon.lower().replace("u", "t")
    if aa == "L":
        return "L1" if anticodon.endswith("ag") else "L2"   # L1 = CUN (tag), L2 = UUR (taa)
    if aa == "S":
        return "S1" if anticodon.endswith("ct") else "S2"   # S1 = AGN (gct/tct), S2 = UCN (tga)
    return aa


def read_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def span(s, e, n):
    """1-based inclusive positions covered by a feature (handles wrap-around)."""
    return set(range(s, e + 1)) if s <= e else set(range(s, n + 1)) | set(range(1, e + 1))


def parse_mitos(bed):
    pcg, trna, rrna = [], [], {}
    for line in open(bed):
        f = line.rstrip("\n").split("\t")
        # BED is 0-based half-open; a gene across the origin has start > end
        s, e, name, strand = int(f[1]) + 1, int(f[2]), f[3], f[5]
        m = re.match(r"trn([A-Z]\d?)\((\w+)\)", name)
        if m:
            aa = m.group(1)[0]
            trna.append(dict(s=s, e=e, strand=strand, label=trna_label(aa, m.group(2)),
                             anticodon=m.group(2).lower(), source="MITOS2"))
        elif re.match(r"rrn[LS]", name):
            g = name[:4]
            lo, hi = rrna.get(g, (s, e))
            rrna[g] = (min(lo, s), max(hi, e))
        elif name in PCG_PRODUCT:
            pcg.append(dict(s=s, e=e, strand=strand, name=name))
    return pcg, trna, rrna


def parse_aragorn(path):
    """ARAGORN -w output lines such as '1   tRNA-Asp   c[1234,1300]   34   (gtc)'."""
    out = []
    for line in open(path):
        m = re.match(r"\s*\d+\s+m?tRNA-(\w+)\s+(c?)\[(-?\d+),(\d+)\]\s+(\S+)\s+\((\w+)\)", line)
        if not m:
            continue
        aa3, comp, s, e, score, ac = m.groups()
        aa = AA1.get(aa3.upper()[:3])
        if not aa:
            continue
        out.append(dict(s=int(s), e=int(e), strand="-" if comp else "+",
                        label=trna_label(aa, ac), anticodon=ac.lower(), source="ARAGORN",
                        score=float(score)))
    return out


def parse_trnascan(path):
    """tRNAscan-SE tabular output (-o)."""
    out = []
    for line in open(path):
        f = line.split()
        if len(f) < 9 or not f[1].isdigit():
            continue
        s, e, aa3, ac = int(f[2]), int(f[3]), f[4], f[5]
        aa = AA1.get(aa3.upper()[:3])
        if not aa or ac.startswith("?") or ac == "NNN":
            continue
        strand = "+" if s <= e else "-"
        out.append(dict(s=min(s, e), e=max(s, e), strand=strand, label=trna_label(aa, ac),
                        anticodon=ac.lower(), source="tRNAscan-SE"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--seqid", required=True)
    ap.add_argument("--mitos-bed", required=True)
    ap.add_argument("--rrna", required=True)
    ap.add_argument("--aragorn", required=True)
    ap.add_argument("--trnascan", required=True)
    ap.add_argument("--aragorn-relaxed", help="ARAGORN -w output at a relaxed score threshold")
    ap.add_argument("--locus-prefix", required=True)
    ap.add_argument("--max-overlap", type=int, default=10)
    ap.add_argument("--max-trna-len", type=int, default=90)
    ap.add_argument("--starts", default="ATG,GTG", help="legal start codons (genetic code 9: ATG, GTG)")
    ap.add_argument("--protein-db", required=True, help="database tag for protein_id (gnl|<tag>|<locus_tag>)")
    ap.add_argument("--tbl", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()

    seq = read_seq(a.fasta)
    n = len(seq)
    pcg, trna, rrna_mitos = parse_mitos(a.mitos_bed)
    report = []
    starts = set(a.starts.split(","))
    for x in pcg:
        codon = seq[x["s"] - 1:x["s"] + 2]
        x["partial5"] = x["strand"] == "+" and codon not in starts
        if x["partial5"]:
            report.append(f"{x['name']}: start codon {codon} at {x['s']} is not a legal start; marked 5' partial")
    aragorn = parse_aragorn(a.aragorn)
    for t in trna:
        length = (t["e"] - t["s"] + 1) % n or n
        if length <= a.max_trna_len:
            continue
        alt = [c for c in aragorn if c["label"] == t["label"] and span(c["s"], c["e"], n) & span(t["s"], t["e"], n)]
        if alt:
            report.append(f"trn{t['label']}: MITOS2 {t['s']}-{t['e']} ({length} bp) replaced by ARAGORN {alt[0]['s']}-{alt[0]['e']}")
            t.update(s=alt[0]["s"], e=alt[0]["e"], source="ARAGORN (bounds)")
        else:
            report.append(f"trn{t['label']}: MITOS2 {t['s']}-{t['e']} is {length} bp; no ARAGORN alternative")

    # rRNAs: reference-projected boundaries where available, else MITOS2
    rrna = {}
    ref = {}
    for line in open(a.rrna):
        f = line.rstrip("\n").split("\t")
        if f[0] in RRNA_PRODUCT:
            ref[f[0]] = (int(f[1]), int(f[2]), int(f[3]))
    for g in ("rrnS", "rrnL"):
        if g in ref and ref[g][2] >= 2:
            rrna[g] = ref[g][:2]
            how = f"median of {ref[g][2]} reference projections"
        elif g in ref:
            lo, hi = rrna_mitos.get(g, ref[g][:2])
            rrna[g] = (min(lo, ref[g][0]), max(hi, ref[g][1]))
            how = "union of MITOS2 and one reference projection"
        else:
            rrna[g] = rrna_mitos[g]
            how = "MITOS2"
        report.append(f"{g}: {rrna[g][0]}-{rrna[g][1]} ({how}; MITOS2 {rrna_mitos.get(g)})")

    # trim rRNAs so they do not overlap protein-coding genes or MITOS2 tRNAs
    others = [(x["s"], x["e"]) for x in pcg + trna]
    for g, (s, e) in list(rrna.items()):
        for os_, oe in others:
            if os_ > oe:
                continue
            if os_ <= s <= oe:
                s = oe + 1
            if os_ <= e <= oe:
                e = os_ - 1
        if (s, e) != rrna[g]:
            report.append(f"{g}: trimmed to {s}-{e} to avoid overlapping neighbouring genes")
        rrna[g] = (s, e)

    # tRNAs MITOS2 missed
    have = {t["label"] for t in trna}
    missing = [x for x in EXPECTED if x not in have]
    report.append(f"tRNAs from MITOS2: {len(trna)}; missing: {', '.join(missing) or 'none'}")
    occupied_strict = set()
    for x in pcg:
        occupied_strict |= span(x["s"], x["e"], n)
    for s, e in rrna.values():
        occupied_strict |= span(s, e, n)
    occupied_trna = set()
    for t in trna:
        occupied_trna |= span(t["s"], t["e"], n)
    candidates = aragorn + parse_trnascan(a.trnascan)
    for lab in missing:
        cands = [c for c in candidates if c["label"] == lab]
        if not cands:
            report.append(f"  {lab}: no ARAGORN or tRNAscan-SE prediction")
            continue
        added = False
        for c in cands:
            cov = span(c["s"], c["e"], n)
            o1, o2 = len(cov & occupied_strict), len(cov & occupied_trna)
            desc = f"{c['source']} {c['s']}-{c['e']} ({c['strand']}) anticodon {c['anticodon']}"
            if o1 > a.max_overlap or o2 > a.max_overlap:
                report.append(f"  {lab}: rejected {desc}: overlaps genes by {o1} bp, tRNAs by {o2} bp")
                continue
            if not added:
                trna.append(c)
                occupied_trna |= cov
                report.append(f"  {lab}: added {desc}")
                added = True
            else:
                report.append(f"  {lab}: also predicted {desc}")

    # last resort for tRNAs still missing: relaxed ARAGORN, expected anticodon only
    expected_ac = {"A": "tgc", "C": "gca", "D": "gtc", "E": "ttc", "F": "gaa", "G": "tcc", "H": "gtg",
                   "I": "gat", "K": "ctt", "L1": "tag", "L2": "taa", "M": "cat", "N": "gtt", "P": "tgg",
                   "Q": "ttg", "R": "tcg", "S1": "tct", "S2": "tga", "T": "tgt", "V": "tac", "W": "tca",
                   "Y": "gta"}
    still = [x for x in EXPECTED if x not in {t["label"] for t in trna}]
    if still and a.aragorn_relaxed:
        main_strand = max("+-", key=lambda s: sum(1 for x in pcg if x["strand"] == s))
        relaxed = parse_aragorn(a.aragorn_relaxed)
        for lab in still:
            cands = [c for c in relaxed if c["label"] == lab and c["anticodon"] == expected_ac[lab]
                     and c["strand"] == main_strand]
            best = None
            for c in sorted(cands, key=lambda c: -c["score"]):
                cov = span(c["s"], c["e"], n)
                if len(cov & occupied_strict) or len(cov & occupied_trna) > a.max_overlap:
                    continue
                best = best or c
            if best:
                best["source"] = "ARAGORN (relaxed, -ps80)"
                trna.append(best)
                occupied_trna |= span(best["s"], best["e"], n)
                report.append(f"  {lab}: added from relaxed ARAGORN {best['s']}-{best['e']} ({best['strand']}) anticodon {best['anticodon']}")
            else:
                report.append(f"  {lab}: no relaxed ARAGORN candidate with anticodon {expected_ac[lab]} in a gap")

    feats = []
    for x in pcg:
        feats.append(("CDS", x["s"], x["e"], x["strand"], x))
    for g, (s, e) in rrna.items():
        feats.append(("rRNA", s, e, "+", g))
    for t in trna:
        feats.append(("tRNA", t["s"], t["e"], t["strand"], t))
    feats.sort(key=lambda f: f[1])

    def intervals(s, e, strand):
        iv = [(s, e)] if s <= e else [(s, n), (1, e)]
        if strand == "-":
            iv = [(b, a_) for a_, b in reversed(iv)]
        return iv

    lines = [f">Feature {a.seqid}"]
    for i, (kind, s, e, strand, obj) in enumerate(feats, 1):
        tag = f"{a.locus_prefix}_M{i:04d}"
        iv = intervals(s, e, strand)
        partial = False
        if kind == "CDS":
            name = obj["name"]
            gene, quals = GENE_NAME.get(name, name), [("product", PCG_PRODUCT[name]),
                                                      ("protein_id", f"gnl|{a.protein_db}|{tag}")]
            if obj["partial5"]:
                partial = True
                quals.append(("note", "start codon not determined"))
        elif kind == "rRNA":
            gene, quals = obj, [("product", RRNA_PRODUCT[obj]),
                                ("note", "boundaries inferred by similarity to other planarian mitogenomes")]
        else:
            aa = obj["label"][0]
            gene = "trn" + obj["label"]
            quals = [("product", f"tRNA-{AA3[aa]}")]
            if obj["source"] == "ARAGORN (bounds)":
                quals.append(("note", "boundaries from ARAGORN"))
            elif obj["source"] != "MITOS2":
                quals.append(("note", f"predicted by {obj['source']}; anticodon {obj['anticodon']}"))
        def first(p, k):
            return f"<{p}" if partial and k == 0 else str(p)
        for k, (p, q) in enumerate(iv):
            lines.append(f"{first(p, k)}\t{q}" + ("\tgene" if k == 0 else ""))
        lines += [f"\t\t\tgene\t{gene}", f"\t\t\tlocus_tag\t{tag}"]
        for k, (p, q) in enumerate(iv):
            lines.append(f"{first(p, k)}\t{q}" + (f"\t{kind}" if k == 0 else ""))
        lines += [f"\t\t\t{k_}\t{v}" for k_, v in quals]

    open(a.tbl, "w").write("\n".join(lines) + "\n")
    counts = {k: sum(1 for f in feats if f[0] == k) for k in ("CDS", "rRNA", "tRNA")}
    report.append(f"Feature table: {counts['CDS']} CDS, {counts['rRNA']} rRNA, {counts['tRNA']} tRNA")
    open(a.report, "w").write("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
