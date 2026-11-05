#!/usr/bin/env python3
"""Assign NCBI-style protein product names, following the International
Protein Nomenclature Guidelines (IPNG) that NCBI uses.

Each gene gets the first name available from, in order:
  1. swissprot  - eggNOG ortholog symbol (orthologous group annotated at the
                  Metazoa level) -> recommended name of that gene's
                  UniProtKB/Swiss-Prot entry (IPNG 3A). Only for symbols given
                  to at most MAX_PER_SYMBOL genes: a name like "TNF
                  receptor-associated factor 6" on dozens of paralogs would
                  claim orthology that eggNOG does not establish.
  2. eggnog     - eggNOG description, when it is a protein name rather than a
                  GO term, phrase or family/domain label (Metazoa level only,
                  same symbol limit)
  3. family     - InterPro Family entry -> "X family protein" (IPNG 3A/3B)
  4. domain     - InterPro Domain/Repeat entry -> "X domain-containing protein"
                  or "X repeat-containing protein" (IPNG 3B: from signature
                  matches, not BLAST hits)
  5. hypothetical protein
Where several InterPro entries qualify, the one covering most of the protein
wins. eggNOG and InterProScan were run on the longest isoform; the name applies
to every isoform of the gene. Every candidate goes through clean(), which
applies the IPNG style rules and rejects what cannot be made into a name.

Usage: annot_ncbi_products.py eggnog.annotations swissprot_names.tsv
                              interproscan.tsv interpro_entry.list braker.gff3 out.tsv
Output columns: gene, product, source, evidence
"""
import re
import sys
from collections import defaultdict

HYPO = "hypothetical protein"
ORTHOLOGY_LEVELS = {"33208"}            # eggNOG max_annot_lvl taxid: Metazoa
MAX_PER_SYMBOL = 3                      # allows a paralog or a leftover haplotig
# Swiss-Prot organisms preferred when a symbol maps to several entries
ORG_RANK = {"9606": 0, "10090": 1, "10116": 2, "7955": 3, "8355": 4,
            "7227": 5, "6239": 6}

# GO-style phrases (functions, processes, locations) that are not protein names
GO_LIKE = re.compile(
    r"\bactivity\b|\bbinding\b(?!.*\b(protein|domain|factor|site)\b)|"
    r"^(positive |negative )?regulation of\b|^response to\b|^cellular\b|"
    r"\bprocess\b|\bpathway\b|^structural constituent\b|"
    r"\b(transport|signaling|signalling|assembly|organization|organisation|"
    r"development|differentiation|localization|metabolism|biosynthesis|"
    r"catabolism|morphogenesis|homeostasis|proliferation|secretion|"
    r"transcription|translation|repair|replication|import|export|adhesion|"
    r"apoptosis|division|segregation|cycle|formation|autophagy|chemotaxis|"
    r"migration|polarity|growth|maturation|silencing|splicing|folding|"
    r"degradation|ubiquitination|phosphorylation|methylation)$|"
    r"\b(complex|membrane|nucleus|cytoplasm|organelle|cytoskeleton|junction)$|"
    r"^(stabilization|maintenance|normal)\b|"
    r"\bcatabolic\b|\bbiosynthetic\b|\bmetabolic\b", re.I)
# Sentences, phrases and organism references (IPNG 2E, 2F)
PHRASE = re.compile(
    r"^(to|and|or|from|in|of|with|for|binds?|mediates?|plays?|required|acts|"
    r"contains|catalyz\w*|involved|repeats?\s+(in|found)|[-'\"])\b|"
    r"\b(involved in|related to|similar to|is|are|was)\b|"
    r"\b(saccharomyces|cerevisiae|e\.\s?coli|escherichia|drosophila|mouse|"
    r"human|worm|fly|yeast|bacteriophage|genome|arabidopsis|elegans)\b|"
    r"\b(and|or|from|of|in|to|with|for|the|the protein)$|\(\s*\)|\bparalog\b|"
    r"^[A-Z0-9]+\.\d+", re.I)
UNNAMED = re.compile(
    r"\b(unknown function|uncharacteri[sz]ed|hypothetical|open reading frame|"
    r"KIAA\d+|RIKEN|zgc|DUF\d+|UPF\d+)\b|^protein \w+$|^(vacuolar|membrane) protein$",
    re.I)
# Labels for domains and families (left to InterPro, IPNG 3B)
DOMAIN_LABEL = re.compile(r"\b(domain|region|motif|repeats?|box|finger|terminal|"
                          r"fold|signature|family|superfamily|subfamily)$", re.I)
# Lower-case text with none of these is a term, not a protein name
PROTEIN_WORD = re.compile(r"(ase|in|ins)\b|\b(protein|factor|receptor|channel|"
                          r"transporter|cytochrome|subunit|domain|family|"
                          r"superfamily|repeat|finger|lectin|histone|toxin|"
                          r"globin|\w+-like)\b|\d", re.I)
COMMON_WORDS = {"family", "protein", "domain", "repair", "recombination",
                "substrate", "diverse", "paralog", "factor", "binding",
                "containing", "subunit", "like", "associated", "related", "type",
                "class", "group", "complex", "chain", "site", "motif", "repeat",
                "plasma", "bovine", "pancreatic", "inhibitor"}
SPELLING = (("isation", "ization"), ("ising", "izing"), ("disulphide", "disulfide"),
            ("sulphur", "sulfur"), ("tumour", "tumor"), ("signalling", "signaling"),
            ("haem", "hem"), ("oestr", "estr"), ("aemia", "emia"))


def clean(name):
    """IPNG style: one concise lowercase name, American spelling, no
    uncertainty words, 'homolog', molecular weights, commas or plurals."""
    d = re.sub(r"\s+", " ", name).strip().rstrip(".").strip()
    if re.search(r"\b[A-Z]\. [a-z]", d):                              # "S. mansoni" (2E)
        return None
    d = d.split(". ")[0]                                             # first sentence
    if d in ("", "-") or re.match(r"[^(]*\)", d):
        return None
    if d.count("(") != d.count(")"):
        d = d[:d.rfind("(")].rstrip(" ,")
    d = re.sub(r"\s*\((EC [^()]*|[^()]*\bkDa?)\)", "", d)          # EC numbers, MW
    d = re.sub(r"[\s,-]*\b\d+(\.\d+)?\s?kDa?\b,?", " ", d).strip()  # MW (2E)
    d = re.sub(r"\s*\([^()]*\)$", "", d).strip(" ,")
    for a, b in SPELLING:
        d = d.replace(a, b)
    d = re.sub(r"^(putative|probable|possible|potential|predicted|likely)\s+",
               "", d, flags=re.I)                                    # 2F, 3B
    d = re.sub(r"\s+homolog(ue)?\b", "", d, flags=re.I)               # 2F
    d = re.sub(r"\s+(metazoa|metazoan|vertebrate|eukaryota)\b", "", d, flags=re.I)  # 2E
    d = re.sub(r"\s*,\s*", " ", d)                                   # 2C
    d = re.sub(r"[;:=%@\\]", " ", d)
    d = re.sub(r"\brepeats\b", "repeat", d)                          # 2A
    d = re.sub(r"\bdomain containing\b", "domain-containing", d, flags=re.I)
    d = re.sub(r"\s+isoform(\s+\w+)?$", "", d)
    d = re.sub(r"(ase)\s+(protein|enzyme)$", r"\1", d)               # 2E
    d = re.sub(r"\s+(precursor|subtype)$", "", d, flags=re.I)
    # Title Case from some InterPro names: lowercase common words, keep
    # proper nouns and symbols ("von Willebrand", "AlkB")
    words = re.split(r"([ /])", d)
    for i in range(2, len(words), 2):
        if re.fullmatch(r"[A-Z][a-z]+", words[i]) and (
                words[i].lower() in COMMON_WORDS
                or re.search(r"(ase|ases|ion|or|er|ing|ive|al|ated)$", words[i])):
            words[i] = words[i].lower()
    d = "".join(words)
    d = re.sub(r"\s+", " ", d).strip().rstrip(".-").strip()
    if (not d or GO_LIKE.search(d) or PHRASE.search(d) or UNNAMED.search(d)
            or len(d) > 100):
        return None
    if d[0].islower() and not PROTEIN_WORD.search(d):
        return None
    if len(d.split()) == 1 and not re.search(r"[a-z]{4,}", d):
        return None                     # bare symbols and acronyms ("TOG")
    if re.match(r"[A-Z][a-z]{2,}", d) and not re.match(r"[A-Z][a-z]+[A-Z]", d):
        d = d[0].lower() + d[1:]        # lowercase except acronyms (2E)
    return d


def recommended_name(field):
    """First name of a UniProt 'Protein names' field: 'Name (Alt) [Cleaved into: ...]'."""
    depth = 0
    for i, c in enumerate(field):
        if c in "([":
            if depth == 0 and i > 0 and field[i - 1] == " ":
                return field[:i - 1]
            depth += 1
        elif c in ")]":
            depth -= 1
    return field


def family_name(desc):
    d = clean(desc)
    if not d:
        return None
    if re.search(r"\b(family|superfamily|subfamily)$", d, re.I):
        return d + " protein"
    if re.search(r"\b(family|superfamily) (protein|member)$|-containing protein$", d, re.I):
        return d
    return d + " family protein"


def domain_name(desc):
    d = clean(desc)
    if not d:
        return None
    if re.search(r"\b(domain-containing protein|repeat-containing protein)$", d):
        return d
    d = re.sub(r"\s+(domain|region)$", "", d, flags=re.I)
    if re.search(r"\b(repeat|motif)$", d, re.I):
        return d + "-containing protein"
    return d + " domain-containing protein"


def merged_length(spans):
    total, end = 0, 0
    for s, e in sorted(spans):
        if e > end:
            total += e - max(s, end + 1) + 1
            end = e
    return total


def main():
    ann, sprot, ipr_tsv, entry_list, gff, out = sys.argv[1:7]

    # Swiss-Prot: symbol (case-insensitive) -> recommended name, by organism rank
    sp = {}
    with open(sprot) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or not f[1]:
                continue
            sym, rank = f[1].upper(), ORG_RANK.get(f[3], 99)
            if sym not in sp or rank < sp[sym][0]:
                sp[sym] = (rank, recommended_name(f[2]))

    egg = {}
    with open(ann) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            egg[f[0].rsplit(".t", 1)[0]] = (f[5].split("|")[0], f[8], f[7])

    etype = {}
    with open(entry_list) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            etype[f[0]] = (f[1], f[2])

    # InterPro entries per gene: accession -> matched spans
    hits = defaultdict(lambda: defaultdict(list))
    with open(ipr_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13 or not f[11].startswith("IPR"):
                continue
            hits[f[0].rsplit(".t", 1)[0]][f[11]].append((int(f[6]), int(f[7])))
    n_genes_with = defaultdict(int)
    for g in hits:
        for acc in hits[g]:
            n_genes_with[acc] += 1

    def best_interpro(g, types):
        # most of the protein covered first; ties go to the rarer (more specific) entry
        cands = [(merged_length(spans), -n_genes_with[acc], acc)
                 for acc, spans in hits.get(g, {}).items()
                 if etype.get(acc, ("",))[0] in types]
        return [c[2] for c in sorted(cands, reverse=True)]

    genes = []
    with open(gff) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                genes.append(re.search(r"ID=([^;]+)", f[8]).group(1))

    per_symbol = defaultdict(int)
    for g in genes:
        level, sym, _ = egg.get(g, ("", "-", "-"))
        if level in ORTHOLOGY_LEVELS and sym != "-":
            per_symbol[sym.upper()] += 1

    counts = defaultdict(int)
    with open(out, "w") as o:
        o.write("gene\tproduct\tsource\tevidence\n")
        for g in genes:
            prod, src, ev = None, "hypothetical", "-"
            level, sym, desc = egg.get(g, ("", "-", "-"))
            if per_symbol.get(sym.upper(), 0) > MAX_PER_SYMBOL:
                level = "too many paralogs"
            if level in ORTHOLOGY_LEVELS and sym != "-" and sym.upper() in sp:
                prod = clean(sp[sym.upper()][1])
                src, ev = "swissprot", f"{sym}: {sp[sym.upper()][1]}"
            if not prod and level in ORTHOLOGY_LEVELS and desc != "-":
                d = clean(desc)
                if d and not DOMAIN_LABEL.search(d) and not desc.lower().startswith("belongs"):
                    prod, src, ev = d, "eggnog", desc
            if not prod:
                for acc in best_interpro(g, {"Family"}):
                    prod = family_name(etype[acc][1])
                    if prod:
                        src, ev = "family", f"{acc}: {etype[acc][1]}"
                        break
            if not prod:
                for acc in best_interpro(g, {"Domain", "Repeat"}):
                    prod = domain_name(etype[acc][1])
                    if prod:
                        src, ev = "domain", f"{acc}: {etype[acc][1]}"
                        break
            if not prod:
                prod, src, ev = HYPO, "hypothetical", "-"
            counts[src] += 1
            o.write(f"{g}\t{prod}\t{src}\t{ev}\n")
    print(f"{len(genes)} genes: " + ", ".join(
        f"{k} {counts[k]}" for k in ("swissprot", "eggnog", "family", "domain",
                                     "hypothetical")), file=sys.stderr)


if __name__ == "__main__":
    main()
