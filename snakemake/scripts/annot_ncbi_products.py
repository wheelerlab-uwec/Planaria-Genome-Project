#!/usr/bin/env python3
"""Turn eggNOG-mapper descriptions into NCBI-style protein product names.

eggNOG's Description column mixes real protein names ("Fatty-acid amide
hydrolase") with Pfam family/domain names ("LicD family") and GO terms
("zinc ion binding", "regulation of brood size"). NCBI wants protein names, so:
  - GO-like terms and empty descriptions become "hypothetical protein"
  - "Belongs to the X family" becomes "X family protein"
  - family/domain/region names become "X family protein" or
    "X domain-containing protein"
  - real protein names are kept, minus "homolog", trailing periods and
    unbalanced parentheses
Names are given per gene (eggNOG was run on the longest isoform) and apply to
every isoform of that gene.

Usage: annot_ncbi_products.py Gdor.emapper.annotations braker.gff3 out.tsv
Output columns: gene, product, source (eggnog|hypothetical), eggNOG description
"""
import re
import sys

HYPO = "hypothetical protein"

# GO-style phrases (functions, processes, locations) that are not protein names
GO_LIKE = re.compile(
    r"\bactivity\b|\bbinding\b(?!.*\b(protein|domain|factor|site)\b)|"
    r"^(positive |negative )?regulation of\b|^response to\b|^cellular\b|"
    r"\bprocess\b|\bis involved\b|^involved in\b|\bpathway\b|"
    r"\b(transport|signaling|signalling|assembly|organization|organisation|"
    r"development|differentiation|localization|metabolism|biosynthesis|"
    r"catabolism|morphogenesis|homeostasis|proliferation|secretion|"
    r"transcription|translation|repair|replication|import|export|"
    r"adhesion|apoptosis|division|segregation|cycle|formation|autophagy|"
    r"macroautophagy|chemotaxis|migration|polarity|growth|maturation|"
    r"silencing|splicing|folding|degradation|ubiquitination|"
    r"phosphorylation|methylation)$|^structural constituent\b|"
    r"\b(complex|membrane|nucleus|cytoplasm|organelle|cytoskeleton)$|"
    r"\bcatabolic\b|\bbiosynthetic\b|\bmetabolic\b",
    re.I)
# Words that mark a family or domain name rather than a protein name
DOMAIN_END = re.compile(r"\b(domain|region|motif|repeat|repeats|box|finger|"
                        r"terminal|fold|signature)$", re.I)
FAMILY_END = re.compile(r"\b(family|superfamily|subfamily)$", re.I)
DUF = re.compile(r"\b(DUF\d+|UPF\d+)\b")
# Lower-case GO-style text with none of these is a term, not a protein name
PROTEIN_WORD = re.compile(r"(ase|in|ins)\b|\b(protein|factor|receptor|channel|"
                          r"transporter|cytochrome|subunit|domain|family|"
                          r"superfamily|repeat|finger|lectin|histone|toxin|"
                          r"globin|\w+-like)\b|\d", re.I)


def clean(desc):
    d = re.sub(r"\s+", " ", desc).strip().rstrip(".").strip()
    d = d.replace("tumour", "tumor")
    if d in ("", "-"):
        return None
    # "Belongs to the X family. Y subfamily" -> most specific family
    if d.lower().startswith("belongs to the"):
        fam = d.split(". ")[0]
        fam = re.sub(r"^belongs to the\s+", "", fam, flags=re.I)
        fam = re.sub(r"\s*\(.*$", "", fam)
        if not FAMILY_END.search(fam):
            fam += " family"
        return fam + " protein"
    # Keep the first sentence only
    d = d.split(". ")[0]
    # Drop an unbalanced or trailing parenthetical
    if re.match(r"[^(]*\)", d):
        return None              # text truncated inside a parenthesis
    if d.count("(") != d.count(")"):
        d = d[:d.rfind("(")].rstrip(" ,")
    m = DUF.search(d)
    if m and re.search(r"unknown function|uncharacteri[sz]ed|predicted", d, re.I):
        return f"{m.group(1)} domain-containing protein"
    d = re.sub(r"\s*\([^()]*\)$", "", d).strip(" ,")
    if not d or GO_LIKE.search(d):
        return None
    if d[0].islower() and not PROTEIN_WORD.search(d):
        return None
    d = re.sub(r"\s+isoform(\s+\w+)?$", "", d)
    # US spelling and NCBI's preferred qualifier, as table2asn would change them
    for a, b in (("isation", "ization"), ("ising", "izing"), ("disulphide", "disulfide"),
                 ("sulphur", "sulfur"), ("tumour", "tumor"), ("signalling", "signaling"), ("haem", "hem"),
                 ("oestr", "estr")):
        d = d.replace(a, b)
    d = re.sub(r"^(possible|probable)\s+", "putative ", d, flags=re.I)
    d = re.sub(r"\bdomain containing\b", "domain-containing", d, flags=re.I)
    d = re.sub(r"\brepeats\b", "repeat", d)
    d = re.sub(r"\s+homolog(ue)?\b", "", d, flags=re.I).strip()
    d = re.sub(r"^putative\s+", "", d, flags=re.I)
    # sentences, organism references and SMART-style descriptions
    if re.search(r"^(to|and|or|from|in|of|with|binds?|mediates?|plays?|required|"
                 r"catalyz\w*|involved|acts|contains|repeats?\s+(in|found)|[-'\"])\b|"
                 r"\b(saccharomyces|cerevisiae|e\.\s?coli|drosophila|mouse|human|"
                 r"worm|fly|yeast|bacteriophage|genome)\b|"
                 r"\b(and|or|from|of|in|to|with)$|\(\)|\b(is|are|was)\b",
                 d, re.I):
        return None
    if re.search(r"\b(unknown function|uncharacteri[sz]ed|hypothetical)\b|"
                 r"\bzgc\b|open reading frame|\bKIAA\d|\bRIKEN\b|"
                 r"genomic stop codon|^protein \w+$|^(vacuolar|membrane) protein$",
                 d, re.I):
        return None
    d = re.sub(r"(?<=\w )Protein\b", "protein", d)
    if d.lower().endswith(" transmembrane"):
        d += " domain-containing protein"
    if FAMILY_END.search(d):
        d += " protein"
    elif DOMAIN_END.search(d) and not d.lower().endswith("containing protein"):
        last = d.split()[-1].lower()
        if last == "finger":
            d += " protein"
        elif last in ("repeat", "motif", "box", "domain"):
            d += "-containing protein"
        elif last == "region":
            d = d[:-len("region")] + "domain-containing protein"
        else:
            d += " domain-containing protein"
    if len(d) > 100 or len(d.split()) == 1 and not re.search(r"[a-z]{4,}", d):
        # bare symbols/acronyms ("TOG") and very long text are not usable names
        return None
    # Lower-case the first letter unless it starts an acronym or symbol
    if re.match(r"[A-Z][a-z]{2,}", d) and not re.match(r"[A-Z][a-z]+[A-Z]", d):
        d = d[0].lower() + d[1:]
    return d


def main():
    ann, gff, out = sys.argv[1:4]
    names = {}
    with open(ann) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            gene = f[0].rsplit(".t", 1)[0]
            names[gene] = (clean(f[7]), f[7])
    genes = []
    with open(gff) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                genes.append(re.search(r"ID=([^;]+)", f[8]).group(1))
    n_named = 0
    with open(out, "w") as o:
        o.write("gene\tproduct\tsource\teggnog_description\n")
        for g in genes:
            prod, desc = names.get(g, (None, "-"))
            src = "eggnog" if prod else "hypothetical"
            n_named += bool(prod)
            o.write(f"{g}\t{prod or HYPO}\t{src}\t{desc}\n")
    print(f"{len(genes)} genes; {n_named} named from eggNOG; "
          f"{len(genes) - n_named} hypothetical protein", file=sys.stderr)


if __name__ == "__main__":
    main()
