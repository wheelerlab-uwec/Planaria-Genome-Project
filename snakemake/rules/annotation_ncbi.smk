# Nuclear genome annotation, part 3: GenBank submission files for the
# annotation of WGS JCFHDI000000000 (same contigs as the submitted FASTA).
# Product names follow the International Protein Nomenclature Guidelines:
# Swiss-Prot names of eggNOG orthologs, then eggNOG protein names, then InterPro
# family and domain names, else "hypothetical protein" (see
# scripts/annot_ncbi_products.py). BRAKER models are given ADAYWY locus tags and
# NCBI ids, and table2asn builds the .sqn with validation (.val) and discrepancy
# (.dr) reports. Everything runs from planaria-env.

NCBI_ANN = f"{ANN}/ncbi"
IPR_RELEASE = config["annotation"]["interproscan_version"].split("-")[1]


rule annot_db_swissprot_names:
    # Reviewed metazoan UniProtKB entries: gene symbol -> protein name. The REST
    # API only serves the current release, so the rule stops if it is not the
    # release in config.
    output:
        f"{ADB}/uniprot/swissprot_metazoa_names.tsv"
    params:
        release = config["annotation"]["swissprot_release"],
        url = ("https://rest.uniprot.org/uniprotkb/stream?query=(reviewed:true)"
               "+AND+(taxonomy_id:33208)&fields=accession,gene_primary,"
               "protein_name,organism_id&format=tsv&compressed=true")
    resources:
        mem_mb = 4000
    log:
        "logs/annot_db_swissprot_names/run.log"
    shell:
        """
        curl -sSfL -D {output}.headers "{params.url}" 2> {log} | gunzip -c > {output}
        grep -i "x-uniprot-release:" {output}.headers | tee -a {log} | grep -q "{params.release}" \
            || {{ echo "UniProt release is not {params.release}" >> {log}; rm {output}; exit 1; }}
        rm {output}.headers
        """


rule annot_db_interpro_entries:
    # InterPro entry types (Family, Domain, Repeat, ...) for the release that
    # InterProScan uses
    output:
        f"{ADB}/interpro/entry.list"
    params:
        url = f"https://ftp.ebi.ac.uk/pub/databases/interpro/releases/{IPR_RELEASE}/entry.list"
    log:
        "logs/annot_db_interpro_entries/run.log"
    shell:
        """
        curl -sSfL -o {output} {params.url} 2> {log}
        """


rule annot_ncbi_products:
    input:
        script = "scripts/annot_ncbi_products.py",
        eggnog = f"{GENES}/eggnog/Gdor.emapper.annotations",
        sprot = f"{ADB}/uniprot/swissprot_metazoa_names.tsv",
        ipr = f"{GENES}/interproscan/Gdor.interproscan.tsv",
        entries = f"{ADB}/interpro/entry.list",
        gff = f"{GENES}/braker/braker.gff3"
    output:
        f"{NCBI_ANN}/products.tsv"
    log:
        "logs/annot_ncbi_products/run.log"
    shell:
        """
        python {input.script} {input.eggnog} {input.sprot} {input.ipr} \
               {input.entries} {input.gff} {output} 2> {log}
        """


rule annot_ncbi_gff:
    # Locus tags, NCBI ids and partial ends; drops proteins < 30 aa and
    # incomplete models that do not reach a contig end
    input:
        script = "scripts/annot_ncbi_gff.py",
        gff = f"{GENES}/braker/braker.gff3",
        products = f"{NCBI_ANN}/products.tsv",
        genome = config["annotation"]["genome"]
    output:
        lengths = f"{NCBI_ANN}/contig_lengths.tsv",
        gff = f"{NCBI_ANN}/Gdor.gff3",
        tags = f"{NCBI_ANN}/locus_tags.tsv"
    params:
        prefix = config["ncbi"]["locus_tag_prefix"],
        db = config["ncbi"]["protein_db"]
    log:
        "logs/annot_ncbi_gff/run.log"
    shell:
        """
        seqkit fx2tab -nl {input.genome} > {output.lengths} 2> {log}
        python {input.script} {input.gff} {input.products} {params.prefix} \
               {params.db} {output.lengths} {output.gff} {output.tags} 30 2>> {log}
        """


rule annot_table2asn:
    # -T looks up the taxonomy remotely (lineage and genetic code); without it
    # the discrepancy report treats the genome as bacterial
    input:
        genome = config["annotation"]["genome"],
        gff = f"{NCBI_ANN}/Gdor.gff3",
        template = config["ncbi"]["template"]
    output:
        fsa = f"{NCBI_ANN}/Gdor.fsa",
        sqn = f"{NCBI_ANN}/Gdor.sqn",
        val = f"{NCBI_ANN}/Gdor.val",
        dr = f"{NCBI_ANN}/Gdor.dr",
        stats = f"{NCBI_ANN}/Gdor.stats"
    params:
        organism = config["ncbi"]["organism"],
        isolate = config["ncbi"]["isolate"]
    threads: 4
    resources:
        mem_mb = 64000
    log:
        "logs/annot_table2asn/run.log"
    shell:
        """
        cp {input.genome} {output.fsa}
        table2asn -M n -J -c w -euk -T -t {input.template} \
            -j "[organism={params.organism}] [isolate={params.isolate}] [gcode=1]" \
            -i {output.fsa} -f {input.gff} -o {output.sqn} -Z > {log} 2>&1
        cat {output.stats} >> {log}
        """
