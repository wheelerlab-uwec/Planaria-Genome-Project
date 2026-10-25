# Nuclear genome annotation, part 3: GenBank submission files for the
# annotation of WGS JCFHDI000000000 (same contigs as the submitted FASTA).
# Products come from eggNOG; BRAKER models are given ADAYWY locus tags and NCBI
# ids, and table2asn builds the .sqn with validation (.val) and discrepancy
# (.dr) reports. Everything runs from planaria-env.

NCBI_ANN = f"{ANN}/ncbi"


rule annot_ncbi_products:
    input:
        eggnog = f"{GENES}/eggnog/Gdor.emapper.annotations",
        gff = f"{GENES}/braker/braker.gff3"
    output:
        f"{NCBI_ANN}/products.tsv"
    log:
        "logs/annot_ncbi_products/run.log"
    shell:
        """
        python scripts/annot_ncbi_products.py {input.eggnog} {input.gff} {output} 2> {log}
        """


rule annot_ncbi_gff:
    # Locus tags, NCBI ids and partial ends; drops proteins < 30 aa and
    # incomplete models that do not reach a contig end
    input:
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
        python scripts/annot_ncbi_gff.py {input.gff} {input.products} {params.prefix} \
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
