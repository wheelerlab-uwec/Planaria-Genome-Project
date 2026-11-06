# Comparison with published triclad genomes: download each assembly and its
# NCBI dataset report, run the same BUSCO setup as the pipeline on it (and on
# our cleaned primary assembly), and tabulate size, contiguity, BUSCO and
# repeat content (published figures from config; ours from RepeatMasker).
# All tools come from planaria-env.

import json

CMP = "results/comparison"
CMP_DB = "resources/comparison"
CMP_ASM = config["comparison"]["assemblies"]     # {accession: [species, family, repeat %, source]}


rule cmp_download:
    output:
        fna = f"{CMP_DB}/{{acc}}.fna",
        report = f"{CMP_DB}/{{acc}}.json"
    wildcard_constraints:
        acc = r"GCA_\d+\.\d+"
    params:
        api = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"
    resources:
        mem_mb = 8000
    log:
        "logs/cmp_download/{acc}.log"
    shell:
        """
        curl -sSfL -o {output.fna}.zip \
            "{params.api}/{wildcards.acc}/download?include_annotation_type=GENOME_FASTA" 2> {log}
        unzip -p {output.fna}.zip "ncbi_dataset/data/{wildcards.acc}/*.fna" > {output.fna} 2>> {log}
        rm {output.fna}.zip
        curl -sSfL -o {output.report} "{params.api}/{wildcards.acc}/dataset_report" 2>> {log}
        """


def cmp_busco_input(wildcards):
    if wildcards.acc == "Gdor":
        return f"{ANN}/genome/genome.fa"
    return f"{CMP_DB}/{wildcards.acc}.fna"


rule cmp_busco:
    # Same lineage and mode as busco_final; "Gdor" is our submitted assembly
    input:
        cmp_busco_input
    output:
        summary = f"{CMP}/busco/{{acc}}/short_summary.txt"
    params:
        lineage = config["busco"]["lineage"],
        outdir = f"{CMP}/busco"
    threads: config["threads"]["busco"]
    log:
        "logs/cmp_busco/{acc}.log"
    shell:
        """
        busco -i {input} -l {params.lineage} -m genome --offline -f \
              --cpu {threads} -o {wildcards.acc} --out_path {params.outdir} > {log} 2>&1
        cp {params.outdir}/{wildcards.acc}/short_summary.specific.*.txt {output.summary}
        """


rule cmp_table:
    input:
        reports = expand(f"{CMP_DB}/{{acc}}.json", acc=CMP_ASM),
        busco = expand(f"{CMP}/busco/{{acc}}/short_summary.txt", acc=list(CMP_ASM) + ["Gdor"]),
        fai = f"{ANN}/genome/genome.fa.fai",
        repeats = f"{ANN}/repeats/genome.fa.tbl"
    output:
        f"{CMP}/planarian_genomes.tsv"
    params:
        meta = lambda wildcards: json.dumps(CMP_ASM)
    log:
        "logs/cmp_table/run.log"
    shell:
        """
        python scripts/cmp_table.py '{params.meta}' {input.fai} {input.repeats} {CMP_DB} {CMP}/busco \\
            > {output} 2> {log}
        """
