# Nuclear genome annotation, part 2: gene prediction, QC and function.
# BRAKER3 (RNA-seq + proteins) runs from braker-env (envs/braker.yaml plus
# GeneMark-ETP cloned into that env); everything else runs from planaria-env.
# InterProScan is EBI's release, unpacked under resources/ and run with the
# planaria-env Java.

GENES = f"{ANN}/genes"
ADB = "resources/annotation"      # downloaded databases


# ---- Databases ---------------------------------------------------------------

rule annot_db_orthodb:
    # OrthoDB Metazoa protein partition (BRAKER's recommended protein evidence)
    output:
        f"{ADB}/orthodb/Metazoa.fa"
    params:
        url = config["annotation"]["orthodb_url"]
    resources:
        mem_mb = 8000
    log:
        "logs/annot_db_orthodb/run.log"
    shell:
        """
        curl -sSfL {params.url} 2> {log} | gunzip -c > {output}
        grep -c ">" {output} >> {log}
        """


rule annot_db_omark:
    output:
        f"{ADB}/omark/LUCA.h5"
    params:
        url = config["annotation"]["omark_db_url"]
    resources:
        mem_mb = 8000
    log:
        "logs/annot_db_omark/run.log"
    shell:
        "curl -sSfL -o {output} {params.url} 2> {log}"


rule annot_db_eggnog:
    # Direct download: download_eggnog_data.py (2.1.15) still points at
    # eggnogdb.embl.de, which no longer resolves
    output:
        db = f"{ADB}/eggnog/eggnog.db",
        taxa = f"{ADB}/eggnog/eggnog.taxa.db",
        dmnd = f"{ADB}/eggnog/eggnog_proteins.dmnd"
    params:
        outdir = f"{ADB}/eggnog",
        url = config["annotation"]["eggnog_url"]
    threads: 4
    resources:
        mem_mb = 8000
    log:
        "logs/annot_db_eggnog/run.log"
    shell:
        """
        cd {params.outdir}
        rm -f eggnog.db.gz eggnog.taxa.tar.gz eggnog_proteins.dmnd.gz
        log=$(readlink -f ../../../{log})
        for f in eggnog.db.gz eggnog.taxa.tar.gz eggnog_proteins.dmnd.gz; do
            curl -sSfL --retry 5 -o $f {params.url}/$f 2>> $log
        done
        pigz -d -p {threads} eggnog.db.gz eggnog_proteins.dmnd.gz
        tar -xzf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz
        ls -l >> $log
        """


rule annot_db_interproscan:
    # EBI release tarball; setup.py builds the HMM indexes once
    output:
        f"{ADB}/interproscan/interproscan.sh"
    params:
        version = config["annotation"]["interproscan_version"],
        outdir = f"{ADB}/interproscan"
    threads: 8
    resources:
        mem_mb = 32000
    log:
        "logs/annot_db_interproscan/run.log"
    shell:
        """
        v={params.version}
        tar=interproscan-$v-64-bit.tar.gz
        log=$(readlink -f {log})
        tmp=$(mktemp -d -p resources/annotation)
        curl -sSfL -o $tmp/$tar https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/$v/$tar 2> $log
        curl -sSfL -o $tmp/$tar.md5 https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/$v/$tar.md5 2>> $log
        (cd $tmp && md5sum -c $tar.md5) >> $log 2>&1
        tar -xzf $tmp/$tar -C $tmp
        rm -rf {params.outdir} && mv $tmp/interproscan-$v {params.outdir} && rm -rf $tmp
        cd {params.outdir} && python3 setup.py -f interproscan.properties >> $log 2>&1
        """


# ---- Gene prediction -----------------------------------------------------------

rule braker3:
    input:
        genome = f"{ANN}/repeats/genome.fa.masked",
        bam = f"{ANN}/rnaseq/rnaseq_merged.bam",
        prot = f"{ADB}/orthodb/Metazoa.fa"
    output:
        gtf = f"{GENES}/braker/braker.gtf",
        gff3 = f"{GENES}/braker/braker.gff3",
        aa = f"{GENES}/braker/braker.aa",
        cds = f"{GENES}/braker/braker.codingseq"
    params:
        workdir = f"{GENES}/braker",
        species = config["annotation"]["braker_species"]
    threads: 48
    resources:
        partition = "month",
        mem_mb = 200000
    log:
        "logs/braker3/run.log"
    shell:
        # A private copy of the AUGUSTUS config keeps training files out of the env
        """
        rm -rf {params.workdir} && mkdir -p {params.workdir}
        genome=$(readlink -f {input.genome}); bam=$(readlink -f {input.bam})
        prot=$(readlink -f {input.prot}); wd=$(readlink -f {params.workdir})
        conda run -n braker-env bash -c "
            set -e
            cp -r \\$CONDA_PREFIX/config $wd/augustus_config
            braker.pl --genome=$genome --bam=$bam --prot_seq=$prot \
                      --species={params.species} --threads={threads} --gff3 \
                      --GENEMARK_PATH=\\$GENEMARK_PATH \
                      --PROTHINT_PATH=\\$PROTHINT_PATH \
                      --AUGUSTUS_CONFIG_PATH=$wd/augustus_config --workingdir=$wd
        " > {log} 2>&1
        """


rule annot_longest_isoform:
    # One protein per gene (longest isoform) for BUSCO, OMArk and eggNOG;
    # stop-codon asterisks removed for InterProScan
    input:
        aa = f"{GENES}/braker/braker.aa"
    output:
        f"{GENES}/proteins.longest.faa"
    log:
        "logs/annot_longest_isoform/run.log"
    shell:
        "python scripts/annot_longest_isoform.py {input.aa} {output} 2> {log}"


rule annot_gene_stats:
    input:
        gtf = f"{GENES}/braker/braker.gtf"
    output:
        f"{GENES}/gene_stats.tsv"
    log:
        "logs/annot_gene_stats/run.log"
    shell:
        "python scripts/annot_gene_stats.py {input.gtf} > {output} 2> {log}"


# ---- QC ----------------------------------------------------------------------

rule annot_busco:
    # Protein-mode BUSCO avoids the genome-mode miniprot issue in planarians
    input:
        f"{GENES}/proteins.longest.faa"
    output:
        f"{GENES}/busco/short_summary.txt"
    params:
        outdir = f"{GENES}/busco",
        lineage = "busco_downloads/lineages/metazoa_odb10"
    threads: 16
    resources:
        mem_mb = 32000
    log:
        "logs/annot_busco/run.log"
    shell:
        """
        busco -i {input} -m proteins -l {params.lineage} --offline -c {threads} \
              --out_path {params.outdir} -o run -f > {log} 2>&1
        cp {params.outdir}/run/short_summary.*.txt {output}
        """


rule annot_omark:
    # Proteome completeness and consistency (contamination, fragments)
    input:
        faa = f"{GENES}/proteins.longest.faa",
        db = f"{ADB}/omark/LUCA.h5"
    output:
        omamer = f"{GENES}/omark/proteins.omamer",
        summary = f"{GENES}/omark/proteins.sum"
    params:
        outdir = f"{GENES}/omark"
    threads: 16
    resources:
        mem_mb = 64000
    log:
        "logs/annot_omark/run.log"
    shell:
        """
        omamer search --db {input.db} --query {input.faa} --out {output.omamer} \
                      --nthreads {threads} > {log} 2>&1
        omark -f {output.omamer} -d {input.db} -o {params.outdir} >> {log} 2>&1
        """


# ---- Function ----------------------------------------------------------------

rule annot_eggnog:
    input:
        faa = f"{GENES}/proteins.longest.faa",
        db = f"{ADB}/eggnog/eggnog.db",
        dmnd = f"{ADB}/eggnog/eggnog_proteins.dmnd"
    output:
        f"{GENES}/eggnog/Gdor.emapper.annotations"
    params:
        outdir = f"{GENES}/eggnog",
        dbdir = f"{ADB}/eggnog"
    threads: 48
    resources:
        mem_mb = 64000
    log:
        "logs/annot_eggnog/run.log"
    shell:
        """
        emapper.py -i {input.faa} --itype proteins -m diamond --cpu {threads} \
                   --data_dir {params.dbdir} --output_dir {params.outdir} --output Gdor \
                   --override > {log} 2>&1
        """


rule annot_interproscan:
    input:
        faa = f"{GENES}/proteins.longest.faa",
        ipr = f"{ADB}/interproscan/interproscan.sh"
    output:
        tsv = f"{GENES}/interproscan/Gdor.interproscan.tsv",
        gff3 = f"{GENES}/interproscan/Gdor.interproscan.gff3"
    params:
        outdir = f"{GENES}/interproscan"
    threads: 48
    resources:
        mem_mb = 128000
    log:
        "logs/annot_interproscan/run.log"
    shell:
        """
        mkdir -p {params.outdir}
        {input.ipr} -i {input.faa} -f tsv,gff3 -goterms -iprlookup -cpu {threads} \
                    -b {params.outdir}/Gdor.interproscan -T {params.outdir}/tmp > {log} 2>&1
        rm -rf {params.outdir}/tmp
        """
