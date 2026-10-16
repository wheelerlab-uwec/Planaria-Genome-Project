# Mitochondrial genome assembly.
# The nuclear Flye run tangles the mitogenome into a ~47 kb high-depth contig,
# and purge_dups drops it as HIGHCOV. Here it is reassembled on its own:
# find the mito contig in the Flye graph, pull the reads that are almost
# entirely mitochondrial, assemble them with Flye, polish with Medaka, rotate
# to cox1 and annotate with MITOS2.

MITO = "results/mito"
MITO_REFS = config["mito"]["references"]


rule mito_references:
    output:
        fasta = "resources/mito/references.fasta",
        gb = "resources/mito/references.gb"
    params:
        ids = ",".join(MITO_REFS)
    log:
        "logs/mito_references/run.log"
    shell:
        """
        E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
        curl -sSf "$E?db=nuccore&id={params.ids}&rettype=fasta&retmode=text" > {output.fasta} 2> {log}
        curl -sSf "$E?db=nuccore&id={params.ids}&rettype=gb&retmode=text" > {output.gb} 2>> {log}
        grep ">" {output.fasta} >> {log}
        """


rule mito_find_seed:
    # Search high-depth nuclear Flye contigs for the mitogenome; keep the
    # contig with the largest total BLAST bitscore against the references.
    input:
        assembly = "results/assembly/flye/assembly.fasta",
        info = "results/assembly/flye/assembly_info.txt",
        refs = "resources/mito/references.fasta"
    output:
        candidates = f"{MITO}/seed/highdepth_contigs.fasta",
        hits = f"{MITO}/seed/blast_refs_vs_highdepth.tsv",
        seed = f"{MITO}/seed/seed.fasta"
    params:
        min_depth = config["mito"]["seed_min_depth"]
    threads: 8
    log:
        "logs/mito_find_seed/run.log"
    shell:
        """
        awk -F'\\t' 'NR>1 && $3>={params.min_depth} {{print $1}}' {input.info} > {output.candidates}.ids
        seqkit grep -f {output.candidates}.ids {input.assembly} > {output.candidates} 2> {log}
        blastn -task dc-megablast -query {input.refs} -subject {output.candidates} \
               -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' \
               > {output.hits} 2>> {log}
        seed=$(awk '{{s[$2]+=$10}} END{{for(k in s) print s[k], k}}' {output.hits} | sort -k1,1nr | head -1 | cut -d' ' -f2)
        echo "seed contig: $seed" >> {log}
        seqkit grep -p "$seed" {output.candidates} > {output.seed} 2>> {log}
        rm {output.candidates}.ids
        """


rule mito_select_reads:
    # Keep reads that are >= min_aligned_frac mitochondrial, which excludes
    # reads from nuclear mito-like segments (NUMTs) with nuclear flanks.
    input:
        seed = f"{MITO}/seed/seed.fasta",
        reads = "results/filtered/filtered_reads.fastq.gz"
    output:
        paf = temp(f"{MITO}/reads/reads_vs_seed.paf"),
        ids = f"{MITO}/reads/mito_read_ids.txt",
        reads = f"{MITO}/reads/mito_reads.fastq.gz"
    params:
        min_frac = config["mito"]["min_aligned_frac"],
        min_len = config["mito"]["min_read_length"]
    threads: 32
    log:
        "logs/mito_select_reads/run.log"
    shell:
        """
        minimap2 -x map-ont -t {threads} {input.seed} {input.reads} > {output.paf} 2> {log}
        # aligned query bases per read (capped at read length) / read length
        awk -v OFS='\\t' '{{a[$1]+=$4-$3; l[$1]=$2}} END{{for(r in a){{x=(a[r]>l[r])?l[r]:a[r]; print r, l[r], x/l[r]}}}}' {output.paf} \
            | awk -v f={params.min_frac} -v m={params.min_len} '$2>=m && $3>=f {{print $1}}' > {output.ids}
        echo "mito reads: $(wc -l < {output.ids})" >> {log}
        seqkit grep -j {threads} -f {output.ids} {input.reads} -o {output.reads} 2>> {log}
        seqkit stats -a {output.reads} >> {log}
        """


rule mito_subsample:
    # Longest, highest-quality reads up to ~target_bases
    input:
        f"{MITO}/reads/mito_reads.fastq.gz"
    output:
        f"{MITO}/reads/mito_reads_subsampled.fastq.gz"
    params:
        target = config["mito"]["target_bases"]
    log:
        "logs/mito_subsample/run.log"
    shell:
        """
        filtlong --target_bases {params.target} {input} 2> {log} | gzip -c > {output}
        """


rule mito_flye:
    input:
        f"{MITO}/reads/mito_reads_subsampled.fastq.gz"
    output:
        assembly = f"{MITO}/flye/assembly.fasta",
        info = f"{MITO}/flye/assembly_info.txt"
    params:
        outdir = f"{MITO}/flye",
        genome_size = config["mito"]["genome_size"]
    threads: 16
    log:
        "logs/mito_flye/run.log"
    shell:
        """
        flye --nano-hq {input} --out-dir {params.outdir} --threads {threads} \
             --genome-size {params.genome_size} --iterations 3 2> {log}
        """


rule mito_pick_contig:
    # Circular contig of plausible mitogenome length with the highest depth
    input:
        assembly = f"{MITO}/flye/assembly.fasta",
        info = f"{MITO}/flye/assembly_info.txt"
    output:
        f"{MITO}/draft/mito_draft.fasta"
    params:
        min_len = config["mito"]["min_length"],
        max_len = config["mito"]["max_length"]
    log:
        "logs/mito_pick_contig/run.log"
    shell:
        """
        cat {input.info} > {log}
        id=$(awk -F'\\t' -v a={params.min_len} -v b={params.max_len} \
              'NR>1 && $4=="Y" && $2>=a && $2<=b' {input.info} | sort -t$'\\t' -k3,3nr | head -1 | cut -f1)
        if [ -z "$id" ]; then
            echo "No circular contig of {params.min_len}-{params.max_len} bp; check {input.info}" >> {log}
            exit 1
        fi
        echo "picked: $id" >> {log}
        seqkit grep -p "$id" {input.assembly} > {output} 2>> {log}
        """


rule mito_polish_draft:
    input:
        draft = f"{MITO}/draft/mito_draft.fasta",
        reads = f"{MITO}/reads/mito_reads_subsampled.fastq.gz"
    output:
        f"{MITO}/polish1/consensus.fasta"
    params:
        outdir = f"{MITO}/polish1",
        model = config["mito"]["medaka_model"]
    threads: 8
    log:
        "logs/mito_polish_draft/run.log"
    shell:
        """
        medaka_consensus -i {input.reads} -d {input.draft} -o {params.outdir} \
                         -m {params.model} -t {threads} -f 2> {log}
        """


rule mito_rotate:
    # Start the circle at cox1 on the forward strand (as in the references)
    input:
        fasta = f"{MITO}/polish1/consensus.fasta",
        gb = "resources/mito/references.gb"
    output:
        f"{MITO}/rotated/mito_rotated.fasta"
    params:
        workdir = f"{MITO}/rotated"
    log:
        "logs/mito_rotate/run.log"
    shell:
        """
        python scripts/mito_rotate.py --fasta {input.fasta} --genbank {input.gb} \
               --workdir {params.workdir} --out {output} 2> {log}
        """


rule mito_polish_final:
    # Second Medaka round on the rotated sequence fixes the old circle seam
    input:
        draft = f"{MITO}/rotated/mito_rotated.fasta",
        reads = f"{MITO}/reads/mito_reads_subsampled.fastq.gz"
    output:
        consensus = f"{MITO}/polish2/consensus.fasta",
        final = f"{MITO}/mito_genome.fasta"
    params:
        outdir = f"{MITO}/polish2",
        model = config["mito"]["medaka_model"],
        name = config["mito"]["seq_name"]
    threads: 8
    log:
        "logs/mito_polish_final/run.log"
    shell:
        """
        medaka_consensus -i {input.reads} -d {input.draft} -o {params.outdir} \
                         -m {params.model} -t {threads} -f 2> {log}
        seqkit replace -p '^(\\S+).*' -r '{params.name} [topology=circular] [location=mitochondrion]' \
            {output.consensus} > {output.final} 2>> {log}
        """


rule mito_mitos_refdata:
    output:
        directory(f"resources/mito/mitos/{config['mito']['mitos_refseq']}")
    params:
        ref = config["mito"]["mitos_refseq"]
    log:
        "logs/mito_mitos_refdata/run.log"
    shell:
        """
        mkdir -p resources/mito/mitos
        curl -sSfL https://zenodo.org/record/4284483/files/{params.ref}.tar.bz2 2> {log} \
            | tar -xj -C resources/mito/mitos
        """


rule mito_annotate:
    input:
        fasta = f"{MITO}/mito_genome.fasta",
        ref = f"resources/mito/mitos/{config['mito']['mitos_refseq']}"
    output:
        gff = f"{MITO}/mitos/result.gff",
        bed = f"{MITO}/mitos/result.bed"
    params:
        outdir = f"{MITO}/mitos",
        code = config["mito"]["genetic_code"],
        ref = config["mito"]["mitos_refseq"]
    log:
        "logs/mito_annotate/run.log"
    shell:
        # MITOS2 lives in its own env (envs/mitos.yaml); the cluster's conda is
        # too old for Snakemake's per-rule conda: environments
        """
        mkdir -p {params.outdir}
        conda run -n mitos-env runmitos.py --input {input.fasta} --code {params.code} --outdir {params.outdir} \
                    --refseqver {params.ref} --refdir resources/mito/mitos 2> {log}
        """


rule mito_summary:
    input:
        fasta = f"{MITO}/mito_genome.fasta",
        reads = f"{MITO}/reads/mito_reads.fastq.gz",
        refs = "resources/mito/references.fasta",
        bed = f"{MITO}/mitos/result.bed",
        info = f"{MITO}/flye/assembly_info.txt"
    output:
        f"{MITO}/mito_summary.txt"
    threads: 8
    log:
        "logs/mito_summary/run.log"
    shell:
        """
        {{
        echo "== Flye contigs"; cat {input.info}
        echo; echo "== Final sequence"; seqkit fx2tab -n -l -g {input.fasta}
        echo; echo "== Mean read depth (all mito reads)"
        minimap2 -ax map-ont -t {threads} {input.fasta} {input.reads} 2>> {log} \
            | samtools sort -@ 4 -o - 2>> {log} | samtools depth -a - | awk '{{s+=$3}} END{{printf "%.0f\\n", s/NR}}'
        echo; echo "== Identity to references (blastn dc-megablast, hits >= 500 bp)"
        blastn -task dc-megablast -query {input.refs} -subject {input.fasta} -evalue 1e-10 \
               -outfmt '6 qseqid pident length qstart qend sstart send' | awk '$3>=500'
        echo; echo "== MITOS2 features"
        cut -f4 {input.bed} | sort | uniq -c | sort -k1,1nr | wc -l
        cut -f4 {input.bed} | tr '\\n' ' '; echo
        }} > {output} 2>> {log}
        """
