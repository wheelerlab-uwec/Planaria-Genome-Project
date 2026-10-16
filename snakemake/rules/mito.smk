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


rule mito_qc:
    # Checks that the assembled length is real and the consensus has settled:
    # depth along the circle, reads spanning (nearly) the whole circle, the
    # cox1-to-cox1 distance inside very long reads, a third Medaka round, and
    # windows where reads disagree with the consensus (repeats, heteroplasmy).
    input:
        fasta = f"{MITO}/mito_genome.fasta",
        reads = f"{MITO}/reads/mito_reads.fastq.gz",
        sub = f"{MITO}/reads/mito_reads_subsampled.fastq.gz",
        bed = f"{MITO}/mitos/result.bed"
    output:
        bam = f"{MITO}/qc/reads_vs_mito.bam",
        depth = f"{MITO}/qc/depth_500bp.tsv",
        repeats = f"{MITO}/qc/cox1_spacing_in_long_reads.tsv",
        disagree = f"{MITO}/qc/low_agreement_1kb.tsv",
        summary = f"{MITO}/qc/mito_qc.txt"
    params:
        outdir = f"{MITO}/qc",
        model = config["mito"]["medaka_model"],
        long = config["mito"]["qc_long_read"],
        span = config["mito"]["qc_span"]
    threads: 16
    log:
        "logs/mito_qc/run.log"
    shell:
        """
        O={params.outdir}
        minimap2 -ax map-ont -t {threads} {input.fasta} {input.reads} 2> {log} \
            | samtools sort -@ 4 -o {output.bam} - 2>> {log}
        samtools index {output.bam}
        samtools depth -a {output.bam} \
            | awk -v OFS='\\t' '{{w=int(($2-1)/500); s[w]+=$3; n[w]++}} END{{for(w in s) print w*500+1, int(s[w]/n[w])}}' \
            | sort -k1,1n > {output.depth}

        # cox1 copies inside reads longer than the genome; spacing = genome length
        awk -v OFS='\\t' '$4=="cox1" {{print $1, $2, $3}}' {input.bed} > $O/cox1.bed
        seqkit subseq --bed $O/cox1.bed {input.fasta} 2>> {log} > $O/cox1.fasta
        seqkit seq -m {params.long} {input.reads} 2>> {log} | seqkit fq2fa > $O/long_reads.fasta
        blastn -query $O/cox1.fasta -subject $O/long_reads.fasta -max_target_seqs 100000 -max_hsps 100 \
               -outfmt '6 sseqid sstart send length pident' 2>> {log} \
            | awk '$4>=1000' | sort -k1,1 -k2,2n \
            | awk -v OFS='\\t' '$1==p {{d=$2-q; print $1, (d<0?-d:d)}} {{p=$1; q=$2}}' > {output.repeats}

        # third Medaka round: how much does the consensus still move?
        medaka_consensus -i {input.sub} -d {input.fasta} -o $O/polish3 -m {params.model} \
                         -t 8 -f >> {log} 2>&1

        # per-base agreement of reads with the consensus, summarised per 1 kb
        samtools mpileup -f {input.fasta} -Q 0 -B -d 0 {output.bam} 2>> {log} \
            | awk -v OFS='\\t' '$4>0 {{s=$5; m=gsub(/[.,]/,"",s); w=int(($2-1)/1000)*1000+1; n[w]++; if(m/$4<0.8) k[w]++}}
                                END{{for(w in n) print w, k[w]+0}}' | sort -k1,1n > {output.disagree}

        {{
        echo "== Length"; seqkit fx2tab -n -l {input.fasta} | cut -f2
        echo; echo "== Depth per 500 bp window (min / median / max)"
        cut -f2 {output.depth} | sort -n | awk '{{a[NR]=$1}} END{{print a[1], a[int(NR/2)+1], a[NR]}}'
        echo; echo "== Primary alignments spanning >= {params.span} bp of the reference"
        samtools view -F 0x904 -q 20 {output.bam} \
            | awk '{{n=0; c=$6; while(match(c,/[0-9]+[MD=X]/)){{n+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}} if(n>={params.span}) k++}} END{{print k+0}}'
        echo; echo "== cox1-to-cox1 spacing in reads >= {params.long} bp (count per 500 bp bin)"
        cut -f2 {output.repeats} | awk '{{print int($1/500)*500}}' | sort -n | uniq -c
        echo; echo "== Medaka round 3 vs final (length; substitutions and indels)"
        seqkit fx2tab -n -l $O/polish3/consensus.fasta | cut -f2
        minimap2 -cx asm5 --cs {input.fasta} $O/polish3/consensus.fasta 2>> {log} \
            | grep -o 'cs:Z:.*' | awk '{{s+=gsub(/\\*/,""); i+=gsub(/[+-]/,"")}} END{{print "substitutions", s+0, "indels", i+0}}'
        echo; echo "== 1 kb windows with >= 20 positions where < 80% of reads agree"
        awk '$2>=20' {output.disagree}
        }} > {output.summary} 2>> {log}
        rm -f $O/long_reads.fasta
        """


rule mito_barcode:
    # Species check: our cox1 against every Girardia COI record in GenBank
    input:
        fasta = f"{MITO}/mito_genome.fasta",
        bed = f"{MITO}/mitos/result.bed"
    output:
        db = "resources/mito/girardia_coi.fasta",
        hits = f"{MITO}/qc/cox1_vs_girardia_coi.tsv",
        summary = f"{MITO}/qc/cox1_barcode_by_species.tsv"
    params:
        query = config["mito"]["barcode_query"],
        outdir = f"{MITO}/qc"
    log:
        "logs/mito_barcode/run.log"
    shell:
        """
        E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
        curl -sSf -G "$E/esearch.fcgi" --data-urlencode db=nuccore --data-urlencode retmax=5000 \
             --data-urlencode "term={params.query}" 2> {log} \
            | grep -o '<Id>[0-9]*' | cut -c5- | paste -sd, > {params.outdir}/coi_ids.txt
        curl -sSf -X POST "$E/efetch.fcgi" --data "db=nuccore&rettype=fasta&retmode=text" \
             --data "id=$(cat {params.outdir}/coi_ids.txt)" > {output.db} 2>> {log}
        echo "COI records: $(grep -c '>' {output.db})" >> {log}
        awk -v OFS='\\t' '$4=="cox1" {{print $1, $2, $3}}' {input.bed} > {params.outdir}/cox1_barcode.bed
        seqkit subseq --bed {params.outdir}/cox1_barcode.bed {input.fasta} 2>> {log} > {params.outdir}/cox1_barcode.fasta
        # accession, identity, aligned length, species (from the FASTA title)
        blastn -query {params.outdir}/cox1_barcode.fasta -subject {output.db} -max_target_seqs 5000 \
               -outfmt '6 sseqid pident length' 2>> {log} \
            | awk '$3>=300' | sort -k2,2nr | awk '!s[$1]++' \
            | awk -v OFS='\\t' 'NR==FNR {{if(/^>/){{a=substr($1,2); sp[a]=$2" "$3}} next}} {{print $1, $2, $3, sp[$1]}}' {output.db} - \
            > {output.hits}
        # best identity and number of records per species
        awk -F'\\t' -v OFS='\\t' '!($4 in b) {{b[$4]=$2}} {{n[$4]++}} END{{for(s in b) print s, b[s], n[s]}}' {output.hits} \
            | sort -t$'\\t' -k2,2nr > {output.summary}
        """
