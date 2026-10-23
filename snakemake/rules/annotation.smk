# Nuclear genome annotation, part 1: repeats and RNA-seq evidence.
# Runs on the cleaned primary assembly submitted to NCBI (WGS JCFHDI000000000),
# so gene models line up with the accessioned sequences. All tools come from
# planaria-env.

ANN = "results/annotation"
RNASEQ = config["annotation"]["rnaseq"]          # {run accession: sample label}


rule annot_genome:
    # Local copy of the submission FASTA (sequence IDs unchanged)
    input:
        config["annotation"]["genome"]
    output:
        fasta = f"{ANN}/genome/genome.fa",
        fai = f"{ANN}/genome/genome.fa.fai"
    log:
        "logs/annot_genome/run.log"
    shell:
        """
        seqkit seq -w 80 {input} > {output.fasta} 2> {log}
        samtools faidx {output.fasta} 2>> {log}
        seqkit stats -a -T {output.fasta} >> {log}
        """


# ---- Repeats ----------------------------------------------------------------

rule repeatmodeler:
    # De novo repeat library, including the LTR structural pipeline
    input:
        f"{ANN}/genome/genome.fa"
    output:
        families = f"{ANN}/repeats/repeatmodeler/Gdor-families.fa",
        stk = f"{ANN}/repeats/repeatmodeler/Gdor-families.stk"
    params:
        workdir = f"{ANN}/repeats/repeatmodeler"
    threads: 48
    resources:
        partition = "month",
        mem_mb = 200000
    log:
        "logs/repeatmodeler/run.log"
    shell:
        """
        genome=$(readlink -f {input})
        log=$(readlink -f {log})
        cd {params.workdir}
        BuildDatabase -name Gdor $genome > $log 2>&1
        RepeatModeler -database Gdor -threads {threads} -LTRStruct >> $log 2>&1
        """


rule repeatmasker:
    # Soft-mask with the species library; keep alignments for the landscape
    input:
        genome = f"{ANN}/genome/genome.fa",
        lib = f"{ANN}/repeats/repeatmodeler/Gdor-families.fa"
    output:
        masked = f"{ANN}/repeats/genome.fa.masked",
        out = f"{ANN}/repeats/genome.fa.out",
        tbl = f"{ANN}/repeats/genome.fa.tbl",
        align = f"{ANN}/repeats/genome.fa.align",
        gff = f"{ANN}/repeats/genome.fa.out.gff"
    params:
        outdir = f"{ANN}/repeats",
        parallel = 12     # RepeatMasker -pa; each rmblast job uses ~4 threads
    threads: 48
    resources:
        partition = "week",
        mem_mb = 128000
    log:
        "logs/repeatmasker/run.log"
    shell:
        """
        RepeatMasker -pa {params.parallel} -lib {input.lib} -xsmall -gff -a \
                     -dir {params.outdir} {input.genome} > {log} 2>&1
        """


rule repeat_landscape:
    # Kimura divergence summary and landscape plot (HTML). The bioconda
    # wrappers for these scripts cannot find RepeatMaskerConfig.pm, so
    # RepeatMasker's share directory is added to PERL5LIB.
    input:
        align = f"{ANN}/repeats/genome.fa.align",
        fai = f"{ANN}/genome/genome.fa.fai"
    output:
        divsum = f"{ANN}/repeats/genome.divsum",
        html = f"{ANN}/repeats/repeat_landscape.html"
    resources:
        mem_mb = 32000
    log:
        "logs/repeat_landscape/run.log"
    shell:
        """
        export PERL5LIB="$(dirname $(readlink -f $(which RepeatMasker)))/../share/RepeatMasker${{PERL5LIB:+:$PERL5LIB}}"
        calcDivergenceFromAlign.pl -s {output.divsum} {input.align} > {log} 2>&1
        size=$(awk '{{s+=$2}} END{{print s}}' {input.fai})
        createRepeatLandscape.pl -div {output.divsum} -g $size > {output.html} 2>> {log}
        """


# ---- RNA-seq (MA-C2 runs from PRJNA317859) -----------------------------------

rule rnaseq_download:
    output:
        r1 = "resources/annotation/rnaseq/{run}_1.fastq.gz",
        r2 = "resources/annotation/rnaseq/{run}_2.fastq.gz"
    params:
        outdir = "resources/annotation/rnaseq"
    threads: 8
    resources:
        mem_mb = 16000
    log:
        "logs/rnaseq_download/{run}.log"
    shell:
        """
        tmp=$(mktemp -d -p {params.outdir})
        prefetch {wildcards.run} --max-size 200G -O $tmp > {log} 2>&1
        fasterq-dump $tmp/{wildcards.run} --split-files -e {threads} -t $tmp -O $tmp >> {log} 2>&1
        pigz -p {threads} -c $tmp/{wildcards.run}_1.fastq > {output.r1}
        pigz -p {threads} -c $tmp/{wildcards.run}_2.fastq > {output.r2}
        rm -rf $tmp
        """


rule rnaseq_fastp:
    input:
        r1 = "resources/annotation/rnaseq/{run}_1.fastq.gz",
        r2 = "resources/annotation/rnaseq/{run}_2.fastq.gz"
    output:
        r1 = temp(f"{ANN}/rnaseq/trimmed/{{run}}_1.fastq.gz"),
        r2 = temp(f"{ANN}/rnaseq/trimmed/{{run}}_2.fastq.gz"),
        json = f"{ANN}/rnaseq/trimmed/{{run}}.fastp.json",
        html = f"{ANN}/rnaseq/trimmed/{{run}}.fastp.html"
    threads: 16
    resources:
        mem_mb = 16000
    log:
        "logs/rnaseq_fastp/{run}.log"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
              --detect_adapter_for_pe -w {threads} -j {output.json} -h {output.html} > {log} 2>&1
        """


rule star_index:
    # Unmasked genome; chrBinNbits sized for ~21,000 contigs over 1.47 Gb
    input:
        f"{ANN}/genome/genome.fa"
    output:
        directory(f"{ANN}/rnaseq/star_index")
    threads: 32
    resources:
        mem_mb = 128000
    log:
        "logs/star_index/run.log"
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} \
             --genomeFastaFiles {input} --genomeSAindexNbases 14 --genomeChrBinNbits 16 \
             --outFileNamePrefix {output}/ > {log} 2>&1
        """


rule star_align:
    # Two-pass mapping; intronMotif strand tags are what BRAKER/GeneMark expect.
    # BAM sorting opens threads x 50 temporary files, more than the default
    # 1,024 open-file limit on the compute nodes, so the soft limit is raised.
    input:
        index = f"{ANN}/rnaseq/star_index",
        r1 = f"{ANN}/rnaseq/trimmed/{{run}}_1.fastq.gz",
        r2 = f"{ANN}/rnaseq/trimmed/{{run}}_2.fastq.gz"
    output:
        bam = f"{ANN}/rnaseq/star/{{run}}.Aligned.sortedByCoord.out.bam",
        final = f"{ANN}/rnaseq/star/{{run}}.Log.final.out"
    params:
        prefix = f"{ANN}/rnaseq/star/{{run}}.",
        max_intron = config["annotation"]["max_intron"]
    threads: 32
    resources:
        mem_mb = 128000
    log:
        "logs/star_align/{run}.log"
    shell:
        """
        ulimit -n 16384
        rm -rf {params.prefix}_STARtmp
        STAR --runThreadN {threads} --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat \
             --twopassMode Basic --outSAMtype BAM SortedByCoordinate \
             --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical \
             --alignIntronMax {params.max_intron} --alignMatesGapMax {params.max_intron} \
             --limitBAMsortRAM 60000000000 --outFileNamePrefix {params.prefix} > {log} 2>&1
        """


rule rnaseq_merge:
    input:
        bams = expand(f"{ANN}/rnaseq/star/{{run}}.Aligned.sortedByCoord.out.bam", run=RNASEQ),
        finals = expand(f"{ANN}/rnaseq/star/{{run}}.Log.final.out", run=RNASEQ)
    output:
        bam = f"{ANN}/rnaseq/rnaseq_merged.bam",
        bai = f"{ANN}/rnaseq/rnaseq_merged.bam.bai",
        summary = f"{ANN}/rnaseq/mapping_summary.tsv"
    params:
        labels = [RNASEQ[r] for r in RNASEQ]
    threads: 16
    resources:
        mem_mb = 32000
    log:
        "logs/rnaseq_merge/run.log"
    shell:
        """
        samtools merge -@ {threads} -f {output.bam} {input.bams} 2> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        {{
        printf "run\\tsample\\tinput_pairs\\tunique_pct\\tmulti_pct\\tsplices\\n"
        labels=({params.labels})
        i=0
        for f in {input.finals}; do
            run=$(basename $f .Log.final.out)
            awk -F'|' -v r=$run -v s=${{labels[$i]}} '
                /Number of input reads/ {{n=$2}}
                /Uniquely mapped reads %/ {{u=$2}}
                /% of reads mapped to multiple loci/ {{m=$2}}
                /Number of splices: Total/ {{sp=$2}}
                END {{gsub(/[ \\t]/,"",n); gsub(/[ \\t]/,"",u); gsub(/[ \\t]/,"",m); gsub(/[ \\t]/,"",sp);
                      print r"\\t"s"\\t"n"\\t"u"\\t"m"\\t"sp}}' $f
            i=$((i+1))
        done
        }} > {output.summary}
        """
