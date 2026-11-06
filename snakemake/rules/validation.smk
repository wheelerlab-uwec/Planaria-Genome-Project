# Assembly validation: map the filtered ONT reads (the Flye input) back to
# the cleaned primary assembly submitted to NCBI. Reports the fraction of
# reads and bases that map and the depth distribution. A 0.2% read subset is
# also aligned to the primary assembly and to the purged haplotigs to measure
# per-read divergence (minimap2 de tag, gap-compressed). The NM-based "error
# rate" from samtools stats is not used: it counts every base of long gaps,
# which inflates it far above read error. All tools come from planaria-env.

VAL = "results/validation"


rule val_map_reads:
    # Primary alignments only; the BAM is large and only feeds the stats below
    input:
        genome = f"{ANN}/genome/genome.fa",
        reads = "results/filtered/filtered_reads.fastq.gz"
    output:
        bam = temp(f"{VAL}/reads/reads.bam"),
        bai = temp(f"{VAL}/reads/reads.bam.bai")
    params:
        tmp = f"{VAL}/reads/sort_tmp"
    threads: 48
    resources:
        partition = "week",
        mem_mb = 160000
    log:
        "logs/val_map_reads/run.log"
    shell:
        """
        minimap2 -ax map-ont --secondary=no -t {threads} {input.genome} {input.reads} 2> {log} | \
            samtools sort -@ 8 -m 4G -T {params.tmp} -o {output.bam} - 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        """


rule val_samtools_stats:
    # Reads/bases mapped and the mismatch ("error rate") of reads vs assembly
    input:
        bam = f"{VAL}/reads/reads.bam"
    output:
        stats = f"{VAL}/reads/samtools_stats.txt"
    threads: 8
    log:
        "logs/val_samtools_stats/run.log"
    shell:
        """
        samtools stats -@ {threads} {input.bam} > {output.stats} 2> {log}
        """


rule val_depth:
    # Per-contig coverage, plus a genome-wide histogram of per-base depth
    # (depths above 1000 are pooled into the 1000 bin)
    input:
        bam = f"{VAL}/reads/reads.bam",
        bai = f"{VAL}/reads/reads.bam.bai"
    output:
        coverage = f"{VAL}/reads/contig_coverage.tsv",
        hist = f"{VAL}/reads/depth_histogram.tsv"
    threads: 8
    log:
        "logs/val_depth/run.log"
    shell:
        """
        samtools coverage {input.bam} > {output.coverage} 2> {log}
        samtools depth -a -@ {threads} {input.bam} 2>> {log} | \
            awk '{{d = $3 > 1000 ? 1000 : $3; h[d]++}} END {{for (d in h) print d "\\t" h[d]}}' | \
            sort -n > {output.hist}
        """


rule val_read_subset:
    input:
        "results/filtered/filtered_reads.fastq.gz"
    output:
        f"{VAL}/divergence/subset.fastq.gz"
    params:
        fraction = 0.002,
        seed = 11
    log:
        "logs/val_read_subset/run.log"
    shell:
        """
        seqkit sample -p {params.fraction} -s {params.seed} {input} -o {output} 2> {log}
        """


rule val_subset_align:
    # Base-level PAF (-c) so each alignment carries de (gap-compressed divergence)
    input:
        reads = f"{VAL}/divergence/subset.fastq.gz",
        primary = f"{ANN}/genome/genome.fa",
        haplotigs = "results/final/haplotigs.fasta"
    output:
        primary = f"{VAL}/divergence/subset_vs_primary.paf.gz",
        haplotigs = f"{VAL}/divergence/subset_vs_haplotigs.paf.gz"
    threads: 32
    log:
        "logs/val_subset_align/run.log"
    shell:
        """
        minimap2 -x map-ont -c --secondary=no -t {threads} {input.primary} {input.reads} 2> {log} | \
            gzip > {output.primary}
        minimap2 -x map-ont -c --secondary=no -t {threads} {input.haplotigs} {input.reads} 2>> {log} | \
            gzip > {output.haplotigs}
        """


rule val_mapping_summary:
    input:
        stats = f"{VAL}/reads/samtools_stats.txt",
        hist = f"{VAL}/reads/depth_histogram.tsv",
        primary = f"{VAL}/divergence/subset_vs_primary.paf.gz",
        haplotigs = f"{VAL}/divergence/subset_vs_haplotigs.paf.gz"
    output:
        f"{VAL}/read_mapping_summary.tsv"
    log:
        "logs/val_mapping_summary/run.log"
    shell:
        """
        python scripts/val_mapping_summary.py {input.stats} {input.hist} \
            {input.primary} {input.haplotigs} > {output} 2> {log}
        """
