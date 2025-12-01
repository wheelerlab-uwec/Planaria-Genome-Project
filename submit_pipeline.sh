# Run Snakemake with cluster submission
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cluster-config cluster.yaml \
    --cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.cpus} --output={cluster.output} --error={cluster.error}" \
    --jobs 20 \
    --latency-wait 60 \
    --restart-times 3 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --reason \
    --use-conda

echo "Pipeline completed!"
echo "End time: $(date)"

# ============================================================================
# Alternative: Dry Run (Test without execution)
# ============================================================================
# Uncomment to do a dry run first:
# snakemake --snakefile Snakefile --configfile config.yaml --dry-run --printshellcmds

# ============================================================================
# Alternative: Local Execution (No Cluster)
# ============================================================================
# snakemake --snakefile Snakefile --configfile config.yaml --cores 64 --printshellcmds