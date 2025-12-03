To run the Snakemake workflow:

1. Open a separate terminal instance (not the terminal of a VS Code sesssion)

2. Activate the module and environment:

    `module load python-libs`

    `conda activate planaria-env`

3. If you haven't recently, update the active environment:

    `mamba env update --file environment.yml`

4. In the Snakefile, update the `racon` command to have the full path to your personal installation, e.g.:

    `/data/users/wheelenj/GitHub/Planaria-Genome-Project/racon/build/bin/racon ...`

5. Dry-run to confirm the correct rules will be run:

    `snakemake -np`

6. After confirming, run the workflow using the custom executor that tailors Slurm resources:

    `snakemake --profile executor/`
    
    Or, to run all jobs in the background:

    `nohup snakemake --profile executor/ &`