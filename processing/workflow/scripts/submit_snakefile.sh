cd /work/users/k/e/kellybce/ovale1r/snakemake/processing
module load anaconda/2021.11
#conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
snakemake --use-conda --cores 16 -s workflow/Snakefile_processing.py --profile slurm
### --profile slurm command tells snakemake to check current directory for a folder caller "slurm" and then for a config.yaml file within
