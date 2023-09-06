cd /work/users/k/e/kellybce/ovale1r/snakemake/processing
module load anaconda/2021.11
#conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
cp workflow/Snakefile_processing.py workflow/Snakefile_submission.py
snakemake --use-conda -s workflow/Snakefile_submission.py --profile slurm --rerun-incomplete
### --profile slurm command tells snakemake to check current directory for a folder caller "slurm" and then for a config.yaml file within
snakemake --use-conda -s workflow/Snakefile_ortholog_prep.py --profile slurm --rerun-incomplete
