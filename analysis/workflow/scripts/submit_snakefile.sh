cd /work/users/k/e/kellybce/ovale1r/snakemake/analysis
module load anaconda/2021.11
#conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
snakemake --use-conda -s workflow/Snakefile_ortholog_prep.py --profile slurm
snakemake --use-conda -s workflow/Snakefile_analysis.py --profile slurm
