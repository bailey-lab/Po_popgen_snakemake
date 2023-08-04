#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=140:00:00
#SBATCH --mem=30G

cd /work/users/k/e/kellybce/ovale1r/snakemake/processing/workflow/scripts
bash submit_snakefile.sh
