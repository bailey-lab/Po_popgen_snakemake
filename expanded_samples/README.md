#README

The associated files contain the needed code for alignment, variant-calling, quality-filtering, and analysis of whole-genome sequencing data
of Plasmodium ovale curtisi and Plasmodium ovale wallikeri. This pipeline is optimized for use of reference genomes PocGH01 and PowCR01, respectively.

This pipeline can be run from an anaconda environment using snakemake "Snakefiles" for pipeline organization and reproducibility. This particular analysis uses three
Snakefiles contained in the workflow directory: one for sequencing alignment to initial variant calling (Snakefile_processing.py), one for determination of orthologues genes between 
P. ovale curtisi and P. ovale wallikeri suitable for comparison (Snakefile_ortholog_prep.py), and one for final variant callset cleaning and analyses (Snakefile_analysis.py).

Each snakefile contains needed coded to load the appropriate anaconda environment and execute the code within. Each Snakefile also contains a list of needed input files for 
the pipeline to run completely. Code and input files may be altered as needed to adjust the exact functioning of this pipeline. The purpose of each snakemake rule is described
in the comments.

All Snakefiles are optimized to run on Anaconda version 2021.11 and Snakemake version 7.24.2. Individual software versions needed for packages employed in specific steps are described in the corresponding ".yaml"
files referenced in the "conda" section of each snakemake rule; these yaml files can be found within the worflow/envs directory. This pipeline was not executed on software versions
other than those described in the yaml files.

This pipeline can be run on the publicly-available data described in the Data Availability section of the corresponding manuscript. The "submit_snakefile.sh" script in the workflow/scripts
directory can be used to organize the three snakefles relative to each other, and the "sbatch_snakefile.sh" script can be used to submit all three Snakefiles to a SLURM resource manager
if applicable. The workflow'scripts directory also includes various scripts that are employed by the Snakemake pipeline itself (files ending in ".py") or separate scripts used for final
analysis and data visulaizations in R version 4.3.3 (files ending in ".R"). 

Installation takes up to 20 minutes on computers that have already installed anaconda. Running of the entire pipeline takes up to 7 days but may vary depending on access to computer clusters resources and use of
the SLURM system. Due to size of corresponding datasets, analysis may not be feasible on a conventional desktop computer without substantial memory enhancements.

After downloading corresponding "fastq.gz" files from the GenBank repository, input files for the pipeline shown in "Snakefile_processing.py" should be constructed using the appropriate reference genome and the
sample metadata, mainly sample names and codes. Input files should be deposited in the input directory, while output files from analysis will appear in the output directory. Sequencing data alignments and 
variant-calling can then be performed by submitting the processing Snakefile manually or to slurm. Of note, the "Rule All" within Snakefiles
can be altered or commented out to determine which final output files will be produced.

After production of a variant callset "vcf.gz" file for each species, the "Snakefile_ortholog_prep.py" pipeline can be submitted (using additional ortholog list input files) to determine a final list of orthologous 
genes with strong sequencing coverage across the included samples.

Finally, the "Snakefile_analysis.py" pipeline can be employed to perform masking and quality filtering of the "vcf.gz" files output from the "Snakefile_processing.py" pipeline. This pipeline will yield various final
variant callsets using prespecified quality filtering thresholds and genomic windows to be excluded (specified in config/config.yaml), as well as output files for analysis. Some of these files represent the final endpoint of a given analysis, while other 
others require additional processing and visualization in the R environment using corresponding ".R" scripts available in the workflow/scripts directory.
