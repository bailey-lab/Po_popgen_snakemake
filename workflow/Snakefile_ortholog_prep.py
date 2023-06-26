### This file is used to generate the necessary ortholog lists for use in comparison of ortholog diversity


#Setting up conda
#navigate to /work/users/k/e/kellybce/ovale1r/snakemake
#Load anaconda through installation or local version
#module load anaconda/2021.11
#set up conda environment
#conda create -c conda-forge -c bioconda -n snakemake snakemake
#conda activate snakemake
#conda config --set channel_priority strict
###then run the snakefile from the snakemake/directory (Snakefile is the in the workflow directory)
#snakemake --use-conda --cores 1 --jobs 1 -s workflow/Snakefile

#needed input files:
#	bed file of chromosome intervals for P. o. curtisi and P. o. wallikeri genomes to include (which excludes extrachromosomal contigs), formatted as "input/beds/{species}_chr.bed"
#	bed file of specific hypervariable gene family intervals for the P. o. curtisi and P. o. wallikeri genomes, formatted as "input/beds/{species}_genemask.bed
#	bed file for P. falciparum core genome to be included, formatted as "input/beds/{species}_core.bed"
#	List of pf-poc one-to-one orthologs in the pf genome, in the same order as the corresponding orthologs from poc
#	List of pf-poc one-to-one orthologs in the poc genome, in the same order as the corresponding orthologs from pf

configfile: "config/config.yaml"





###Determines the final files to be output by the snakemake pipeline
rule all:
#Previous rules should use wildcards for the project and species, but this final rule should employ the actual names and terms for the final file to be created
	input:
		pf_orthos_masked = config["input_beds"]+ config["pf_orthos_masked"],
		poc_orthos_masked = config["input_beds"]+ config["poc_orthos_masked"]

rule ortholog_masker:
# this pipeline requires, as input, a set of files containing the gene ID, chromosome, start position, and end position of 1-to-1 orthologs between the PocGH01 and Pf3D7 genome
# these two files must be in the same order (so that corresponding lines are the paired orthologs from each species)
### This rule uses a python script to remove any ortholog pairs in which one (or both) orthologs are partially/fully masked during the earlier vcf masking steps
	input:
		pf_orthos = config["input_beds"]+ config["pf_orthos"],
		poc_orthos = config["input_beds"]+ config["poc_orthos"],
		poc_chr = config["input_beds"]+"curtisigh01_chr.bed",
		poc_gene_mask = config["input_beds"]+"curtisigh01_genemask.bed",
		pf_core = config["input_beds"]+"pf3d7_core.bed"	
	output:
		pf_orthos_masked = config["input_beds"]+config["pf_orthos_masked"],
		poc_orthos_masked = config["input_beds"]+config["poc_orthos_masked"]
	script:
	###scripts are run from the location of the Snakefile, not the directory from which the pipeline is run
		"scripts/ortholog_masker.py"


