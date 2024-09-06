### This file is used to generate the necessary ortholog lists for use in comparison of ortholog diversity 
#by checking for masking and insufficient coverage of any orthologues in each set
# initial lists of one-to-one orthologue sets drawn from OrthoMCL


#Setting up conda
#navigate tocorrect directory
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
#	List of pf-poc-pow one-to-one-to-one orthologs in each genome, all in the same order
#	List of poc-pow one-to-one orthologs in each genome, all in the same order

configfile: "config/config.yaml"

#create list of individual poc samples
pocsamples = []
for f in open(config["input_lists"]+"dbimport_curtisigh01_speciescall_samplemap.txt").readlines():
        pocsamples.append(f.split("\t")[0])

#create list of individual pow samples
powsamples = []
for f in open(config["input_lists"]+"dbimport_wallikericr01_speciescall_samplemap.txt").readlines():
        powsamples.append(f.split("\t")[0])

#create list of individual pf samples
pfsamples = []
for f in open(config["input_lists"]+"dbimport_pf3d7_samplemap.txt").readlines():
		pfsamples.append(f.split("\t")[0].rstrip("\n"))

###Determines the final files to be output by the snakemake pipeline
rule all:
#Previous rules should use wildcards for the project and species, but this final rule should employ the actual names and terms for the final file to be created
	input:
		pf_triple_orthos_masked = config["output"]+ config["pf_triple_orthos_masked"],
		poc_triple_orthos_masked = config["output"]+ config["poc_triple_orthos_masked"],
		pow_triple_orthos_masked = config["output"]+ config["pow_triple_orthos_masked"],
		poc_ovale_orthos_masked = config["output"]+ config["poc_ovale_orthos_masked"],
		pow_ovale_orthos_masked = config["output"]+ config["pow_ovale_orthos_masked"],
		triple_beds = expand(config["output"]+"ortholog_beds/{species}_pf-poc-pow_orthologs_masked.bed", species = ["curtisigh01","wallikericr01","pf3d7"]),
		ovale_beds = expand(config["output"]+"ortholog_beds/{species}_poc-pow_orthologs_masked.bed", species = ["curtisigh01","wallikericr01"]),
		

rule ovale_coverage_bed_files:
### generate bed files for each sample by species call that show the intervals of the genome with over 5X covergae
	input:
		config["output"]+"ovale_alignments_cleaned/{samplename}_{species}.bam"
	output:
		config["output"]+"statistics/coverage_beds/{samplename}_{species}_cov{depth}.bed"
	conda:
		"envs/bedtools.yaml"
	resources:
		mem_mb = 20000
	shell:
		"bedtools genomecov -ibam {input} -bg | awk '$4>={wildcards.depth}' > {output}"
		
rule curtisigh01_coverage_overall:
### generate a bed file reflecting intervals where {prop} of samples had at least {cov_filter} coverage (see config.yaml)
#use len(open(dbimportspeciescall).readlines())*{prop} to find threshold for including a window
	input:
		beds = expand(config["output"]+"statistics/coverage_beds/{samplename}_curtisigh01_cov{depth}.bed", samplename = pocsamples, allow_missing = True),
	output:
		config["output"]+"statistics/coverage_beds/curtisigh01_{prop}cov{depth}_unmerged.bed"
	params:
		integer = lambda wildcards, output: (float(wildcards.prop)*len(pocsamples)) if (float(wildcards.prop) * len(pocsamples)).is_integer() else int(float(wildcards.prop)*len(pocsamples))+1
	conda:
		"envs/bedtools.yaml"
	resources:
		mem_mb = 20000
	shell:
		"""echo {params.integer}
		bedtools multiinter -i {input.beds} | awk '$4>={params.integer}' > {output}"""

rule wallikericr01_coverage_overall:
### generate a bed file reflecting intervals where {prop} of samples had at least {cov_filter} coverage (see config.yaml)
#use len(open(dbimportspeciescall).readlines())*{prop} to find threshold for including a window
	input:
		beds = expand(config["output"]+"statistics/coverage_beds/{samplename}_wallikericr01_cov{depth}.bed", samplename = powsamples, allow_missing = True),
	output:
		config["output"]+"statistics/coverage_beds/wallikericr01_{prop}cov{depth}_unmerged.bed"
	params:
		integer = lambda wildcards, output: (float(wildcards.prop)*len(powsamples)) if (float(wildcards.prop) * len(powsamples)).is_integer() else int(float(wildcards.prop)*len(powsamples))+1
	conda:
		"envs/bedtools.yaml"
	resources:
		mem_mb = 20000
	shell:
		"bedtools multiinter -i {input.beds} | awk '$4>={params.integer}' > {output}"

###P. falciparum sample coverage rules

rule pf_coverage_bed_files:
### generate bed files for each sample by species call that show the intervals of the genome with over 5X coverage
	input:
		config["input"]+"pf_alignments/{samplename}_{species}.bam"
	output:
		config["output"]+"statistics/coverage_beds/{samplename}_{species}_cov{depth}.bed"
	conda:
		"envs/bedtools.yaml"
	resources:
		mem_mb = 20000
	shell:
		"bedtools genomecov -ibam {input} -bg | awk '$4>={wildcards.depth}' > {output}"
		
rule pf3d7_coverage_overall:
### generate a bed file reflecting intervals where {prop} of samples had at least {cov_filter} coverage (see config.yaml)
#use len(open(dbimportspeciescall).readlines())*{prop} to find threshold for including a window
	input:
		beds = expand(config["output"]+"statistics/coverage_beds/{samplename}_pf3d7_cov{depth}.bed", samplename = pfsamples, allow_missing = True),
	output:
		config["output"]+"statistics/coverage_beds/pf3d7_{prop}cov{depth}_unmerged.bed"
	params:
		integer = lambda wildcards, output: (float(wildcards.prop)*len(pfsamples)) if (float(wildcards.prop) * len(pfsamples)).is_integer() else int(float(wildcards.prop)*len(pfsamples))+1,
		pf_count = len(pfsamples)
	conda:
		"envs/bedtools.yaml"
	resources:
		mem_mb = 20000
	shell:
		"""echo {params.integer};
		echo {params.pf_count};
		bedtools multiinter -i {input.beds} | awk '$4>={params.integer}' > {output}"""

rule merge_coverage_intervals:
### take the coverage window bed files generated in the previous step and merge any directely adjacent windows into single windows with a buffer of 10bp
# This means that continuous regions with at least 60% of samples giving 5X coverage will be combined into a single interval for ortholog masking
	input:
		config["output"]+"statistics/coverage_beds/{species}_{prop}cov{depth}_unmerged.bed"
	output:
		config["output"]+"statistics/coverage_beds/{species}_{prop}cov{depth}_merged.bed"
	conda:
		"envs/bedtools.yaml"
	resources:
		mem_mb = 20000
	shell:
		"bedtools merge -i {input} -d 10 > {output}"


rule triple_ortholog_masker:
# this requires as input a file in which each line has a group ID followed by gene ID, chrom, start, and stop for each ortholog from pf, poc, and pow in that order
### This rule uses a python script to remove any ortholog groups where any orthologs within would be partially or fully masked by our chromosome, gene, and core genome filters
### It also relies on the previous rules to limit only to ovale orthologs within windows with at least some level of coverage (defined in config/config.yaml)
#outputs a single file structured as the input but without ortholog groups where one or multiple would be masked
#unfiltered output does not consider coverage when determining whether to keep a given ortholog set. This collection of orthologs was not used for analysis
	input:
		all_orthos = config["input_beds"]+ config["all_triple_orthos"],
		poc_chr = config["input_beds"]+"curtisigh01_chr.bed",
		poc_gene_mask = config["input_beds"]+"curtisigh01_genemask.bed",
		poc_coverage_mask = config["output"]+"statistics/coverage_beds/curtisigh01_"+config["prop"]+"cov"+config["cov_filter"]+"_merged.bed",
		pow_chr = config["input_beds"]+"wallikericr01_chr.bed",
		pow_gene_mask = config["input_beds"]+"wallikericr01_genemask.bed",
		pow_coverage_mask = config["output"]+"statistics/coverage_beds/wallikericr01_"+config["prop"]+"cov"+config["cov_filter"]+"_merged.bed",
		pf_core = config["input_beds"]+"pf3d7_core.bed"	,
		pf_coverage_mask = config["output"]+"statistics/coverage_beds/pf3d7_"+config["prop"]+"cov"+config["cov_filter"]+"_merged.bed",
	output:
		pf_orthos_masked = config["output"]+config["pf_triple_orthos_masked"],
		poc_orthos_masked = config["output"]+config["poc_triple_orthos_masked"],
		pow_orthos_masked = config["output"]+config["pow_triple_orthos_masked"],
		pf_orthos_unfiltered = config["output"]+config["pf_triple_orthos_unfiltered"],	
		poc_orthos_unfiltered = config["output"]+config["poc_triple_orthos_unfiltered"],
		pow_orthos_unfiltered = config["output"]+config["pow_triple_orthos_unfiltered"],		
	resources:
		mem_mb = 20000	
	script:
	###scripts are run from the location of the Snakefile, not the directory from which the pipeline is run
		"scripts/pf-poc-pow_ortholog_masker.py"


rule ovale_ortholog_masker:
# this requires as input a file in which each line has a group ID followed by gene ID, chrom, start, and stop for each ortholog from poc and pow in that order
### This rule uses a python script to remove any ortholog groups where any orthologs within would be partially or fully masked by our chromosome or gene filters
### It also relies on the previous rules to limit only to ovale orthologs within windows with at least some level of coverage (defined in config/config.yaml)
#outputs a single file structured as the input but without ortholog groups where one or multiple would be masked
#unfiltered output does not consider coverage when determining whether to keep a given ortholog set. This collection of orthologs was not used for analysis
	input:
		ovale_orthos = config["input_beds"]+ config["all_ovale_orthos"],
		poc_chr = config["input_beds"]+"curtisigh01_chr.bed",
		poc_gene_mask = config["input_beds"]+"curtisigh01_genemask.bed",
		poc_coverage_mask = config["output"]+"statistics/coverage_beds/curtisigh01_"+config["prop"]+"cov"+config["cov_filter"]+"_merged.bed",
		pow_chr = config["input_beds"]+"wallikericr01_chr.bed",
		pow_gene_mask = config["input_beds"]+"wallikericr01_genemask.bed",
		pow_coverage_mask = config["output"]+"statistics/coverage_beds/wallikericr01_"+config["prop"]+"cov"+config["cov_filter"]+"_merged.bed",
	output:
		poc_orthos_masked = config["output"]+config["poc_ovale_orthos_masked"],
		pow_orthos_masked = config["output"]+config["pow_ovale_orthos_masked"],
		poc_orthos_unfiltered = config["output"]+config["poc_ovale_orthos_unfiltered"],
		pow_orthos_unfiltered = config["output"]+config["pow_ovale_orthos_unfiltered"],		
	resources:
		mem_mb = 20000	
	script:
	###scripts are run from the location of the Snakefile, not the directory from which the pipeline is run
		"scripts/poc-pow_ortholog_masker.py"

rule ortholog_bed_maker:
### generate bed files from output table files created in previous two rules
	input:
		config["output"]+"ortholog_beds/{species}_{set}_orthologs_{coverage}.tsv"
	output:
		config["output"]+"ortholog_beds/{species}_{set}_orthologs_{coverage}.bed"
	shell:
		'''cat {input}| while read ID CHROM START STOP; do echo "$CHROM\t$START\t$STOP" >> {output}; done'''
