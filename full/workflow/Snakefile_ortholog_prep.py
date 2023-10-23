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

pocsamples = []
for f in open(config["input_lists"]+"dbimport_curtisigh01_speciescall_samplemap.txt").readlines():
        pocsamples.append(f.split("\t")[0])

print(len(pocsamples))
print(len(pocsamples)* 0.8)

powsamples = []
for f in open(config["input_lists"]+"dbimport_wallikericr01_speciescall_samplemap.txt").readlines():
        powsamples.append(f.split("\t")[0])

pfsamples = []
for f in open(config["input_lists"]+"dbimport_pf3d7_samplemap.txt").readlines():
		pfsamples.append(f.split("\t")[0].rstrip("\n"))

###Determines the final files to be output by the snakemake pipeline
rule all:
#Previous rules should use wildcards for the project and species, but this final rule should employ the actual names and terms for the final file to be created
	input:
		pf_orthos_masked = config["input_beds"]+ config["pf_orthos_masked"],
		poc_orthos_masked = config["input_beds"]+ config["poc_orthos_masked"],
		pow_orthos_masked = config["input_beds"]+ config["pow_orthos_masked"],

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
#use len(open(dbimportspeciescall).readlines())*0.8 to find threshold for including a window
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
#use len(open(dbimportspeciescall).readlines())*0.8 to find threshold for including a window
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
### generate bed files for each sample by species call that show the intervals of the genome with over 5X covergae
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
#use len(open(dbimportspeciescall).readlines())*0.8 to find threshold for including a window
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
# This means that continuous regions with at least 80% of samples giving 5X coverage will be combined into a single interval for ortholog masking
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


rule ortholog_masker:
# this requires as input a file in which each line has a group ID followed by gene ID, chrom, start, and stop for each ortholog from pf, poc, and pow in that order
### This rule uses a python script to remove any ortholog groups where any orthologs within would be partially or fully masked by our chromosome, gene, and core genome filters
### It also relies on the previous rules to limit only to ovale orthologs within windows with at least some level of coverage (defined in config/config.yaml)
#outputs a single file structured as the input but without ortholog groups where one or multiple would be masked
	input:
		all_orthos = config["input_beds"]+ config["all_orthos"],
		poc_chr = config["input_beds"]+"curtisigh01_chr.bed",
		poc_gene_mask = config["input_beds"]+"curtisigh01_genemask.bed",
		poc_coverage_mask = config["output"]+"statistics/coverage_beds/curtisigh01_"+config["prop"]+"cov"+config["cov_filter"]+"_merged.bed",
		pow_chr = config["input_beds"]+"wallikericr01_chr.bed",
		pow_gene_mask = config["input_beds"]+"wallikericr01_genemask.bed",
		pow_coverage_mask = config["output"]+"statistics/coverage_beds/wallikericr01_"+config["prop"]+"cov"+config["cov_filter"]+"_merged.bed",
		pf_core = config["input_beds"]+"pf3d7_core.bed"	,
		pf_coverage_mask = config["output"]+"statistics/coverage_beds/pf3d7_"+config["prop"]+"cov"+config["cov_filter"]+"_merged.bed",
	output:
		pf_orthos_masked = config["input_beds"]+config["pf_orthos_masked"],
		poc_orthos_masked = config["input_beds"]+config["poc_orthos_masked"],
		pow_orthos_masked = config["input_beds"]+config["pow_orthos_masked"],
		pf_orthos_unfiltered = config["input_beds"]+config["pf_orthos_unfiltered"],	
		poc_orthos_unfiltered = config["input_beds"]+config["poc_orthos_unfiltered"],
		pow_orthos_unfiltered = config["input_beds"]+config["pow_orthos_unfiltered"],		
	resources:
		mem_mb = 20000	
	script:
	###scripts are run from the location of the Snakefile, not the directory from which the pipeline is run
		"scripts/ortholog_masker.py"

