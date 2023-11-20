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

#TODO:



#needed input files:
#   final analysis vcfs for each species
#   lists of monoclonal sample names for each species

configfile: "config/config.yaml"

###create list of sample names
curtisigh01_names = []
for i in open(config["input_lists"]+"sample_names_monoclonal_curtisigh01.txt").readlines():
        curtisigh01_names.append(i.rstrip("\n"))

wallikericr01_names = []
for i in open(config["input_lists"]+"sample_names_monoclonal_wallikericr01.txt").readlines():
        wallikeri_names.append(i.rstrip("\n"))

###for nsl calculation, we also need a list of chromosome names for each species
###nSL is calculated for haplotypes on each chromosome, so we must subset vcfs by chromosome. This function derives a list of chromosome names from the chr.bed files
def	chrom_names(bed_file):
	chr_dict={}
	count = 0
	for line in open(bed_file):
		if count > 0:
			chrom = line.split("\t")[0]
			start = line.split("\t")[1]
			chr_dict[chrom]= [start]
		count += 1
	return chr_dict

pf3d7_chr = chrom_names(config["input_beds"]+ "pf3d7_chr.bed")
#print(pf3d7_chr)
curtisigh01_chr = chrom_names(config["input_beds"]+ "curtisigh01_chr.bed")
#print(curtisigh01_chr)
wallikericr01_chr = chrom_names(config["input_beds"]+ "wallikericr01_chr.bed")
#print(wallikericr01_chr)


###Determines the final files to be output by the snakemake pipeline
rule all:
#Previous rules should use wildcards for the project and species, but this final rule should employ the actual names and terms for the final file to be created
	input:
		###generate consensus sequences per sample
		poc_consensus = expand(config["output"]+"consensus_seq/ov1_curtisigh01_speciescall-{sample}.fasta", sample = curtisigh01_names),
		pow_consensus = expand(config["output"]+"consensus_seq/ov1_wallikericr01_speciescall-{sample}.fasta", sample = wallikericr01_names),
		###Snakemake administrative###
		config = config["output"]+"pipeline/config.yaml",
		snakefile = config["output"]+"pipeline/Snakefile_phylogeny.py"


###Transfers a copy of the config file and Snakefile to the output for future reference of those results
rule copy_snakefileandconfig:
	input:
		config = "config/config.yaml",
		snakefile = "workflow/Snakefile_analysis.py"
	output:
		config = config["output"]+"pipeline/config.yaml",
		snakefile = config["output"]+"pipeline/Snakefile_analysis.py"
	shell:
		"""cp {input.config} {output.config}
		cp {input.snakefile} {output.snakefile}"""


rule vcf_2_fasta:
	input:
		vcf = config["output"]+"sample_sets/{project}_{species}_{mixed}_monoclonal_allmaf.vcf.gz",
		index = config["output"]+"sample_sets/{project}_{species}_{mixed}_monoclonal_allmaf.vcf.gz.tbi",
		ref = config["input_genomes"]+"{species}.fasta"
	output: 
		config["output"]+"consensus_seq/{project}_{species}_{mixed}-{sample}.fasta"
	conda:
		"envs/bcftools.yaml"
	shell:
		#-I argument causes consensus sequences to use IUPAC codes for ambiguous base calls
		"bcftools consensus -I -f {input.ref} -s {wildcards.sample} -o {output} {input.vcf}"



