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
curtisigh01_names = {}
for i in open(config["input_lists"]+"sample_names_monoclonal_curtisigh01.txt").readlines():
		name = i.split("\t")[0]
		curtisigh01_names[name] = [i.split("\t")[1].rstrip("\n")]



wallikericr01_names = {}
for i in open(config["input_lists"]+"sample_names_monoclonal_wallikericr01.txt").readlines():
		name = i.split("\t")[0]
		wallikericr01_names[name] = [i.split("\t")[1].rstrip("\n")]

all_ovale_names = {}
for i in open(config["input_lists"]+"sample_names_monoclonal_ovale.txt").readlines():
		name = i.split("\t")[0]
		all_ovale_names[name] = [i.split("\t")[1].rstrip("\n")]

#samplenames = [curtisigh01_names.keys(),wallikericr01_names.keys()]

#print(samplenames)

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
		###test mafft
		test = config["output"]+"consensus_seq/multi_alignment/test_chr.fasta",
		###generate consensus sequences per sample
		poc_consensus = expand(config["output"]+"consensus_seq/no-header_fasta/ov1_curtisigh01_speciescall-{sample}.fasta", sample = curtisigh01_names.keys()),
		pow_consensus = expand(config["output"]+"consensus_seq/no-header_fasta/ov1_wallikericr01_speciescall-{sample}.fasta", sample = wallikericr01_names.keys()),
		poc_headers = expand(config["output"]+"consensus_seq/formatted_headers/ov1_curtisigh01_speciescall-{sample}.fasta", sample = curtisigh01_names.keys()),
		pow_headers = expand(config["output"]+"consensus_seq/formatted_headers/ov1_wallikericr01_speciescall-{sample}.fasta", sample = wallikericr01_names.keys()),
		poc_fasta = expand(config["output"]+"consensus_seq/formatted_sample_fasta/ov1_curtisigh01_speciescall-{sample}.fasta", sample = curtisigh01_names.keys()),
		pow_fasta = expand(config["output"]+"consensus_seq/formatted_sample_fasta/ov1_wallikericr01_speciescall-{sample}.fasta", sample = wallikericr01_names.keys()),
		multi_fasta = expand(config["output"]+"consensus_seq/multi_sample_fasta/ov1_{species}_speciescall.fasta", species = ["curtisigh01","wallikericr01","ovale"]),
		#multi_align = expand(config["output"]+"consensus_seq/multi_alignment/ov1_{species}_speciescall.fasta", species = ["curtisigh01","wallikericr01","ovale"]),
		nexus = expand(config["output"]+"consensus_seq/nexus/ov1_{species}_speciescall.nex", species = ["curtisigh01","wallikericr01"]),
		test_nexus = config["output"]+"consensus_seq/nexus/test_chr.nex",
		tree = expand(config["output"]+"consensus_seq/trees/ov1_{ovale}01_speciescall.nex1.mcmc", ovale = ["curtisigh","wallikericr"]),
		all_ovale_tree = config["output"]+"consensus_seq/trees/ov1_wallikericr01_all.nex1.mcmc",
		###Snakemake administrative###
		config = config["output"]+"pipeline/config.yaml",
		snakefile = config["output"]+"pipeline/Snakefile_phylogeny.py"


###Transfers a copy of the config file and Snakefile to the output for future reference of those results
rule copy_snakefileandconfig:
	input:
		config = "config/config.yaml",
		snakefile = "workflow/Snakefile_phylogeny.py"
	output:
		config = config["output"]+"pipeline/config.yaml",
		snakefile = config["output"]+"pipeline/Snakefile_phylogeny.py"
	shell:
		"""cp {input.config} {output.config}
		cp {input.snakefile} {output.snakefile}"""


rule vcf_2_fasta:
	input:
		vcf = config["output"]+"sample_sets/{project}_{species}_{mixed}_monoclonal_allmaf.vcf.gz",
		index = config["output"]+"sample_sets/{project}_{species}_{mixed}_monoclonal_allmaf.vcf.gz.tbi",
		ref = config["input_genomes"]+"{species}.fasta"
	output: 
		config["output"]+"consensus_seq/sample_fasta/{project}_{species}_{mixed}-{sample}.fasta"
	conda:
		"envs/bcftools.yaml"
	shell:
		#-I argument causes consensus sequences to use IUPAC codes for ambiguous base calls
		"bcftools consensus -I -f {input.ref} -s {wildcards.sample} -o {output} {input.vcf}"

rule fasta_masker:
	input:
		fasta = config["output"]+"consensus_seq/sample_fasta/{project}_{species}_{mixed}-{sample}.fasta",
		bed = config["input"]+"beds/{species}_genemask.bed"
	output:
		config["output"]+"consensus_seq/masked_fasta/{project}_{species}_{mixed}_genemasked-{sample}.fasta"
	conda:
		"envs/bedtools.yaml"
	shell:
		"bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output}"

###for introns, exons, cds, genes, and intergenic region beds, need to remove extrachromosomal contigs and have sorted
rule remove_extrachromosomal:
	input:
		fasta = config["output"]+"consensus_seq/masked_fasta/{project}_{species}_{mixed}_genemasked-{sample}.fasta",
		bed = config["input_beds"]+"{species}_chr.bed",
	output:
		config["output"]+ "consensus_seq/masked_fasta/{project}_{species}_{mixed}_masked-{sample}.fasta"
	conda:
		"envs/bedtools.yaml"
	shell:
		"bedtools getfasta -fi {input.fasta} -bed {input.bed} -fo {output}"



rule remove_fasta_headers:
	input:
		config["output"]+"consensus_seq/masked_fasta/{project}_{species}_{mixed}_masked-{sample}.fasta"
	output:
		config["output"]+"consensus_seq/no-header_fasta/{project}_{species}_{mixed}-{sample}.fasta"	
	shell:
		"""awk '$1 !~ /^>/ {{print $0}}' {input} > {output}"""

rule new_poc_headers:
# creates text file containing only the first line (header) of the sample fasta
	input:
		config["output"]+"consensus_seq/sample_fasta/{project}_curtisigh01_{mixed}-{sample}.fasta"
	output:
		config["output"]+"consensus_seq/formatted_headers/{project}_curtisigh01_{mixed}-{sample}.fasta"
	params:
		sample_code = lambda wildcards, output: curtisigh01_names.get(wildcards.sample)		
	shell:
		"echo '>poc_{params.sample_code}' > {output}" 

rule new_pow_headers:
# creates text file containing only the first line (header) of the sample fasta
	input:
		config["output"]+"consensus_seq/sample_fasta/{project}_wallikericr01_speciescall-{sample}.fasta"
	output:
		config["output"]+"consensus_seq/formatted_headers/{project}_wallikericr01_speciescall-{sample}.fasta"
	params:
		sample_code = lambda wildcards, output: wallikericr01_names.get(wildcards.sample)		
	shell:
		"echo '>pow_{params.sample_code}' > {output}" 

rule new_ovale_headers_pow_ref:
# creates text file containing only the first line (header) of the sample fasta
	input:
		config["output"]+"consensus_seq/sample_fasta/{project}_wallikericr01_all-{sample}.fasta"
	output:
		config["output"]+"consensus_seq/formatted_headers/{project}_wallikericr01_all-{sample}.fasta"
	params:
		sample_code = lambda wildcards, output: all_ovale_names.get(wildcards.sample)		
	shell:
		"echo '>pow_{params.sample_code}' > {output}" 


rule add_new_sample_header:
	input:
		fasta = config["output"]+"consensus_seq/no-header_fasta/{project}_{species}_{mixed}-{sample}.fasta",
		header = config["output"]+"consensus_seq/formatted_headers/{project}_{species}_{mixed}-{sample}.fasta"
	output:
		config["output"]+"consensus_seq/formatted_sample_fasta/{project}_{species}_{mixed}-{sample}.fasta"	
	shell:
		"cat {input.header} {input.fasta} > {output}"	

rule concat_ovale_fastas:
	input:
		ovale_fasta = expand(config["output"]+"consensus_seq/formatted_sample_fasta/{project}_wallikericr01_all-{sample}.fasta", sample = all_ovale_names.keys(), allow_missing = True),
	output:
		config["output"]+"consensus_seq/multi_sample_fasta/{project}_wallikericr01_all.fasta"
	shell:
		"cat {input.ovale_fasta} > {output}"

#rule ovale_multi_aligner:
#	input:
#		config["output"]+"consensus_seq/multi_sample_fasta/{project}_ovale_{mixed}.fasta"
#	output:
#		config["output"]+"consensus_seq/multi_alignment/{project}_ovale_{mixed}.fasta"
#	conda:
#		"envs/mafft.yaml"
#	resources:
#		mem_mb = 450000
#	shell:
#		#employs L-INS-i approach
#		"mafft --auto --memsave {input} > {output}"


rule concat_poc_fastas:
	input:
		poc_fasta = expand(config["output"]+"consensus_seq/formatted_sample_fasta/{project}_curtisigh01_{mixed}-{sample}.fasta", sample = curtisigh01_names.keys(), allow_missing = True),
		#pow_fasta = expand(config["output"]+"consensus_seq/formatted_sample_fasta/{project}_wallikericr01_{mixed}-{sample}.fasta", sample = wallikericr01_names.keys(), allow_missing = True),
	output:
		config["output"]+"consensus_seq/multi_sample_fasta/{project}_curtisigh01_{mixed}.fasta"
	shell:
		"cat {input.poc_fasta} > {output}"

#rule poc_multi_aligner:
#	input:
#		config["output"]+"consensus_seq/multi_sample_fasta/{project}_curtisigh01_{mixed}.fasta"
#	output:
#		config["output"]+"consensus_seq/multi_alignment/{project}_curtisigh01_{mixed}.fasta"
#	conda:
#		"envs/mafft.yaml"
#	resources:
#		mem_mb = 450000
#	shell:
		#employs G-INS-i approach
		#"mafft --globalpair --maxiterate 1000 --memsave --thread {threads} {input} > {output}"
#		"mafft --auto --memsave {input} > {output}"



rule concat_pow_fastas:
	input:
		#poc_fasta = expand(config["output"]+"consensus_seq/formatted_sample_fasta/{project}_curtisigh01_{mixed}-{sample}.fasta", sample = curtisigh01_names.keys(), allow_missing = True),
		pow_fasta = expand(config["output"]+"consensus_seq/formatted_sample_fasta/{project}_wallikericr01_speciescall-{sample}.fasta", sample = wallikericr01_names.keys(), allow_missing = True),
	output:
		config["output"]+"consensus_seq/multi_sample_fasta/{project}_wallikericr01_speciescall.fasta"
	shell:
		"cat {input.pow_fasta} > {output}"





#rule pow_multi_aligner:
#	input:
#		config["output"]+"consensus_seq/multi_sample_fasta/{project}_wallikericr01_{mixed}.fasta"
#	output:
#		config["output"]+"consensus_seq/multi_alignment/{project}_wallikericr01_{mixed}.fasta"
#	conda:
#		"envs/mafft.yaml"
#	resources:
#		mem_mb = 450000
#	shell:
		#employs L-INS-i approach
		#"mafft --localpair --maxiterate 20 --memsave --thread {threads} {input} > {output}"
#		"mafft --auto --memsave {input} > {output}"

#rule test_mafft:
#	input:
#		config["output"]+"consensus_seq/sample_fasta/test_chr.fasta"
#	output:
#		config["output"]+"consensus_seq/multi_alignment/test_chr.fasta"
#	conda:
#		"envs/mafft.yaml"
#	resources:
#		mem_mb = 450000,
#		nodes = 15
	#threads: 30
#	shell:
		#employs L-INS-i approach
		#"mafft --localpair --maxiterate 1000 --memsave --thread {threads} {input} > {output}"
		#"mafft --localpair --maxiterate 2 --memsave --thread {threads} {input} > {output}"
#		"mafft --retree 2 --memsave {input} > {output}"
		#"mafft --retree 2 --memsave --thread {threads} {input} > {output}"


rule fasta_to_nexus:
	input:
		config["output"]+"consensus_seq/multi_sample_fasta/{project}_{species}_{mixed}.fasta"
	output:
		config["output"]+"consensus_seq/nexus/{project}_{species}_{mixed}.nex"
	conda:
		"envs/seqmagick.yaml"
	resources:
		mem_mb = 450000
	shell:
		"seqmagick convert --alphabet dna {input} {output}"	

rule test_fasta_to_nexus:
	input:
		config["output"]+"consensus_seq/sample_fasta/test_chr.fasta"
	output:
		config["output"]+"consensus_seq/nexus/test_chr.nex"
	conda:
		"envs/seqmagick.yaml"
	resources:
		mem_mb = 450000
	shell:
		"seqmagick convert --alphabet dna {input} {output}"	

rule write_mrbayes_nexus_execute_file:
	input:
		config["output"]+"consensus_seq/nexus/{project}_{species}_{mixed}.nex"
	output:
		config["output"]+"consensus_seq/execute_files/{project}_{species}_{mixed}_execute.nex"
	params:
		outfile = config["output"]+"consensus_seq/trees/{project}_{species}_{mixed}.nex1"
	shell:
		"""echo 'begin mrbayes;
		set autoclose=yes nowarnings=yes;
		execute {input};
		lset nst=6 rates=gamma;
		mcmc nruns=1 ngen=10000 samplefreq=10 file={params.outfile};
		end;' > {output}"""

rule mrbayes_treemaker:
	input:
		execute = config["output"]+"consensus_seq/execute_files/{project}_{ovale}01_{mixed}_execute.nex",		
	output:
		outfile = config["output"]+"consensus_seq/trees/{project}_{ovale}01_{mixed}.nex1.mcmc",
		log = config["output"]+"consensus_seq/mrbayes_logs/{project}_{ovale}01_{mixed}.nex1_log.txt"
	conda:
		"envs/mrbayes.yaml"
	resources:
		mem_mb = 450000
	shell:
		"mb {input.execute} > {output.log}"
