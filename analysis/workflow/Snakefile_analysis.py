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
#	determine origin of the Pf3D7 core genome bed file I got from karamoko
#	add config file system. define input and output directory at top of snakefile, then reference for all inputs (and give specific name for output directory)
#	Add steps to call variants and develop vcf
# 	include those valentin joste samples in our ovale dataset
#	add bcftools normalize step before limiting to biallelic snps? Would left-align variants but otherwise might not be important if we then just exclude indels
#	nSl calculations screwed up by PowCR01 chromosomes including "|" in the name.
#     Make it clear that we manually renamed, uses symlinks to avoid problematic characters, or simply realign once we have updated Pow reference
#   pf3d7 chr3 vcf from pf3k release 6 is wrong. either find this original source online or revert to release 5?
#  should tajima's D calculation include missing sites? currently not using sites with any missingness

#needed input files:
#	vcf of variants called for sample(s) after aligning to a specific reference genome, formatted as "input/vcfs/{project}_{species}.vcf.gz"
#	the reference genome to which sample reads were aligned, formatted as "input/genomes/{species}.fasta". Of note, PowCR01 needed to have its 14 chromosomes renamed manually
#   index files for each reference genome, generated by "bwa-mem2 index", "samtools faidx", and "gatk CreateSequenceDictionary -R"
#	bed file of chromosome intervals for P. o. curtisi and P. o. wallikeri genomes to include (which excludes extrachromosomal contigs), formatted as "input/beds/{species}_chr.bed", also excludes Poc chr10 because of low quality of the contig
#	bed file of specific hypervariable gene family intervals for the P. o. curtisi and P. o. wallikeri genomes, formatted as "input/beds/{species}_genemask.bed
#	bed file for P. falciparum core genome to be included, formatted as "input/beds/{species}_core.bed"
#	Determined hard quality filtering thresholds for ovale variants, given lack of literature on genomics of this organism
#	Use soft quality filter or preapplied quality filter for falciparum variants. Samples from Pf6k use the VQSLOD score filter
#	List of pf-poc one-to-one orthologs in the pf genome, in the same order as the corresponding orthologs from poc, after subsetting and filtering using Snakefile_ortholog_prep.py
#	List of pf-poc one-to-one orthologs in the poc genome, in the same order as the corresponding orthologs from pf, after subsetting and filtering using Snakefile_ortholog_prep.py
#	A map file of old and new names of the PowCR01 chromosomes; the new names cannot contain "|" characters

configfile: "config/config.yaml"


###for ortholog filtering later, we must create sets of pf and poc orthologs
def parse_bed_file(bed_file):
	gene_dict={}
	for line in open(bed_file):
		name = line.split("\t")[0]
		chrom = line.split("\t")[1]
		start = line.split("\t")[2]
		stop = line.split("\t")[3]
		window = int(stop) - int(start) + 1
		gene_dict[name]= [chrom, start, stop, window]
	return gene_dict

pf3d7_orthos = parse_bed_file(config["input_beds"]+ config["pf_orthos_masked"])
curtisigh01_orthos = parse_bed_file(config["input_beds"]+config["poc_orthos_masked"])

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
		###vcf filtering###
		#files = config["input_vcfs"]+"ov1_wallikericr01.vcf.gz",
		#pofiles = expand(config["output"]+"masked_vcfs/ov1_{species}_chrmasked.recode.vcf", species = ["wallikericr01", "curtisigh01"]),
		#pffiles = config["output"]+"masked_vcfs/ov1_pf3d7_masked.recode.vcf",
		#files = expand(config["output"]+"masked_biallelic_vcfs/ov1_{species}_masked_biallelic.vcf.gz", species = ["wallikericr01", "curtisigh01", "pf3d7"]),
		#files = expand(config["output"]+"masked_vcfs/ov1_{species}_masked.vcf.gz", species = ["wallikericr01", "curtisigh01", "pf3d7"]),
		#files = expand(config["output"]+"filtered_vcfs/ov1_{species}_masked_biallelic_filtertag.vcf.gz", species = ["wallikericr01", "curtisigh01"]),
		#files = config["output"]+"filtered_vcfs/ov1_wallikericr01_masked_biallelic_filtertag.vcf.gz",
		#files = expand(config["output"]+"filtered_vcfs/ov1_{species}_masked_biallelic_qfiltered.vcf.gz", species = ["wallikericr01", "curtisigh01", "pf3d7"]),
		#files = expand(config["output"]+"analysis_vcfs/ov1_{species}.vcf.gz", species = ["wallikericr01", "curtisigh01", "pf3d7"]),
		#files = config["output"]+"analysis_vcfs/ov1_pf3d7.vcf.gz",
		#files = expand(config["output"]+"sample_sets/ov1_{species}_monoclonal.vcf.gz", species = ["wallikericr01", "curtisigh01", "pf3d7"]),
		###ortholog subsetting, processing, and management
		#vcfspo = expand(config["output"]+"orthologs/vcfs/ov1_curtisigh01_{geneid}.vcf.gz", geneid=curtisigh01_orthos.keys()),
		#vcfspf = expand(config["output"]+"orthologs/vcfs/ov1_pf3d7_{geneid}.vcf.gz", geneid=pf3d7_orthos.keys()),
		#popi = expand(config["output"]+"orthologs/stats/ov1_curtisigh01_{geneid}.windowed.pi", geneid=curtisigh01_orthos.keys()),
		#pfpi = expand(config["output"]+"orthologs/stats/ov1_pf3d7_{geneid}.windowed.pi", geneid=pf3d7_orthos.keys()),
		#pfpi = expand(config["output"]+"orthologs/stats/ov1_pf3d7_{geneid}_pi.txt", geneid=pf3d7_orthos.keys()),
		#popi = expand(config["output"]+"orthologs/stats/ov1_curtisigh01_{geneid}_pi.txt", geneid=curtisigh01_orthos.keys()),
		#pfpi = config["output"]+"orthologs/ov1_pf3d7_ortholog_pi.txt",
		#popi = config["output"]+"orthologs/ov1_pocgh01_ortholog_pi.txt",
		###Selection signals###
		#nomissvcf = expand(config["output"]+"sample_sets/ov1_{species}_monoclonal_nomissing.vcf.gz", species = ["wallikericr01", "curtisigh01", "pf3d7"]),
		#pfnsl = expand(config["output"]+"selection/ov1_pf3d7_nsl_{chromosome}.nsl.out", chromosome = pf3d7_chr.keys()),
		#pocnsl = expand(config["output"]+"selection/ov1_curtisigh01_nsl_{chromosome}.nsl.out", chromosome = curtisigh01_chr.keys()),
		#pownsl = expand(config["output"]+"selection/ov1_wallikericr01_nsl_{chromosome}.nsl.out", chromosome = wallikericr01_chr.keys()),
		#nsl_compile_poc = config["output"]+"selection/ov1_curtisigh01_nsl_total.txt",
		#nsl_compile = expand(config["output"]+"selection/ov1_{species}_nsl_total.txt", species = ["wallikericr01", "curtisigh01"]),
		### Tajima's D calculations
		#tajimad_poc = config["output"]+"selection/ov1_curtisigh01.Tajima.D",
		tajima_poc_genes = config["output"]+"selection/ov1_curtisigh01_tajimad_genes.txt",
		###Snakemake administrative###
		config = config["output"]+"pipeline/config.yaml",
		snakefile = config["output"]+"pipeline/Snakefile_analysis.py"

### Step _: rename Po wallikeri CR01 chromsomes so they no longer includes "|" characters
rule rename_pow_chr:
 	input: 
 		vcf = config["input_vcfs"]+"{project}_wallikericr01_original.vcf.gz",
 		chrnames = config["input_beds"]+ "wallikericr01_rename_chr.txt"
 	output: 
 		config["input_vcfs"]+"{project}_wallikericr01.vcf.gz"
 	conda:
 		"envs/bcftools.yaml"
	shell:
		"bcftools annotate --rename-chrs {input.chrnames} {input.vcf} | bgzip > {output}"
 		

### Step _: filter vcfs to core genome only
###
###


###For samples aligned to the PocGH01 or PowCR01 genomes, this limits the SNPs to chromosomal contigs only
rule mask_Po_chr:
	input:
	#This step requires that the ovale species names both end in "01", while the falciparum species name does not. This causes this step to only be run on ovale samples
		vcf = config["input_vcfs"]+"{project}_{ovale}01.vcf.gz",
		bed = config["input_beds"]+"{ovale}01_chr.bed"
	params:
		prefix = config["output"]+"masked_vcfs/{project}_{ovale}01_chrmasked"
	output:
		config["output"]+"masked_vcfs/{project}_{ovale}01_chrmasked.recode.vcf"
	conda:
		"envs/masker.yaml"
	shell:
		"vcftools --gzvcf {input.vcf} --out {params.prefix} --bed {input.bed} --recode --recode-INFO-all"

###Recodes default output of mask_Po_chr as gzipped vcf files
rule recode_Po_chrmask:
	input:
		config["output"]+"masked_vcfs/{project}_{ovale}01_chrmasked.recode.vcf"
	output:
		config["output"]+"masked_vcfs/{project}_{ovale}01_chrmasked.vcf.gz"
	
	shell:
		"bgzip -c {input} > {output}"

###For samples aligned to PocGH01 or PowCR01, maks SNPs in specific hypervariable gene families shown in Rutledge 2017 (with PowCR01 orthologs identified by blast)
rule mask_Po_variablegenes:
	input:
		vcf = config["output"]+"masked_vcfs/{project}_{ovale}01_chrmasked.vcf.gz",
		bed = config["input_beds"]+"{ovale}01_genemask.bed"
	params:
		prefix = config["output"]+"masked_vcfs/{project}_{ovale}01_masked"
	output:
		config["output"]+"masked_vcfs/{project}_{ovale}01_masked.recode.vcf"
	conda:
		"envs/masker.yaml"
	shell:
		"""vcftools --gzvcf {input.vcf} --out {params.prefix} --exclude-bed {input.bed} --recode --recode-INFO-all"""


###For samples aligned to Pf3D7, this limits SNPs to the standard Pf3D7 core genome (masking telomeres and some specific gene families/regions)
rule mask_Pf_acessory_genome:
	input: 
		vcf = config["input_vcfs"]+"{project}_pf3d7.vcf.gz",
		bed = config["input_beds"]+"pf3d7_core.bed"
	params:
		prefix = config["output"]+"masked_vcfs/{project}_pf3d7_masked"
	output: config["output"]+"masked_vcfs/{project}_pf3d7_masked.recode.vcf"
		#this command outputs "masked_vcfs/{project}_pf3d7_masked.recode.vcf"
	conda:
		"envs/masker.yaml"
	shell:
		"""vcftools --gzvcf {input.vcf} --out {params.prefix} --bed {input.bed} --recode --recode-INFO-all"""


###final recoding after masking for all species
rule recode_masked:
	input:
		config["output"]+"masked_vcfs/{project}_{species}_masked.recode.vcf"
	output:
		config["output"]+"masked_vcfs/{project}_{species}_masked.vcf.gz"
	shell:
		"bgzip -c {input} > {output}"

### produces index (.tbi) files for all masked vcfs
rule index_masked_vcfs:
	input:
		config["output"]+"masked_vcfs/{project}_{species}_masked.vcf.gz"
	output:
		config["output"]+"masked_vcfs/{project}_{species}_masked.vcf.gz.tbi"
	conda:
		"envs/masker.yaml"
	shell:
		"""bcftools index -t {input}"""



### Step _: quality filter sites based on QC scores, minor allele frequency, and individual-level missingness
###
###


###limits to biallelic SNPs
rule biallelic:
	input:
		vcf = config["output"]+"masked_vcfs/{project}_{species}_masked.vcf.gz",
		index = config["output"]+"masked_vcfs/{project}_{species}_masked.vcf.gz.tbi"
	output:
		config["output"]+"masked_biallelic_vcfs/{project}_{species}_masked_biallelic.vcf.gz"
	params:
		reference = config["input_genomes"]+"{species}.fasta"
	conda:
		"envs/biallelic.yaml"
	shell:
		"gatk SelectVariants -R {params.reference}  -V {input.vcf} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O {output}"

###adds a filter-tag to all SNPs who fail to pass the following hard quality thresholds
rule po_qualfiltertagger:
	input:
		vcf = config["output"]+"masked_biallelic_vcfs/{project}_{ovale}01_masked_biallelic.vcf.gz",
		reference = config["input_genomes"]+"{ovale}01.fasta"
	output:
		config["output"]+"filtered_vcfs/{project}_{ovale}01_masked_biallelic_filtertag.vcf"
	conda:
		"envs/biallelic.yaml"
	shell:
		"""gatk VariantFiltration -R {input.reference} -V {input.vcf} -O {output} \
		--filter-name 'lowQD' \
		--filter-expression """+config["qd_filter"]+""" \
		--filter-name 'highFS' \
		--filter-expression """+config["fs_filter"]+""" \
		--filter-name 'lowMQ' \
		--filter-expression """+config["mq_filter"]+""" \
		--filter-name 'lowMQRankSum' \
		--filter-expression """+config["mqranksum_filter"]+""" \
		--filter-name 'lowReadPosRankSum' \
		--filter-expression """+config["readposranksum_filter"]

rule po_qual_gzip:
	input:
		config["output"]+"filtered_vcfs/{project}_{ovale}01_masked_biallelic_filtertag.vcf"
	output:
		config["output"]+"filtered_vcfs/{project}_{ovale}01_masked_biallelic_filtertag.vcf.gz"
	shell:
		"bgzip -c {input} > {output}"		

###adds a filter-tag to all SNPs who fail the VQSLOD quality score, a default metric provided by Pf6k
###This score and filter-tag are already included, so this rule just copies and indexes the same vcf file
rule pf_qualfiltertagger:
	input:
		vcf = config["output"]+"masked_biallelic_vcfs/{project}_pf3d7_masked_biallelic.vcf.gz",
		reference = config["input_genomes"]+"pf3d7.fasta"
	output:
		config["output"]+"filtered_vcfs/{project}_pf3d7_masked_biallelic_filtertag.vcf.gz"
	conda:
		"envs/biallelic.yaml"
	shell:
		"""cp {input.vcf} {output}
		gatk IndexFeatureFile -I {output}"""

###produces index (.tbi) files for all qc-filter-tagged vcfs
rule index_tagged_vcfs:
		input:
			config["output"]+"filtered_vcfs/{project}_{species}_masked_biallelic_filtertag.vcf.gz"
		output:
			config["output"]+"filtered_vcfs/{project}_{species}_masked_biallelic_filtertag.vcf.gz.tbi"
		conda:
			"envs/masker.yaml"
		shell:
			"""bcftools index -t {input}"""

			
###Filter out all SNPs that fail QC (by hard filtering threshold for Po or the VQSLOD score for Pf)
rule filter_by_tag:
	input: 
		vcf = config["output"]+"filtered_vcfs/{project}_{species}_masked_biallelic_filtertag.vcf.gz",
		index = config["output"]+"filtered_vcfs/{project}_{species}_masked_biallelic_filtertag.vcf.gz.tbi",
		reference = config["input_genomes"]+"{species}.fasta"
	output:
		config["output"]+"filtered_vcfs/{project}_{species}_masked_biallelic_qfiltered.vcf.gz"
	conda:
		"envs/biallelic.yaml"
	shell:
		"gatk SelectVariants -R {input.reference} -V {input.vcf} -O {output} --exclude-filtered"


###remove sites with a MAF lower than a threshold
rule filter_maf:
	input:
		config["output"]+"filtered_vcfs/{project}_{species}_masked_biallelic_qfiltered.vcf.gz"
	output:
		config["output"]+"filtered_vcfs/{project}_{species}_masked_biallelic_qmaffiltered.vcf.gz"
	conda:
		"envs/filter.yaml"
	shell:
#the following command removes any sites with a minor allele frequency less than 5%
		"vcftools --gzvcf {input} --out {output} --recode --maf 0.05 --stdout | bgzip > {output}"

###only keep sites which are present (nonmissing) in XX% of samples
rule filter_missing:
	input:
		config["output"]+"filtered_vcfs/{project}_{species}_masked_biallelic_qmaffiltered.vcf.gz"
	output:
		config["output"]+"filtered_vcfs/{project}_{species}_masked_biallelic_qmafmissfiltered.vcf.gz"
	conda:
		"envs/filter.yaml"
	shell:
#This command limits the vcf to sites which were present in at least 80% of individuals
		"vcftools --gzvcf {input} --out {output} --max-missing 0.8 --recode --stdout | bgzip > {output}"

###Creates final cleaned vcf for downstream analysis
rule final_set:
#moves annotated processed vcfs to a new directory for further analysis
#vcfs in the analysis_vcfs directory have been limited to the core genome (in ovale, this means chromosomal contigs with a specific list 
# of hypervariable genes manually masked based on Rutledge 2017; limited to biallelic snps; and quality filtered by hard thresholds (in ovale) or the 
# default VQSLOD filter applied in the Pf6k dataset (for falciparum), minor allele frequency of 5%, and presence in >80% of samples 
	input:
		config["output"]+"filtered_vcfs/{project}_{species}_masked_biallelic_qmafmissfiltered.vcf.gz"
	output:
		config["output"]+"analysis_vcfs/{project}_{species}.vcf.gz"
	shell:
		"cp {input} {output}"


rule index_final_vcfs:
	input:
		config["output"]+"analysis_vcfs/{project}_{species}.vcf.gz"
	output:
		config["output"]+"analysis_vcfs/{project}_{species}.vcf.gz.tbi"
	conda:
		"envs/masker.yaml"
	shell:
		"""bcftools index -t {input}"""


rule subset_monoclonal:
###excludes polyclonal samples, as determined using RealMcCOIl
	input:
		vcf = config["output"]+"analysis_vcfs/{project}_{species}.vcf.gz",
		index = config["output"]+"analysis_vcfs/{project}_{species}.vcf.gz.tbi",
		samplelist = config["input_lists"] + "{project}_{species}_polyclonalsamples.args",
		reference = config["input_genomes"]+"{species}.fasta"
	output:
		config["output"]+"sample_sets/{project}_{species}_monoclonal.vcf.gz"
	conda:
		"envs/biallelic.yaml"
	shell:
		"gatk SelectVariants --exclude-sample-name {input.samplelist} -R {input.reference} -V {input.vcf} -O {output}"

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


### Section ___: Calculate nucleotide diversity (pi) for 1-to-1 orthologs of PocGH01 and Pf3D7


		
###create distinct vcf files for each one-to-one ortholog between PocGH01 and Pf3D7

rule subset_po_orthologs:
###Accepts vcfs and a modified bed file (with gene IDs, chromosome IDs, and start/stop positions of 1-to-1 orthologs)
###and subsets to those loci, outputting a unique vcf for SNPs within each gene
	input:
		povcf = config["output"]+"sample_sets/{project}_curtisigh01_monoclonal.vcf.gz",
		poc_orthos_masked = config["input_beds"]+ config["poc_orthos_masked"]
	output: 
		config["output"]+"orthologs/vcfs/{project}_curtisigh01_{pogeneid}.recode.vcf"
	params:
		pooutfile = config["output"]+"orthologs/vcfs/{project}_curtisigh01_{pogeneid}",
		chrom = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[0],
		start = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[1],
		stop = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[2]
	conda:
		"envs/filter.yaml"
	shell:
		"""vcftools --gzvcf {input.povcf} --out {params.pooutfile} --chr {params.chrom} --from-bp {params.start} --to-bp {params.stop} --recode --recode-INFO-all"""


rule subset_pf_orthologs:
###Accepts vcfs and a modified bed file (with gene IDs, chromosome IDs, and start/stop positions of 1-to-1 orthologs)
###and subsets to those loci, outputting a unique vcf for SNPs within each gene
	input:
		pfvcf = config["output"]+"sample_sets/{project}_pf3d7_monoclonal.vcf.gz",
		pf_orthos_masked = config["input_beds"]+ config["pf_orthos_masked"]
	output:
		config["output"]+"orthologs/vcfs/{project}_pf3d7_{pfgeneid}.recode.vcf"
	params:
		pfoutfile = config["output"]+"orthologs/vcfs/{project}_pf3d7_{pfgeneid}",
		chrom = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[0],
		start = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[1],
		stop = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[2]
	conda:
		"envs/filter.yaml"
	shell:
		"""vcftools --gzvcf {input.pfvcf} --out {params.pfoutfile} --chr {params.chrom} --from-bp {params.start} --to-bp {params.stop} --recode --recode-INFO-all"""

rule gzip_orthologs:
	input:
		genes = config["output"]+"orthologs/vcfs/{project}_{species}_{geneid}.recode.vcf"
	output:
		outfile = config["output"]+"orthologs/vcfs/{project}_{species}_{geneid}.vcf.gz"
	shell:
		"bgzip -c {input.genes} > {output.outfile}"

rule calculate_po_ortho_pi:
###outputs individual summary files for reach set of orthologs containing the pi of each ortholog
	input:
		povcf = config["output"]+"sample_sets/{project}_curtisigh01_monoclonal.vcf.gz",
		poc_orthos_masked = config["input_beds"]+config["poc_orthos_masked"]
	output:
		popi = config["output"]+"orthologs/stats/{project}_curtisigh01_{pogeneid}.windowed.pi"
	params:
		pooutfile = config["output"]+"orthologs/stats/{project}_curtisigh01_{pogeneid}",
		chrom = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[0],
		start = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[1],
		stop = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[2],
		window = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[3]
	conda:
		"envs/filter.yaml"
	shell:
		"""vcftools --gzvcf {input.povcf} --out {params.pooutfile}  --chr {params.chrom} --from-bp {params.start} --to-bp {params.stop} --window-pi {params.window} --window-pi-step 1"""

		
rule calculate_pf_ortho_pi:
	input:
		pfvcf = config["output"]+"sample_sets/{project}_pf3d7_monoclonal.vcf.gz",
		pf_orthos_masked = config["input_beds"]+config["pf_orthos_masked"]
	output:
		pfpi = config["output"]+"orthologs/stats/{project}_pf3d7_{pfgeneid}.windowed.pi"
	params:
		pfoutfile = config["output"]+"orthologs/stats/{project}_pf3d7_{pfgeneid}",
		chrom = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[0],
		start = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[1],
		stop = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[2],
		window = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[3]
	conda:
		"envs/filter.yaml"
	shell:
		"""vcftools --gzvcf {input.pfvcf} --out {params.pfoutfile}  --chr {params.chrom} --from-bp {params.start} --to-bp {params.stop} --window-pi {params.window} --window-pi-step 1"""
	
rule select_pfpi:
	input:
		pfpi = config["output"]+"orthologs/stats/{project}_pf3d7_{pfgeneid}.windowed.pi"
	output:
		pfpitable = config["output"]+"orthologs/stats/{project}_pf3d7_{pfgeneid}_pi.txt"
	params:
		chrom = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[0],
		start = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[1],
		stop = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[2],
		window = lambda wildcards, output: list(pf3d7_orthos[wildcards.pfgeneid])[3]
	shell:
		'''if cat {input.pfpi} | while read CHROM START END N_VARIANTS PI; do echo "{wildcards.pfgeneid} $CHROM $START $END $N_VARIANTS $PI"; done | grep "{params.start} {params.stop}"; then
		cat {input.pfpi} | while read CHROM START END N_VARIANTS PI; do echo "{wildcards.pfgeneid} $CHROM $START $END $N_VARIANTS $PI"; done | grep "{params.start} {params.stop}" >> {output.pfpitable};
		else echo "{wildcards.pfgeneid} {params.chrom} {params.start} {params.stop} 0 0" >> {output.pfpitable};fi'''

rule select_popi:
	input:
		popi = config["output"]+"orthologs/stats/{project}_curtisigh01_{pogeneid}.windowed.pi"
	output:
		popitable = config["output"]+"orthologs/stats/{project}_curtisigh01_{pogeneid}_pi.txt"
	params:
		chrom = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[0],
		start = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[1],
		stop = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[2],
		window = lambda wildcards, output: list(curtisigh01_orthos[wildcards.pogeneid])[3]
	shell:
		'''if cat {input.popi} | while read CHROM START END N_VARIANTS PI; do echo "{wildcards.pogeneid} $CHROM $START $END $N_VARIANTS $PI"; done | grep "{params.start} {params.stop}"; then
				cat {input.popi} | while read CHROM START END N_VARIANTS PI; do echo "{wildcards.pogeneid} $CHROM $START $END $N_VARIANTS $PI"; done | grep "{params.start} {params.stop}" >> {output.popitable};
				else echo "{wildcards.pogeneid} {params.chrom} {params.start} {params.stop} 0 0" >> {output.popitable};fi'''
				
rule compile_pi:
	input:
		pfpitable = expand(config["output"]+"orthologs/stats/{project}_pf3d7_{pfgeneid}_pi.txt", pfgeneid=pf3d7_orthos.keys(), allow_missing=True),
		popitable = expand(config["output"]+"orthologs/stats/{project}_curtisigh01_{pogeneid}_pi.txt", pogeneid=curtisigh01_orthos.keys(), allow_missing=True),
	output:
		pfpi = config["output"]+"orthologs/{project}_pf3d7_ortholog_pi.txt",
		popi = config["output"]+"orthologs/{project}_pocgh01_ortholog_pi.txt"
	shell:
		"""for i in {input.pfpitable}; do cat $i >> {output.pfpi}; done; for i in {input.popitable}; do cat $i >> {output.popi}; done"""


### Section _: Signatures of Selection

rule gunzip_vcfs:
	input:
		vcf = config["output"]+"sample_sets/{project}_{species}_monoclonal.vcf.gz"
	output:
		config["output"]+"sample_sets/{project}_{species}_monoclonal.vcf"
	shell:
		"gunzip -c {input.vcf} > {output}"

rule remove_missing:
	###selscan requires no missing genotypes in order to impute selection, so we will filter out all sites that are missing any genotypes
	input: 
		config["output"]+"sample_sets/{project}_{species}_monoclonal.vcf"
	output: 
		config["output"]+"sample_sets/{project}_{species}_monoclonal_nomissing.recode.vcf"
	params:
		outfile = config["output"]+"sample_sets/{project}_{species}_monoclonal_nomissing"
	conda: 
		"envs/filter.yaml"
	shell: 
		"vcftools --gzvcf {input} --out {params.outfile} --max-missing 1 --recode --recode-INFO-all"

rule gzip_nsl:
	input:
		config["output"]+"sample_sets/{project}_{species}_monoclonal_nomissing.recode.vcf"
	output:
		config["output"]+"sample_sets/{project}_{species}_monoclonal_nomissing.vcf.gz"
	shell:
		"bgzip -c {input} > {output}"

rule index_nomissing_vcfs:
	input:
		config["output"]+"sample_sets/{project}_{species}_monoclonal_nomissing.vcf.gz"
	output:
		config["output"]+"sample_sets/{project}_{species}_monoclonal_nomissing.vcf.gz.tbi"
	conda:
		"envs/masker.yaml"
	shell:
		"""bcftools index -t {input}"""

rule subset_chr:
	input:
		vcf = config["output"]+"sample_sets/{project}_{species}_monoclonal_nomissing.vcf.gz",
		index = config["output"]+"sample_sets/{project}_{species}_monoclonal_nomissing.vcf.gz.tbi",
		bed = config["input_beds"]+"{species}_chr.bed"
	output:
		config["output"]+"chr_sets/{project}_{species}_monoclonal_nomissing_{chromosome}.vcf.gz"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools view -r {wildcards.chromosome} {input.vcf} | gzip > {output}"
		
		
		
rule calc_n_sl:
	input:
		vcf = config["output"]+"chr_sets/{project}_{species}_monoclonal_nomissing_{chromosome}.vcf.gz"
	output:
		config["output"]+"selection/{project}_{species}_nsl_{chromosome}.nsl.out"
	conda:
		"envs/selscan.yaml"
	params:
		outfile = config["output"]+"selection/{project}_{species}_nsl_{chromosome}"
	shell:
		"selscan --nsl --vcf {input.vcf} --out {params.outfile}"

rule compile_poc_n_sl:
	input:
		poc_nsl = expand(config["output"]+"selection/{project}_curtisigh01_nsl_{chromosome}.nsl.out", chromosome = curtisigh01_chr.keys(), allow_missing=True),
	output:
		poc_total = config["output"]+"selection/{project}_curtisigh01_nsl_total.txt"
	shell:
#		"""for i in {input.poc_nsl}; do cat $i | while read ID POS AF L1 L2 NSL; do echo "${{i%.nsl.out}} $POS $NSL" >> {output.poc_total}; done; done"""
		"""for i in {input.poc_nsl}; do cat $i | while read ID POS AF L1 L2 NSL; do echo "${{i:45:2}} $POS $NSL" >> {output.poc_total}; done; done;"""

rule compile_pow_n_sl:
	input:
		pow_nsl = expand(config["output"]+"selection/{project}_wallikericr01_nsl_{chromosome}.nsl.out", chromosome = wallikericr01_chr.keys(), allow_missing=True),
	output:
		pow_total = config["output"]+"selection/{project}_wallikericr01_nsl_total.txt"
	shell:
		"""for i in {input.pow_nsl}; do cat $i | while read ID POS AF L1 L2 NSL; do echo "${{i:45:2}} $POS $NSL" >> {output.pow_total}; done; done;"""



rule tajiima_d:
	#currently using monoclonal samples only and with only non-missing sites
	input:
		config["output"]+"sample_sets/{project}_{species}_monoclonal_nomissing.vcf.gz"
	output:
		config["output"]+"selection/{project}_{species}.Tajima.D"
	params:
		outfile = config["output"]+"selection/{project}_{species}",
		window = config["tajima_window"]
	conda:
		"envs/filter.yaml"
	shell:
		"vcftools --gzvcf {input} --out {params.outfile} --TajimaD {params.window}"


rule select_tajima_d:
	input:
		stats = config["output"]+"selection/{project}_{species}.Tajima.D",
		genes = config["input_beds"]+ "{species}_genes.bed",
		exons = config["input_beds"]+ "{species}_exons.bed"	
	output:
		genes = config["output"]+"selection/{project}_{species}_tajimad_genes.txt",
		exons = config["output"]+"selection/{project}_{species}_tajimad_exons.txt"
	params:
		window = config["tajima_window"]
	script:
		"scripts/gene_exon_tajima.py"
		
		
