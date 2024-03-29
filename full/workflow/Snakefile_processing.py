## This file allows for alignment and variant calling of whole-genome sequencing data of P. ovale parasite DNA isolates

###Setting up conda
#navigate to directory
#Load anaconda through installation or local version
#module load anaconda/2021.11
#set up conda environment
#conda create -c conda-forge -c bioconda -n snakemake snakemake
#conda activate snakemake
#conda config --set channel_priority strict
###then run the snakefile from the snakemake/directory (Snakefile is the in the workflow directory)
#snakemake --use-conda --cores 1 --jobs 1 -s workflow/Snakefile

###needed input files:
#  Fastq.gz files containing sequencing reads for all samples. This pipeline used paired-end reads (R1 and R2 files) and, for some samples, multiple sequencing lanes that will be merged
#  List of sample names
#  List of sample file names
#  List of sample file names that must be merged (run in parallel on separate sequencing lanes)
#  List of sample file names that do not need to be merged (only run on a single sequencing lane)
#  Table containing readgroup information for all samples
#  Bed file containing desired genomic contigs to be analyzed. In this case, chromosomal contigs for both P. ovale species
#  Reference genome for alignment and processing in .fasta format. Must align with bed files
#  sample name map for DBimportation step, showing which gvcfs for which samples and species should be compiled



configfile: "config/config.yaml"


#Collects names of samples and sample files into python lists. Also determines separate lists of sample names for those with two lanes of sequencing to be merged vs only one
samplenames = []
for i in open(config["input_lists"]+"sample_names.txt").readlines():
	samplenames.append(i.rstrip("\n"))


samplefiles = []
for f in open(config["input_lists"]+"sample_files.txt").readlines():
	samplefiles.append(f.rstrip("\n"))

samples_to_merge = []
for g in open(config["input_lists"]+"samples_to_merge.txt").readlines():
	samples_to_merge.append(g.rstrip("\n"))

#print(samples_to_merge)

samples_no_merge = []
for h in open(config["input_lists"]+"samples_no_merge.txt").readlines():
	samples_no_merge.append(h.rstrip("\n"))

#print(samples_no_merge)
#pull readgroups from comma-delimited file containing sample codes and metadata about sequencing run
def pull_readgroups(readgroup_csv):
        readgroup_dict={}
        for line in open(readgroup_csv):
                samplename = line.split(",")[0]
                RGID = line.split(",")[1]
                RGLB = line.split(",")[2]
                RGPU = line.split(",")[3]
                RGPL = line.split(",")[4]
                RGSM = line.split(",")[5].rstrip("\n")
                readgroup_dict[samplename]= [RGID,RGLB,RGPU,RGPL,RGSM]
        return readgroup_dict

readgroups = pull_readgroups(config["input_lists"]+"readgroups.csv")

def     chrom_names(bed_file):
        chr_dict={}
        count = 0
        for line in open(bed_file):
                if count > 0:
                        chrom = line.split("\t")[0]
                        start = line.split("\t")[1]
                        stop = line.split("\t")[2].rstrip("\n")
                        chr_dict[chrom]= [start,stop]
                count += 1
        return chr_dict

curtisigh01_chr = chrom_names(config["input_beds"]+ "curtisigh01_chr.bed")
wallikericr01_chr = chrom_names(config["input_beds"]+ "wallikericr01_chr.bed")


###Determines the final files to be output by the snakemake pipeline
rule all:
#Previous rules should use wildcards for the project and species, but this final rule should employ the actual names and terms for the final file to be created
#Comment in/out specific file endpoints to adjust the portion of the pipeline that will be run
	input:
		###confirm the presence of the dual reference genomes and needed indices
		dualovale_pf = expand(config["input"]+"genomes/{species}-pf3d7.fasta", species = ["curtisigh01","wallikericr01"]),
		fai = expand(config["input"]+"genomes/{species}-pf3d7.fasta.fai", species = ["curtisigh01","wallikericr01"]),
		bwa_index = expand(config["input"]+"genomes/{species}-pf3d7.fasta.0123", species = ["curtisigh01","wallikericr01"]),
		gatk_dict = expand(config["input"]+"genomes/{species}-pf3d7.dict", species = ["curtisigh01","wallikericr01"]),
		###Check initial quality of reads
		fastqc = expand(config["output"]+"fastqc/untrimmed/{samplefile}_fastqc.zip", samplefile = samplefiles),
		###align to both reference genomes
		alignments_bothlanes = expand(config["output"]+"alignments/{samplename}_{lane}_{species}-pf3d7.bam", samplename = samples_to_merge, lane = ["L001","L002"], species = ["curtisigh01","wallikericr01"]),
		alignments_onelane = expand(config["output"]+"alignments/{samplename}_L001_{species}-pf3d7.bam", samplename = samples_no_merge, species = ["curtisigh01","wallikericr01"]),
		###generate sorted sam files
		sort_bothlanes = expand(config["output"]+"alignments/sorted/{samplename}_{lane}_{species}-pf3d7_sorted.bam", samplename = samples_to_merge, lane = ["L001","L002"], species = ["curtisigh01","wallikericr01"]),
		sort_onelane = expand(config["output"]+"alignments/sorted/{samplename}_L001_{species}-pf3d7_sorted.bam", samplename = samples_no_merge, species = ["curtisigh01","wallikericr01"]),
		###merge across sequencing lanes for samples run on multiple lanes
		merged = expand(config["output"]+"alignments/sorted/merged/{samplename}_{species}-pf3d7_sorted_merged.bam", samplename = samples_to_merge, species = ["curtisigh01","wallikericr01"]),
		unmerged = expand(config["output"]+"alignments/sorted/merged/{samplename}_{species}-pf3d7_sorted_merged.bam", samplename = samples_no_merge, species = ["curtisigh01","wallikericr01"]),
		totalmerge = expand(config["output"]+"merged_alignments/{samplename}_{species}-pf3d7.bam", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		###Add readgroup informaton
		readgrouped = expand(config["output"]+"merged_alignments/readgrouped/{samplename}_{species}-pf3d7.bam", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		###Perform deduplication step
		dedupped = expand(config["output"]+"merged_alignments/readgrouped/dedupped/{samplename}_{species}-pf3d7.bam", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		###Select reads that align better to P. ovale reference genome than P. falciparum
		ovaleselected = expand(config["output"]+"ovale_alignments/{samplename}_{species}.bam", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		###Perform samtools clean step
		cleaned = expand(config["output"]+"ovale_alignments/cleaned/{samplename}_{species}.bam", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		###Generate various alignment statistics
		ovale_cov = expand(config["output"]+"statistics_alignments/genomecov/ovale/{samplename}_{species}_genomecov{depth}.txt", samplename = samplenames, species = ["curtisigh01","wallikericr01"], depth = ["1","5","10"]),
		competitive_cov = expand(config["output"]+"statistics_alignments/genomecov/competitive/{samplename}_{dual}_genomecov{depth}.txt", samplename = samplenames, dual = ["curtisigh01-pf3d7","wallikericr01-pf3d7"], depth = ["1","5","10"]),
		ovale_flagstats = expand(config["output"]+"statistics_alignments/flagstats/ovale/{samplename}_{species}_flagstats.txt", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		competitive_flagstats = expand(config["output"]+"statistics_alignments/flagstats/competitive/{samplename}_{dual}_flagstats.txt", samplename = samplenames, dual = ["curtisigh01-pf3d7","wallikericr01-pf3d7"]),
		coverage = expand(config["output"]+"statistics_alignments/coverage/{samplename}_{species}_coverage.txt", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		overall_coverage = expand(config["output"]+"statistics_alignments/overall_coverage/{species}_coverage-at-{depth}.txt", species = ["curtisigh01","wallikericr01"], depth = ["1","5","10"]),
		###Call variants across all samples using both reference genomes, then combine gvcfs only among samples of the appropriate species as determined by separate species-specific qPCR
		gvcf = expand(config["output"]+"gvcfs/samples/{samplename}_{species}.g.vcf.gz", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		db_species = expand(config["output"]+"gvcfs/grouped/{species}_{mixed}", species = ["curtisigh01","wallikericr01"], mixed = ["andmixed","only"]),
		vcf = expand(config["output"]+"vcfs/ov1/ov1_{species}_{mixed}.vcf.gz", species = ["curtisigh01","wallikericr01"], mixed = ["speciescall"]),
		###Generate administrative documents showing the pipeline and config inputs
		config = config["output"]+"pipeline/config.yaml",
		snakefile = config["output"]+"pipeline/Snakefile_processing.py"


rule copy_snakefileandconfig:
###Transfers a copy of the config file and Snakefile to the output for future reference of those results
	input:
		config = "config/config.yaml",
		snakefile = "workflow/Snakefile_processing.py"
	output:
		config = config["output"]+"pipeline/config.yaml",
		snakefile = config["output"]+"pipeline/Snakefile_processing.py"
	shell:
		"""cp {input.config} {output.config}
		cp {input.snakefile} {output.snakefile}"""

rule dual_reference:
###prepares concatenated reference genomes of curtisi and wallikeri
	input:
		ovale = config["input_genomes"]+"{ovale}01.fasta",
		pf3d7 = config["input_genomes"]+"pf3d7.fasta"
	output:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta"
	shell:
		"cat {input.ovale} {input.pf3d7} > {output}"

rule samtools_index:
#generate samtools index file of reference genome
	input:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta"
	output:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta.fai"
	conda:
		"envs/samtools.yaml"
	resources:
		mem_mb = 60000
	shell:
		"samtools faidx {input}"

rule bwa_index:
#generate bwa-mem index file of reference genome
	input:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta"
	output:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta.0123"
	conda:
		"envs/bwa-mem2.yaml"
	shell:
		"bwa-mem2 index {input}"

rule gatk_index:
#generate gatk index file of reference genome
	input:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta"
	output:
		config["input_genomes"]+"{ovale}01-pf3d7.dict"
	conda:
		"envs/gatk.yaml"
	shell:
		"gatk CreateSequenceDictionary -R {input}"
		

rule fastqc_raw:
# generate fastqc files representing sequencing quality metrics
	input:
		config["input"]+"reads/{samplefile}.fastq.gz"
	output:
		config["output"]+"fastqc/untrimmed/{samplefile}_fastqc.zip"
	params:
		outdirectory = config["output"]+"fastqc/untrimmed/"
	conda:
		"envs/fastqc.yaml"
	shell:
		"fastqc {input} -o {params.outdirectory}"
		
rule trimmomatic:
# trim sequencing adapters
	input:
		R1 = config["input_reads"]+"{samplename}_{lane}_R1.fastq.gz",
		R2 = config["input_reads"]+"{samplename}_{lane}_R2.fastq.gz"
	output:
		P1 = config["output"]+"trimmed_reads/{samplename}_{lane}_1P.fastq.gz",
		#U1 = config["output"]+"trimmed_reads/{samplename}_{lane}_1U.fastq.gz",
		P2 = config["output"]+"trimmed_reads/{samplename}_{lane}_2P.fastq.gz",
		#U2 = config["output"]+"trimmed_reads/{samplename}_{lane}_2U.fastq.gz",
	params:
		outfile = config["output"]+"trimmed_reads/{samplename}_{lane}.fastq.gz"
	conda:
		"envs/trimmomatic.yaml"
	shell:
		"trimmomatic PE -threads 8 {input.R1} {input.R2} -baseout {params.outfile} ILLUMINACLIP:/nas/longleaf/apps/trimmomatic/0.36/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads"

rule fastqc_trimmed:
# generate fastqc files of sequencing data following trimming
	input:
		config["output"]+"trimmed_reads/{samplename}_{lane}_{paired}.fastq.gz"
	output:
		config["output"]+"fastqc/trimmed/{samplename}_{lane}_{paired}_fastqc.zip"
	params:
		outdirectory = config["output"]+"fastqc/trimmed/"
	conda:
		"envs/fastqc.yaml"
	shell:
		"fastqc {input} -o {params.outdirectory}"


rule align:
# align sample reads to P. ovale reference genomes
	input:
		P1 = config["output"]+"trimmed_reads/{samplename}_{lane}_1P.fastq.gz",
		P2 = config["output"]+"trimmed_reads/{samplename}_{lane}_2P.fastq.gz",
		reference = config["input_genomes"]+"{dual}.fasta",
		fai_index = config["input_genomes"]+"{dual}.fasta.fai",
		gatk_index = config["input_genomes"]+"{dual}.dict",
		bwa_index = config["input_genomes"]+"{dual}.fasta.0123",
	output:
		config["output"]+"alignments/{samplename}_{lane}_{dual}.bam"
	conda:
		"envs/bwa-mem2.yaml"
	resources:
		mem_mb = 60000
	shell:
		"bwa-mem2 mem -t 8 {input.reference} {input.P1} {input.P2} -o {output}"
		
rule sort:
# sort aligned reads
	input:
		config["output"]+"alignments/{samplename}_{lane}_{species}.bam"
	output:
		config["output"]+"alignments_sorted/{samplename}_{lane}_{species}.bam"
	conda:
		"envs/picard.yaml"
	resources:
		mem_mb = 60000
	shell:
		"picard SortSam I={input} O={output} sort_order=coordinate create_index=true"

rule prep_merge:
# prep samples that were sequenced in parallel lanes to be merged across lanes
	input:
		expand(config["output"]+"alignments_sorted/{samplename}_{lane}_{species}-pf3d7.bam", samplename = samples_to_merge, species = ["curtisigh01","wallikericr01"], lane = ["L001","L002"])
	output:
		expand(config["output"]+"to_merge/{samplename}_{lane}_{species}-pf3d7.bam", samplename = samples_to_merge, species = ["curtisigh01","wallikericr01"], lane = ["L001","L002"])
	params:
		outdir = config["output"]+"to_merge/"
	shell:
		"for i in {input}; do cp ${{i}} {params.outdir}; done"

rule prep_nomerge:
# copy sample files for samples which were run on a single sequencing lane and do not need to be merge
	input:
		expand(config["output"]+"alignments_sorted/{samplename}_L001_{species}-pf3d7.bam", samplename = samples_no_merge, species = ["curtisigh01","wallikericr01"])
	output:
		expand(config["output"]+"to_merge/{samplename}_L001_{species}-pf3d7.bam", samplename = samples_no_merge, species = ["curtisigh01","wallikericr01"])
	params:
		outdir = config["output"]+"to_merge/"
	shell:
		"for i in {input}; do cp ${{i}} {params.outdir}; done"

rule merge:
# perform merge of samples that were sequenced in parallel lanes
	input:
		lane1 = config["output"]+"to_merge/{samplename}_L001_{species}-pf3d7.bam"
	params:
		lane2 = config["output"]+"to_merge/{samplename}_L002_{species}-pf3d7.bam"
	output:
		config["output"]+"merged_alignments/{samplename}_{species}-pf3d7.bam"
	conda:
		"envs/samtools.yaml"
	resources:
		mem_mb = 60000
	script:
		"scripts/merger.py"

rule add_readgroup:
# add readgroup information for all samples
	input:
		config["output"]+"merged_alignments/{samplename}_{dual}.bam"
	output:
		config["output"]+"readgrouped/{samplename}_{dual}.bam"
	params:
		RGID = lambda wildcards, output: list(readgroups[wildcards.samplename])[0],
		RGLB = lambda wildcards, output: list(readgroups[wildcards.samplename])[1],
		RGPU = lambda wildcards, output: list(readgroups[wildcards.samplename])[2],
		RGPL = lambda wildcards, output: list(readgroups[wildcards.samplename])[3],
		RGSM = lambda wildcards, output: list(readgroups[wildcards.samplename])[4]
	conda:
		"envs/picard.yaml"
	resources:
		mem_mb = 60000
	shell:
		"picard AddOrReplaceReadGroups I={input} O={output} SORT_ORDER=coordinate RGID={params.RGID} RGLB={params.RGLB} RGPL={params.RGPL} RGPU={params.RGPU} RGSM={params.RGSM}"
		
rule deduplicate:
# remove duplicate reads
	input:
		config["output"]+"readgrouped/{samplename}_{dual}.bam"
	output:
		out = config["output"]+"dedupped/{samplename}_{dual}.bam",
		metrics = config["output"]+"dedup_metrics/{samplename}_{dual}_dupmetrics.txt"
	conda:
		"envs/gatk.yaml"
	resources:
		mem_mb = 60000
	shell:
		"gatk MarkDuplicatesSpark -I {input} -O {output.out} -M {output.metrics}"

rule ovale_selecter:
# select reads that align better to P. ovale reference genome portion of dual reference genome than the P. falciparum section.
	input:
		bam = config["output"]+"dedupped/{samplename}_{species}-pf3d7.bam",
		ref =  config["input_genomes"]+"{species}.fasta",
		bed = config["input_beds"]+"{species}.bed"
	output:
		config["output"]+"ovale_alignments/{samplename}_{species}.bam"
	conda:
		"envs/samtools.yaml"
	resources:
		mem_mb = 60000
	shell:
		"samtools view -b -h {input.bam} -T {input.ref} -L {input.bed} > {output}"

rule sam_cleaner:
# perform samtools cleaning step
	input:
		config["output"]+"ovale_alignments/{samplename}_{species}.bam"
	output:
		config["output"]+"ovale_alignments_cleaned/{samplename}_{species}.bam"
	conda:
		"envs/gatk.yaml"
	resources:
		mem_mb = 60000
	shell:
		"gatk CleanSam -I {input} -O {output}"

rule genome_coverage_ovale:
# calculate coverage of reference genome for each sample at a given read depth after selecton of ovale reads
	input:
		config["output"]+"ovale_alignments_cleaned/{samplename}_{species}.bam"
	output:
		config["output"]+"statistics_alignments/genomecov/ovale/{samplename}_{species}_genomecov{depth}.txt"
	params:
		depth = lambda wildcards, output: wildcards.depth
	conda:
		"envs/bedtools.yaml"
	resources:
		mem_mb = 60000
	shell:
		"bedtools genomecov -max {params.depth} -ibam {input} > {output}"

rule genome_coverage_competitive:
# calculate coverage of reference genome for each sample at a given read depth prior to selecton of ovale reads
	input:
		config["output"]+"dedupped/{samplename}_{dual}.bam"
	output:
		config["output"]+"statistics_alignments/genomecov/competitive/{samplename}_{dual}_genomecov{depth}.txt"
	params:
		depth = lambda wildcards, output: wildcards.depth
	conda:
		"envs/bedtools.yaml"
	resources:
		mem_mb = 60000
	shell:
		"bedtools genomecov -max {params.depth} -ibam {input} > {output}"

rule coverage_calc:
###calculates weighted coverage of the fourteen chromosomes for each species after selection of ovale reads
	input:
		cov = expand(config["output"]+"statistics_alignments/genomecov/ovale/{samplename}_{species}_genomecov{depth}.txt", samplename = samplenames, allow_missing = True),
		bed = config["input_beds"]+"{species}_chr.bed"
	output:
		config["output"]+"statistics_alignments/overall_coverage/{species}_coverage-at-{depth}.txt"
	script:
		"scripts/coverager.py"


rule flagstats_ovale:
# calculate mapping percentages and total read alignment after selection of ovale reads
	input:
		config["output"]+"ovale_alignments_cleaned/{samplename}_{species}.bam"
	output:
		config["output"]+"statistics_alignments/flagstats/ovale/{samplename}_{species}_flagstats.txt"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools flagstat {input} > {output}"

rule flagstats_competitive:
# calculate mapping percentages and total read alignment prior to selection of ovale reads
	input:
		config["output"]+"dedupped/{samplename}_{dual}.bam",
	output:
		config["output"]+"statistics_alignments/flagstats/competitive/{samplename}_{dual}_flagstats.txt"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools flagstat {input} > {output}"

rule coveragestats:
# calculate coverage per chromosome for each sample
	input:
		config["output"]+"ovale_alignments_cleaned/{samplename}_{species}.bam"
	output:
		config["output"]+"statistics_alignments/coverage/{samplename}_{species}_coverage.txt"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools coverage {input} > {output}"

###Variant Calling
rule index_bams:
# generate index of cleaned, final alignment .bam files
	input:
		config["output"]+"ovale_alignments_cleaned/{samplename}_{species}.bam"
	output:
		config["output"]+"ovale_alignments_cleaned/{samplename}_{species}.bam.bai"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools index {input}"

rule haplotypecaller_curtisigh01:
# call variants across P. ovale curtisi GH01 reference genome for all samples, on a chromosome-by-chromosome basis for speed
	input:
		bam = config["output"]+"ovale_alignments_cleaned/{samplename}_curtisigh01.bam",
		index = config["output"]+"ovale_alignments_cleaned/{samplename}_curtisigh01.bam.bai",
		ref = config["input_genomes"]+"curtisigh01.fasta",
		#bed = config["input_beds"]+"{species}.bed"
	output:
		config["output"]+"gvcfs/chromosomes/{samplename}_curtisigh01_{chromosome}.g.vcf.gz"
	params:
		java_opts=config["gatk_haplotypecaller"]["java"],
		extra=config["gatk_haplotypecaller"]["extra"]
	conda:
		"envs/gatk.yaml"
	resources:
		mem_mb = 20000,
		nodes = 2
	threads: 2
	shell:
		'''gatk --java-options "{params.java_opts}" HaplotypeCaller --native-pair-hmm-threads 2 {params.extra} -R {input.ref} -L {wildcards.chromosome} -I {input.bam} -O {output}'''


rule combine_chrom_gvcfs_curtisigh01:
# compile chromosomes from each sample's PocGH01 g.vcf files
	input:
		gvcfs = expand(config["output"]+"gvcfs/chromosomes/{samplename}_curtisigh01_{chromosome}.g.vcf.gz", chromosome = curtisigh01_chr.keys(), allow_missing = True),
		ref = config["input_genomes"]+"curtisigh01.fasta",
	output:
		config["output"]+"gvcfs/samples/{samplename}_curtisigh01.g.vcf.gz"
	params:
		var_command = expand("--variant "+config["output"]+"gvcfs/chromosomes/{samplename}_curtisigh01_{chromosome}.g.vcf.gz", chromosome = curtisigh01_chr.keys(), allow_missing = True),
	conda:
		"envs/gatk.yaml"
	resources:
		mem_mb = 20000
	shell:
		"gatk CombineGVCFs -R {input.ref} {params.var_command} -O {output}"

rule haplotypecaller_wallikericr01:
# call variants across P. ovale wallikeri CR01 reference genome for all samples, on a chromosome-by-chromosome basis for speed
	input:
		bam = config["output"]+"ovale_alignments_cleaned/{samplename}_wallikericr01.bam",
		index = config["output"]+"ovale_alignments_cleaned/{samplename}_wallikericr01.bam.bai",
		ref = config["input_genomes"]+"wallikericr01.fasta",
		#bed = config["input_beds"]+"{species}.bed"
	output:
		config["output"]+"gvcfs/chromosomes/{samplename}_wallikericr01_{chromosome}.g.vcf.gz"
	params:
		java_opts=config["gatk_haplotypecaller"]["java"],
		extra=config["gatk_haplotypecaller"]["extra"]
	conda:
		"envs/gatk.yaml"
	resources:
		mem_mb = 20000,
		nodes = 2
	threads: 2
	shell:
		'''gatk --java-options "{params.java_opts}" HaplotypeCaller --native-pair-hmm-threads 2 {params.extra} -R {input.ref} -L {wildcards.chromosome} -I {input.bam} -O {output}'''


rule combine_chrom_gvcfs_wallikericr01:
# compile chromosomes from each sample's PowCR01 g.vcf files
	input:
		gvcfs = expand(config["output"]+"gvcfs/chromosomes/{samplename}_wallikericr01_{chromosome}.g.vcf.gz", chromosome = wallikericr01_chr.keys(), allow_missing = True),
		ref = config["input_genomes"]+"wallikericr01.fasta",
	output:
		config["output"]+"gvcfs/samples/{samplename}_wallikericr01.g.vcf.gz"
	params:
		var_command = expand("--variant "+config["output"]+"gvcfs/chromosomes/{samplename}_wallikericr01_{chromosome}.g.vcf.gz", chromosome = wallikericr01_chr.keys(), allow_missing = True),
	conda:
		"envs/gatk.yaml"
	resources:
		mem_mb = 20000
	shell:
		"gatk CombineGVCFs -R {input.ref} {params.var_command} -O {output}"


rule dbimport_species:
# generate dbimport directory and associated files, compiling g.vcfs for the samples that were identified as one species
	input:
		map = config["input_lists"]+"dbimport_{species}_{mixed}_samplemap.txt",
		gvcfs = expand(config["output"]+"gvcfs/samples/{samplename}_{species}.g.vcf.gz", samplename = samplenames, allow_missing = True),
		bed = config["input_beds"]+"{species}.bed"
	output:
		directory(config["output"]+"gvcfs/grouped/{species}_{mixed}")
	params:
		dir = config["output"]+"gvcfs/grouped/{species}_{mixed}"
	conda:
		"envs/gatk.yaml"
	resources:
		mem_mb = 20000
	shell:
		"gatk GenomicsDBImport --sample-name-map {input.map} --genomicsdb-workspace-path {params.dir} -L {input.bed}"


rule genotype_sample_gvcfs:
# For each species' dbimport space, generate a g-zipped vcf file containing all called variants for all samples of that species
        input:
              	dir = config["output"]+"gvcfs/grouped/{species}_{mixed}",
                ref = config["input_genomes"]+"{species}.fasta"
        output:
               	config["output"]+"vcfs/{project}/{project}_{species}_{mixed}.vcf.gz"
        params:
               	java_opts = config["gatk_genotypegvcfs"]["java"]
        conda:
              	"envs/gatk.yaml"
        resources:
                mem_mb = 2000,
                nodes = 2
        threads: 2
        shell:
              	'''gatk --java-options "{params.java_opts}" GenotypeGVCFs -R {input.ref} -V gendb://{input.dir} -O {output}'''
