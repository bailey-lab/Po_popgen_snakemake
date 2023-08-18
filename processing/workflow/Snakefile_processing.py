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

###TODO:
#	update wallikeru reference genome for analysis pipelin
#   when to check with flagstats? after selecting to ovale genome only
#   add a basequalityscorerecalibrator step, which will involde taking in the unrecalibrated vcf?


###needed input files:
#	sample name map for DBimportation step, showing which gvcfs for which samples and species should be compiled


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

print(samples_to_merge)

samples_no_merge = []
for h in open(config["input_lists"]+"samples_no_merge.txt").readlines():
	samples_no_merge.append(h.rstrip("\n"))

print(samples_no_merge)
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



###Determines the final files to be output by the snakemake pipeline
rule all:
#Previous rules should use wildcards for the project and species, but this final rule should employ the actual names and terms for the final file to be created
	input:
		#confirm the presence of the dual reference genomes
		#dualovale_pf = expand(config["input"]+"genomes/{species}-pf3d7.fasta", species = ["curtisigh01","wallikericr01"]),
		#fai = expand(config["input"]+"genomes/{species}-pf3d7.fasta.fai", species = ["curtisigh01","wallikericr01"]),
		#bwa_index = expand(config["input"]+"genomes/{species}-pf3d7.fasta.0123", species = ["curtisigh01","wallikericr01"]),
		#gatk_dict = expand(config["input"]+"genomes/{species}-pf3d7.dict", species = ["curtisigh01","wallikericr01"]),
		###Check initial quality of reads
		#fastqc = expand(config["output"]+"fastqc/untrimmed/{samplefile}_fastqc.zip", samplefile = samplefiles),
		###confirm trimmomatic step works, but only checks for lane 1 and paired 1; all other files should also be produced
		#trimmed_bothlanes = expand(config["output"]+"trimmed_reads/{samplename}_{lane}_1P.fastq.gz", samplename = samples_to_merge, lane =["L001","L002"]),
		#trimmed_onelane = expand(config["output"]+"trimmed_reads/{samplename}_L001_1P.fastq.gz", samplename = samples_no_merge),
		###fastqc trimmed files, must manually check all forms
		#fastqc_bothlanes = expand(config["output"]+"fastqc/trimmed/{samplename}_{lane}_{pair}_fastqc.zip", samplename = samples_to_merge, lane = ["L001","L002"], pair = ["1P","2P"]),
		#fastqc_onelane = expand(config["output"]+"fastqc/trimmed/{samplename}_L001_{pair}_fastqc.zip", samplename = samples_no_merge, pair = ["1P", "2P"]),
		###align to both reference genomes
		alignments_bothlanes = expand(config["output"]+"alignments/{samplename}_{lane}_{species}-pf3d7.bam", samplename = samples_to_merge, lane = ["L001","L002"], species = ["curtisigh01","wallikericr01"]),
		alignments_onelane = expand(config["output"]+"alignments/{samplename}_L001_{species}-pf3d7.bam", samplename = samples_no_merge, species = ["curtisigh01","wallikericr01"]),
		###generate sorted sam files
		sort_bothlanes = expand(config["output"]+"alignments/sorted/{samplename}_{lane}_{species}-pf3d7_sorted.bam", samplename = samples_to_merge, lane = ["L001","L002"], species = ["curtisigh01","wallikericr01"]),
		sort_onelane = expand(config["output"]+"alignments/sorted/{samplename}_L001_{species}-pf3d7_sorted.bam", samplename = samples_no_merge, species = ["curtisigh01","wallikericr01"]),
		###merge across sequencing lanes, let fail for files with only 
		#merged = expand(config["output"]+"alignments/sorted/merged/{samplename}_{species}-pf3d7_sorted_merged.bam", samplename = samples_to_merge, species = ["curtisigh01","wallikericr01"]),
		#unmerged = expand(config["output"]+"alignments/sorted/merged/{samplename}_{species}-pf3d7_sorted_merged.bam", samplename = samples_no_merge, species = ["curtisigh01","wallikericr01"]),
		totalmerge = expand(config["output"]+"merged_alignments/{samplename}_{species}-pf3d7.bam", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		###readgrouping
		#readgrouped = expand(config["output"]+"merged_alignments/readgrouped/{samplename}_{species}-pf3d7.bam", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		###Deduplication
		#dedupped = expand(config["output"]+"merged_alignments/readgrouped/dedupped/{samplename}_{species}-pf3d7.bam", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		###Ovale selection
		ovaleselected = expand(config["output"]+"ovale_alignments/{samplename}_{species}.bam", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		###alignment cleaning
		cleaned = expand(config["output"]+"ovale_alignments/cleaned/{samplename}_{species}.bam", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		###Alignment statistics
		cov = expand(config["output"]+"statistics_alignments/genomecov/{samplename}_{species}_genomecov{depth}.txt", samplename = samplenames, species = ["curtisigh01","wallikericr01"], depth = ["1","5","10"]),
		flagstats = expand(config["output"]+"statistics_alignments/flagstats/{samplename}_{species}_flagstats.txt", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		coverage = expand(config["output"]+"statistics_alignments/coverage/{samplename}_{species}_coverage.txt", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		overall_coverage = expand(config["output"]+"statistics_alignments/overall_coverage/{samplename}_{species}_coverage-at-{depth}.txt", samplename = samplenames, species = ["curtisigh01","wallikericr01"], depth = ["1","5","10"]),
		###Variant calling
		gvcf = expand(config["output"]+"gvcfs/{samplename}_{species}.g.vcf.gz", samplename = samplenames, species = ["curtisigh01","wallikericr01"]),
		db_species = expand(directory(config["output"]+"gvcfs/grouped/{species}"), species = ["curtisigh01","wallikericr01"]),
		db_speciesandmixed = expand(directory(config["output"]+"gvcfs/grouped/{species}andmixed"), species = ["curtisigh01","wallikericr01"]),
		###administrative documents showing the pipeline and config inputs
		config = config["output"]+"pipeline/config.yaml",
		snakefile = config["output"]+"pipeline/Snakefile_processing.py"

###Transfers a copy of the config file and Snakefile to the output for future reference of those results
rule copy_snakefileandconfig:
	input:
		config = "config/config.yaml",
		snakefile = "workflow/Snakefile_processing.py"
	output:
		config = config["output"]+"pipeline/config.yaml",
		snakefile = config["output"]+"pipeline/Snakefile_processing.py"
	shell:
		"""cp {input.config} {output.config}
		cp {input.snakefile} {output.snakefile}"""

###prepares concatenated reference genomes of curtisi and wallikeri
rule dual_reference:
	input:
		ovale = config["input_genomes"]+"{ovale}01.fasta",
		pf3d7 = config["input_genomes"]+"pf3d7.fasta"
	output:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta"
	shell:
		"cat {input.ovale} {input.pf3d7} > {output}"

rule samtools_index:
	input:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta"
	output:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta.fai"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools faidx {input}"

rule bwa_index:
	input:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta"
	output:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta.0123"
	conda:
		"envs/bwa-mem2.yaml"
	shell:
		"bwa-mem2 index {input}"

rule gatk_index:
	input:
		config["input_genomes"]+"{ovale}01-pf3d7.fasta"
	output:
		config["input_genomes"]+"{ovale}01-pf3d7.dict"
	conda:
		"envs/gatk.yaml"
	shell:
		"gatk CreateSequenceDictionary -R {input}"
		
###process sample reads
rule fastqc_raw:
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
	input:
		R1 = config["input_reads"]+"{samplename}_{lane}_R1.fastq.gz",
		R2 = config["input_reads"]+"{samplename}_{lane}_R2.fastq.gz"
	output:
		P1 = config["output"]+"trimmed_reads/{samplename}_{lane}_1P.fastq.gz",
		#U1 = config["output"]+"trimmed_reads/{samplename}_{lane}_1U.fastq.gz",
		#P2 = config["output"]+"trimmed_reads/{samplename}_{lane}_2P.fastq.gz",
		#U2 = config["output"]+"trimmed_reads/{samplename}_{lane}_2U.fastq.gz",
	params:
		outfile = config["output"]+"trimmed_reads/{samplename}_{lane}.fastq.gz"
	conda:
		"envs/trimmomatic.yaml"
	shell:
		"trimmomatic PE -threads 8 {input.R1} {input.R2} -baseout {params.outfile} ILLUMINACLIP:/nas/longleaf/apps/trimmomatic/0.36/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads"

rule fastqc_trimmed:
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
	input:
		P1 = config["output"]+"trimmed_reads/{samplename}_{lane}_1P.fastq.gz",
		P2 = config["output"]+"trimmed_reads/{samplename}_{lane}_2P.fastq.gz",
		reference = config["input_genomes"]+"{dual}.fasta"
	output:
		config["output"]+"alignments/{samplename}_{lane}_{dual}.bam"
	conda:
		"envs/bwa-mem2.yaml"
	resources:
		mem_mb = 150000
	shell:
		"bwa-mem2 mem -t 8 {input.reference} {input.P1} {input.P2} -o {output}"
		
rule sort:
	input:
		config["output"]+"alignments/{samplename}_{lane}_{species}.bam"
	output:
		config["output"]+"alignments/sorted/{samplename}_{lane}_{species}_sorted.bam"
	conda:
		"envs/picard.yaml"
	resources:
		mem_mb = 150000
	shell:
		"picard SortSam I={input} O={output} sort_order=coordinate create_index=true"

rule merge:
	input:
		lane1 = config["output"]+"alignments/sorted/{samplename}_L001_{species}-pf3d7_sorted.bam"
	params:
		lane2 = config["output"]+"alignments/sorted/{samplename}_L002_{species}-pf3d7_sorted.bam"
	output:
		config["output"]+"merged_alignments/{samplename}_{species}-pf3d7.bam"
	conda:
		"envs/samtools.yaml"
	resources:
		mem_mb = 150000
	script:
		"scripts/merger.py"


# rule merge:
	# input:
		# L001 = expand(config["output"]+"alignments/sorted/{samplename}_L001_{species}-pf3d7_sorted.bam", samplename = samples_to_merge, allow_missing=True),#, species = ["curtisigh01","wallikericr01"]),
		# L002 = expand(config["output"]+"alignments/sorted/{samplename}_L002_{species}-pf3d7_sorted.bam", samplename = samples_to_merge, allow_missing = True)#, species = ["curtisigh01","wallikericr01"])
	# output:
		# expand(config["output"]+"alignments/sorted/merged/{samplename}_{species}-pf3d7_sorted_merged.bam", samplename = samples_to_merge, allow_missing=True)#, species = ["curtisigh01","wallikericr01"])
	# params:
		# index = len(samples_to_merge)
	# conda:
		# "envs/samtools.yaml"
	# resources:
		# mem_mb = 150000
	# script:
		# "scripts/merger.py"
	# #shell:
	# #	"for i in {{0..{params.index}}}; do samtools merge -o {output[${{i}}]} {input.L001[${{i}}]} {input.L002[${{i}}]}; done"
	# #run:
	# #	for i in range(0,params.index):
	# #		subprocess.call(["samtools","merge","-o",output[i],input.L001[i],input.L002[i]])
# rule nomerge:
	# input:
		# expand(config["output"]+"alignments/sorted/{samplename}_L001_{species}-pf3d7_sorted.bam", samplename = samples_no_merge, allow_missing=True)#, species = ["curtisigh01","wallikericr01"])
	# output:
		# expand(config["output"]+"alignments/sorted/merged/{samplename}_{species}-pf3d7_sorted_merged.bam", samplename = samples_no_merge, allow_missing=True)#, species = ["curtisigh01","wallikericr01"])
	# params:
		# index = len(samples_no_merge)
	# resources:
		# mem_mb = 150000
	# #shell:
	# #	"for i in ((1..(params.index}}}; do cp {input} {output}; done"
	# script:
		# "scripts/nomerger.py"
# 
# rule clean_merge_names:
	# input:
		# config["output"]+"alignments/sorted/merged/{samplename}_{dual}_sorted_merged.bam"
	# output:
		# config["output"]+"merged_alignments/{samplename}_{dual}.bam"
	# shell:
		# "cp {input} {output}"

#rule clean_unmerge_names:
#	input:
#		config["output"]+"alignments/sorted/merged/{samplename}_{dual}_sorted_unmerged.bam"
#	output:
#		config["output"]+"merged_alignments/{samplename}_{dual}.bam"
#	shell:
#		"cp {input} {output}"

rule add_readgroup:
	input:
		config["output"]+"merged_alignments/{samplename}_{dual}.bam"
	output:
		config["output"]+"merged_alignments/readgrouped/{samplename}_{dual}.bam"
	params:
		RGID = lambda wildcards, output: list(readgroups[wildcards.samplename])[0],
		RGLB = lambda wildcards, output: list(readgroups[wildcards.samplename])[1],
		RGPU = lambda wildcards, output: list(readgroups[wildcards.samplename])[2],
		RGPL = lambda wildcards, output: list(readgroups[wildcards.samplename])[3],
		RGSM = lambda wildcards, output: list(readgroups[wildcards.samplename])[4]
	conda:
		"envs/picard.yaml"
	resources:
		mem_mb = 300000
	shell:
		"picard AddOrReplaceReadGroups I={input} O={output} SORT_ORDER=coordinate RGID={params.RGID} RGLB={params.RGLB} RGPL={params.RGPL} RGPU={params.RGPU} RGSM={params.RGSM}"
		
rule deduplicate:
	input:
		config["output"]+"merged_alignments/readgrouped/{samplename}_{dual}.bam"
	output:
		out = config["output"]+"merged_alignments/readgrouped/dedupped/{samplename}_{dual}.bam",
		metrics = config["output"]+"merged_alignments/readgrouped/dedupped/metrics/{samplename}_{dual}_dupmetrics.txt"
	conda:
		"envs/gatk.yaml"
	resources:
		mem_mb = 300000
	shell:
		"gatk MarkDuplicatesSpark -I {input} -O {output.out} -M {output.metrics}"

rule ovale_selecter:
	input:
		bam = config["output"]+"merged_alignments/readgrouped/{samplename}_{species}-pf3d7.bam",
		ref =  config["input_genomes"]+"{species}.fasta",
		bed = config["input_beds"]+"{species}.bed"
	output:
		config["output"]+"ovale_alignments/{samplename}_{species}.bam"
	conda:
		"envs/samtools.yaml"
	resources:
		mem_mb = 300000
	shell:
		"samtools view -b -h {input.bam} -T {input.ref} -L {input.bed} > {output}"

rule sam_cleaner:
	input:
		config["output"]+"ovale_alignments/{samplename}_{species}.bam"
	output:
		config["output"]+"ovale_alignments/cleaned/{samplename}_{species}.bam"
	conda:
		"envs/gatk.yaml"
	resources:
			mem_mb = 300000
	shell:
		"gatk CleanSam -I {input} -O {output}"

rule genome_coverage:
	input:
		config["output"]+"ovale_alignments/cleaned/{samplename}_{species}.bam"
	output:
		config["output"]+"statistics_alignments/genomecov/{samplename}_{species}_genomecov{depth}.txt"
	params:
		depth = lambda wildcards, output: wildcards.depth
	conda:
		"envs/bedtools.yaml"
	shell:
		"bedtools genomecov -max {params.depth} -ibam {input} > {output}"

rule coverage_calc:
	input:
		cov = config["output"]+"statistics_alignments/genomecov/{samplename}_{species}_genomecov{depth}.txt",
		bed = config["input_beds"]+"{species}_chr.bed"
	output:
		config["output"]+"statistics_alignments/overall_coverage/{samplename}_{species}_coverage-at-{depth}.txt"
	script:
		"scripts/coverager.py"


rule flagstats:
	input:
		config["output"]+"ovale_alignments/cleaned/{samplename}_{species}.bam"
	output:
		config["output"]+"statistics_alignments/flagstats/{samplename}_{species}_flagstats.txt"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools flagstat {input} > {output}"

rule coveragestats:
	input:
		config["output"]+"ovale_alignments/cleaned/{samplename}_{species}.bam"
	output:
		config["output"]+"statistics_alignments/coverage/{samplename}_{species}_coverage.txt"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools coverage {input} > {output}"

###Variant Calling
rule index_bams:
	input:
		config["output"]+"ovale_alignments/cleaned/{samplename}_{species}.bam"
	output:
		config["output"]+"ovale_alignments/cleaned/{samplename}_{species}.bam.bai"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools index {input}"
rule haplotypecaller:
	input:
		bam = config["output"]+"ovale_alignments/cleaned/{samplename}_{species}.bam",
		index = config["output"]+"ovale_alignments/cleaned/{samplename}_{species}.bam.bai",
		ref = config["input_genomes"]+"{species}.fasta",
		bed = config["input_beds"]+"{species}.bed"
	output:
		config["output"]+"gvcfs/{samplename}_{species}.g.vcf.gz"
	conda:
		"envs/gatk.yaml"
	resources:
		mem_mb = 300000
	shell:
		"gatk HaplotypeCaller -R {input.ref} -L {input.bed} -I {input.bam} -O {output} -ERC GVCF"

rule dbimport_species:
	input:
		map = config["input_lists"]+"dbimport_{species}_samplemap.txt",
		gvcfs = expand(config["output"]+"gvcfs/{samplename}_{species}.g.vcf.gz", samplename = samplenames, allow_missing = True),
		bed = config["input_beds"]+"{species}.bed"
	output:
		directory(config["output"]+"gvcfs/grouped/{species}")
	params:
		dir = config["output"]+"gvcfs/grouped/{species}
	conda:
		"envs/gatk.yaml"
	shell:
		"gatk GenomicsDBImport --sample-name-map {input.map} --genomicsdb-workspace-path {params.dir} -L {input.bed}"


rule dbimport_speciesandmixed:
	input:
		map = config["input_lists"]+"dbimport_{species}andmixed_samplemap.txt",
		gvcfs = expand(config["output"]+"gvcfs/{samplename}_{species}.g.vcf.gz", samplename = samplenames, allow_missing = True),
		bed = config["input_beds"]+"{species}.bed"
	output:
		directory(config["output"]+"gvcfs/grouped/{species}andmixed")
	params:
		dir = config["output"]+"gvcfs/grouped/{species}andmixed
	conda:
		"envs/gatk.yaml"
	shell:
		"gatk GenomicsDBImport --sample-name-map {input.map} --genomicsdb-workspace-path {params.dir} -L {input.bed}"

