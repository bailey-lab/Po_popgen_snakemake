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
#	update wallikeru reference genome for analysis pipelin
#  add in joste pow samples
# fix profile slurm; where is the profile file?
#needed input files:


configfile: "config/config.yaml"

samplenames = []
for i in open(config["input"]+"lists/sample_names.txt").readlines():
	samplenames.append(i.rstrip("\n"))


samplefiles = []
for f in open(config["input"]+"lists/sample_files.txt").readlines():
	samplefiles.append(f.rstrip("\n"))

print(samplenames)
print(samplefiles)
###Determines the final files to be output by the snakemake pipeline
rule all:
#Previous rules should use wildcards for the project and species, but this final rule should employ the actual names and terms for the final file to be created
	input:
		fastqc = expand(config["output"]+"fastqc/untrimmed/{samplefile}.fastqc", samplefile = samplefiles),
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

rule fastqc_raw:
	input:
		config["input"]+"reads/{samplefile}.fastq.gz"
	output:
		config["output"]+"fastqc/untrimmed/{samplefile}.fastqc"
	conda:
		"envs/fastqc.yaml"
	shell:
		"fastqc {input} -o {output}"
		
