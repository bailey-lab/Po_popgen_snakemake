import yaml
import subprocess

with open("config/config.yaml") as file:
        config = yaml.safe_load(file)

#read in input files
species = snakemake.wildcards.species
genes = open(config["output"]+ config[species+"_genes"]).readlines()
exons = open(config["output"]+ config[species+"_exons"]).readlines()
project = snakemake.wildcards.project
mixed = snakemake.wildcards.mixed
window = snakemake.params.window
stats = open(config["output"]+ "selection/"+project+"_"+species+"_"+mixed+"-"+window+".Tajima.D").readlines()

#create output files
tajima_genes = open(config["output"]+"selection/"+project+"_"+species+"_"+mixed+"_tajimad-"+window+"_genes.bed", "w")
tajima_exons = open(config["output"]+"selection/"+project+"_"+species+"_"+mixed+"_tajimad-"+window+"_exons.bed", "w")
#right header line to each output
tajima_genes.writelines(["#",stats[0].split("\t")[0],"\t",stats[0].split("\t")[1],"\t",stats[0].split("\t")[2],"\t",stats[0].split("\t")[5]])
tajima_exons.writelines(["#",stats[0].split("\t")[0],"\t",stats[0].split("\t")[1],"\t",stats[0].split("\t")[2],"\t",stats[0].split("\t")[5]])

#for each line of the input Tajima D, parse out the pieces and confirmed that they fall within genes or exons and are not missing
for i in range(1,(len(stats)-1)):
	#define chrom, start position, number of variants, and tajima's D for the line being read in
	chrom = stats[i].split("\t")[0]
	start = int(stats[i].split("\t")[1])
	stop = int(stats[i].split("\t")[2])
	nsites = int(stats[i].split("\t")[3])
	nsnps = int(stats[i].split("\t")[4])
	tajima = stats[i].split("\t")[5].rstrip("\n")
	geneflag = 0
	exonflag = 0
	window = int(snakemake.params.window)
	#If tajima is defined for this window, we will perform two for loops to check whether it falls within a gene or an exon
	if tajima.lstrip("-").replace(".","").rstrip("\n").isnumeric():
		#iterate through all genes in the gene bed file
		for j in range(0,(len(genes)-1)):
			#if a given gene matches the chromosome of the tajima calculation, we'll check the interval
			if chrom == genes[j].split("\t")[0]:
				print(genes[j].split("\t")[0])
				#if the tajima window starts after this gene window, we'll continue considering
				if int(start) >= int(genes[j].split("\t")[1]):
					#if the end of the tajima window (stop) is before the end of the gene, we'll use this tajima window
					if (int(stop)) <= int(genes[j].split("\t")[2].rstrip("\n")):
						#only if the tajima window matches chromosome and falls within the start and stop of the gene will we keep it
						geneflag = 1
		#iterate through all genes in the gene bed file
		for k in range(0,(len(exons)-1)):
			#if a given exon matches the chromosome of the tajima calculation, we'll check the interval
			if chrom == exons[k].split("\t")[0]:
				#if the tajima window starts after this exon window, we'll continue considering
				if int(start) >= int(exons[k].split("\t")[1]):
					#if the end of the tajima window (start + window length) is before the end of the exon, we'll use this tajima window
					if (int(stop)) <= int(exons[k].split("\t")[2].rstrip("\n")):
						exonflag = 1
	#tajima windows will be written to gene list if the window falls within one gene window in the genome
	if geneflag == 1:
		tajima_genes.writelines([str(chrom),"\t",str(start),"\t",str(stop),"\t",str(tajima),"\n"])
	#tajima windows will be written to exon list if the window falls within one exon window in the genome
	if exonflag == 1:
		tajima_exons.writelines([str(chrom),"\t",str(start),"\t",str(stop),"\t",str(tajima),"\n"])
		
			
		
