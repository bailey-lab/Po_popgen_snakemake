import yaml
import subprocess

with open("config/config.yaml") as file:
        config = yaml.safe_load(file)

#read in input files
species = snakemake.wildcards.species
#species = "curtisigh01"
print(species)
genes = open(config["input_beds"]+ config[species+"_genes"]).readlines()
#print(genes)
exons = open(config["input_beds"]+ config[species+"_exons"]).readlines()
project = snakemake.wildcards.project
#project = "ov1"
print(project)
stats = open(config["output"]+ "selection/"+project+"_"+species+".Tajima.D").readlines()

#create output files
tajima_genes = open(config["output"]+"selection/"+project+"_"+species+"_tajimad_genes.txt", "w")
tajima_exons = open(config["output"]+"selection/"+project+"_"+species+"_tajimad_exons.txt", "w")
#right header line to each output
tajima_genes.writelines(stats[0])
tajima_exons.writelines(stats[0])

#for each line of the input Tajima D, parse out the pieces and confirmed that they fall within genes or exons and are not missing
for i in range(1,(len(stats)-1)):
	chrom = stats[i].split("\t")[0]
	#print(chrom)
	start = stats[i].split("\t")[1]
	#print(start)
	n = stats[i].split("\t")[2]
	tajima = stats[i].split("\t")[3].rstrip("\n")
	#print(tajima.lstrip("-").replace(".",""))
	#print(tajima.lstrip("-").replace(".","").rstrip("\n").isnumeric())
	geneflag = 0
	exonflag = 0
	window = snakemake.params.window
	#window = 50
	if tajima.lstrip("-").replace(".","").rstrip("\n").isnumeric():
		#print(tajima)
		for j in range(1,(len(genes)-1)):
			#print(genes[j].split("\t")[1])
			if chrom == genes[j].split("\t")[1]:
				if int(start) >= int(genes[j].split("\t")[2]):
					#print(genes[j].split("\t")[3].rstrip("\n"))
					if (int(start)+int(window)) <= int(genes[j].split("\t")[3].rstrip("\n")):
						geneflag = 1
						print(chrom)
						print(tajima)
		for k in range(1,(len(exons)-1)):
			if chrom == exons[k].split("\t")[1]:
				if int(start) >= int(exons[k].split("\t")[2]):
					if (int(start)+int(window)) <= int(exons[k].split("\t")[3].rstrip("\n")):
						exonflag = 1
	if geneflag == 1:
		tajima_genes.writelines(stats[i])
	if exonflag == 1:
		tajima_exons.writelines(stats[i])
		
			
		
