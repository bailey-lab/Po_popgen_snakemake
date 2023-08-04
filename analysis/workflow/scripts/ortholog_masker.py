import yaml
import subprocess

with open("config/config.yaml") as file:
        config = yaml.safe_load(file)

all_orthos = []
masked_orthos = []
f = open(config["input_beds"]+ config["pf_orthos"])
c = open(config["input_beds"]+ config["poc_orthos"])
flines = f.readlines()
clines = c.readlines()
#print(flines)
for i in range(0,(len(flines)-1)):
	all_orthos.append(clines[i].split("\t")[0]+"\t"+clines[i].split("\t")[1]+"\t"+clines[i].split("\t")[2]+"\t"+clines[i].split("\t")[3]+"\t"+flines[i].split("\t")[0]+"\t"+flines[i].split("\t")[1]+"\t"+flines[i].split("\t")[2]+"\t"+flines[i].split("\t")[3])

poc_chr = open(config["input_beds"]+"curtisigh01_chr.bed").readlines()
poc_gene_mask = open(config["input_beds"]+"curtisigh01_genemask.bed").readlines()
pf_core = open(config["input_beds"]+"pf3d7_core.bed").readlines()

for i in all_orthos:
	poc_chrom_list = []
	pocchr_flag = 1
	for chrom in poc_chr:
		if i.split("\t")[1] == chrom.split("\t")[0]:
			pocchr_flag = 0
	pocgene_flag = 0
	for gene in poc_gene_mask:
		#if the ortholog and hypervariable gene are on the same chromosome, we'll check if they overlap using an if...elif...elif statement
		if gene.split("\t")[0] == i.split("\t")[1]:
			#if start of ortholog is within the hypervariavle gene window, we'll drop
			if gene.split("\t")[1] < i.split("\t")[2] <= gene.split("\t")[2]:
				pocgene_flag = 1
			#if end of ortholog is within the hypervariable gene window, we'll drop
			elif gene.split("\t")[1] < i.split("\t")[3] <= gene.split("\t")[2]:
				pocgene_flag = 1
			#if start of ortholog precedes the hypervariable gene and the end comes after, we'll drop
			elif ( i.split("\t")[2] < gene.split("\t")[1] and gene.split("\t")[2] < i.split("\t")[3]):
				pocgene_flag = 1
	#to determine that the ortholog falls within the core of the pf genome, we will use a separate flag that excludes a given ortholog unless we show it is in the pf core genome
	pf_flag = 1
	for core in pf_core:		
		if core.split("\t")[0] == i.split("\t")[5]:
			if i.split("\t")[6] > core.split("\t")[1] and i.split("\t")[7] < core.split("\t")[2]:
				pfflag = 0
	if pfflag == 0 and pocchr_flag == 0 and pocgene_flag == 0:
		masked_orthos.append(i)

poc_masked_orthos = []
pf_masked_orthos = []				
for line in masked_orthos:
	poc_masked_orthos.append(line.split("\t")[0]+"\t"+line.split("\t")[1]+"\t"+line.split("\t")[2]+"\t"+line.split("\t")[3]+"\t"+" ")
	pf_masked_orthos.append(line.split("\t")[4]+"\t"+line.split("\t")[5]+"\t"+line.split("\t")[6]+"\t"+line.split("\t")[7]+"\t"+" ")

pomo = open(config["input_beds"]+config["poc_orthos_masked"], "w")
for i in poc_masked_orthos:
	pomo.writelines(i+"\n")
pfmo = open(config["input_beds"]+config["pf_orthos_masked"], "w")
for j in pf_masked_orthos:
	pfmo.writelines(j+"\n")

