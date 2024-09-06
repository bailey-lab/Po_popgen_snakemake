import yaml
import subprocess

with open("config/config.yaml") as file:
        config = yaml.safe_load(file)

all_orthos = []
masked_orthos = []
unfiltered_orthos = []
#f = open(config["input_beds"]+ config["pf_orthos"])
#c = open(config["input_beds"]+ config["poc_orthos"])
#flines = f.readlines()
#clines = c.readlines()
#print(flines)

####Replace with single file of all orthologs
all_orthos = open(config["input_beds"]+config["all_ovale_orthos"]).readlines()[1:]
#for i in range(0,(len(flines)-1)):
#	all_orthos.append(clines[i].split("\t")[0]+"\t"+clines[i].split("\t")[1]+"\t"+clines[i].split("\t")[2]+"\t"+clines[i].split("\t")[3]+"\t"+flines[i].split("\t")[0]+"\t"+flines[i].split("\t")[1]+"\t"+flines[i].split("\t")[2]+"\t"+flines[i].split("\t")[3])

poc_chr = open(config["input_beds"]+"curtisigh01_chr.bed").readlines()
poc_gene_mask = open(config["input_beds"]+"curtisigh01_genemask.bed").readlines()
poc_cov_mask =  open(snakemake.input.poc_coverage_mask).readlines()
pow_chr = open(config["input_beds"]+"wallikericr01_chr.bed").readlines()
pow_gene_mask = open(config["input_beds"]+"wallikericr01_genemask.bed").readlines()
pow_cov_mask = open(snakemake.input.pow_coverage_mask).readlines()
for i in all_orthos:
	###Poc masking
	#Remove any orthologs groups where Poc ortholog falls outside of Poc chromosome or overlaps with masked region
	poc_chrom_list = []
	pocchr_flag = 1
	for chrom in poc_chr:
		if i.split("\t")[2] == chrom.split("\t")[0]:
			pocchr_flag = 0
	pocgene_flag = 0
	for gene in poc_gene_mask:
		#if the ortholog and hypervariable gene are on the same chromosome, we'll check if they overlap using an if...elif...elif statement
		if gene.split("\t")[0] == i.split("\t")[2]:
			#if start of ortholog is within the hypervariable gene window, we'll drop
			if (int(gene.split("\t")[1]) < int(i.split("\t")[3]) and int(i.split("\t")[3]) <= int(gene.split("\t")[2])):
				pocgene_flag = 1
			#if end of ortholog is within the hypervariable gene window, we'll drop
			elif (int(gene.split("\t")[1]) < int(i.split("\t")[4]) and int(i.split("\t")[4]) <= int(gene.split("\t")[2])):
				pocgene_flag = 1
			#if start of ortholog precedes the hypervariable gene and the end comes after, we'll drop
			elif (int(i.split("\t")[3]) < int(gene.split("\t")[1]) and int(gene.split("\t")[2]) < int(i.split("\t")[4])):
				pocgene_flag = 1
	###limiting to poc regions with sufficient coverage
	poc_cov_flag = 1
	for interval in poc_cov_mask:		
		if interval.split("\t")[0] == i.split("\t")[2]:
			#If the poc ortholog in question falls entirely within one of the sufficient coverage windows, we'll keep
			if (int(interval.split("\t")[1]) < int(i.split("\t")[3]) and int(i.split("\t")[4]) < int(interval.split("\t")[2])):
				poc_cov_flag = 0
	###Pow masking
	#Remove any orthologs groups where Pow ortholog falls outside of Pow chromosome or overlaps with masked region
	pow_chrom_list = []
	powchr_flag = 1
	for chrom in pow_chr:
		if i.split("\t")[6] == chrom.split("\t")[0]:
			powchr_flag = 0
	powgene_flag = 0
	for gene in pow_gene_mask:
		#if the ortholog and hypervariable gene are on the same chromosome, we'll check if they overlap using an if...elif...elif statement
		if gene.split("\t")[0] == i.split("\t")[6]:
			#if start of ortholog is within the hypervariable gene window, we'll drop
			if (int(gene.split("\t")[1]) < int(i.split("\t")[7]) and int(i.split("\t")[7]) <= int(gene.split("\t")[2])):
				powgene_flag = 1
			#if end of ortholog is within the hypervariable gene window, we'll drop
			elif (int(gene.split("\t")[1]) < int(i.split("\t")[8]) and int(i.split("\t")[8]) <= int(gene.split("\t")[2])):
				powgene_flag = 1
			#if start of ortholog precedes the hypervariable gene and the end comes after, we'll drop
			elif (int(i.split("\t")[7]) < int(gene.split("\t")[1]) and int(gene.split("\t")[2]) < int(i.split("\t")[8])):
				powgene_flag = 1
	###limiting to pow regions with sufficient coverage
	pow_cov_flag = 1
	for interval in pow_cov_mask:		
		if interval.split("\t")[0] == i.split("\t")[6]:
			#If the pow ortholog in question falls entirely within one of the sufficient coverage windows, we'll keep
			if (int(interval.split("\t")[1]) < int(i.split("\t")[7]) and int(i.split("\t")[8]) < int(interval.split("\t")[2])):
				pow_cov_flag = 0
	#### add cov filter flags
	if (pocchr_flag == 0 and pocgene_flag == 0 and poc_cov_flag == 0 and powchr_flag == 0 and powgene_flag == 0 and pow_cov_flag == 0):
		masked_orthos.append(i)
	if (pocchr_flag == 0 and pocgene_flag == 0 and powchr_flag == 0 and powgene_flag == 0):
		unfiltered_orthos.append(i)
	

poc_masked_orthos = []
pow_masked_orthos = []
				
for line in masked_orthos:
	poc_masked_orthos.append(line.split("\t")[1]+"\t"+line.split("\t")[2]+"\t"+line.split("\t")[3]+"\t"+line.split("\t")[4]+"\t"+" ")
	pow_masked_orthos.append(line.split("\t")[5]+"\t"+line.split("\t")[6]+"\t"+line.split("\t")[7]+"\t"+line.split("\t")[8].rstrip("\n")+"\t"+" ")

pocmo = open(config["output"]+config["poc_ovale_orthos_masked"], "w")
for k in poc_masked_orthos:
	pocmo.writelines(k+"\n")

powmo = open(config["output"]+config["pow_ovale_orthos_masked"], "w")
for l in pow_masked_orthos:
	powmo.writelines(l+"\n")



poc_unfiltered_orthos = []
pow_unfiltered_orthos = []


for line in unfiltered_orthos:
	poc_unfiltered_orthos.append(line.split("\t")[1]+"\t"+line.split("\t")[2]+"\t"+line.split("\t")[3]+"\t"+line.split("\t")[4]+"\t"+" ")
	pow_unfiltered_orthos.append(line.split("\t")[5]+"\t"+line.split("\t")[6]+"\t"+line.split("\t")[7]+"\t"+line.split("\t")[8]+"\t"+" ")

pocuo = open(config["output"]+config["poc_ovale_orthos_unfiltered"], "w")
for k in poc_unfiltered_orthos:
	pocuo.writelines(k+"\n")

powuo = open(config["output"]+config["pow_ovale_orthos_unfiltered"], "w")
for l in pow_unfiltered_orthos:
	powuo.writelines(l+"\n")
