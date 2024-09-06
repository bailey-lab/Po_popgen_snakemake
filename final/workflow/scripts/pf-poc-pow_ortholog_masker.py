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
all_orthos = open(config["input_beds"]+config["all_triple_orthos"]).readlines()[1:]
#for i in range(0,(len(flines)-1)):
#	all_orthos.append(clines[i].split("\t")[0]+"\t"+clines[i].split("\t")[1]+"\t"+clines[i].split("\t")[2]+"\t"+clines[i].split("\t")[3]+"\t"+flines[i].split("\t")[0]+"\t"+flines[i].split("\t")[1]+"\t"+flines[i].split("\t")[2]+"\t"+flines[i].split("\t")[3])

poc_chr = open(config["input_beds"]+"curtisigh01_chr.bed").readlines()
poc_gene_mask = open(config["input_beds"]+"curtisigh01_genemask.bed").readlines()
poc_cov_mask =  open(snakemake.input.poc_coverage_mask).readlines()
pow_chr = open(config["input_beds"]+"wallikericr01_chr.bed").readlines()
pow_gene_mask = open(config["input_beds"]+"wallikericr01_genemask.bed").readlines()
pow_cov_mask = open(snakemake.input.pow_coverage_mask).readlines()
pf_core = open(config["input_beds"]+"pf3d7_core.bed").readlines()
pf_cov_mask = open(snakemake.input.pf_coverage_mask).readlines()

for i in all_orthos:
	###Poc masking
	#Remove any orthologs groups where Poc ortholog falls outside of Poc chromosome or overlaps with masked region
	poc_chrom_list = []
	pocchr_flag = 1
	for chrom in poc_chr:
		if i.split("\t")[6] == chrom.split("\t")[0]:
			pocchr_flag = 0
	pocgene_flag = 0
	for gene in poc_gene_mask:
		#if the ortholog and hypervariable gene are on the same chromosome, we'll check if they overlap using an if...elif...elif statement
		if gene.split("\t")[0] == i.split("\t")[6]:
			#if start of ortholog is within the hypervariable gene window, we'll drop
			if (int(gene.split("\t")[1]) < int(i.split("\t")[7]) and int(i.split("\t")[7]) <= int(gene.split("\t")[2])):
				pocgene_flag = 1
			#if end of ortholog is within the hypervariable gene window, we'll drop
			elif (int(gene.split("\t")[1]) < int(i.split("\t")[8]) and int(i.split("\t")[8]) <= int(gene.split("\t")[2])):
				pocgene_flag = 1
			#if start of ortholog precedes the hypervariable gene and the end comes after, we'll drop
			elif (int(i.split("\t")[7]) < int(gene.split("\t")[1]) and int(gene.split("\t")[2]) < int(i.split("\t")[8])):
				pocgene_flag = 1
	###limiting to poc regions with sufficient coverage
	poc_cov_flag = 1
	for interval in poc_cov_mask:
		if interval.split("\t")[0] == i.split("\t")[6]:
			#If the poc ortholog in question falls entirely within one of the sufficient coverage windows, we'll keep
			if (int(interval.split("\t")[1]) < int(i.split("\t")[7]) and int(i.split("\t")[8]) < int(interval.split("\t")[2])):
				poc_cov_flag = 0
	###Pow masking
	#Remove any orthologs groups where Pow ortholog falls outside of Pow chromosome or overlaps with masked region
	pow_chrom_list = []
	powchr_flag = 1
	for chrom in pow_chr:
		if i.split("\t")[10] == chrom.split("\t")[0]:
			powchr_flag = 0
	powgene_flag = 0
	for gene in pow_gene_mask:
		#if the ortholog and hypervariable gene are on the same chromosome, we'll check if they overlap using an if...elif...elif statement
		if gene.split("\t")[0] == i.split("\t")[10]:
			#if start of ortholog is within the hypervariable gene window, we'll drop
			if (int(gene.split("\t")[1]) < int(i.split("\t")[11]) and int(i.split("\t")[11]) <= int(gene.split("\t")[2])):
				powgene_flag = 1
			#if end of ortholog is within the hypervariable gene window, we'll drop
			elif (int(gene.split("\t")[1]) < int(i.split("\t")[12]) and int(i.split("\t")[12]) <= int(gene.split("\t")[2])):
				powgene_flag = 1
			#if start of ortholog precedes the hypervariable gene and the end comes after, we'll drop
			elif (int(i.split("\t")[11]) < int(gene.split("\t")[1]) and int(gene.split("\t")[2]) < int(i.split("\t")[12])):
				powgene_flag = 1
	###limiting to pow regions with sufficient coverage
	pow_cov_flag = 1
	for interval in pow_cov_mask:
		if interval.split("\t")[0] == i.split("\t")[10]:
			#If the pow ortholog in question falls entirely within one of the sufficient coverage windows, we'll keep
			if (int(interval.split("\t")[1]) < int(i.split("\t")[11]) and int(i.split("\t")[12]) < int(interval.split("\t")[2])):
				pow_cov_flag = 0
	#to determine that the ortholog falls within the core of the pf genome, we will use a separate flag that excludes a given ortholog unless we show it is in the pf core genome
	pf_flag = 1
	for core in pf_core:
		if core.split("\t")[0] == i.split("\t")[2]:
			#If the pf ortholog in question falls entirely within one of the core windows, we'll keep
			if (int(core.split("\t")[1]) < int(i.split("\t")[3]) and int(i.split("\t")[4]) < int(core.split("\t")[2])):
				pf_flag = 0
	pf_cov_flag = 1
	for interval in pf_cov_mask:
		if interval.split("\t")[0] == i.split("\t")[2]:
			#If the pf ortholog in question falls entirely within one of the sufficient coverage windows, we'll keep
			if (int(interval.split("\t")[1]) < int(i.split("\t")[3]) and int(i.split("\t")[4]) < int(interval.split("\t")[2])):
				pf_cov_flag = 0
	#### add cov filter flags
	if (pf_flag == 0 and pf_cov_flag == 0 and pocchr_flag == 0 and pocgene_flag == 0 and poc_cov_flag == 0 and powchr_flag == 0 and powgene_flag == 0 and pow_cov_flag == 0):
		masked_orthos.append(i)
	if (pf_flag == 0 and pocchr_flag == 0 and pocgene_flag == 0 and powchr_flag == 0 and powgene_flag == 0):
		unfiltered_orthos.append(i)
	

poc_masked_orthos = []
pow_masked_orthos = []
pf_masked_orthos = []
				
for line in masked_orthos:
	pf_masked_orthos.append(line.split("\t")[1]+"\t"+line.split("\t")[2]+"\t"+line.split("\t")[3]+"\t"+line.split("\t")[4]+"\t"+" ")
	poc_masked_orthos.append(line.split("\t")[5]+"\t"+line.split("\t")[6]+"\t"+line.split("\t")[7]+"\t"+line.split("\t")[8]+"\t"+" ")
	pow_masked_orthos.append(line.split("\t")[9]+"\t"+line.split("\t")[10]+"\t"+line.split("\t")[11]+"\t"+line.split("\t")[12]+"\t"+" ")


pfmo = open(config["output"]+config["pf_triple_orthos_masked"], "w")
for j in pf_masked_orthos:
	pfmo.writelines(j+"\n")

pocmo = open(config["output"]+config["poc_triple_orthos_masked"], "w")
for k in poc_masked_orthos:
	pocmo.writelines(k+"\n")

powmo = open(config["output"]+config["pow_triple_orthos_masked"], "w")
for l in pow_masked_orthos:
	powmo.writelines(l+"\n")



poc_unfiltered_orthos = []
pow_unfiltered_orthos = []
pf_unfiltered_orthos = []


for line in unfiltered_orthos:
	pf_unfiltered_orthos.append(line.split("\t")[1]+"\t"+line.split("\t")[2]+"\t"+line.split("\t")[3]+"\t"+line.split("\t")[4]+"\t"+" ")
	poc_unfiltered_orthos.append(line.split("\t")[5]+"\t"+line.split("\t")[6]+"\t"+line.split("\t")[7]+"\t"+line.split("\t")[8]+"\t"+" ")
	pow_unfiltered_orthos.append(line.split("\t")[9]+"\t"+line.split("\t")[10]+"\t"+line.split("\t")[11]+"\t"+line.split("\t")[12]+"\t"+" ")

pfuo = open(config["output"]+config["pf_triple_orthos_unfiltered"], "w")
for j in pf_unfiltered_orthos:
	pfuo.writelines(j+"\n")

pocuo = open(config["output"]+config["poc_triple_orthos_unfiltered"], "w")
for k in poc_unfiltered_orthos:
	pocuo.writelines(k+"\n")

powuo = open(config["output"]+config["pow_triple_orthos_unfiltered"], "w")
for l in pow_unfiltered_orthos:
	powuo.writelines(l+"\n")
