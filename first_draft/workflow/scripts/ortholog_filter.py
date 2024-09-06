import yaml
import subprocess

with open("config/config.yaml") as file:
	config = yaml.safe_load(file)
	
#print(config[0])
#print(config["input_beds"])

def parse_bed_file(bed_file):
	gene_dict={}
	for line in open(bed_file):
		name = line.split("\t")[0]
		chrom = line.split("\t")[1]
		start = line.split("\t")[2]
		stop = line.split("\t")[3]
		gene_dict[name] =  [chrom, start, stop]
		#print(name)
		#print(gene_dict[name])
	return gene_dict

pf3d7_orthos = parse_bed_file(config["input_beds"]+ "pf3d7_poc-pf-orthologs.tsv")
curtisigh01_orthos = parse_bed_file(config["input_beds"]+ "pocgh01_poc-pf-orthologs.tsv")


#pfinput = input.pfvcf
print(snakemake.input.pfvcf)
#pfoutfile = params.pfoutfile
print(snakemake.params.pfoutfile)
pfgeneid = snakemake.wildcards.pfgeneid
#print(pfgeneid)
chrom = list(pf3d7_orthos[pfgeneid])[0]
#print(chrom)
#print(list(pf3d7_orthos[pfgeneid])[0])
start = list(pf3d7_orthos[pfgeneid])[1]
#print(start)
stop = list(pf3d7_orthos[pfgeneid])[2]
#print(stop)

subprocess.call(['vcftools', '--gzvcf', snakemake.input.pfvcf, '--out', snakemake.params.pfoutfile, '--chr', chrom, '--from-bp', start, '--to-bp', stop, '--recode', '--recode-INFO-all'])
