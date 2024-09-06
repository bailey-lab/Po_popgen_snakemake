import yaml
import subprocess

with open("config/config.yaml") as file:
        config = yaml.safe_load(file)

#load wildcard strings
species = snakemake.wildcards.species
project = snakemake.wildcards.project
mixed = snakemake.wildcards.mixed

#read in bed files as lists
genome_bed = open(snakemake.input.genome_bed).readlines()
gene_bed = open(snakemake.input.gene_bed).readlines()
exon_bed = open(snakemake.input.exon_bed).readlines()
intron_bed = open(snakemake.input.intron_bed).readlines()
intergenic_bed = open(snakemake.input.intergenic_bed).readlines()
cds_bed = open(snakemake.input.cds_bed).readlines()

#count number of variants in each table (-1 for the header line)
genome_snp_count = int(len(open(snakemake.input.genome_table).readlines())-1)
gene_snp_count = int(len(open(snakemake.input.gene_table).readlines())-1)
exon_snp_count = int(len(open(snakemake.input.exon_table).readlines())-1)
cds_snp_count = int(len(open(snakemake.input.cds_table).readlines())-1)
intron_snp_count = int(len(open(snakemake.input.intron_table).readlines())-1)
intergenic_snp_count = int(len(open(snakemake.input.intergenic_table).readlines())-1)

#create output files
outfile = open(str(snakemake.output), "w")

#calculate SNP density for genome
genome_interval = 0
for i in range(1,(len(genome_bed)-1)):
	start = int(genome_bed[i].split("\t")[1])
	stop = int(genome_bed[i].split("\t")[2].rstrip("\n"))
	window = stop - start
	genome_interval += window
genome_dens = (genome_snp_count)/(genome_interval/1000)
#calculate SNP density for genes
gene_interval = 0
for i in range(1,(len(gene_bed)-1)):
	start = int(gene_bed[i].split("\t")[1])
	stop = int(gene_bed[i].split("\t")[2].rstrip("\n"))
	window = stop - start
	gene_interval += window
gene_dens = (gene_snp_count)/(gene_interval/1000)
#calculate SNP density for exons
exon_interval = 0
for i in range(1,(len(exon_bed)-1)):
	start = int(exon_bed[i].split("\t")[1])
	stop = int(exon_bed[i].split("\t")[2].rstrip("\n"))
	window = stop - start
	exon_interval += window
exon_dens = (exon_snp_count)/(exon_interval/1000)

#calculate SNP density for introns
intron_interval = 0
for i in range(1,(len(intron_bed)-1)):
	start = int(intron_bed[i].split("\t")[1])
	stop = int(intron_bed[i].split("\t")[2].rstrip("\n"))
	window = stop - start
	intron_interval += window
intron_dens = (intron_snp_count)/(intron_interval/1000)

#calculate SNP density for intergenic regions
intergenic_interval = 0
for i in range(1,(len(intergenic_bed)-1)):
	start = int(intergenic_bed[i].split("\t")[1])
	stop = int(intergenic_bed[i].split("\t")[2].rstrip("\n"))
	window = stop - start
	intergenic_interval += window
intergenic_dens = (intergenic_snp_count)/(intergenic_interval/1000)

#calculate SNP density for cds
cds_interval = 0
for i in range(1,(len(cds_bed)-1)):
	start = int(cds_bed[i].split("\t")[1])
	stop = int(cds_bed[i].split("\t")[2].rstrip("\n"))
	window = stop - start
	cds_interval += window
cds_dens = (cds_snp_count)/(cds_interval/1000)


#write densities to outfile
outfile.writelines(["Genome SNP density: ", str(genome_dens), "; SNPs: ", str(genome_snp_count), "; interval: ",str(genome_interval)," bases\n"])
outfile.writelines(["Gene SNP density: ", str(gene_dens), "; SNPs: ", str(gene_snp_count), "; interval: ",str(gene_interval)," bases\n"])
outfile.writelines(["Exon SNP density: ", str(exon_dens), "; SNPs: ", str(exon_snp_count), "; interval: ",str(exon_interval)," bases\n"])
outfile.writelines(["Intron SNP density: ", str(intron_dens), "; SNPs: ", str(intron_snp_count), "; interval: ",str(intron_interval)," bases\n"])
outfile.writelines(["CDS SNP density: ", str(cds_dens), "; SNPs: ", str(cds_snp_count), "; interval: ",str(cds_interval)," bases\n"])
outfile.writelines(["Intergenic SNP density: ", str(intergenic_dens), "; SNPs: ", str(intergenic_snp_count), "; interval: ",str(intergenic_interval)," bases\n"])







