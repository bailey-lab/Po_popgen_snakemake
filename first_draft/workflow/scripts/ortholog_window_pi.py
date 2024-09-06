import yaml
import subprocess

with open("config/config.yaml") as file:
	config = yaml.safe_load(file)

print(len(snakemake.output.pi)-1)
print(snakemake.params.ortholist[1200][1])
print(snakemake.output.pi[len(snakemake.output.pi)-1][0])
for n in range(0,len(snakemake.output.pi)):
	subprocess.call(['vcftools', '--gzvcf',snakemake.input.vcf,'--out',snakemake.params.outfile[n],'--chr',snakemake.params.ortholist[n][1],'--from-bp',snakemake.params.ortholist[n][2],'--to-bp',snakemake.params.ortholist[n][3],'--window-pi',str(snakemake.params.ortholist[n][4]),'--window-pi-step',' 1'])
	
#subprocess.call(['vcftools', '--gzvcf', snakemake.input.pfvcf, '--out', snakemake.params.pfoutfile, '--chr', chrom, '--from-bp', start, '--to-bp', stop, '--recode', '--recode-INFO-all'])
