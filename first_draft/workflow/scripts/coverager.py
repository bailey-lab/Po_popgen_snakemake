import yaml

with open("config/config.yaml") as file:
	config = yaml.safe_load(file)
outfile = open(str(snakemake.output), "w")
for i in snakemake.input.cov:
	total = 0
	weighted_cov = 0
	with open(snakemake.input.bed) as l:
		next(l)
		for line in l:
			chrom = line.split("\t")[0]
			length = int(line.split("\t")[2].rstrip("\n"))
			total += length
			for row in open(i):
				spot = row.split("\t")[0]
				depth = row.split("\t")[1]
				coverage = float(row.split("\t")[4].rstrip("\n"))
				if (spot == chrom) and (depth == snakemake.wildcards.depth):
					weighted_cov += length*coverage
	print(weighted_cov)
	print(total)
	total_coverage = weighted_cov/total
	print(total_coverage)
	outfile.write(i+"\t"+str(total_coverage)+"\n")
