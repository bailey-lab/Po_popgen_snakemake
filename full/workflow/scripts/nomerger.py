###Probably unnecessary for running
##consider deleting



import yaml
import subprocess

with open("config/config.yaml") as file:
	config = yaml.safe_load(file)
	
for i in range(0,snakemake.params.index):
	subprocess.call(["cp",snakemake.input[i],snakemake.output[i]])
