import yaml
import subprocess
import os

with open("config/config.yaml") as file:
	config = yaml.safe_load(file)

lane1 = str(snakemake.input.lane1)
lane2 = str(snakemake.params.lane2)
output = str(snakemake.output)
print(lane2)
print(output)
if os.path.exists(lane2):
	subprocess.call(["samtools","merge","-o",output,lane1,lane2])
else:
	subprocess.call(["cp",lane1,output])
