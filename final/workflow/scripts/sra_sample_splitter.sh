#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=10-23:00:00
#SBATCH --mem=100G

module load sratoolkit
module load parallel-fastq-dump

#from list of sample names and metadata, download split fastq files from SRA and deposit in reads/ directory
cd /work/users/k/e/kellybce/ovale1r/snakemake/expanded_samples/input/original_reads/
cat ../lists/expanded_samples.txt | while read sraid sampleid; 
do 
parallel-fastq-dump -s $sraid --split-3
gzip $sraid"_1.fastq"
gzip $sraid"_2.fastq"
cp $sraid"_1.fastq.gz" "../reads/"$sraid$sampleid"_L001_R1.fastq.gz"
cp $sraid"_2.fastq.gz" "../reads/"$sraid$sampleid"_L001_R2.fastq.gz"
done

