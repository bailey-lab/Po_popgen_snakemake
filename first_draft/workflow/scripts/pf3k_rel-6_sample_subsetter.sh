#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=10-23:00:00
#SBATCH --mem=100G

module load vcftools

#Filter Pf6k Pf3D7 vcf down to 24 specific samples selected to match distribution of ovale WGS samples across Sub-Saharan Africa
cd /work/users/k/e/kellybce/ovale1r/snakemake/full/input/pf6k_vcfs
for i in Pf_60*v3.final.vcf.gz;
do vcftools --gzvcf ${i} --indv QG0007-C --indv QG0015-C --indv QG0140-C --indv QG0189-C --indv QG0270-C --indv QS0109-C --indv QS0110-C --indv PE0538-C --indv QS0154-C  --indv QS0159-C --indv QS0168-C --indv QS0169-C --indv QS0170-C --indv PE0129-C --indv PE0130-C --indv PE0109-C --indv QJ0159-C --indv QQ0111-C --indv QJ0008-C --indv QV0030-C --indv QJ0172-C  --non-ref-ac-any 1 --recode --stdout | gzip -c > ../pf_subset_vcfs/${i%.final.vcf.gz}ov1_no-invariant_sites.vcf.gz;
done

cd /work/users/k/e/kellybce/ovale1r/snakemake/full/input/pf_subset_vcfs

vcf-concat Pf_60*ov1_no-invariant_sites.vcf.gz | bgzip -c > ../../output/vcfs/ov1/ov1_pf3d7_only.vcf.gz


#module load bcftools
#bcftools index -t ../../output/vcfs/ov1/ov1_pf3d7_only.vcf.gz
