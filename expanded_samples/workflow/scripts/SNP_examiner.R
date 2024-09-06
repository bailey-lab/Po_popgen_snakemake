####Examine overall SNP calls for density in genomic regions and SNPs in key genes of interest
setwd("~/Documents/P. Ovale Genomic Analysis")

poc <- read.table("ov1_curtisigh01_speciescall_genome.table", header = TRUE)
pow <- read.table("ov1_wallikericr01_speciescall_genome.table", header = TRUE)


###Curtisi: Determine whether there are snps in certain key drug resistance genes
poc$CHR <- substr(poc$CHROM,9,10)
#dhfr: last SNP has abs(nSl) of 2.8 (in top 1%)
length(print(poc$POS[(poc$CHR == "05") & (762457 < poc$POS) & (poc$POS < 764370)]))
#Chloroquine resistance transporter
length(print(poc$POS[(poc$CHR == "01") & (323118 < poc$POS) & (poc$POS < 326399)]))
#chloroquine resistance associated protein
length(print(poc$POS[(poc$CHR == "01") & (343585 < poc$POS) & (poc$POS < 347363)]))
#mdr1
length(print(poc$POS[(poc$CHR == "10") & (323034 < poc$POS) & (poc$POS < 327332)]))
#k13: abs(nSl) of 2.3
length(print(poc$POS[(poc$CHR == "12") & (404824 < poc$POS) & (poc$POS < 407001)]))
#dhps
length(print(poc$POS[(poc$CHR == "14") & (1119845 < poc$POS) & (poc$POS < 1122576)]))
#cytochrome b: first SNP has nSL of 2.25
length(print(poc$POS[(poc$CHR == "14") & (2059533 < poc$POS) & (poc$POS < 2060640)]))


###wallikeri: Determine whether there are snps in certain key drug resistance genes
pow$chr <- as.numeric(substr(pow$CHROM,7,8)) - 4

#dhfr LT594509:843,815..845,740, LT594516:1,183,447..1,185,319
length(print(pow$POS[(pow$chr == 5) & (843815 < pow$POS) & (pow$POS < 845740)]))
length(print(pow$POS[(pow$chr == 12) & (1183447 < pow$POS) & (pow$POS < 1183447)]))
#Chloroquine resistance transporter LT594505:331,639..334,593. One SNP had nSl of 2.805 (in top 1%)
length(print(pow$POS[(pow$chr == 1) & (331639 < pow$POS) & (pow$POS < 334593)]))
#chloroquine resistance associated protein LT594505:352,724..356,765
length(print(pow$POS[(pow$chr == 1) & (352724 < pow$POS) & (pow$POS < 356765)]))
#mdr1 LT594514:355,729..360,030
length(print(pow$POS[(pow$chr == 10) & (355729 < pow$POS) & (pow$POS < 360030)]))
#k13 LT594516:455,780..457,957
length(print(pow$POS[(pow$chr == 12) & (455780 < pow$POS) & (pow$POS < 457957)]))
#dhps LT594518:1,121,839..1,124,323
length(print(pow$POS[(pow$chr == 14) & (1121839 < pow$POS) & (pow$POS < 1124323)]))
#cytocrhome b LT594518:2,083,698..2,084,367
length(print(pow$POS[(pow$chr == 14) & (2083698 < pow$POS) & (pow$POS < 2084367)]))


#SNP Density line graphs
densities <- read.table("Revised - P. o. pop genomics data tables - SNP Density by Species.tsv", header = FALSE)
library("ggplot2")
library("plyr")
colnames(densities) <- c("Density", "Region", "Species")
region_order <- c("Genome","Genes","Exons","CDS","Introns","Intergenic")
species_order <- c("Poc","Pow","Pf")
ggplot(densities, aes(fill=factor(Region, region_order), y=Density, x=factor(Species, species_order))) + 
  geom_bar(position="dodge", stat="identity", order()) +
  scale_fill_manual(values = c("Genome"= "#666666","Genes" = "#66A61E","Exons"= "skyblue3","CDS" = "#7570B3","Introns" = "#E7298A","Intergenic" = "#E6AB02"), labels = c("Whole Genome","Genes","Exons","Coding Sequences","Introns","Intergenic Regions")) +
  labs(fill = "Functional Region", y = "SNP Density (SNPs/kb)", x = element_blank()) + 
  scale_x_discrete(labels = c("Poc (n=21)", "Pow (n=24)","Pf (n=33)")) + 
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 14),axis.text.x= element_text(size = 14), axis.text.y = element_text(size=14),legend.title = element_text(size = 16), legend.text = element_text(size = 14))
ggsave("sequencing/snp_density.png", plot = last_plot(), device = "png", dpi = 600, width = 6.5, height = 3)


#Examine region assignment by powdhfr mutants (Phe57Leu+Ser58Arg and Ile169Leu)
powdhfr <- read.table("Revised - P. o. pop genomics data tables - Powdhfr analysis.tsv", header = FALSE, skip = 8, sep = "\t")
#limit only to columns used for this graph
powdhfr <- powdhfr[,1:3]
library("plyr")
colnames(powdhfr) <- c("Sample","Region","Genotype")
#create new dataframe with table structure
powdhfr.dataframe <- data.frame(table(powdhfr$Region,powdhfr$Genotype))
names(powdhfr.dataframe) <- c("Region_raw","Genotype","Count")
#Manually adjust totals
powdhfr.dataframe$total <- c(rep(c(11,11,3),3))
powdhfr.dataframe$Region_factor <- factor(powdhfr.dataframe$Region_raw, levels = c("East Africa", "Central Africa", "West Africa"))
powdhfr.dataframe$Region <- revalue(powdhfr.dataframe$Region_factor, c("East Africa"="East Africa\n(n=11)", "Central Africa"="Central Africa\n(n=11)","West Africa"="West Africa\n(n=3)"))
powdhfr.dataframe$Proportion <- powdhfr.dataframe$Count / powdhfr.dataframe$total
library("ggplot2")
ggplot(powdhfr.dataframe, aes(fill=Genotype, x=Region, y=Proportion)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  labs(fill = "powdhfr Genotype") +
  theme(plot.title = element_text(hjust =0.5)) + theme_bw()+ theme(axis.title.y = element_text(size = 15),axis.title.x = element_text(size=15), title = element_text(size=20),axis.text= element_text(size = 10),legend.text = element_text(size = 15))
library(ggpubr)
ggsave("PCA/powdhfr_genotypesbyregion.png", plot = last_plot(), device = "png", dpi = 600, width = 6, height = 3)

#Examine region assignment by pocdhfr mutants (Ala15Ser, Ala15Ser-Ser58Arg, Ala15Ser-Ser113Thr, Pro58His, and Ile169Thr)
pocdhfr <- read.table("Revised - P. o. pop genomics data tables - Pocdhfr analysis.tsv", header = FALSE, skip = 14, sep = "\t")
#limit only to columns used for this graph
pocdhfr <- pocdhfr[,1:3]
library("plyr")
colnames(pocdhfr) <- c("Sample","Region","Genotype")
#create new dataframe with table structure
pocdhfr.dataframe <- data.frame(table(pocdhfr$Region,pocdhfr$Genotype))
names(pocdhfr.dataframe) <- c("Region_raw","Genotype","Count")
#Manually adjust totals
pocdhfr.dataframe$total <- c(rep(c(10,6,6),3))
pocdhfr.dataframe$Region_factor <- factor(pocdhfr.dataframe$Region_raw, levels = c("East Africa", "Central Africa", "West Africa"))
pocdhfr.dataframe$Region <- revalue(pocdhfr.dataframe$Region_factor, c("East Africa"="East Africa\n(n=6)", "Central Africa"="Central Africa\n(n=10)","West Africa"="West Africa\n(n=6)"))
pocdhfr.dataframe$Proportion <- pocdhfr.dataframe$Count / pocdhfr.dataframe$total
library("ggplot2")
ggplot(pocdhfr.dataframe, aes(fill=Genotype, x=Region, y=Proportion)) + 
  geom_bar(stat="identity") +
  labs(fill = "pocdhfr Genotype") +
  theme_bw() + theme(plot.title = element_text(hjust =0.5)) + theme_bw()+ theme(axis.title.y = element_text(size = 15),axis.title.x = element_text(size=15), title = element_text(size=20),axis.text= element_text(size = 10),legend.text = element_text(size = 15))
library(ggpubr)
ggsave("PCA/pocdhfr_genotypesbyregion.png", plot = last_plot(), device = "png", dpi = 600, width = 6, height = 3)

