###Plot genome-wide coverage and mapping percentage among all samples
setwd("~/Documents/P. Ovale Genomic Analysis")
library(ggplot2)
cov <- read.table("P. o. pop genomics data tables - Po coverage and depth.tsv", header = T, sep = "\t")
cov <- cov[1:29,]
cov$group <- as.factor(cov$species)
cov$tech <- as.factor(cov$sample.type)
cov$Coverage10X <- cov$Coverage..excluding.chr10.
cov$Mapping <- cov$Mapping.proportion

ggplot(cov, aes(x=Coverage10X, fill = group)) + geom_histogram(binwidth = 0.05) + theme_bw() + scale_fill_manual(values = c("curtisi" = "skyblue3","mixed" = "mediumpurple3","wallikeri" = "green4")) + labs(x="Proportion of genome with >10 reads mapped", y = "Frequency",fill="Species") + ggtitle("Histogram of 10X genome coverage per sample") + theme(plot.title = element_text(hjust =0.5)) + theme_bw()+ theme(axis.title.y = element_text(size = 22),axis.title.x = element_text(size=22),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12)) 
ggplot(cov, aes(x=Mapping, fill = group)) + geom_histogram(binwidth = 0.1) + theme_bw() + scale_fill_manual(values = c("curtisi" = "skyblue3","mixed" = "mediumpurple3","wallikeri" = "green4")) + labs(x="% Mapping to P. ovale reference genome after competitive alignment", y = "Frequency",fill="Species") + ggtitle("Histogram of ovale Mapping Percentage") + theme(plot.title = element_text(hjust =0.5)) + theme_bw()+ theme(axis.title.y = element_text(size = 22),axis.title.x = element_text(size=15), title = element_text(size=20),axis.text= element_text(size = 12),legend.text = element_text(size = 15))
ggplot(cov, aes(x=Coverage10X, fill = tech)) + geom_histogram(binwidth = 0.05) + theme_bw() + scale_fill_manual(values = c("HC" = "#CC6666","LDB" = "#9999CC","sWGA" = "#66CC99")) + labs(x="Proportion of genome with >10 reads mapped", y = "Frequency",fill="Technique") + ggtitle("Histogram of 10X genome coverage per sample") + theme(plot.title = element_text(hjust =0.5)) + theme_bw()+ theme(axis.title.y = element_text(size = 22),axis.title.x = element_text(size=22),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12))
ggplot(cov, aes(x=Mapping, fill = tech)) + geom_histogram(binwidth = 0.1) + theme_bw() + scale_fill_manual(values = c("HC" = "#CC6666","LDB" = "#9999CC","sWGA" = "#66CC99")) + labs(x="% Mapping to P. ovale reference genome after competitive alignment", y = "Frequency",fill="Technique") + ggtitle("Histogram of ovale Mapping Percentage") + theme(plot.title = element_text(hjust =0.5)) + theme_bw()+ theme(axis.title.y = element_text(size = 22),axis.title.x = element_text(size=15), title = element_text(size=20),axis.text= element_text(size = 12),legend.text = element_text(size = 15))

species_cov <- ggplot(cov, aes(x=Coverage10X, fill = group)) + 
  geom_histogram(binwidth = 0.05) + 
  theme_bw() + 
  scale_fill_manual(values = c("curtisi" = "green4","mixed" = "mediumpurple3","wallikeri" = "skyblue3")) + 
  labs(x="Proportion with 10x coverage", y = "Frequency",fill="Species") + 
  theme(plot.title = element_text(hjust =0.5)) + 
  theme_bw()+ 
  theme(axis.title.y = element_text(size = 15),axis.title.x = element_text(size=15),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12))
  #geom_text(data=data.frame(), aes(label = 'A', x = 0, y = 100), hjust = 0, vjust = 1)
species_map <- ggplot(cov, aes(x=Mapping, fill = group)) + geom_histogram(binwidth = 0.1) + theme_bw() + scale_fill_manual(values = c("curtisi" = "green4","mixed" = "mediumpurple3","wallikeri" = "skyblue3")) + labs(x="Proportion mapping", y = element_blank(),fill="Species") + theme(plot.title = element_text(hjust =0.5)) + theme_bw()+ theme(axis.title.y = element_text(size = 15),axis.title.x = element_text(size=15), title = element_text(size=15),axis.text= element_text(size = 12),legend.text = element_text(size = 15))
tech_cov <- ggplot(cov, aes(x=Coverage10X, fill = tech)) + geom_histogram(binwidth = 0.05) + theme_bw() + scale_fill_manual(values = c("HC" = "#CC6666","LDB" = "#9999CC","sWGA" = "#66CC99")) + labs(x="Proportion with 10x coverage", y = "Frequency",fill="Technique")  + theme(plot.title = element_text(hjust =0.5)) + theme_bw()+ theme(axis.title.y = element_text(size = 15),axis.title.x = element_text(size=15),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12))
tech_map <- ggplot(cov, aes(x=Mapping, fill = tech)) + geom_histogram(binwidth = 0.1) + theme_bw() + scale_fill_manual(values = c("HC" = "#CC6666","LDB" = "#9999CC","sWGA" = "#66CC99")) + labs(x="Proportion mapping", y = element_blank(),fill="Technique")  + theme(plot.title = element_text(hjust =0.5)) + theme_bw()+ theme(axis.title.y = element_text(size = 15),axis.title.x = element_text(size=15), title = element_text(size=20),axis.text= element_text(size = 12),legend.text = element_text(size = 15))
library("ggpubr")

ggarrange(species_cov, species_map, ncol=2, nrow=1, common.legend = TRUE, legend="right", widths = c(5,5), labels = c("B","D"), hjust = 0)
ggsave("sequencing/coverage_species.png", plot = last_plot(), device = "png", dpi = 600, width = 8, height = 3)

ggarrange(tech_cov, tech_map, ncol=2, nrow=1, common.legend = TRUE, legend="right", widths = c(5,5),labels = c("A","C"), hjust = 0)
ggsave("sequencing/coverage_technique.png", plot = last_plot(), device = "png", dpi = 600, width = 8, height = 3)
