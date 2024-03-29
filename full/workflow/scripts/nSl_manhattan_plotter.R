setwd("~/Documents/P. Ovale Genomic Analysis")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
poc <- read.table("ov1_curtisigh01_speciescall_nsl_total.txt")
pow <- read.table("ov1_wallikericr01_speciescall_nsl_total.txt")
pocchr <- read.table("curtisigh01_chr.bed")
powchr <- read.table("wallikericr01_chr.bed")
#Plot mapping quality by position across every chromosome
library("stringr")
#install.packages("ggplot2")
library("ggplot2")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("ggtext")
library("ggtext")
## set up manhattan plot



### Poc nSl for all variants

#rename variables from table
colnames(poc) <- c("chrom","pos","nsl")
poc$chr <- as.numeric(substr(poc$chrom,9,10))
colnames(pocchr) <- c("chrom","start","end")
pocchr$chr <- as.numeric(substr(pocchr$chrom,9,10))

#determine absolute nSl values
poc$nsl_abs <- abs(poc$nsl)

#determine cumulative position of each chromosome in genome
chrom_cum <- pocchr %>% 
  #determine the length of all preceding chromosomes
  mutate(pos_add = lag(cumsum(end), default = 0)) %>% 
  #then add half the length of the current chromosome to find the overall center in the genome
  mutate(center = pos_add+(end/2)) %>%
  select(chr, pos_add, center)


#add pos_cum variable, reflecting positing along the entire genome
poc <- poc %>% 
  inner_join(chrom_cum, by = "chr") %>% 
  mutate(pos_cum = pos + pos_add)


#determine high and low values for nSl, top and bottom 0.5%
low_poc_nsl <- quantile(poc$nsl,0.005)
high_poc_nsl <- quantile(poc$nsl,0.995)
high_abs_poc_nsl <- quantile(poc$nsl_abs,0.995)

#generate plot
manhplot <- ggplot(poc, aes(x = pos_cum, y = nsl_abs,
                    color = factor(nsl_abs), shape = factor(chr), size = 3)) +
                    geom_hline(yintercept = high_abs_poc_nsl, color = "grey40", linetype = "dashed") + 
                    geom_point(alpha = 0.75, aes(color = cut(nsl_abs, c(-Inf, high_abs_poc_nsl, Inf), labels = c("normal","high")))) +#, shape = factor(chr))) +
                    scale_x_continuous(label = chrom_cum$chr, breaks = chrom_cum$center) +
                    scale_y_continuous(expand = c(0,0), limits = c(0, 3.5)) +
                    scale_color_manual(values = c("normal" = "#276FBF", "high" = "red")) +
                    scale_shape_manual(values = rep(c(20,3), length(unique(poc$chr)))) +
                    scale_size_continuous(range = c(0.5,3)) +
                    labs(x = "Chromosome", 
                          y = "Absolute nSl") + 
                    theme_minimal() +
         theme(plot.title = element_text(hjust = 0.5, size = 20),
           legend.position = "none",
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank(),
               axis.title.y = element_text(size=18),
              axis.title.x = element_text(size=18),
              axis.text.y = element_text(size = 16),
              axis.text.x = element_text(angle = 60, size = 16, vjust = 0.5)) +
            ggtitle("Absolute nSL for variants in Poc samples")
    
print(manhplot)

#print loci with nSl above high_poc_nsl
print(poc[poc$nsl_abs >= high_abs_poc_nsl, c("chrom","chr","pos", "nsl")])
print(poc[poc$nsl_abs >= high_abs_poc_nsl, c("chrom")])
print(poc[poc$chr == 9 & poc$pos > 1196000, c("chrom","chr","pos", "nsl")])
#Investigate particular genes of interest
#dhfr
print(poc$nsl_abs[(poc$chr == 5) & (762457 < poc$pos) & (poc$pos < 764370)])
#chloroquone resistance transporter
print(poc$nsl_abs[(poc$chr == 1) & (323118 < poc$pos) & (poc$pos < 326399)])
#chloroquine resistance associated protein
print(poc$nsl_abs[(poc$chr == 1) & (343585 < poc$pos) & (poc$pos < 347363)])
#k13
print(poc$nsl_abs[(poc$chr == 12) & (404824 < poc$pos) & (poc$pos < 407001)])
#dhps
print(poc$nsl_abs[(poc$chr == 14) & (1119845 < poc$pos) & (poc$pos < 1122576)])
#cytochrome b
print(poc$nsl_abs[(poc$chr == 14) & (2059533 < poc$pos) & (poc$pos < 2060640)])


### Pow nSl plot by variant

colnames(pow) <- c("chrom","pos","nsl")
pow$chr <- as.numeric(substr(pow$chrom,7,8))-4
colnames(powchr) <- c("chrom","start","end")
powchr$chr <- as.numeric(substr(powchr$chrom,7,8))-4

#Determine absolute nSl values
pow$nsl_abs <- abs(pow$nsl)

#determine cumulative position of each chromosome in genome
chrom_cum <- powchr %>% 
  #determine the length of all preceding chromosomes
  mutate(pos_add = lag(cumsum(end), default = 0)) %>% 
  #then add half the length of the current chromosome to find the overall center in the genome
  mutate(center = pos_add+(end/2)) %>%
  select(chr, pos_add, center)

#add pos_cum variable, reflecting positing along the entire genome
pow <- pow %>% 
  inner_join(chrom_cum, by = "chr") %>% 
  mutate(pos_cum = pos + pos_add)


#determine high and low values for nSl, top and bottom 0.5%
low_pow_nsl <- quantile(pow$nsl,0.005)
high_pow_nsl <- quantile(pow$nsl,0.995)
high_abs_pow_nsl <- quantile(pow$nsl_abs,0.995)

#generate plot
manhplot <- ggplot(pow, aes(x = pos_cum, y = nsl_abs,
                            color = factor(nsl_abs), shape = factor(chr), size = 3)) +
  geom_hline(yintercept = high_abs_pow_nsl, color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75, aes(color = cut(nsl_abs, c(-Inf, high_abs_pow_nsl, Inf), labels = c("normal","high")))) +
  scale_x_continuous(label = chrom_cum$chr, breaks = chrom_cum$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 3.5)) +
  scale_color_manual(values = c("normal" = "#276FBF", "high" = "red")) +
  scale_size_continuous(range = c(0.5,3)) +
  scale_shape_manual(values = rep(c(20,3), length(unique(pow$chr)))) +
  labs(x = "Chromosome", 
       y = "Absolute nSl") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(angle = 60, size = 16, vjust = 0.5)) +
  ggtitle("Absolute nSl for non-missing variants in monoclonal Pow samples")

print(manhplot)

#print loci with nSl above high_pow_nsl
print(pow[pow$nsl_abs >= high_abs_pow_nsl, c("chrom","chr","pos", "nsl")])
print(pow[pow$chr == 9 & pow$pos > 1220000, c("chrom","chr","pos", "nsl")])

#Investigate particular SNPs of interest
#dhfr LT594509:843,815..845,740, LT594516:1,183,447..1,185,319
print(pow$nsl_abs[(pow$chr == 5) & (843815 < pow$pos) & (pow$pos < 845740)])
print(pow$nsl_abs[(pow$chr == 12) & (1183447 < pow$pos) & (pow$pos < 1183447)])
#Chloroquine resistance transporter LT594505:331,639..334,593
print(pow$nsl_abs[(pow$chr == 1) & (331639 < pow$pos) & (pow$pos < 334593)])
#chloroquine resistance associated protein LT594505:352,724..356,765
print(pow$nsl_abs[(pow$chr == 1) & (352724 < pow$pos) & (pow$pos < 356765)])
#mdr1 LT594514:355,729..360,030
print(pow$nsl_abs[(pow$chr == 10) & (355729 < pow$pos) & (pow$pos < 360030)])
#k13 LT594516:455,780..457,957
print(pow$nsl_abs[(pow$chr == 12) & (455780 < pow$pos) & (pow$pos < 457957)])
#dhps LT594518:1,121,839..1,124,323
print(pow$nsl_abs[(pow$chr == 14) & (1121839 < pow$pos) & (pow$pos < 1124323)])
#cytocrhome b LT594518:2,083,698..2,084,367
print(pow$nsl_abs[(pow$chr == 14) & (2083698 < pow$pos) & (pow$pos < 2084367)])

#prepare Tajima D plots
pocexons <- read.table("ov1_curtisigh01_speciescall_tajimad_exons.txt", header = TRUE)
pocgenes <- read.table("ov1_curtisigh01_speciescall_tajimad_genes.txt", header = TRUE)
powexons <- read.table("ov1_wallikericr01_speciescall_tajimad_exons.txt", header = TRUE)
powgenes <- read.table("ov1_wallikericr01_speciescall_tajimad_genes.txt", header = TRUE)

pocchr <- read.table("curtisigh01_chr.bed")
powchr <- read.table("wallikericr01_chr.bed")


library("stringr")
#install.packages("ggplot2")
library("ggplot2")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("ggtext")
library("ggtext")
## set up manhattan plot

###Curtisi exons
pocexons$chr <- as.numeric(substr(pocexons$CHROM,9,10))
pocexons$pos <- pocexons$BIN_START
colnames(pocchr) <- c("chrom","start","end")
pocchr$chr <- as.numeric(substr(pocchr$chrom,9,10))

#determine cumulative position of each chromosome in genome
chrom_cum <- pocchr %>% 
  #determine the length of all preceding chromosomes
  mutate(pos_add = lag(cumsum(end), default = 0)) %>% 
  #then add half the length of the current chromosome to find the overall center in the genome
  mutate(center = pos_add+(end/2)) %>%
  select(chr, pos_add, center)


#add pos_cum variable, reflecting positing along the entire genome
pocexons <- pocexons %>% 
  inner_join(chrom_cum, by = "chr") %>% 
  mutate(pos_cum = pos + pos_add)


#determine high and low values for Tajima's D, top and bottom 0.5%
low_pocexons_d <- quantile(pocexons$TajimaD,0.005)
high_pocexons_d <- quantile(pocexons$TajimaD,0.995)

#Graph Tajima D across genomes
manhplot <- ggplot(pocexons, aes(x = pos_cum, y = TajimaD,
                            color = factor(TajimaD), size = 3)) +
  geom_hline(yintercept = high_pocexons_d, color = "grey40", linetype = "dashed") + 
  #geom_hline(yintercept = low_pocexons_d, color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75, aes(color = cut(TajimaD, c(-Inf, low_pocexons_d, high_pocexons_d, Inf), labels = c("low","normal","high")))) +
  scale_x_continuous(label = chrom_cum$chr, breaks = chrom_cum$center) +
  scale_y_continuous(expand = c(0,0), limits = c(-3, 3)) +
  scale_color_manual(values = c("low" = "red", "normal" = "#276FBF", "high" = "red")) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome", 
       y = "Tajima's D") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size= 20),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(angle = 60, size = 15, vjust = 0.5)) +
  ggtitle("Tajima's D for exons of monoclonal Poc samples")

print(manhplot)


###Curtisi genes
pocgenes$chr <- as.numeric(substr(pocgenes$CHROM,9,10))
pocgenes$pos <- pocgenes$BIN_START
colnames(pocchr) <- c("chrom","start","end")
pocchr$chr <- as.numeric(substr(pocchr$chrom,9,10))

#determine cumulative position of each chromosome in genome
chrom_cum <- pocchr %>% 
  #determine the length of all preceding chromosomes
  mutate(pos_add = lag(cumsum(end), default = 0)) %>% 
  #then add half the length of the current chromosome to find the overall center in the genome
  mutate(center = pos_add+(end/2)) %>%
  select(chr, pos_add, center)


#add pos_cum variable, reflecting positing along the entire genome
pocgenes <- pocgenes %>% 
  inner_join(chrom_cum, by = "chr") %>% 
  mutate(pos_cum = pos + pos_add)


#determine high and low values for Tajima's D, top and bottom 0.5%
low_pocgenes_d <- quantile(pocgenes$TajimaD,0.005)
high_pocgenes_d <- quantile(pocgenes$TajimaD,0.995)

#Graph Tajima D across genomes
manhplot <- ggplot(pocgenes, aes(x = pos_cum, y = TajimaD,
                                 color = factor(TajimaD), size = 3)) +
  geom_hline(yintercept = high_pocgenes_d, color = "grey40", linetype = "dashed") + 
  #geom_hline(yintercept = low_pocgenes_d, color = "grey40", linetype = "dashed") + 
  #geom_point(alpha = 0.75, aes(color = cut(TajimaD, c(-Inf, low_pocgenes_d, high_pocgenes_d, Inf), labels = c("low","normal","high")))) +
  geom_point(alpha = 0.75, aes(color = cut(TajimaD, c(-Inf,high_pocgenes_d, Inf), labels = c("normal","high")))) +
  scale_x_continuous(label = chrom_cum$chr, breaks = chrom_cum$center) +
  scale_y_continuous(expand = c(0,0), limits = c(-2.5, 2.5)) +
  scale_color_manual(values = c("low" = "red", "normal" = "#276FBF", "high" = "red")) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome", 
       y = "Tajima's D") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size= 20),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(angle = 60, size = 15, vjust = 0.5)) +
  ggtitle("Tajima's D for genes of monoclonal Poc samples")

print(manhplot)

#Output windows with high Tajima's D
print(pocgenes[pocgenes$TajimaD > high_pocgenes_d, c("CHROM","BIN_START","TajimaD")])
print(pocgenes[pocgenes$TajimaD < low_pocgenes_d, c("CHROM","BIN_START","TajimaD")])
###Examine curtisi windows of interest
poc <- pocgenes
#dhfr
print(poc$TajimaD[(poc$chr == 5) & (762457 < poc$pos) & (poc$pos < 764370)])
#chloroquone resistance transporter
print(poc$TajimaD[(poc$chr == 1) & (323118 < poc$pos) & (poc$pos < 326399)])
#chloroquine resistance associated protein
print(poc$TajimaD[(poc$chr == 1) & (343585 < poc$pos) & (poc$pos < 347363)])
#k13
print(poc$TajimaD[(poc$chr == "12") & (404824 < poc$pos) & (poc$pos < 407001)])
#dhps
print(poc$TajimaD[(poc$chr == "14") & (1119845 < poc$pos) & (poc$pos < 1122576)])
#cytochrome b
print(poc$TajimaD[(poc$chr == "14") & (2059533 < poc$pos) & (poc$pos < 2060640)])

###wallikeri exons
powexons$chr <- as.numeric(substr(powexons$CHROM,7,8))-4
powexons$pos <- powexons$BIN_START
colnames(powchr) <- c("chrom","start","end")
powchr$chr <- as.numeric(substr(powchr$chrom,7,8))-4

#determine cumulative position of each chromosome in genome
chrom_cum <- powchr %>% 
  #determine the length of all preceding chromosomes
  mutate(pos_add = lag(cumsum(end), default = 0)) %>% 
  #then add half the length of the current chromosome to find the overall center in the genome
  mutate(center = pos_add+(end/2)) %>%
  select(chr, pos_add, center)


#add pos_cum variable, reflecting positing along the entire genome
powexons <- powexons %>% 
  inner_join(chrom_cum, by = "chr") %>% 
  mutate(pos_cum = pos + pos_add)


#determine high and low values for Tajima's D, top and bottom 0.5%
low_powexons_d <- quantile(powexons$TajimaD,0.005)
high_powexons_d <- quantile(powexons$TajimaD,0.995)

#Graph Tajima D across genomes
manhplot <- ggplot(powexons, aes(x = pos_cum, y = TajimaD,
                                 color = factor(TajimaD), size = 3)) +
  geom_hline(yintercept = high_powexons_d, color = "grey40", linetype = "dashed") + 
  #geom_hline(yintercept = low_powexons_d, color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75, aes(color = cut(TajimaD, c(-Inf, low_powexons_d, high_powexons_d, Inf), labels = c("low","normal","high")))) +
  scale_x_continuous(label = chrom_cum$chr, breaks = chrom_cum$center) +
  scale_y_continuous(expand = c(0,0), limits = c(-3, 3)) +
  scale_color_manual(values = c("low" = "red", "normal" = "#276FBF", "high" = "red")) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome", 
       y = "Tajima's D") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size= 20),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(angle = 60, size = 15, vjust = 0.5)) +
  ggtitle("Tajima's D for exons of monoclonal Pow samples")

print(manhplot)


###wallikeri genes
powgenes$chr <- as.numeric(substr(powgenes$CHROM,7,8))-4
powgenes$pos <- powgenes$BIN_START
colnames(powchr) <- c("chrom","start","end")
powchr$chr <- as.numeric(substr(powchr$chrom,7,8))-4

#determine cumulative position of each chromosome in genome
chrom_cum <- powchr %>% 
  #determine the length of all preceding chromosomes
  mutate(pos_add = lag(cumsum(end), default = 0)) %>% 
  #then add half the length of the current chromosome to find the overall center in the genome
  mutate(center = pos_add+(end/2)) %>%
  select(chr, pos_add, center)


#add pos_cum variable, reflecting positing along the entire genome
powgenes <- powgenes %>% 
  inner_join(chrom_cum, by = "chr") %>% 
  mutate(pos_cum = pos + pos_add)


#determine high and low values for tajima's, top and bottom 0.5%
low_powgenes_d <- quantile(powgenes$TajimaD,0.005)
high_powgenes_d <- quantile(powgenes$TajimaD,0.995)

#Graph Tajima D across genomes
manhplot <- ggplot(powgenes, aes(x = pos_cum, y = TajimaD,
                                 color = factor(TajimaD), size = 3)) +
  geom_hline(yintercept = high_powgenes_d, color = "grey40", linetype = "dashed") + 
  #geom_hline(yintercept = low_powgenes_d, color = "grey40", linetype = "dashed") + 
  #geom_point(alpha = 0.75, aes(color = cut(TajimaD, c(-Inf, low_powgenes_d, high_powgenes_d, Inf), labels = c("low","normal","high")))) +
  geom_point(alpha = 0.75, aes(color = cut(TajimaD, c(-Inf,high_powgenes_d, Inf), labels = c("normal","high")))) +
  scale_x_continuous(label = chrom_cum$chr, breaks = chrom_cum$center) +
  scale_y_continuous(expand = c(0,0), limits = c(-2.5, 2.75)) +
  scale_color_manual(values = c("low" = "red", "normal" = "#276FBF", "high" = "red")) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome", 
       y = "Tajima's D") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size= 20),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(angle = 60, size = 15, vjust = 0.5)) +
  ggtitle("Tajima's D for genes of monoclonal Pow samples")

print(manhplot)

#Print windows above the 0.5% cut-off
print(powgenes[powgenes$TajimaD > high_powgenes_d, c("CHROM","BIN_START","TajimaD", "N_SNPS")])

#Investigate particular SNPs of interest
pow <- powgenes
#dhfr LT594509:843,815..845,740, LT594516:1,183,447..1,185,319
print(pow$TajimaD[(pow$chr == 5) & (843815 < pow$pos) & (pow$pos < 845740)])
print(pow$TajimaD[(pow$chr == 12) & (1183447 < pow$pos) & (pow$pos < 1183447)])
#Chloroquine resistance transporter LT594505:331,639..334,593
print(pow$TajimaD[(pow$chr == 1) & (331639 < pow$pos) & (pow$pos < 334593)])
#chloroquine resistance associated protein LT594505:352,724..356,765
print(pow$TajimaD[(pow$chr == 1) & (352724 < pow$pos) & (pow$pos < 356765)])
#mdr1 LT594514:355,729..360,030
print(pow$TajimaD[(pow$chr == 10) & (355729 < pow$pos) & (pow$pos < 360030)])
#k13 LT594516:455,780..457,957
print(pow$TajimaD[(pow$chr == 12) & (455780 < pow$pos) & (pow$pos < 457957)])
#dhps LT594518:1,121,839..1,124,323
print(pow$TajimaD[(pow$chr == 14) & (1121839 < pow$pos) & (pow$pos < 1124323)])
#cytocrhome b LT594518:2,083,698..2,084,367
print(pow$TajimaD[(pow$chr == 14) & (2083698 < pow$pos) & (pow$pos < 2084367)])


#if need to adjust y-axis to contain all values
#ylim <- gwas_data %>% 
#  filter(p == min(p)) %>% 
#  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
#  pull(ylim)


