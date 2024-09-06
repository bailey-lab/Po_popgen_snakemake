####Examine heterozygous SNP calls among polyclonal samples
setwd("~/Documents/P. Ovale Genomic Analysis")

library("stringr")
#install.packages("ggplot2")
library("ggplot2")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("ggtext")
library("ggtext")
## set up manhattan plot

pocchr <- read.table("curtisigh01_chr.bed")
powchr <- read.table("wallikericr01_chr.bed")

s10 <- read.table("ov1_curtisigh01_speciescall_366152-S10.table", header = TRUE)
s5 <- read.table("ov1_wallikericr01_speciescall_113135-S5.table", header = TRUE)

s10$chr <- as.numeric(substr(s10$CHROM,9,10))
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
s10 <- s10 %>% 
  inner_join(chrom_cum, by = "chr") %>% 
  mutate(pos_cum = POS + pos_add)

#convert position (POS) to kilobases
s10$position <- s10$POS/1000

s10$call <- as.factor(s10$HET)

manhplot <- ggplot(s10, aes(x = pos_cum, y = HET, size = 1)) +
  #geom_hline(yintercept = high_abs_poc_nsl, color = "grey40", linetype = "dashed") + 
  geom_point() +#, shape = factor(chr))) +
  scale_x_continuous(label = chrom_cum$chr, breaks = chrom_cum$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  #scale_color_manual(values = c("normal" = "#276FBF", "high" = "red")) +
  #scale_shape_manual(values = rep(c(20,3), length(unique(poc$chr)))) +
  #scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome", 
       y = "Het?") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(angle = 60, size = 16, vjust = 0.5)) +
  ggtitle("Het calls among sample S10")

print(manhplot)
max(s5chr6$POS)


s10chr1 <- subset(s10, chr==1, select=c(chr,POS,HET,call))
ggplot(s10chr1, aes(x = POS, y = HET, size = 1)) + geom_point(size=0.1)
s10chr1_plot <- ggplot(s10chr1, aes(x=POS, fill = factor(HET, levels = c("1","0")))) + 
  geom_histogram(binwidth = 100000) + theme_bw() + 
  scale_fill_manual(values = c("1" = "green4","0" = "skyblue3"), labels = c("1" = "heterozygous", "0" = "homozygous alternate")) + 
  labs(x="Position, Poc Chromosome 1", y = "Frequency",fill="Base call") + 
  scale_x_continuous(breaks = seq(0, 781771, by = 100000)) +
  #ggtitle("Histogram of base calls across Poc Chromosome 1 in sample 366152") + 
  theme(plot.title = element_text(hjust =0.5),axis.title.y = element_text(size = 18),axis.title.x = element_text(size=18),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12)) 

s10chr1_plot
s10chr2 <- subset(s10, chr==2, select=c(chr,POS,HET,call))
ggplot(s10chr2, aes(x = POS, y = HET, size = 1)) + geom_point(size=0.1)
s10chr2_plot <- ggplot(s10chr2, aes(x=POS, fill = factor(HET, levels = c("1","0")))) + 
  geom_histogram(binwidth = 100000) + theme_bw() + 
  scale_fill_manual(values = c("1" = "gray","0" = "gray7"), labels = c("1" = "heterozygous", "0" = "homozygous alternate")) + 
  labs(x="Position, Poc Chromosome 2", y = element_blank(),fill="Base call") + 
  scale_x_continuous(breaks = seq(0, 600000, by = 100000)) +
  #ggtitle("Histogram of base calls across Poc Chromosome 2 in sample 366152") + 
  theme(plot.title = element_text(hjust =0.5),axis.title.y = element_text(size = 22),axis.title.x = element_text(size=22),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12)) 


s10chr3 <- subset(s10, chr==3, select=c(chr,POS,HET,call))
ggplot(s10chr3, aes(x = POS, y = HET, size = 1)) + geom_point(size=0.1)
s10chr3_plot <- ggplot(s10chr3, aes(x=POS, fill = factor(HET, levels = c("1","0")))) + 
  geom_histogram(binwidth = 100000) + theme_bw() + 
  scale_fill_manual(values = c("1" = "gray","0" = "gray7"), labels = c("1" = "heterozygous", "0" = "homozygous alternate")) + 
  labs(x="Position, Poc Chromosome 3", y = element_blank(),fill="Base call") + 
  #ggtitle("Histogram of base calls across Poc Chromosome 1 in sample 366152") + 
  theme(plot.title = element_text(hjust =0.5),axis.title.y = element_text(size = 22),axis.title.x = element_text(size=22),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12)) 



s10chr4 <- subset(s10, chr==4, select=c(chr,POS,HET))
ggplot(s10chr4, aes(x = POS, y = HET, size = 1)) + geom_point(size=0.1)
s10chr4_plot <- ggplot(s10chr4, aes(x=POS, fill = factor(HET, levels = c("1","0")))) + 
  geom_histogram(binwidth = 100000) + theme_bw() + 
  scale_fill_manual(values = c("1" = "green4","0" = "skyblue3"), labels = c("1" = "heterozygous", "0" = "homozygous alternate")) + 
  labs(x="Position, Poc Chromosome 4", y = element_blank(),fill="Base call") + 
  #ggtitle("Histogram of base calls across Poc Chromosome 1 in sample 366152") + 
  theme(plot.title = element_text(hjust =0.5),axis.title.y = element_text(size = 22),axis.title.x = element_text(size=22),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12)) 

s10chr5 <- subset(s10, chr==5, select=c(chr,POS,HET))
ggplot(s10chr5, aes(x = POS, y = HET, size = 1)) + geom_point(size=0.1)
s10chr5_plot <- ggplot(s10chr5, aes(x=POS, fill = factor(HET, levels = c("1","0")))) + 
  geom_histogram(binwidth = 100000) + theme_bw() + 
  scale_fill_manual(values = c("1" = "green4","0" = "skyblue3"), labels = c("1" = "heterozygous", "0" = "homozygous alternate")) + 
  labs(x="Position, Poc Chromosome 5", y = element_blank(),fill="Base call") + 
  #ggtitle("Histogram of base calls across Poc Chromosome 1 in sample 366152") + 
  theme(plot.title = element_text(hjust =0.5),axis.title.y = element_text(size = 22),axis.title.x = element_text(size=22),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12)) 

s10chr6 <- subset(s10, chr==6, select=c(chr,POS,position,HET))
s10chr6_plot <- ggplot(s10chr6, aes(x=position)) + 
  geom_histogram(binwidth = 100, fill = "green4") + theme_bw() + 
  labs(x="Poc Chr. 6 (Kb)", y = "Frequency",fill="Variant call") + 
  scale_x_continuous(breaks = seq(0, 850, by = 400)) +
  #ggtitle("Histogram of base calls across Poc Chromosome 1 in sample 366152") + 
  theme(plot.title = element_text(hjust =0.5),axis.title.y = element_text(size = 16),axis.title.x = element_text(size=16),title = element_text(size = 20), legend.text = element_text(size = 12),axis.text= element_text(size = 12), legend.title = element_text(size=16)) 
s10chr6_plot


s5 <- read.table("ov1_wallikericr01_speciescall_113135-S5.table", header = TRUE)
s5$chr <- as.numeric(substr(s5$CHROM,7,8))-4
#convert position (POS) to kilobases
s5$position <- s5$POS/1000

s5chr1 <- subset(s5, chr==1, select=c(chr,POS,HET))


s5chr2 <- subset(s5, chr==2, select=c(chr,POS,HET))
ggplot(s5chr2, aes(x = POS, y = HET, size = 1)) + geom_point(size=0.1)

s5chr3 <- subset(s5, chr==3, select=c(chr,POS,HET))
ggplot(s5chr3, aes(x = POS, y = HET, size = 1)) + geom_point(size=0.1)

s5chr4 <- subset(s5, chr==4, select=c(chr,POS,HET))
ggplot(s5chr4, aes(x = POS, y = HET, size = 1)) + geom_point(size=0.1)

s5chr5 <- subset(s5, chr==5, select=c(chr,POS,HET))
ggplot(s5chr5, aes(x = POS, y = HET, size = 1)) + geom_point(size=0.1)
ggplot(s5chr1, aes(x = POS, y = HET, size = 1)) + geom_point(size=0.1)
s5chr5_plot <- ggplot(s5chr5, aes(x=POS, fill = factor(HET, levels = c("1","0")))) + 
  geom_histogram(binwidth = 100000) + theme_bw() + 
  scale_fill_manual(values = c("1" = "green4","0" = "skyblue3"), labels = c("1" = "heterozygous", "0" = "homozygous alternate")) + 
  labs(x="Position, Pow Chromosome 5", y = element_blank(),fill="Base call") + 
  #ggtitle("Histogram of base calls across Poc Chromosome 1 in sample 366152") + 
  theme(plot.title = element_text(hjust =0.5),axis.title.y = element_text(size = 22),axis.title.x = element_text(size=22),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12)) 

s5chr6 <- subset(s5, chr==6, select=c(chr,POS,HET,position))
s5chr6_plot <- ggplot(s5chr6, aes(x=position)) + 
  geom_histogram(binwidth = 100, fill = "skyblue3") + theme_bw() + 
  labs(x="Pow Chr. 6 (Kb)", y = element_blank()) + 
  scale_x_continuous(breaks = seq(0, 1000, by = 500)) +
  #ggtitle("Histogram of base calls across Poc Chromosome 1 in sample 366152") + 
  theme(plot.title = element_text(hjust =0.5),axis.title.y = element_text(size = 16),axis.title.x = element_text(size=16),title = element_text(size = 20), legend.text = element_text(size = 12),axis.text= element_text(size = 12), legend.title = element_text(size=16)) 


ggarrange(s10chr6_plot, s5chr6_plot, ncol=2, nrow=1, common.legend = TRUE, legend="right", widths = c(5,5), labels = c("A","B"), hjust = 0)
ggsave("sequencing/polyclonal_chr6.png", plot = last_plot(), device = "png", dpi = 600, width = 8, height = 4)

