## Data visualization for output of Complexity of Infection estimates from THE REAL McCOIL
setwd("~/Documents/P. Ovale Genomic Analysis")
library("dplyr")
library("ggplot2")
library("RColorBrewer")

combined_COI <- read.table("P. o. pop genomics data tables - COI_outputtable.tsv", header = T, sep = "\t")

combined_COI$abbrev <- 

combined_COI$species_ordered <- factor(combined_COI$species, levels = c("P. ovale curtisi (n=21)", "P. ovale wallikeri (n=24)", "P. falciparum (n=2,077)"))

COI <- combined_COI %>% 
  ggplot() +
  geom_violin(aes(x = species_ordered, y = Median.COI, color = species)) +
  scale_color_manual(values = c("P. ovale curtisi (n=21)"="green4","P. ovale wallikeri (n=24)"="skyblue3","P. falciparum (n=2,077)"="orange3")) +
  geom_jitter(aes(x = species_ordered, y = Median.COI, color = species), 
              alpha = 0.7, size = 0.7, height=0.1) +
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(hjust =0.5, size =18),axis.title = element_text(size = 24),axis.text= element_text(size = 18, face = 'italic'), title = element_text(size = 24)) + 
  ylim(0.8,5) +
  #ggtitle("Median Estimated Complexity of Infection") +
  labs(y = "Complexity of Infection", x = element_blank())
COI
ggsave("COI/COI_ov1.png", plot = last_plot(), device = "png", dpi = 600, width = 10, height = 5)
ggsave("COI/COI_ov1_thin.png", plot = last_plot(), device = "png", dpi = 600, width = 5, height = 5)

