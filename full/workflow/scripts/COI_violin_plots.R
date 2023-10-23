setwd("~/Documents/P. Ovale Genomic Analysis")
rawcurtisi <- read.table("ov1_curtisigh01_speciescall_coi.txt", header = F)
limitcurtisi <- rawcurtisi[1:999,2:14]
rawwallikeri <- read.table("ov1_wallikericr01_speciescall_coi.txt", header = F)
limitwallikeri <- rawwallikeri[1:999,2:17]
rawfalciparum <- read.table("ov1_pf3d7_only_coi.txt", header = F)
limitfalciparum <- rawfalciparum[1:999,2:25]
limitcurtisit <- t(limitcurtisi)
limitwallikerit <- t(limitwallikeri)
limitfalciparumt <- t(limitfalciparum)
coi <- data.frame(rbind(limitcurtisit,limitwallikerit,limitfalciparumt))
coit <- data.frame(t(coi))
species <- c(rep("curtisi",13),rep("wallikeri",16),rep("falciparum",24))
colors <- c(rep("skyblue3",13),rep("green4",16), rep("red4",24))
colnames(coit) <- c("X01","X02","X03","X04","X05","X06",
                                          "X07","X08","X09","X10","X11","X12",
                                          "X13","X14","X23","X15","X16","X24",
                                          "X17","X25","X26","X27","X28","X18",
                                          "X19","X20","X21","X22","X29","X30",
                                          "X31","X32","X33","X34","X35","X36",
                                          "X37","X38","X39","X40","X41","X42",
                                          "X43","X44","X45","X46","X47","X48",
                                          "X49","X50","X51","X52","X53")
                      #coit$Xtest <- rep(1:10,100)
            library(tidyverse)
            coit %>% 
            mutate(id = row_number()) %>% 
            pivot_longer(
            cols = starts_with("X"),
            names_to = "names",
            values_to = "values"
            ) %>%
            ggplot(aes(x=names, y=values, fill=names)) +
            geom_violin() +
            theme_bw() +
            scale_fill_manual(values = colors) +
            labs(y="Esimated number of clones per iteration",x="Sample") + theme(legend.position="none", axis.text.x =element_blank(),axis.ticks.x = element_blank()) +
            theme(axis.title = element_text(size = 18),axis.text= element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 15)) 
                      

                      
                      
                      
 