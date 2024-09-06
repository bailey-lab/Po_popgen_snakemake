###Principal Components Analysis of P. o. genomic data
setwd("~/Documents/P. Ovale Genomic Analysis")

#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n=8,"Dark2")
brewer.pal(n=8,"Dark2")

#Load in Poc datatables (including country of origin and eigenvalues)
poc <- read.table("P. o. pop genomics data tables - Poc PCA r2-0.3.tsv", header = TRUE, sep = "\t")
poceigenval <- read.table("P. o. pop genomics data tables - Poc Eigenvalues r2-0.3.tsv", header = FALSE, sep = "\t")
colnames(poceigenval) <- c("eigenvalues")

### P. o. curtisi
###Determine contribution of each principal component
poc_eigen_total <- sum(poceigenval) 
poceigenval$prop <- poceigenval$eigenvalues/poc_eigen_total
poceigenval$percent <- round(poceigenval$prop*100,1)
poceigenval$pc <- c(1:length(poceigenval$eigenvalues))
plot(poceigenval$pc, poceigenval$prop, type="o", col="red", xlab="PC", ylab="Proportion of variance explained", main="Scree plot I")
###Rename principal component variables to include amount of total variation explained
for (i in 1:length(poceigenval$eigenvalues)){
  colnames(poc)[i+3] <- paste(colnames(poc)[i+3]," (",poceigenval$percent[i],"%)", sep = "")
}

#PCA plots for geographical clustering from plink
library(ggplot2)

pca_colors <- c("DRC (West)" = "#D95F02", "DRC (East)"= "#7570B3","Tanzania (East)"="#E7298A", "Tanzania (West)" = "skyblue3", "Ethiopia"="#66A61E", "Cameroon"="#E6AB02", "Senegal"="#A6761D","Ivory Coast"="#666666")

poc$format <- "format"
ggplot(poc, aes(`PC1 (12.3%)`,`PC2 (12%)`, color=as.factor(Country))) + geom_point(size = 5) + scale_color_manual(values = pca_colors)
ggplot(poc, aes(`PC1 (12.3%)`,`PC2 (12%)`, color=as.factor(Country), label = IID)) + geom_text() + geom_point(size = 5) + scale_color_manual(values = pca_colors)
poc$`-PC1 (12.3%)` <- -1*poc$`PC1 (12.3%)`
poc$`-PC2 (12%)` <- -1*poc$`PC2 (12%)`
ggplot(poc, aes(`PC1 (12.3%)`,`PC2 (12%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. ovale curtisi samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = pca_colors)
ggplot(poc, aes(`-PC1 (12.3%)`,`-PC2 (12%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("P. ovale curtisi") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = pca_colors)
ggsave("PCA/poc_pc1-pc2.png", plot = last_plot(), device = "png", dpi = 600, width = 6, height = 5)
ggplot(poc, aes(`-PC1 (12.3%)`,`PC3 (10%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. ovale curtisi samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = pca_colors)
#PC3 shows differences between the Kinshasa samples
ggplot(poc, aes(`-PC1 (12.3%)`,`PC4 (9.9%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. ovale curtisi samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = pca_colors)
#PC4 shows differences between Tanzania, Cameroon, and Bas-Uele in DRC

###P. o. wallikeri
#Load in Pow datatables (including country of origin and eigenvalues)
pow <- read.table("P. o. pop genomics data tables - Pow PCA r2-0.3.tsv", header = TRUE, sep = "\t")
poweigenval <- read.table("P. o. pop genomics data tables - Pow Eigenvalues r2-0.3.tsv", header = FALSE, sep = "\t")
colnames(poweigenval) <- c("eigenvalues")
###Determine contribution of each principal component
pow_eigen_total <- sum(poweigenval) 
poweigenval$prop <- poweigenval$eigenvalues/pow_eigen_total
poweigenval$percent <- round(poweigenval$prop*100,1)
poweigenval$pc <- c(1:length(poweigenval$eigenvalues))
plot(poweigenval$pc, poweigenval$prop, type="o", col="red", xlab="PC", ylab="Proportion of variance explained", main="Scree plot I")
###Rename principal component variables to include amount of total variation explained
for (i in 1:length(poweigenval$eigenvalues)){
  colnames(pow)[i+3] <- paste(colnames(pow)[i+3]," (",poweigenval$percent[i],"%)", sep = "")
}

#PCA plots for geographical clustering from plink
library(ggplot2)


pow$format <- "format"
ggplot(pow, aes(`PC1 (12%)`,`PC2 (9.8%)`, color=as.factor(Country))) + geom_point(size = 5) + scale_color_manual(values = c("DRC (West)" = "#D95F02", "DRC (East)"= "#7570B3","Tanzania (East)"="#E7298A", "Tanzania (West)" = "skyblue3", "Ethiopia"="#66A61E", "Cameroon"="#E6AB02", "Senegal"="#A6761D","Ivory Coast"="#666666"))+ labs(color = "Country of origin") + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)
pow$`-PC2 (9.8%)` <- -1*pow$`PC2 (9.8%)`
ggplot(pow, aes(`PC1 (12%)`,`PC2 (9.8%)`, color=as.factor(Country),label=IID)) + geom_point(size = 5) + geom_text() + scale_color_manual(values = pca_colors) + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)
ggplot(pow, aes(`PC1 (12%)`,`PC2 (9.8%)`, color=as.factor(Country))) + geom_point(size = 5) + scale_color_manual(values = pca_colors)+ labs(color = "Country of origin") + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)+ scale_y_continuous(breaks = seq(-0.5, 0.75, by = 0.25))
ggsave("PCA/pow_pc1-pc2.png", plot = last_plot(), device = "png", dpi = 600, width = 6, height = 4.75)
#PC3 splits the ivory coast sample from the cameroon and senegal especially
ggplot(pow, aes(`PC1 (12%)`,`PC3 (9%)`, color=as.factor(Country))) + geom_point(size = 5) + scale_color_manual(values = pca_colors)+ labs(color = "Country of origin") + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)
ggplot(pow, aes(`PC1 (12%)`,`PC3 (9%)`, label=IID)) + geom_point(size = 5) + geom_text() + scale_color_manual(values = pca_colors) + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)
#PC4 splits Tanzania towards the top and some DRC (Kinshasa and Sud-Kivu) to the bottom
ggplot(pow, aes(`PC1 (12%)`,`PC4 (7.9%)`, color=as.factor(Country))) + geom_point(size = 5)+ scale_color_manual(values = pca_colors)+ labs(color = "Country of origin") + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)
ggplot(pow, aes(`PC1 (12%)`,`PC4 (7.9%)`, label=IID)) + geom_point(size = 5) + geom_text() + scale_color_manual(values = pca_colors) + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)

###arrange ovale plots next to each other
#install.packages("ggpubr")
library("ggpubr")

ggplot(poc, aes(`-PC1 (12.3%)`,`-PC2 (12%)`, color=factor(Country, levels = c("Ethiopia","Tanzania (West)", "Tanzania (East)","DRC (East)","DRC (West)","Cameroon")))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle(expression(paste(italic("P. ovale curtisi")))) + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = pca_colors, drop = FALSE)
ggsave("PCA/poc_pc1-pc2.png", plot = last_plot(), device = "png", dpi = 600, width = 6, height = 5)
ggplot(pow, aes(`PC1 (12%)`,`-PC2 (9.8%)`, color=factor(Country, levels = c("Ethiopia","Tanzania (West)", "Tanzania (East)","DRC (East)","DRC (West)","Cameroon","Ivory Coast","Senegal")))) + geom_point(size = 5) + scale_color_manual(values = pca_colors, drop = FALSE)+ labs(color = "Country of origin") + ggtitle(expression(paste(italic("P. ovale wallikeri")))) +theme_bw(base_size = 15)+ scale_y_continuous(breaks = seq(-0.5, 0.75, by = 0.25)) + scale_x_continuous(breaks = seq(-0.5, 0.75, by = 0.25))
ggsave("PCA/pow_pc1-pc2.png", plot = last_plot(), device = "png", dpi = 600, width = 6, height = 5)


poc_pc34_plot <- ggplot(poc, aes(`PC3 (10%)`,`PC4 (9.9%)`, color=factor(Country, levels = c("Ethiopia","Tanzania (East)", "Tanzania (West)","DRC (East)","DRC (West)","Cameroon","Ivory Coast","Senegal")))) + 
            geom_point(size = 5) + 
            theme_bw(base_size = 15) +
            labs(color = "Country of origin") + 
            theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 24),axis.text= element_text(size = 20), title = element_text(size = 20), legend.text = element_text(size = 20)) + 
            scale_color_manual(values = pca_colors, drop = FALSE) +
            scale_x_continuous(breaks = seq(-0.5, 0.75, by = 0.25))
            #+ ggtitle(expression(paste(italic("P. ovale curtisi")))) 
pow_pc34_plot <- ggplot(pow, aes(`PC3 (9%)`,`PC4 (7.9%)`, color=factor(Country, levels = c("Ethiopia","Tanzania (East)", "Tanzania (West)","DRC (East)","DRC (West)","Cameroon","Ivory Coast","Senegal")))) + 
            geom_point(size = 5) + 
            theme_bw(base_size = 15) + 
            theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 24),axis.text= element_text(size = 20), title = element_text(size = 20), legend.text = element_text(size = 20)) + 
            scale_color_manual(values = pca_colors, drop = FALSE) + 
            labs(color = "Country of origin") + 
            scale_y_continuous(breaks = seq(-0.5, 0.75, by = 0.25)) + 
            scale_x_continuous(breaks = seq(-0.5, 0.75, by = 0.25)) #+ ggtitle(expression(paste(italic("P. ovale wallikeri")))) +theme_bw(base_size = 15)
#grid.arrange(poc_plot, pow_plot, ncol=2)
ggarrange(poc_pc34_plot, pow_pc34_plot, ncol=2, nrow=1, common.legend = TRUE, legend="right", widths = c(5.1,4.75))
ggsave("PCA/poc-pow_pc3-pc4.png", plot = last_plot(), device = "png", dpi = 600, width = 11.5, height = 5)
poc_plot

poc_plot <- ggplot(poc, aes(`-PC1 (12.3%)`,`-PC2 (12%)`, color=factor(Country, levels = c("Ethiopia","Tanzania (East)", "Tanzania (West)","DRC (East)","DRC (West)","Cameroon","Ivory Coast","Senegal")))) + 
  geom_point(size = 5) + 
  theme_bw(base_size = 15) +
  labs(color = "Country of origin") + 
  theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 24),axis.text= element_text(size = 20), title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  scale_color_manual(values = pca_colors, drop = FALSE) +
  scale_x_continuous(breaks = seq(-0.5, 0.75, by = 0.25))

pow_plot <- ggplot(pow, aes(`PC1 (12%)`,`-PC2 (9.8%)`, color=factor(Country, levels = c("Ethiopia","Tanzania (East)", "Tanzania (West)","DRC (East)","DRC (West)","Cameroon","Ivory Coast","Senegal")))) + 
  geom_point(size = 5) + 
  theme_bw(base_size = 15) + 
  theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 24),axis.text= element_text(size = 20), title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  scale_color_manual(values = pca_colors, drop = FALSE) + 
  labs(color = "Country of origin") + 
  scale_y_continuous(breaks = seq(-0.5, 0.75, by = 0.25)) + 
  scale_x_continuous(breaks = seq(-0.5, 0.75, by = 0.25))
#Load in Pf datatables (including country of origin and eigenvalues)
####pdate with r2 0.3 pruned data
###pf <- read.table("P. o. pop genomics data tables - Pf PCA.tsv", header = TRUE, sep = "\t")
###pfeigenval <- read.table("P. o. pop genomics data tables - Pf Eigenvalues.tsv", header = FALSE, sep = "\t")
colnames(pfeigenval) <- c("eigenvalues")

### P. falciparum
###Determine contribution of each principal component
pf_eigen_total <- sum(pfeigenval) 
pfeigenval$prop <- pfeigenval$eigenvalues/pf_eigen_total
pfeigenval$percent <- round(pfeigenval$prop*100,1)
pfeigenval$pc <- c(1:length(pfeigenval$eigenvalues))
plot(pfeigenval$pc, pfeigenval$prop, type="o", col="red", xlab="PC", ylab="Proportion of variance explained", main="Scree plot I")
###Rename principal component variables to include amount of total variation explained
for (i in 1:length(pfeigenval$eigenvalues)){
  colnames(pf)[i+3] <- paste(colnames(pf)[i+3]," (",pfeigenval$percent[i],"%)", sep = "")
}

#PCA plots for geographical clustering from plink
library(ggplot2)

pf$format <- "format"
ggplot(pf, aes(`PC1 (14.5%)`,`PC2 (10.2%)`, color=as.factor(Country))) + geom_point(size = 5) + scale_color_manual(values = pca_colors)
ggplot(pf, aes(`PC1 (14.5%)`,`PC2 (10.2%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. falciparum samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = pca_colors)
ggplot(pf, aes(`PC1 (14.5%)`,`PC3 (7.5%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. falciparum samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = pca_colors)
#PC3 shows more differences in ethiopian samples, splits non-Ethiopian samples some
ggplot(pf, aes(`PC1 (14.5%)`,`PC4 (6.7%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. falciparum samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = pca_colors)
#PC4 shows more splitting on non-Ethiopian samples, though still not very starkly

###Plot contribution to each PC by SNP
###curtisi PC1
install.packages("qqman")
library("qqman")
#load table
####pdate with r2 0.3 pruned data
curtisi_eigenloci <- read.table("ov1_curtisigh01_speciescall.eigenvec.var",header = T, sep = "")
#isolate BP from VAR code
curtisi_eigenloci$BP <- as.numeric(substr(curtisi_eigenloci$VAR, 15, 24))
#isolate chromosome from var code
curtisi_eigenloci$CHROM <- substr(curtisi_eigenloci$VAR,1,13)
curtisi_eigenloci$CHR <- as.numeric(substr(curtisi_eigenloci$CHROM,9,10))
#Rename ID variable
curtisi_eigenloci$SNP <- curtisi_eigenloci$VAR
#Select principal component
curtisi_eigenloci$P <- curtisi_eigenloci$PC1
manhattan(curtisi_eigenloci, logp=FALSE, ylim = c(-4,4),col = c("blue4", "orange3"))
#determine cut-off for highlighting
pcval = 2.5
high_abs_pc1 <- quantile(abs(curtisi_eigenloci$PC1),0.995)
print(curtisi_eigenloci$VAR[abs(curtisi_eigenloci$PC1) > high_abs_pc1])
highlight <- curtisi_eigenloci$VAR[abs(curtisi_eigenloci$PC1) > high_abs_pc1]
manhattan(curtisi_eigenloci, logp=FALSE, ylim = c(-3,3),col = c("blue4", "orange3"),highlight = highlight,main = "Contribution of each SNP to PC1 (East-West) of P. o. curtisi")
print(curtisi_eigenloci[abs(curtisi_eigenloci$PC1) > high_abs_pc1, c("VAR","PC1")])

high_abs_pc2 <- quantile(abs(curtisi_eigenloci$PC2),0.995)
print(curtisi_eigenloci[abs(curtisi_eigenloci$PC2) >= high_abs_pc2, c("VAR","PC2")])

high_abs_pc3 <- quantile(abs(curtisi_eigenloci$PC3),0.995)
print(curtisi_eigenloci[abs(curtisi_eigenloci$PC3) >= high_abs_pc2, c("VAR","PC3")])

high_abs_pc4 <- quantile(abs(curtisi_eigenloci$PC4),0.995)
print(curtisi_eigenloci[abs(curtisi_eigenloci$PC4) >= high_abs_pc4, c("VAR","PC4")])

###wallikeri PC1
library("qqman")
#load table
####pdate with r2 0.3 pruned data
wallikeri_eigenloci <- read.table("ov1_wallikericr01_speciescall.eigenvec.var",header = T, sep = "")
#isolate BP from VAR code
wallikeri_eigenloci$BP <- as.numeric(substr(wallikeri_eigenloci$VAR, 10, 25))
#isolate chromosome from var code
wallikeri_eigenloci$CHROM <- substr(wallikeri_eigenloci$VAR,1,8)
wallikeri_eigenloci$CHR <- as.numeric(substr(wallikeri_eigenloci$CHROM,7,8))-4
#Rename ID variable
wallikeri_eigenloci$SNP <- wallikeri_eigenloci$VAR
#Select principal component
wallikeri_eigenloci$P <- wallikeri_eigenloci$PC1
manhattan(wallikeri_eigenloci,col = c("blue4", "orange3"),logp=FALSE, ylim = c(-3,3))
pcval = 2.75
high_abs_pc1 <- quantile(abs(wallikeri_eigenloci$PC1),0.995)
highlight <- wallikeri_eigenloci$VAR[abs(wallikeri_eigenloci$PC1) > high_abs_pc1]
manhattan(wallikeri_eigenloci,col = c("blue4", "orange3"),logp=FALSE, ylim = c(-3,3),highlight = highlight,main = "Contribution of each SNP to PC1 (East-West) of P. o. wallikeri")
print(wallikeri_eigenloci[abs(wallikeri_eigenloci$PC1) > high_abs_pc1, c("VAR","PC1")])

high_abs_pc2 <- quantile(abs(wallikeri_eigenloci$PC2),0.995)
print(wallikeri_eigenloci[abs(wallikeri_eigenloci$PC2) > high_abs_pc2, c("VAR","PC2")])

high_abs_pc3 <- quantile(abs(wallikeri_eigenloci$PC3),0.995)
print(wallikeri_eigenloci[abs(wallikeri_eigenloci$PC3) > high_abs_pc3, c("VAR","PC3")])

high_abs_pc4 <- quantile(abs(wallikeri_eigenloci$PC4),0.995)
print(wallikeri_eigenloci[abs(wallikeri_eigenloci$PC4) > high_abs_pc4, c("VAR","PC4")])

###Adegenet for dAPC
install.packages("adegenet")
library("adegenet")
