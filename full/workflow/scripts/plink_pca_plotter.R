###Principal Components Analysis of P. o. genomic data
setwd("~/Documents/P. Ovale Genomic Analysis")

#Load in Poc datatables (including country of origin and eigenvalues)
poc <- read.table("P. o. pop genomics data tables - Poc PCA.tsv", header = TRUE, sep = "\t")
poceigenval <- read.table("P. o. pop genomics data tables - Poc Eigenvalues.tsv", header = FALSE, sep = "\t")
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

poc$format <- "format"
ggplot(poc, aes(`PC1 (12.9%)`,`PC2 (11.9%)`, color=as.factor(Country))) + geom_point(size = 5) + scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))
poc$`-PC1 (12.9%)` <- -1*poc$`PC1 (12.9%)`
ggplot(poc, aes(`PC1 (12.9%)`,`PC2 (11.9%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. ovale curtisi samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))
ggplot(poc, aes(`-PC1 (12.9%)`,`PC2 (11.9%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. ovale curtisi samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))
ggplot(poc, aes(`-PC1 (12.9%)`,`PC3 (10%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. ovale curtisi samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))
#PC3 shows differences between the Kinshasa samples
ggplot(poc, aes(`-PC1 (12.9%)`,`PC4 (9.7%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. ovale curtisi samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))
#PC4 shows differences between Tanzania, Cameroon, and Bas-Uele in DRC

###P. o. wallikeri
#Load in Pow datatables (including country of origin and eigenvalues)
pow <- read.table("P. o. pop genomics data tables - Pow PCA.tsv", header = TRUE, sep = "\t")
poweigenval <- read.table("P. o. pop genomics data tables - Pow Eigenvalues.tsv", header = FALSE, sep = "\t")
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
ggplot(pow, aes(`PC1 (14.1%)`,`PC2 (9.4%)`, color=as.factor(Country))) + geom_point(size = 5) + scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))+ labs(color = "Country of origin") + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)
ggplot(pow, aes(`PC1 (14.1%)`,`PC2 (9.4%)`, color=as.factor(Country),label=IID)) + geom_point(size = 5) + geom_text() + scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7")) + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)
#PC3 splits the ivory coast sample from the cameroon and senegal especially
ggplot(pow, aes(`PC1 (14.1%)`,`PC3 (8.1%)`, color=as.factor(Country))) + geom_point(size = 5) + scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))+ labs(color = "Country of origin") + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)
ggplot(pow, aes(`PC1 (14.1%)`,`PC3 (8.1%)`, label=IID)) + geom_point(size = 5) + geom_text() + scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7")) + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)
#PC4 splits Tanzania towards the top and some DRC (Kinshasa and Sud-Kivu) to the bottom
ggplot(pow, aes(`PC1 (14.1%)`,`PC4 (7.6%)`, color=as.factor(Country))) + geom_point(size = 5)+ scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))+ labs(color = "Country of origin") + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)
ggplot(pow, aes(`PC1 (14.1%)`,`PC4 (7.6%)`, label=IID)) + geom_point(size = 5) + geom_text() + scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7")) + ggtitle("PCA of P. ovale wallikeri samples") +theme_bw(base_size = 15)


#Load in Poc datatables (including country of origin and eigenvalues)
pf <- read.table("P. o. pop genomics data tables - Pf PCA.tsv", header = TRUE, sep = "\t")
pfeigenval <- read.table("P. o. pop genomics data tables - Pf Eigenvalues.tsv", header = FALSE, sep = "\t")
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
ggplot(pf, aes(`PC1 (14.5%)`,`PC2 (10.2%)`, color=as.factor(Country))) + geom_point(size = 5) + scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))
ggplot(pf, aes(`PC1 (14.5%)`,`PC2 (10.2%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. falciparum samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))
ggplot(pf, aes(`PC1 (14.5%)`,`PC3 (7.5%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. falciparum samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))
#PC3 shows more differences in ethiopian samples, splits non-Ethiopian samples some
ggplot(pf, aes(`PC1 (14.5%)`,`PC4 (6.7%)`, color=as.factor(Country))) + geom_point(size = 5) + labs(color = "Country of origin") + ggtitle("PCA of P. falciparum samples") + theme(plot.title = element_text(hjust =0.5, size =30),axis.title = element_text(size = 30),axis.text= element_text(size = 14), title = element_text(size = 20), legend.text = element_text(size = 15)) + theme_bw(base_size = 15)+ scale_color_manual(values = c("DRC (West)" = "#F8766D", "DRC (East)"= "#C49A00","Tanzania"="#00C094", "Ethiopia"="#53B400", "Cameroon"="#00B6EB", "Senegal"="#A58AFF","Ivory Coast"="#FB61D7"))
#PC4 shows more splitting on non-Ethiopian samples, though still not very starkly

###Plot contribution to each PC by SNP
###curtisi PC1
library("qqman")
#load table
curtisi_eigenloci <- read.table("ov1_curtisigh01_speciescall.eigenvec.var",header = T, sep = "")
#isolate BP from VAR code
curtisi_eigenloci$BP <- as.numeric(substr(curtisi_eigenloci$VAR, 15, 24))
#isolate chromosome from var code
curtisi_eigenloci$CHROM <- curtisi_eigenloci$CHR
curtisi_eigenloci$CHR <- as.numeric(substr(curtisi_eigenloci$CHROM,9,10))
#Rename ID variable
curtisi_eigenloci$SNP <- curtisi_eigenloci$VAR
#Select principal component
curtisi_eigenloci$P <- curtisi_eigenloci$PC1
manhattan(curtisi_eigenloci, logp=FALSE, ylim = c(-4,4),col = c("blue4", "orange3"))
#determine cut-off for highlighting
pcval = 2.5
print(curtisi_eigenloci$VAR[abs(curtisi_eigenloci$PC1) > pcval])
highlight <- curtisi_eigenloci$VAR[abs(curtisi_eigenloci$PC1) > pcval]
manhattan(curtisi_eigenloci, logp=FALSE, ylim = c(-3,3),col = c("blue4", "orange3"),highlight = highlight,main = "Contribution of each SNP to PC1 (East-West) of P. o. curtisi")
print(curtisi_eigenloci$VAR[abs(curtisi_eigenloci$PC1) > 2.5])
print(curtisi_eigenloci$PC1[abs(curtisi_eigenloci$PC1) > 2.5])


###wallikeri PC1
library("qqman")
#load table
wallikeri_eigenloci <- read.table("ov1_wallikericr01_speciescall.eigenvec.var",header = T, sep = "")
#isolate BP from VAR code
wallikeri_eigenloci$BP <- as.numeric(substr(wallikeri_eigenloci$VAR, 10, 25))
#isolate chromosome from var code
wallikeri_eigenloci$CHROM <- wallikeri_eigenloci$CHR
wallikeri_eigenloci$CHR <- as.numeric(substr(wallikeri_eigenloci$CHROM,7,8))-4
#Rename ID variable
wallikeri_eigenloci$SNP <- wallikeri_eigenloci$VAR
#Select principal component
wallikeri_eigenloci$P <- wallikeri_eigenloci$PC1
manhattan(wallikeri_eigenloci,col = c("blue4", "orange3"),logp=FALSE, ylim = c(-3,3))
pcval = 2.75
highlight <- wallikeri_eigenloci$VAR[abs(wallikeri_eigenloci$PC1) > pcval]
manhattan(wallikeri_eigenloci,col = c("blue4", "orange3"),logp=FALSE, ylim = c(-3,3),highlight = highlight,main = "Contribution of each SNP to PC1 (East-West) of P. o. wallikeri")
pow_sig_loci <- wallikeri_eigenloci$VAR[abs(wallikeri_eigenloci$PC1) > pcval]
print(wallikeri_eigenloci$PC1[abs(wallikeri_eigenloci$PC1) > pcval])
