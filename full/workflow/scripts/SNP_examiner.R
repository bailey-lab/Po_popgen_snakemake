####Examine overall SNP calls
setwd("~/Documents/P. Ovale Genomic Analysis")

poc <- read.table("ov1_curtisigh01_speciescall.table", header = TRUE)
pow <- read.table("ov1_wallikericr01_speciescall.table", header = TRUE)

###Curtisi: Determine whether there are snps in certain key drug resistance genes
poc$CHR <- substr(poc$CHROM,9,10)
#dhfr: last SNP has abs(nSl) of 2.8 (in top 1%)
print(poc$POS[(poc$CHR == "05") & (762457 < poc$POS) & (poc$POS < 764370)])
#Chloroquine resistance transporter
print(poc$POS[(poc$CHR == "01") & (323118 < poc$POS) & (poc$POS < 326399)])
#chloroquine resistance associated protein
print(poc$POS[(poc$CHR == "01") & (343585 < poc$POS) & (poc$POS < 347363)])
#mdr1
print(poc$POS[(poc$CHR == "10") & (323034 < poc$POS) & (poc$POS < 327332)])
#k13: abs(nSl) of 2.3
print(poc$POS[(poc$CHR == "12") & (404824 < poc$POS) & (poc$POS < 407001)])
#dhps
print(poc$POS[(poc$CHR == "14") & (1119845 < poc$POS) & (poc$POS < 1122576)])
#cytochrome b: first SNP has nSL of 2.25
print(poc$POS[(poc$CHR == "14") & (2059533 < poc$POS) & (poc$POS < 2060640)])


###wallikeri: Determine whether there are snps in certain key drug resistance genes
pow$chr <- as.numeric(substr(pow$CHROM,7,8)) - 4

#dhfr LT594509:843,815..845,740, LT594516:1,183,447..1,185,319
print(pow$POS[(pow$chr == 5) & (843815 < pow$POS) & (pow$POS < 845740)])
print(pow$POS[(pow$chr == 12) & (1183447 < pow$POS) & (pow$POS < 1183447)])
#Chloroquine resistance transporter LT594505:331,639..334,593. One SNP had nSl of 2.805 (in top 1%)
print(pow$POS[(pow$chr == 1) & (331639 < pow$POS) & (pow$POS < 334593)])
#chloroquine resistance associated protein LT594505:352,724..356,765
print(pow$POS[(pow$chr == 1) & (352724 < pow$POS) & (pow$POS < 356765)])
#mdr1 LT594514:355,729..360,030
print(pow$POS[(pow$chr == 10) & (355729 < pow$POS) & (pow$POS < 360030)])
#k13 LT594516:455,780..457,957
print(pow$POS[(pow$chr == 12) & (455780 < pow$POS) & (pow$POS < 457957)])
#dhps LT594518:1,121,839..1,124,323
print(pow$POS[(pow$chr == 14) & (1121839 < pow$POS) & (pow$POS < 1124323)])
#cytocrhome b LT594518:2,083,698..2,084,367
print(pow$POS[(pow$chr == 14) & (2083698 < pow$POS) & (pow$POS < 2084367)])

