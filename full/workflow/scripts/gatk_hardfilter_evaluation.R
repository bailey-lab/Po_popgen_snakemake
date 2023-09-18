setwd("~/Downloads")

#set table name here
curt_only <- read.table("ov1_curtisigh01_only.table", header = TRUE)
wall_only <- read.table("ov1_wallikericr01_only.table", header = TRUE)
curt_andmixed <- read.table("ov1_curtisigh01_andmixed.table", header = TRUE)
wall_andmixed <- read.table("ov1_wallikericr01_andmixed.table", header = TRUE)
curt_all <- read.table("ov1_curtisigh01_all.table", header = TRUE)
wall_all <- read.table("ov1_wallikericr01_all.table", header = TRUE)


###Variant filtering

###Evaluate QD cutoff (quality and depth composite score)
#Set dataset and candidate filter threshold
x <- curt_only
qdthreshold <- 3
#default recommendation from GATK is 2.0

#Show initial distiribution
dd <- density(x$QD, na.rm = TRUE)
plot(dd, t = "l")
plot(dd)
polygon(dd, col = rgb(1, 0, 1, 0.5))

#Apply filter
qdfilter <- subset(x, QD <= qdthreshold)

ddfilter <- density(qdfilter$QD)
polygon(ddfilter, col = rgb(1, 0, 0, 0.5))

#Proportion of baseline variants remaining after filter
(dd$n-ddfilter$n)/dd$n


###Evaluate FS cutoff (fisher strand bias)
#Set dataset and candidate filter threshold
x <- wall_andmixed
fsthreshold <- 50
#GATK recommends 60

#Show initial distiribution
dd <- density(x$FS, na.rm = TRUE)
plot(dd, t = "l")
plot(dd)
polygon(dd, col = rgb(1, 0, 1, 0.5))

#Apply filter
fsfilter <- subset(x, FS >= fsthreshold)

ddfilter <- density(fsfilter$FS)
polygon(ddfilter, col = rgb(1, 0, 0, 0.5))

#Proportion of baseline variants remaining after filter
(dd$n-ddfilter$n)/dd$n


###Evaluate SOR cutoff (Strand Odds Ratio)
#Set dataset and candidate filter threshold
x <- wall_andmixed
sorthreshold <- 3
#GATK recommends 3.0

#Show initial distiribution
dd <- density(x$SOR, na.rm = TRUE)
plot(dd, t = "l")
plot(dd)
polygon(dd, col = rgb(1, 0, 1, 0.5))

#Apply filter
sorfilter <- subset(x, SOR > sorthreshold)

ddfilter <- density(sorfilter$SOR)
polygon(ddfilter, col = rgb(1, 0, 0, 0.5))

#Proportion of baseline variants remaining after filter
(dd$n-ddfilter$n)/dd$n

###Evaluate MQ cutoff (RMS Mapping Quality)
#Set dataset and candidate filter threshold
x <- wall_andmixed
mqthreshold <- 50
#GATK recommends 40

#Show initial distiribution
dd <- density(x$MQ, na.rm = TRUE)
plot(dd, t = "l")
plot(dd)
polygon(dd, col = rgb(1, 0, 1, 0.5))

#Apply filter
mqfilter <- subset(x, MQ < mqthreshold)

ddfilter <- density(mqfilter$MQ)
polygon(ddfilter, col = rgb(1, 0, 0, 0.5))

#Proportion of baseline variants remaining after filter
(dd$n-ddfilter$n)/dd$n

###Evaluate MQRankSum cutoff (Mapping Quality Rank Sum)
#Set dataset and candidate filter threshold
x <- wall_andmixed
mqrsthreshold <- -2.5
#GATK recommends -12.5

#Show initial distiribution
dd <- density(x$MQRankSum, na.rm = TRUE)
plot(dd, t = "l")
plot(dd)
polygon(dd, col = rgb(1, 0, 1, 0.5))

#Apply filter
mqrsfilter <- subset(x, MQRankSum < mqrsthreshold)

ddfilter <- density(mqrsfilter$MQRankSum)
polygon(ddfilter, col = rgb(1, 0, 0, 0.5))

#Proportion of baseline variants remaining after filter
(dd$n-ddfilter$n)/dd$n

###Evaluate ReadPosRankSum cutoff (Read Position Rank Sum)
#Set dataset and candidate filter threshold
x <- curt_only
rdrsthreshold <- -3
#GATK recommends -8.0

#Show initial distiribution
dd <- density(x$ReadPosRankSum, na.rm = TRUE)
plot(dd, t = "l")
plot(dd)
polygon(dd, col = rgb(1, 0, 1, 0.5))

#Apply filter
rdrsfilter <- subset(x, ReadPosRankSum < rdrsthreshold)

ddfilter <- density(rdrsfilter$ReadPosRankSum)
polygon(ddfilter, col = rgb(1, 0, 0, 0.5))

#Proportion of baseline variants remaining after filter
(dd$n-ddfilter$n)/dd$n

