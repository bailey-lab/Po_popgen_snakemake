### Examine extended haplotype homozygosity at specified genomic regions deteremined by selection analysis

library(rehh)
library(vcfR)
setwd("~/Documents/P. Ovale Genomic Analysis")

#examine EHH on Poc chromosome 5
poc_chr5 <- data2haplohh(hap_file = "ov1_curtisigh01_speciescall_monoclonal_allmaf_nomissing_PocGH01_05_v2.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
#search poc_chr5@positions to find the number of the right marker
#calculate EHH for hit near TS-DHFR
poc_chr5_dhfr <- calc_ehh(poc_chr5, mrk = 2059)
png(file ="ehh/poc_dhfr.png", res = 300, width = 2000, height = 1500)
par(mar = c(4.2,4.5,1.5,1.5))
plot(poc_chr5_dhfr, main = "", 
     xlim = c(760000,790000), 
     type = "l", 
     col = c("red","#276FBF"), 
     xlab = "Position on Poc Chromosome 5 (Kb)", 
     ylab = "EHH", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     lwd = 3, legend = "")
legend("topright", c("selected allele","unselected allele"), col = c("red","#276FBF"), pch = 20, cex = 1.5)
dev.off()

poc_chr5_dhfr_furcation <- calc_furcation(poc_chr5, mrk = 2059)
png(file ="ehh/poc_dhfr_furcation.png", res = 300, width = 2000, height = 1500)
par(mar = c(4.2,0,0,0))
plot(poc_chr5_dhfr_furcation, main = "", 
     xlim = c(763000,790000), 
     type = "l", 
     col = c("red","#276FBF"), 
     xlab = "Position on Poc Chromosome 5 (Kb)",
     cex.lab = 4, 
     cex.axis = 1.5, 
     lwd = 1, legend = NULL)
#legend("top", c("selected allele","unselected allele"), col = c("red","#276FBF"), pch = 20, cex = 1.5)
dev.off()

pow_chr5 <- data2haplohh(hap_file = "ov1_wallikericr01_speciescall_monoclonal_allmaf_nomissing_LT594509.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
#search pow_chr5@positions to find the number of the right marker
#Calculate EHH for hit near TS-DHFR
pow_chr5_dhfr <- calc_ehh(pow_chr5, mrk = 1057)
png(file ="ehh/pow_dhfr.png", res = 300, width = 2000, height = 1500)
par(mar = c(4.2,4.5,1.5,1.5))
plot(pow_chr5_dhfr, main = "", 
     xlim = c(823000,870000), 
     type = "l", 
     col = c("red","#276FBF"), 
     xlab = "Position on Pow Chromosome 5 (Kb)", 
     ylab = "EHH", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     lwd = 3, legend = "")
legend("topright", c("selected allele","unselected allele"), col = c("red","#276FBF"), pch = 20, cex = 1.5)
dev.off()

pow_chr5_dhfr_furcation <- calc_furcation(pow_chr5, mrk = 1057)
png(file ="ehh/pow_dhfr_furcation.png", res = 300, width = 2000, height = 1500)
par(mar = c(4.2,0,0,0))
plot(pow_chr5_dhfr_furcation, main = "", 
     xlim = c(823000,870000), 
     type = "l", 
     col = c("red","#276FBF"), 
     xlab = "Position on Pow Chromosome 5 (Kb)",
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     lwd = 1, legend = NULL)
#legend("top", c("selected allele","unselected allele"), col = c("red","#276FBF"), pch = 20, cex = 1.5)
dev.off()

#examine EHH on Poc chromosome 14
poc_chr14 <- data2haplohh(hap_file = "ov1_curtisigh01_speciescall_monoclonal_allmaf_nomissing_PocGH01_14_v2.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
#search poc_chr14@positions to find the number of the right marker
poc_chr14_mrp2 <- calc_ehh(poc_chr14, mrk = 5495)
png(file ="ehh/poc_mrp2.png", res = 300, width = 2000, height = 1500)
par(mar = c(4.2,4.5,1.5,1.5))
plot(poc_chr14_mrp2, main = "", 
     #xlim = c(1880000,1920000), 
     type = "l", 
     col = c("red","#276FBF"), 
     xlab = "Position on Poc Chromosome 14 (Mb)", 
     ylab = "EHH", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     lwd = 3, legend = "")
legend("topright", c("selected allele","unselected allele"), col = c("red","#276FBF"), pch = 20, cex = 1.5)
dev.off()

poc_chr14_mrp2_furcation <- calc_furcation(poc_chr14, mrk = 5495)
png(file ="ehh/poc_mrp2_furcation.png", res = 300, width = 2000, height = 1500)
par(mar = c(4.2,0,0,0))
plot(poc_chr14_mrp2_furcation, main = "", 
     xlim = c(1880000,1930000),
     type = "l", 
     col = c("red","#276FBF"), 
     xlab = "Position on Poc Chromosome 14 (Mb)",
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     lwd = 1, legend = NULL)
#legend("topleft", c("selected allele","unselected allele"), col = c("red","#276FBF"), pch = 20, cex = 1.5)
dev.off()

#calculate EHH for hit near FD1
poc_chr14_fd1 <- calc_ehh(poc_chr14, mrk = 6980)
plot(poc_chr14_fd1, main = "chr05 - FD1 - Poc")
  
#examine EHH on Pow chromosome 13
pow_chr13 <- data2haplohh(hap_file = "ov1_wallikericr01_speciescall_monoclonal_allmaf_nomissing_LT594517.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
pow_chr13_pir <- calc_ehh(pow_chr13, mrk = 128)
#mrks 128
pow_chr13_pir <- calc_ehh(pow_chr13, mrk = 143)
png(file ="ehh/pow_pir.png", res = 300, width = 2000, height = 1500)
par(mar = c(4.2,4.5,1.5,1.5))
plot(pow_chr13_pir, main = "", 
     xlim = c(138000,183000), 
     type = "l", 
     col = c("red","#276FBF"), 
     xlab = "Position on Poc Chromosome 14 (Kb)", 
     ylab = "EHH", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     lwd = 3, legend = "")
#drawBox(x = 170000, y=0, width = 1000, height= 1)
rect(170.74,0,171.95,1,density = 15, col = "gray", border = "gray")
rect(166.967,0,168.228,1,density = 15, col = "gray", border = "gray")
rect(157.979,0,159.240,1,density = 15, col = "gray", border = "gray")
rect(152.862,0,154.132,1,density = 15, col = "gray", border = "gray")
rect(149.062,0,150.275,1,density = 15, col = "gray", border = "gray")
rect(144.651,0,145.928,1,density = 15, col = "gray", border = "gray")
rect(140.395,0,141.604,1,density = 15, col = "gray", border = "gray")
rect(182.722,0,183.075,1,density = 15, col = "orange3", border = "orange3")
legend("topleft", c("selected allele","unselected allele", "masked PIR gene", "masked ETRAMP gene"), col = c("red","#276FBF","grey","orange3"), pch = 20, cex = 1.5)
dev.off()


pow_chr13_furcation <- calc_furcation(pow_chr13, mrk = 143)
png(file ="ehh/pow_pir_furcation.png", res = 300, width = 2000, height = 1500)
par(mar = c(4.2,0,0,0))
plot(pow_chr13_furcation, main = "", 
     xlim = c(150000,190000),
     type = "l", 
     col = c("red","#276FBF"), 
     xlab = "Position on Pow Chromosome 13 (Kb)",
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     lwd = 1, legend = NULL)
# rect(170.74,0,171.95,1,density = 15, col = "gray", border = "gray")
# rect(166.967,0,168.228,1,density = 15, col = "gray", border = "gray")
# rect(157.979,0,159.240,1,density = 15, col = "gray", border = "gray")
# rect(152.862,0,154.132,1,density = 15, col = "gray", border = "gray")
# rect(149.062,0,150.275,1,density = 15, col = "gray", border = "gray")
# rect(144.651,0,145.928,1,density = 15, col = "gray", border = "gray")
# rect(140.395,0,141.604,1,density = 15, col = "gray", border = "gray")
# rect(182.722,0,183.075,1,density = 15, col = "orange3", border = "orange3")

#legend("topleft", c("selected allele","unselected allele"), col = c("red","#276FBF"), pch = 20, cex = 1.5)
dev.off()

View(poc)
chr02_relapse_scan <- scan_hh(chr02_relapse)
chr02_background_scan <- scan_hh(chr02_background)
chr02_xpehh <- ies2xpehh(scan_pop1 = chr02_relapse_scan,
                         scan_pop2 = chr02_background_scan,
                         popname1 = "Relapse",
                         popname2 = "Background")
CHR <- c("PvP01_02_v1")
POSITION <- c("686070")
chr02_markers <- data.frame(CHR = c("PvP01_02_v1"),POSITION = c("686070"), ID = c("ApiAP2"))
manhattanplot(chr02_xpehh, mrk = chr02_markers, mrk.col = "RED", mrk.pch = 19, mrk.cex = 3)
chr02_ApiAP2_res <- calc_ehh(chr02_relapse, mrk = 12067)
plot(chr02_ApiAP2_res, main = "chr02 - ApiAP2 - Relapse")


ama1plot <- ggplot(ama1§Tajima, aes(x = BIN§START, y = TajimaD)) + 
  geom§point(aes(color = TajimaD > abs(2)), show.legend = FALSE) +
  scale§color§manual(values = c("black", "red")) +
  theme§classic() +
  ylim(-3,3) +
  labs(y= "Tajima's D", x = "Position") +
  geom§hline(yintercept = 0) +
  geom§hline(yintercept = 2, linetype = "dashed", color = "blue") +
  geom§hline(yintercept = -2, linetype = "dashed", color = "blue") +
  annotate("rect", xmin=AMA1§missing§coords$Start, xmax = AMA1§missing§coords$Stop, ymin = -3, ymax = 3, alpha = 0.2) +
  scale§x§continuous(label = comma) +
  ggtitle(expression(paste(italic("Pvama1")))) +
  theme(axis.text.x = element§text(angle = 45, vjust = 1, hjust = 1))

chr02_ApiAP2_furcation <- calc_furcation(chr02_relapse, mrk = 12067)
plot(chr02_ApiAP2_furcation, main = "chr02 - ApiAP2 - Relapse", xlim = c(685000,687000), hap.names = hap.names(chr02_relapse), cex.lab = 0.3)



chr02_ApiAP2_res2 <- calc_ehh(chr02_background, mrk = 12067)
plot(chr02_ApiAP2_res2, main = "chr02 - ApiAP2 - Background")
chr02_ApiAP2_furcation2 <- calc_furcation(chr02_background, mrk = 12067)
plot(chr02_ApiAP2_furcation2, main = "chr02 - ApiAP2 - Background", xlim = c(685000,687000), hap.names = hap.names(chr02_background), cex.lab = 0.3)

chr07_relapse <- data2haplohh(hap_file = "Pvivax_PvP01_07_v1_new_mask.recode_relapses_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr07_background <- data2haplohh(hap_file = "Pvivax_PvP01_07_v1_new_mask.recode_background_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr07_relapse_scan <- scan_hh(chr07_relapse)
chr07_background_scan <- scan_hh(chr07_background)
chr07_xpehh <- ies2xpehh(scan_pop1 = chr07_relapse_scan,
                         scan_pop2 = chr07_background_scan,
                         popname1 = "Relapse",
                         popname2 = "Background")
chr07_markers <- data.frame(CHR = c("PvP01_07_v1", "PvP01_07_v1", "PvP01_07_v1", "PvP01_07_v1"), POSITION = c("1395879", "1395885", "1395891", "1395897"), ID = c("GDV1", "GDV1", "GDV1", "GDV1"))
manhattanplot(chr07_xpehh, mrk = chr07_markers, mrk.col = "RED", mrk.pch = 19, mrk.cex = 3)
GDVI_1_res <- calc_ehh(chr07_relapse, mrk = 27363)
plot(GDVI_1_res, main = "GDVI SNP1 Relapse")
GDVI_1_furcation <- calc_furcation(chr07_relapse, mrk = 27363)
plot(GDVI_1_furcation, main = "GDVI SNP1 Relapse", hap.names = hap.names(chr07_relapse), cex.lab = 0.3)
GDVI_1_res2 <- calc_ehh(chr07_background, mrk = 27363)
plot(GDVI_1_res2, main = "GDVI SNP1 Background")
GDVI_1_furcation2 <- calc_furcation(chr07_background, mrk = 27363)
plot(GDVI_1_furcation2, main = "GDVI SNP1 Background", hap.names = hap.names(chr07_background), xlim = c(1394000, 1397000), cex.lab = 0.3)
GDVI_2_res <- calc_ehh(chr07_relapse, mrk = 27364)
plot(GDVI_2_res, main = "GDVI SNP2 Relapse")
GDVI_2_furcation <- calc_furcation(chr07_relapse, mrk = 27364)
plot(GDVI_2_furcation, main = "GDVI SNP2 Relapse", hap.names = hap.names(chr07_relapse), cex.lab = 0.3)
GDVI_2_res2 <- calc_ehh(chr07_background, mrk = 27364)
plot(GDVI_2_res2, main = "GDVI SNP2 Background")
GDVI_2_furcation2 <- calc_furcation(chr07_background, mrk = 27364)
plot(GDVI_2_furcation2, main = "GDVI SNP2 Background", hap.names = hap.names(chr07_background), xlim = c(1394000, 1397000), cex.lab = 0.3)
GDVI_3_res <- calc_ehh(chr07_relapse, mrk = 27366)
plot(GDVI_3_res, main = "GDVI SNP3 Relapse")
GDVI_3_furcation <- calc_furcation(chr07_relapse, mrk = 27366)
plot(GDVI_3_furcation, main = "GDVI SNP3 Relapse", hap.names = hap.names(chr07_relapse), cex.lab = 0.3)
GDVI_3_res2 <- calc_ehh(chr07_background, mrk = 27366)
plot(GDVI_3_res2, main = "GDVI SNP3 Background")
GDVI_3_furcation2 <- calc_furcation(chr07_background, mrk = 27366)
plot(GDVI_3_furcation2, main = "GDVI SNP3 Background", hap.names = hap.names(chr07_background), xlim = c(1394000, 1397000), cex.lab = 0.3)
GDVI_4_res <- calc_ehh(chr07_relapse, mrk = 27369)
plot(GDVI_4_res, main = "GDVI SNP4 Relapse")
GDVI_4_furcation <- calc_furcation(chr07_relapse, mrk = 27369)
plot(GDVI_4_furcation, main = "GDVI SNP4 Relapse", hap.names = hap.names(chr07_relapse), cex.lab = 0.3)
GDVI_4_res2 <- calc_ehh(chr07_background, mrk = 27369)
plot(GDVI_4_res2, main = "GDVI SNP4 Background")
GDVI_4_furcation2 <- calc_furcation(chr07_background, mrk = 27369)
plot(GDVI_4_furcation2, main = "GDVI SNP4 Background", hap.names = hap.names(chr07_background), xlim = c(1394000, 1397000), cex.lab = 0.3)

chr08_relapse <- data2haplohh(hap_file = "Pvivax_PvP01_08_v1_new_mask.recode_relapses_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr08_background <- data2haplohh(hap_file = "Pvivax_PvP01_08_v1_new_mask.recode_background_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr08_relapse_scan <- scan_hh(chr08_relapse)
chr08_background_scan <- scan_hh(chr08_background)
chr08_xpehh <- ies2xpehh(scan_pop1 = chr08_relapse_scan,
                         scan_pop2 = chr08_background_scan,
                         popname1 = "Relapse",
                         popname2 = "Background")
chr08_markers <- data.frame(CHR = c("PvP01_08_v1", "PvP01_08_v1", "PvP01_08_v1"), POSITION = c("980636", "980638", "1046302"), ID = c("STOML", "STOML", "YOP1"))
manhattanplot(chr08_xpehh, mrk = chr08_markers, mrk.col = "RED", mrk.pch = 19, mrk.cex = 3)
STOML_1_res <- calc_ehh(chr08_relapse, mrk = 18536)
plot(STOML_1_res, main = "STOML SNP1 Relapse")
STOML_1_furcation <- calc_furcation(chr08_relapse, mrk = 18536)
plot(STOML_1_furcation, main = "STOML SNP1 Relapse", hap.names = hap.names(chr08_relapse), cex.lab = 0.3)
STOML_1_res2 <- calc_ehh(chr08_background, mrk = 18536)
plot(STOML_1_res2, main = "STOML SNP1 Background")
STOML_1_furcation2 <- calc_furcation(chr08_background, mrk = 18536)
plot(STOML_1_furcation2, main = "STOML SNP1 Background", hap.names = hap.names(chr08_background), xlim = c(975000, 1000000), cex.lab = 0.3)
STOML_2_res <- calc_ehh(chr08_relapse, mrk = 18537)
plot(STOML_2_res, main = "STOML SNP2 Relapse")
STOML_2_furcation <- calc_furcation(chr08_relapse, mrk = 18537)
plot(STOML_2_furcation, main = "STOML SNP2 Relapse", hap.names = hap.names(chr08_relapse), cex.lab = 0.3)
STOML_2_res2 <- calc_ehh(chr08_background, mrk = 18537)
plot(STOML_2_res2, main = "STOML SNP2 Background")
STOML_2_furcation2 <- calc_furcation(chr08_background, mrk = 18537)
plot(STOML_2_furcation2, main = "STOML SNP2 Background", hap.names = hap.names(chr08_background), xlim = c(975000, 1000000), cex.lab = 0.3)
YOP1_res <- calc_ehh(chr08_relapse, mrk = 20010)
plot(YOP1_res, main = "YOP1 Relapse")
YOP1_furcation <- calc_furcation(chr08_relapse, mrk = 20010)
plot(YOP1_furcation, main = "YOP1 Relapse", hap.names = hap.names(chr08_relapse), cex.lab = 0.3)
YOP1_res2 <- calc_ehh(chr08_background, mrk = 20010)
plot(YOP1_res2, main = "YOP1 Background")
YOP1_furcation2 <- calc_furcation(chr08_background, mrk = 20010)
plot(YOP1_furcation2, main = "YOP1 Background", hap.names = hap.names(chr08_background), xlim = c(1040000, 1120000), cex.lab = 0.3)

chr10_relapse <- data2haplohh(hap_file = "Pvivax_PvP01_10_v1_new_mask.recode_relapses_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr10_background <- data2haplohh(hap_file = "Pvivax_PvP01_10_v1_new_mask.recode_background_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr10_relapse_scan <- scan_hh(chr10_relapse)
chr10_background_scan <- scan_hh(chr10_background)
chr10_xpehh <- ies2xpehh(scan_pop1 = chr10_relapse_scan,
                         scan_pop2 = chr10_background_scan,
                         popname1 = "Relapse",
                         popname2 = "Background")
chr10_markers <- data.frame(CHR = c("PvP01_10_v1", "PvP01_10_v1", "PvP01_10_v1", "PvP01_10_v1","PvP01_10_v1", "PvP01_10_v1", "PvP01_10_v1", "PvP01_10_v1","PvP01_10_v1"), POSITION = c("1355309",  "1355321",  "1355322",  "1355338",  "1355343", "1355355",  "1355683", "1355795", "1359950"), ID = c("MSP3.3", "MSP3.3", "MSP3.3", "MSP3.3", "MSP3.3", "MSP3.3", "MSP3.3", "MSP3.3", "MSP3.3"))
manhattanplot(chr10_xpehh, mrk = chr10_markers, mrk.col = "RED", mrk.pch = 19, mrk.cex = 3)
MSP3.3_1_res <- calc_ehh(chr10_relapse, mrk = 59732)
plot(MSP3.3_1_res, main = "MSP3.3 SNP1 Relapse")
MSP3.3_1_furcation <- calc_furcation(chr10_relapse, mrk = 59732)
plot(MSP3.3_1_furcation, main = "MSP3.3 SNP1 Relapse", hap.names = hap.names(chr10_relapse), cex.lab = 0.3)
MSP3.3_1_res2 <- calc_ehh(chr10_background, mrk = 59732)
plot(MSP3.3_1_res2, main = "MSP3.3 SNP1 Background")
MSP3.3_1_furcation2 <- calc_furcation(chr10_background, mrk = 59732)
plot(MSP3.3_1_furcation2, main = "MSP3.3 SNP1 Background", hap.names = hap.names(chr10_background), xlim = c(1345000, 1365000), cex.lab = 0.3)
MSP3.3_2_res <- calc_ehh(chr10_relapse, mrk = 59734)
plot(MSP3.3_2_res, main = "MSP3.3 SNP2 Relapse")
MSP3.3_2_furcation <- calc_furcation(chr10_relapse, mrk = 59734)
plot(MSP3.3_2_furcation, main = "MSP3.3 SNP2 Relapse", hap.names = hap.names(chr10_relapse), cex.lab = 0.3)
MSP3.3_2_res2 <- calc_ehh(chr10_background, mrk = 59734)
plot(MSP3.3_2_res2, main = "MSP3.3 SNP2 Background")
MSP3.3_2_furcation2 <- calc_furcation(chr10_background, mrk = 59734)
plot(MSP3.3_2_furcation2, main = "MSP3.3 SNP2 Background", hap.names = hap.names(chr10_background), xlim = c(1345000, 1365000), cex.lab = 0.3)
MSP3.3_3_res <- calc_ehh(chr10_relapse, mrk = 59735)
plot(MSP3.3_3_res, main = "MSP3.3 SNP3 Relapse")
MSP3.3_3_furcation <- calc_furcation(chr10_relapse, mrk = 59735)
plot(MSP3.3_3_furcation, main = "MSP3.3 SNP3 Relapse", hap.names = hap.names(chr10_relapse), cex.lab = 0.3)
MSP3.3_3_res2 <- calc_ehh(chr10_background, mrk = 59735)
plot(MSP3.3_3_res2, main = "MSP3.3 SNP3 Background")
MSP3.3_3_furcation2 <- calc_furcation(chr10_background, mrk = 59735)
plot(MSP3.3_3_furcation2, main = "MSP3.3 SNP3 Background", hap.names = hap.names(chr10_background), xlim = c(1345000, 1365000), cex.lab = 0.3)
MSP3.3_4_res <- calc_ehh(chr10_relapse, mrk = 59738)
plot(MSP3.3_4_res, main = "MSP3.3 SNP4 Relapse")
MSP3.3_4_furcation <- calc_furcation(chr10_relapse, mrk = 59738)
plot(MSP3.3_4_furcation, main = "MSP3.3 SNP4 Relapse", hap.names = hap.names(chr10_relapse), cex.lab = 0.3)
MSP3.3_4_res2 <- calc_ehh(chr10_background, mrk = 59738)
plot(MSP3.3_4_res2, main = "MSP3.3 SNP4 Background")
MSP3.3_4_furcation2 <- calc_furcation(chr10_background, mrk = 59738)
plot(MSP3.3_4_furcation2, main = "MSP3.3 SNP4 Background", hap.names = hap.names(chr10_background), xlim = c(1345000, 1365000), cex.lab = 0.3)
MSP3.3_5_res <- calc_ehh(chr10_relapse, mrk = 59739)
plot(MSP3.3_5_res, main = "MSP3.3 SNP5 Relapse")
MSP3.3_5_furcation <- calc_furcation(chr10_relapse, mrk = 59739)
plot(MSP3.3_5_furcation, main = "MSP3.3 SNP5 Relapse", hap.names = hap.names(chr10_relapse), cex.lab = 0.3)
MSP3.3_5_res2 <- calc_ehh(chr10_background, mrk = 59739)
plot(MSP3.3_5_res2, main = "MSP3.3 SNP5 Background")
MSP3.3_5_furcation2 <- calc_furcation(chr10_background, mrk = 59739)
plot(MSP3.3_5_furcation2, main = "MSP3.3 SNP5 Background", hap.names = hap.names(chr10_background), xlim = c(1345000, 1365000), cex.lab = 0.3)
MSP3.3_6_res <- calc_ehh(chr10_relapse, mrk = 59740)
plot(MSP3.3_6_res, main = "MSP3.3 SNP6 Relapse")
MSP3.3_6_furcation <- calc_furcation(chr10_relapse, mrk = 59740)
plot(MSP3.3_6_furcation, main = "MSP3.3 SNP6 Relapse", hap.names = hap.names(chr10_relapse), cex.lab = 0.3)
MSP3.3_6_res2 <- calc_ehh(chr10_background, mrk = 59740)
plot(MSP3.3_6_res2, main = "MSP3.3 SNP6 Background")
MSP3.3_6_furcation2 <- calc_furcation(chr10_background, mrk = 59740)
plot(MSP3.3_6_furcation2, main = "MSP3.3 SNP6 Background", hap.names = hap.names(chr10_background), xlim = c(1345000, 1365000), cex.lab = 0.3)
MSP3.3_7_res <- calc_ehh(chr10_relapse, mrk = 59784)
plot(MSP3.3_7_res, main = "MSP3.3 SNP7 Relapse")
MSP3.3_7_furcation <- calc_furcation(chr10_relapse, mrk = 59784)
plot(MSP3.3_7_furcation, main = "MSP3.3 SNP7 Relapse", hap.names = hap.names(chr10_relapse), cex.lab = 0.3)
MSP3.3_7_res2 <- calc_ehh(chr10_background, mrk = 59784)
plot(MSP3.3_7_res2, main = "MSP3.3 SNP7 Background")
MSP3.3_7_furcation2 <- calc_furcation(chr10_background, mrk = 59784)
plot(MSP3.3_7_furcation2, main = "MSP3.3 SNP7 Background", hap.names = hap.names(chr10_background), xlim = c(1345000, 1365000), cex.lab = 0.3)
MSP3.3_8_res <- calc_ehh(chr10_relapse, mrk = 59819)
plot(MSP3.3_8_res, main = "MSP3.3 SNP8 Relapse")
MSP3.3_8_furcation <- calc_furcation(chr10_relapse, mrk = 59819)
plot(MSP3.3_8_furcation, main = "MSP3.3 SNP8 Relapse", hap.names = hap.names(chr10_relapse), cex.lab = 0.3)
MSP3.3_8_res2 <- calc_ehh(chr10_background, mrk = 59819)
plot(MSP3.3_8_res2, main = "MSP3.3 SNP8 Background")
MSP3.3_8_furcation2 <- calc_furcation(chr10_background, mrk = 59819)
plot(MSP3.3_8_furcation2, main = "MSP3.3 SNP8 Background", hap.names = hap.names(chr10_background), xlim = c(1345000, 1365000), cex.lab = 0.3)
MSP3.3_9_res <- calc_ehh(chr10_relapse, mrk = 60486)
plot(MSP3.3_9_res, main = "MSP3.3 SNP9 Relapse")
MSP3.3_9_furcation <- calc_furcation(chr10_relapse, mrk = 60486)
plot(MSP3.3_9_furcation, main = "MSP3.3 SNP9 Relapse", hap.names = hap.names(chr10_relapse), cex.lab = 0.3)
MSP3.3_9_res2 <- calc_ehh(chr10_background, mrk = 60486)
plot(MSP3.3_9_res2, main = "MSP3.3 SNP9 Background")
MSP3.3_9_furcation2 <- calc_furcation(chr10_background, mrk = 60486)
plot(MSP3.3_9_furcation2, main = "MSP3.3 SNP9 Background", hap.names = hap.names(chr10_background), xlim = c(1345000, 1365000), cex.lab = 0.3)

chr11_relapse <- data2haplohh(hap_file = "Pvivax_PvP01_11_v1_new_mask.recode_relapses_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr11_background <- data2haplohh(hap_file = "Pvivax_PvP01_11_v1_new_mask.recode_background_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr11_relapse_scan <- scan_hh(chr11_relapse)
chr11_background_scan <- scan_hh(chr11_background)
chr11_xpehh <- ies2xpehh(scan_pop1 = chr11_relapse_scan,
                         scan_pop2 = chr11_background_scan,
                         popname1 = "Relapse",
                         popname2 = "Background")
chr11_markers <- data.frame(CHR = c("PvP01_11_v1", "PvP01_11_v1", "PvP01_11_v1", "PvP01_11_v1"), POSITION = c("1266053",  "1686756",  "1687370",  "1926839"), ID = c("PSOP25", "ZIP1", "ZIP1", "CEP76"))
manhattanplot(chr11_xpehh, mrk = chr11_markers, mrk.col = "RED", mrk.pch = 19, mrk.cex = 3)
PSOP25_res <- calc_ehh(chr11_relapse, mrk = 21261)
plot(PSOP25_res, main = "PSOP25 Relapse")
PSOP25_furcation <- calc_furcation(chr11_relapse, mrk = 21261)
plot(PSOP25_furcation, main = "PSOP25 Relapse", hap.names = hap.names(chr11_relapse), cex.lab = 0.3)
PSOP25_res2 <- calc_ehh(chr11_background, mrk = 21261)
plot(PSOP25_res2, main = "PSOP25 Background")
PSOP25_furcation2 <- calc_furcation(chr11_background, mrk = 21261)
plot(PSOP25_furcation2, main = "PSOP25 Background", hap.names = hap.names(chr11_background), xlim = c(1180000, 1320000), cex.lab = 0.3)
ZIP1_1_res <- calc_ehh(chr11_relapse, mrk = 29990)
plot(ZIP1_1_res, main = "ZIP1 SNP 1 Relapse")
ZIP1_1_furcation <- calc_furcation(chr11_relapse, mrk = 29990)
plot(ZIP1_1_furcation, main = "ZIP1 SNP 1 Relapse", hap.names = hap.names(chr11_relapse), cex.lab = 0.3)
ZIP1_1_res2 <- calc_ehh(chr11_background, mrk = 29990)
plot(ZIP1_1_res2, main = "ZIP1 SNP 1 Background")
ZIP1_1_furcation2 <- calc_furcation(chr11_background, mrk = 29990)
plot(ZIP1_1_furcation2, main = "ZIP1 SNP 1 Background", hap.names = hap.names(chr11_background), xlim = c(1684000, 1690000), cex.lab = 0.3)
ZIP1_2_res <- calc_ehh(chr11_relapse, mrk = 30004)
plot(ZIP1_2_res, main = "ZIP1 SNP2 Relapse")
ZIP1_2_furcation <- calc_furcation(chr11_relapse, mrk = 30004)
plot(ZIP1_2_furcation, main = "ZIP1 SNP2 Relapse", hap.names = hap.names(chr11_relapse), cex.lab = 0.3)
ZIP1_2_res2 <- calc_ehh(chr11_background, mrk = 30004)
plot(ZIP1_2_res2, main = "ZIP1 SNP2 Background")
ZIP1_2_furcation2 <- calc_furcation(chr11_background, mrk = 30004)
plot(ZIP1_2_furcation2, main = "ZIP1 SNP2 Background", hap.names = hap.names(chr11_background), xlim = c(1682000, 1694000), cex.lab = 0.3)
CEP76_res <- calc_ehh(chr11_relapse, mrk = 35162)
plot(CEP76_res, main = "CEP76 Relapse")
CEP76_furcation <- calc_furcation(chr11_relapse, mrk = 35162)
plot(CEP76_furcation, main = "CEP76 Relapse", hap.names = hap.names(chr11_relapse), cex.lab = 0.3)
CEP76_res2 <- calc_ehh(chr11_background, mrk = 35162)
plot(CEP76_res2, main = "CEP76 Background")
CEP76_furcation2 <- calc_furcation(chr11_background, mrk = 35162)
plot(CEP76_furcation2, main = "CEP76 Background", hap.names = hap.names(chr11_background), xlim = c(1915000, 1930000), cex.lab = 0.3)

chr12_relapse <- data2haplohh(hap_file = "Pvivax_PvP01_12_v1_new_mask.recode_relapses_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr12_background <- data2haplohh(hap_file = "Pvivax_PvP01_12_v1_new_mask.recode_background_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr12_relapse_scan <- scan_hh(chr12_relapse)
chr12_background_scan <- scan_hh(chr12_background)
chr12_xpehh <- ies2xpehh(scan_pop1 = chr12_relapse_scan,
                         scan_pop2 = chr12_background_scan,
                         popname1 = "Relapse",
                         popname2 = "Background")
chr12_markers <- data.frame(CHR = c("PvP01_12_v1", "PvP01_12_v1", "PvP01_12_v1", "PvP01_12_v1"), POSITION = c("2102751",  "2102752",  "2102759",  "2358739"), ID = c("FER2", "FER2", "FER2", "ApiAP2"))
manhattanplot(chr12_xpehh, mrk = chr12_markers, mrk.col = "RED", mrk.pch = 19, mrk.cex = 3)
FER2_1_res <- calc_ehh(chr12_relapse, mrk = 40666)
plot(FER2_1_res, main = "FER2 SNP1 Relapse")
FER2_1_furcation <- calc_furcation(chr12_relapse, mrk = 40666)
plot(FER2_1_furcation, main = "FER2 SNP1 Relapse", hap.names = hap.names(chr12_relapse), cex.lab = 0.3)
FER2_1_res2 <- calc_ehh(chr12_background, mrk = 40666)
plot(FER2_1_res2, main = "FER2 SNP1 Background")
FER2_1_furcation2 <- calc_furcation(chr12_background, mrk = 40666)
plot(FER2_1_furcation2, main = "FER2 SNP1 Background", hap.names = hap.names(chr12_background), xlim = c(2090000, 2105000), cex.lab = 0.3)
FER2_2_res <- calc_ehh(chr12_relapse, mrk = 40667)
plot(FER2_2_res, main = "FER2 SNP2 Relapse")
FER2_2_furcation <- calc_furcation(chr12_relapse, mrk = 40667)
plot(FER2_2_furcation, main = "FER2 SNP2 Relapse", hap.names = hap.names(chr12_relapse), cex.lab = 0.3)
FER2_2_res2 <- calc_ehh(chr12_background, mrk = 40667)
plot(FER2_2_res2, main = "FER2 SNP2 Background")
FER2_2_furcation2 <- calc_furcation(chr12_background, mrk = 40667)
plot(FER2_2_furcation2, main = "FER2 SNP2 Background", hap.names = hap.names(chr12_background), xlim = c(2090000, 2105000), cex.lab = 0.3)
FER2_3_res <- calc_ehh(chr12_relapse, mrk = 40668)
plot(FER2_3_res, main = "FER2 SNP3 Relapse")
FER2_3_furcation <- calc_furcation(chr12_relapse, mrk = 40668)
plot(FER2_3_furcation, main = "FER2 SNP3 Relapse", hap.names = hap.names(chr12_relapse), cex.lab = 0.3)
FER2_3_res2 <- calc_ehh(chr12_background, mrk = 40668)
plot(FER2_3_res2, main = "FER2 SNP3 Background")
FER2_3_furcation2 <- calc_furcation(chr12_background, mrk = 40668)
plot(FER2_3_furcation2, main = "FER2 SNP3 Background", hap.names = hap.names(chr12_background), xlim = c(2090000, 2105000), cex.lab = 0.3)
chr12_ApiAP2_res <- calc_ehh(chr12_relapse, mrk = 45922)
plot(chr12_ApiAP2_res, main = "chr12 - ApiAP2 Relapse")
chr12_ApiAP2_furcation <- calc_furcation(chr12_relapse, mrk = 45922)
plot(chr12_ApiAP2_furcation, main = "chr12 - ApiAP2 Relapse", hap.names = hap.names(chr12_relapse), cex.lab = 0.3)
chr12_ApiAP2_res2 <- calc_ehh(chr12_background, mrk = 45922)
plot(chr12_ApiAP2_res2, main = "chr12 - ApiAP2 Background")
chr12_ApiAP2_furcation2 <- calc_furcation(chr12_background, mrk = 45922)
plot(chr12_ApiAP2_furcation2, main = "chr12 - ApiAP2 Background", hap.names = hap.names(chr12_background), xlim = c(2345000, 2365000), cex.lab = 0.3)


chr13_relapse <- data2haplohh(hap_file = "Pvivax_PvP01_13_v1_new_mask.recode_relapses_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr13_background <- data2haplohh(hap_file = "Pvivax_PvP01_13_v1_new_mask.recode_background_only.vcf.gz", polarize_vcf = FALSE, vcf_reader = "vcfR")
chr13_relapse_scan <- scan_hh(chr13_relapse)
chr13_background_scan <- scan_hh(chr13_background)
chr13_xpehh <- ies2xpehh(scan_pop1 = chr13_relapse_scan,
                         scan_pop2 = chr13_background_scan,
                         popname1 = "Relapse",
                         popname2 = "Background")
chr13_markers <- data.frame(CHR = c("PvP01_13_v1"), POSITION = c("879014"), ID = c("PMIX"))
manhattanplot(chr13_xpehh, mrk = chr13_markers, mrk.col = "RED", mrk.pch = 19, mrk.cex = 3)
PMIX_res <- calc_ehh(chr13_relapse, mrk = 17428)
plot(PMIX_res, main = "PMIX Relapse")
PMIX_furcation <- calc_furcation(chr13_relapse, mrk = 17428)
plot(PMIX_furcation, main = "PMIX Relapse", hap.names = hap.names(chr13_relapse), cex.lab = 0.3)
PMIX_res2 <- calc_ehh(chr13_background, mrk = 17428)
plot(PMIX_res2, main = "PMIX Background")
PMIX_furcation2 <- calc_furcation(chr13_background, mrk = 17428)
plot(PMIX_furcation2, main = "PMIX Background", hap.names = hap.names(chr13_background), xlim = c(874000, 886000), cex.lab = 0.3)
