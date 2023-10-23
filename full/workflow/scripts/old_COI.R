####Older version; check longleaf for correct version

#Complexity of Infection using Real McCoil
#prepare data
library("devtools")
#install.packages("vcfR")
library("vcfR")
library("tidyverse")
setwd("~/Downloads")

file <-	"ov1_pf3d7-pf6k_masked_biallelicsnps_vqslodfiltermaf5miss80"
test <- vcfR::read.vcfR(paste(c(file,".vcf.gz"),collapse=""))
test <- vcfR::read.vcfR("ov1_curtisi_masked_biallelicsnps_hardfiltereddefmaf5miss80ld05.vcf")
test <- vcfR::read.vcfR("test_S11.vcf")
test <- vcfR::read.vcfR("ov1_wallikeriandmixed_masked_biallelicsnps_hardfiltereddefmaf5miss80.vcf")
ad_drc <- vcfR::extract.gt(test, element = "AD", as.numeric = T)
vcfR::heatmap.bp(ad_drc) 

loci <- test@fix[,1:2] %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(POS = as.numeric(POS),
                Key = 1:dplyr::n()) %>% 
  dplyr::select(c("CHROM", "POS", "Key"))

test_long <- test %>% 
  vcfR::extract_gt_tidy() %>% 
  dplyr::mutate(country = "test")

combined_long <- dplyr::bind_rows(test_long) %>% 
  # now lets merge the loci information with the individual level information
  dplyr::full_join(x = loci, y = ., by  = "Key") %>% 
  # don't need Key anymore
  dplyr::select(-c("Key"))

combined_long <- combined_long %>% 
  # use `mutate()` to add a new column
  dplyr::mutate(
    gt = dplyr::case_when(gt_GT == "0/0" ~ 0,
                          gt_GT == "0/1" ~ 0.5,
                          gt_GT == "1/1" ~ 1)  # change `gt_GT` column to `gt`
  )

n_hets <- combined_long %>% 
  dplyr::group_by(Indiv, country) %>% 
  dplyr::summarise(
    # multiple by 2 to assume diploid genome
    n_sites = length(gt),
    n_hets = sum(gt == 0.5, na.rm=T)
  )

n_hets %>% 
  ggplot() +
  geom_boxplot(aes(x = country, y = n_hets, color = country), outlier.shape = NA) +
  geom_jitter(aes(x = country, y = n_hets, color = country),
              alpha = 0.3,size = 0.5) +
  theme_linedraw() +
  theme(legend.position = "none") + 
  labs(y = "No. Heterozygous SNPs", x = "Country") 

RMCL <- combined_long %>% 
  dplyr::mutate(loci = paste(CHROM, POS, sep = "|")) %>% 
  dplyr::select(c("loci", "Indiv", "gt")) %>% 
  # liftover missing for RMCL 
  dplyr::mutate(gt = ifelse(is.na(gt), -1, gt)) %>% 
  tidyr::pivot_wider(names_from = "Indiv",
                     values_from = "gt")

orig_wd <- getwd()
setwd("~/Downloads/THEREALMcCOIL-master/categorical_method/")
source("McCOIL_categorical.R")



McCOIL_categorical(RMCL,
                   maxCOI=25, threshold_ind=2,
                   threshold_site=1,
                   totalrun=1000, burnin=100, M0=15,
                   e1=0.05, e2=0.05,
                   err_method=3, 
                   path=getwd(),
                   output=orig_wd # Up from THEREALMcCOIL/categorical/
)

####For testing

testsmpls <- colnames(test@gt)[2:ncol(test@gt)]
test_RMCL <- RMCL[,c("loci", testsmpls)]

test_RMCLmat <- as.matrix(test_RMCL[2:ncol(as.matrix(test_RMCL))])
rownames(test_RMCLmat) <- test_RMCL[["loci"]]
test_RMCLmat <- t(test_RMCLmat)
doubletest_RMCLmat <- matrix(rep(t(test_RMCLmat),2),ncol=ncol(test_RMCLmat),byrow=TRUE)
doubletest_RMCLmat <- doubletest_RMCLmat[1:2,1:5000]


#example dataset
data0= read.table("~/Downloads/THEREALMcCOIL-master/categorical_method/input_test.txt", head=T)
data = data0[1,2:33]
data=data0[,-1]
rownames(data)=data0[,1]
data = t(data)
McCOIL_categorical(data,maxCOI=25, threshold_ind=20, threshold_site=20, totalrun=1000, burnin=100, M0=15, e1=0.05, e2=0.05, err_method=3, path=getwd(), output="output_test.txt" )

orig_wd <- getwd()
setwd("~/Downloads/THEREALMcCOIL-master/categorical_method/")
source("McCOIL_categorical.R")
