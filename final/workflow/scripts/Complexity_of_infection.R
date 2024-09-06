#Complexity of Infection using Real McCoil
#prepare data
library("devtools")
library("vcfR")
library("tidyverse")
#setwd("/pine/scr/k/e/kellybce/ovale1r/novaseq/ov1variants/")
#setwd("~/Documents/P. Ovale Genomic Analysis")


test <- vcfR::read.vcfR(snakemake@input)
#test <- vcfR::read.vcfR("pow.vcf")
ad_drc <- vcfR::extract.gt(test, element = "AD", as.numeric = T)
#vcfR::heatmap.bp(ad_drc) 

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

testsmpls <- colnames(test@gt)[2:ncol(test@gt)]
test_RMCL <- RMCL[,c("loci", testsmpls)]

test_RMCLmat <- as.matrix(test_RMCL[2:ncol(as.matrix(test_RMCL))])
rownames(test_RMCLmat) <- test_RMCL[["loci"]]
test_RMCLmat <- t(test_RMCLmat)
#doubletest_RMCLmat <- matrix(rep(t(test_RMCLmat),2),ncol=ncol(test_RMCLmat),byrow=TRUE)
#doubletest_RMCLmat <- doubletest_RMCLmat[1:2,1:5000]

orig_wd <- getwd()
#setwd("/pine/scr/k/e/kellybce/ovale1r/novaseq/ov1variants/statistics/")
devtools::install_github("OJWatson/McCOILR")
library("McCOILR")

invisible(McCOIL_categorical(test_RMCLmat,
                   maxCOI=snakemake@params[['maxcoi']], threshold_ind=snakemake@params[['threshold_ind']],
                   threshold_site=snakemake@params[['threshold_site']],
                   totalrun=snakemake@params[['totalrun']], burnin=snakemake@params[['burnin']], M0=15,
                   e1=0.05, e2=0.05,
                   err_method=3, 
                   path=getwd(),
                   output=snakemake@output # Up from THEREALMcCOIL/categorical/
))
