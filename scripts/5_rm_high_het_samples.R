if(interactive()){.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
                  setwd("~/re_gecip/cancer_pan/malvarez/gwas_cosmic/scripts")}
library(tidyverse)
library(plinkQC)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

# rm samples with heterozygosity rate +/-3 standard deviations from the mean

args = commandArgs(trailingOnly=TRUE)

het_QC2 = ifelse(interactive(),
                 yes = Sys.glob("../work/*/*/trimmed_geno_maf_hwe_mind.het")[1],
                 no = args[1]) %>% 
  # keep prefix name only
  gsub(".*\\/|.het", "", .)

qc_het_SD = ifelse(interactive(),
                   yes = "3",
                   no = args[2]) %>%
  as.numeric

res = ifelse(interactive(),
             yes = Sys.glob("../work/*/*/trimmed_geno_maf_hwe_mind.het")[1] %>% 
                      gsub("trimmed_geno_maf_hwe_mind.het", "", .),
             no = ".") %>% 
  evaluate_check_het_and_miss(.,
                              het_QC2,
                              imissTh = 0.03,
                              hetTh = qc_het_SD)

write_tsv(res$fail_het %>% 
            select(FID, IID),
          "het_outliers", col_names = F)
