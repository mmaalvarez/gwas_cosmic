if(interactive()){.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
                  setwd("~/re_gecip/cancer_pan/malvarez/gwas_cosmic/scripts")}
library(tidyverse)
library(RNOmni)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
            
args = commandArgs(trailingOnly=TRUE)

PCs_covariates = ifelse(interactive(),
                          yes = "20",
                          no = args[1]) %>% 
  as.numeric

eigenvec = ifelse(interactive(),
                  yes = Sys.glob("../work/*/*/*_het_afreq_pruned.eigenvec")[1],
                  no = args[2]) %>% 
  read_delim(col_names = F) %>% 
  `colnames<-`(c("FID", "IID", paste0("PC", seq(1,length(names(.)))))) %>% 
  # keep only PCs up to PCs_covariates
  select(FID, IID, matches(paste0("PC", seq(1, PCs_covariates))))

sample_metadata = ifelse(interactive(),
                         yes = "../data/sample_metadata/cancer_analysis_metadata.tsv",
                         no = args[3]) %>%
  read_tsv(., guess_max = 10000) %>% 
  select(IID, age, sex, tumor_type) %>% 
  distinct %>% 
  left_join(eigenvec, by = 'IID') %>% 
  # some samples from samples_metadata are not in eigenvec, probably were removed in QC
  drop_na %>% 
  relocate(FID, .before = IID) %>% 
  ## WARNING in case some samples have duplicated rows, where e.g. 1 row is age 31 and primary tumor, and the other is age 32 and metastases, with all the remaining columns being the same -- keep the first row
  group_by(FID, IID, sex, across(starts_with("PC"))) %>% 
  slice_head(n = 1) %>% 
  ungroup

write_tsv(sample_metadata,
          "age_sex_tumor_pcs.tsv")


## get original_SBS_exposures for the final samples
original_continuous_phenotype = ifelse(interactive(),
                                       yes = "../data/degasperi_cosmic_exposures/RefSig_SBS_Exposures_v2.03__mutburden_norm.tsv",
                                       no = args[4]) %>% 
  read_tsv %>% 
  # a couple of samples from degasperi are not in sample_metadata, remove them
  right_join(sample_metadata) %>%
  drop_na %>% 
  # back to only SBS columns, for RINT
  select(IID, starts_with("SBS")) %>%
  column_to_rownames("IID")

## take each SBS column (all samples) and return its RINT (rank-based inverse normal transformation)
RINT_continuous_phenotype = apply(original_continuous_phenotype, 2, RankNorm) %>% 
  data.frame %>% 
  rownames_to_column("FID") %>% 
  mutate(IID = FID) %>% 
  relocate(IID, .after = FID)

write_tsv(RINT_continuous_phenotype,
          "RINT_continuous_phenotype.tsv")
