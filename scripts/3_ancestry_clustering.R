if(interactive()){setwd("/g/strcombio/fsupek_data/users/malvarez/projects/gwas_cosmic/scripts")}
library(tidyverse)
library(tclust)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

args = commandArgs(trailingOnly=TRUE)

# first PCs used for clustering using the R package tclust (version 1.4.2)
# trim some % of the outlying samples

# samples grouped into k clusters (k independent gwas cohorts)

samples_ancestry_labels = ifelse(interactive(),
                                 yes = "../data/sample_metadata/cancer_analysis_metadata.tsv",
                                 no = args[1]) %>%
  read_tsv(., guess_max = 10000) %>% 
  select(IID, ancestry) %>% 
  distinct

PCs_ancestry_clustering = ifelse(interactive(),
								 yes = "10",
								 no = args[2]) %>%
  as.numeric

LD_pruned_eigenvec = ifelse(interactive(),
                            yes = Sys.glob("../work/*/*/LD_pruned.eigenvec")[1],
                            no = args[3]) %>%
  read_delim(., col_names = F) %>% 
  select(-X2) %>% # this is just the same sample name as in X1
  column_to_rownames("X1") %>% 
  `colnames<-`(paste0("PC", seq(1, length(names(.))))) %>% 
  select(paste0("PC", seq(1,PCs_ancestry_clustering))) %>% 
  rownames_to_column("IID") %>% 
  left_join(samples_ancestry_labels) %>% 
  arrange(ancestry)

# using these labels to get max k, as it's more detailed than `Genetically Inferred Ancestry (No)Thr`, and other columns have many NAs
max_n_ancestries = length(unique(LD_pruned_eigenvec$ancestry))

# The proportion of observations to be trimmed
clustering_outliers_trimmed = ifelse(interactive(),
                                     yes = "0.01", # Mischan used 1%
                                     no = args[4]) %>% 
  as.numeric

# restriction factor >=1, 1 is the most strict (less cohorts)
clustering_restriction_factor = ifelse(interactive(),
				                       yes = "1",
				                       no = args[5]) %>% 
  as.numeric


## run tclust
tclust_res = tclust(as.matrix(LD_pruned_eigenvec %>% 
                                select(starts_with("PC"))),
                    k = max_n_ancestries,
                    alpha = clustering_outliers_trimmed,
                    restr.fact = clustering_restriction_factor)

## keep only samples with other samples in its same cluster, and exclude also cohort==0 as these are the clustering_outliers_trimmed trimmed samples
samples_after_trimming = tclust_res$cluster %>% 
  data.frame %>% 
  ## WARNING this assumes that the samples order kept being the same after running tclust
  cbind(LD_pruned_eigenvec$IID, LD_pruned_eigenvec$ancestry) %>% 
  `colnames<-`(c("cohort", "IID", "ancestry")) %>% 
  arrange(cohort) %>% 
  ## remove trimmed samples (cohort=='0')
  filter(cohort != 0) %>%
  mutate(FID = IID) %>% 
  select(FID, IID)

write_tsv(samples_after_trimming, "samples_after_trimming.tsv", col_names = F)
