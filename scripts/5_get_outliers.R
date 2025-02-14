if(interactive()){setwd("/g/strcombio/fsupek_data/users/malvarez/projects/gwas_cosmic/scripts")}
library(tidyverse)
library(plinkQC)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

args = commandArgs(trailingOnly=TRUE)

# read QC2_trimmed_pruned.het prefix --> for sample heterozygosity outliers
het_QC2 = ifelse(interactive(),
                 yes = "QC2_trimmed_pruned.het",
                 no = args[1]) %>% 
  # keep prefix name only
  gsub(".*\\/|.het", "", .)

qc_het_SD = ifelse(interactive(),
                   yes = "3",
                   no = args[2]) %>%
  as.numeric

# read QC2_trimmed.bed/bim/fam prefix --> for snp allele frequency outliers
QC2_unpruned = ifelse(interactive(),
                     yes = "QC2_trimmed.bed",
                     no = args[3]) %>% 
  # keep prefix name only
  gsub(".*\\/|.bed", "", .)
                     
gold_standard_allele_freqs = ifelse(interactive(),
                                    yes = Sys.glob("../data/UKBB_SNP_freqs/ukbb_SNP_freqs_UK.tsv"),
                                    no = args[4]) %>% 
  read_tsv(col_names = F) %>% 
  `colnames<-`(c("CHR", "pos", "A1", "A2", "MAF")) %>% 
  unite("SNP", CHR, pos, sep = ":") %>% 
  # if MAF > 0.5, swap alleles and get 1-MAF
  mutate(A1_gold_standard = ifelse(MAF > 0.5, A2, A1),
         A2_gold_standard = ifelse(MAF > 0.5, A1, A2),
         MAF_gold_standard = ifelse(MAF > 0.5, 1-MAF, MAF)) %>% 
  select(SNP, A1_gold_standard, A2_gold_standard, MAF_gold_standard)
gc()

allele_freq_deviation = ifelse(interactive(),
                               yes = "0.1",
                               no = args[5]) %>%
  as.numeric

plink_path = ifelse(interactive(),
                    yes = "/tools/aws-workspace-apps/re_admin/source_code/plink/1.9/plink",
                    no = args[6])


### get samples with heterozygosity rate +/-X standard deviations from the mean, to remove in next process
het_outliers = ifelse(interactive(),
                      yes = Sys.glob("../work/*/*/QC2_trimmed_pruned.het")[1] %>% 
                        gsub("/QC2_trimmed_pruned.het", "", .),
                      no = ".") %>% 
  evaluate_check_het_and_miss(.,
                              het_QC2,
                              imissTh = 0.03,
                              hetTh = qc_het_SD) %>% 
  pluck("fail_het") %>% 
  select(FID, IID)
write_tsv(het_outliers, "het_outliers", col_names = F)
gc()


### get SNPs whose MAF diverges from golden standard, e.g. UK biobank 
allele_freq_deviated = ifelse(interactive(),
                              yes = Sys.glob("../work/00/206757e823244c52e8123f05a533ae/QC2_trimmed.bed")[1] %>% 
                                gsub("/QC2_trimmed.bed", "", .),
                              no = '.') %>% 
  check_maf(.,
            name = QC2_unpruned,
            # include all mafs
            macTh = NULL, mafTh = 0.51,
            path2plink = plink_path,
            # remove het_outliers from QC2_unpruned, before calculating the maf
            remove_individuals = './het_outliers') %>% 
  # get all MAFs
  pluck("fail_maf") %>% 
  # swap A1/A2, as in plink they are the other way around
  rename("A1_tmp" = "A1",
         "A1" = "A2") %>% 
  rename("A2" = "A1_tmp") %>% 
  select(SNP, A1, A2, MAF) %>% 
  # merge gold standard MAFs
  left_join(gold_standard_allele_freqs) %>% 
  # remove rows where there is no gold standard
  filter(!is.na(MAF_gold_standard)) %>% 
  # remove rows where alleles or allele order do not match (typically because MAF=~0.5)
  filter(A1==A1_gold_standard & A2==A2_gold_standard) %>% 
  # keep rows where MAF deviate more than allele_freq_deviation (to remove these SNPs in next process)
  filter(abs(MAF_gold_standard - MAF) > allele_freq_deviation) %>% 
  select(SNP)
write_tsv(allele_freq_deviated, "allele_freq_deviated", col_names = F)
gc()
