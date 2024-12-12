if(interactive()){.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
  setwd("~/re_gecip/cancer_pan/malvarez/gwas_cosmic/scripts")}
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

args = commandArgs(trailingOnly=TRUE)

# p-values of all SNPs (no pruning, pre-clumping)
preclump_flattened_output = ifelse(interactive(),
                                   yes = Sys.glob("../work/*/*/linear_preclump")[1],
                                   no = args[1]) %>% 
  read_delim

gwas_method = unique(preclump_flattened_output$method)

# harmonize column names between GWAS methods
if(gwas_method == "regenie"){
  
  preclump_flattened_output = preclump_flattened_output %>% 
    rename("#CHROM" = "CHROM",
           "POS" = "GENPOS") %>% 
    mutate(P = 1*10^-LOG10P)
}

# SNPs kept after pruning (originally obtained to do the final PCA)
bim_Q3_pruned = ifelse(interactive(),
                       yes = Sys.glob("../work/*/*/QC3_het_pruned.bim")[1],
                       no = args[2]) %>% 
  read_tsv(col_names = F) %>%
  rename("#CHROM" = "X1",
         "POS" = "X4") %>% 
  select(`#CHROM`, POS)

## keep only the p-values of the SNPs that would have been kept after pruning
preclump_flattened_output_pruned = preclump_flattened_output %>%
  right_join(bim_Q3_pruned)

# expected uniform distribution
n = length(preclump_flattened_output$P)
expected = seq(1/(2*n), 1 - 1/(2*n), length.out = n)

# plot
qqplot = ggplot(tibble(expected = expected,
                       observed = sort(preclump_flattened_output$P)),
                aes(x = expected,
                    y = observed)) + 
  coord_fixed() +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightblue") +
  labs(x = "EXP p-value", 
       y = "OBS p-value") +
  theme_classic()

ggsave(filename = paste0("qqplot_", gwas_method, ".jpg"),
       plot = qqplot,
       width = 12.5,
       height = 7,
       dpi = 300,
       bg = "white")
