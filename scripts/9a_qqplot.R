if(interactive()){.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
  setwd("~/re_gecip/cancer_pan/malvarez/gwas_cosmic/scripts")}
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

args = commandArgs(trailingOnly=TRUE)

# p-values of all SNPs (pruned, pre-clumping)
preclump_pruned_flattened_output = ifelse(interactive(),
                                   yes = Sys.glob("../res/regenie_preclump_pruned_10percent")[1],
                                   no = args[1]) %>% 
  read_delim
gc()

gwas_method = unique(preclump_pruned_flattened_output$method)

# harmonize column names between GWAS methods
if(gwas_method == "regenie"){
  
  preclump_pruned_flattened_output = preclump_pruned_flattened_output %>% 
    rename("#CHROM" = "CHROM",
           "POS" = "GENPOS") %>% 
    mutate(P = 1*10^-LOG10P)
}

# remove NA pvalues
preclump_pruned_flattened_output = preclump_pruned_flattened_output %>% 
  filter(!is.na(P))


# expected uniform distribution
n = length(preclump_pruned_flattened_output$P)
expected = seq(1/(2*n), 1 - 1/(2*n), length.out = n)


## calculate FDR thresholds
fdr_levels = c(0.25, 0.2, 0.15, 0.1, 0.05, 0.01)

fdr_thresholds = sapply(fdr_levels, function(q) {
    # Sort p-values in ascending order
    p_sorted = sort(preclump_pruned_flattened_output$P)
    # Calculate critical values
    critical_vals = (1:n) * q / n
    # Find largest p-value meeting BH criterion
    passing = p_sorted <= critical_vals
    if(any(passing)) {
      max(p_sorted[passing])
    } else {
      0
    }
  }) %>% 
  as_tibble() %>% 
  cbind(fdr_levels) %>% 
  rename("p-value threshold" = "value",
         "FDR level" = "fdr_levels") %>% 
  select(`FDR level`, `p-value threshold`)

write_tsv(fdr_thresholds, 
          paste0("fdr_thresholds_", gwas_method, ".tsv"))


# plot
qqplot = ggplot(tibble(expected = expected,
                       observed = sort(preclump_pruned_flattened_output$P)),
                aes(x = expected,
                    y = observed)) + 
  coord_fixed() +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "lightblue") +
  labs(x = "EXP p-value", 
       y = "OBS p-value") +
  ggtitle(paste0(gsub("linear", "PLINK", gwas_method), " QQ plot")) + 
  theme_classic()

ggsave(filename = paste0("qqplot_", gwas_method, ".jpg"),
       plot = qqplot,
       width = 12.5,
       height = 7,
       dpi = 300,
       bg = "white")
