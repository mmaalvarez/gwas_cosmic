if(interactive()){.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
  setwd("~/re_gecip/cancer_pan/malvarez/gwas_cosmic/scripts")}
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

args = commandArgs(trailingOnly=TRUE)

clumped_flattened_output = ifelse(interactive(),
                                  yes = Sys.glob("../res/linear_clumped")[1],
                                  no = args[1]) %>% 
  read_delim %>% 
  filter(!is.na(SBS))

gwas_method = unique(clumped_flattened_output$method)


thr_signif_pval = ifelse(interactive(),
                         yes = "1e-25",
                         no = args[2]) %>% 
  as.numeric


jet_colors = colorRampPalette(c("gray", "red", "yellow", "green", "cyan", "blue", "magenta", "black"))

clumped_flattened_output = clumped_flattened_output %>% 
  mutate(`#CHROM` = factor(`#CHROM`, levels = seq(1, 25), ordered = T)) %>% 
  arrange(`#CHROM`, POS)

manhattan = ggplot(clumped_flattened_output,
                   aes(x = POS,
                       y = -log10(P))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_point(size = 0.8,
             alpha = 0.8,
             aes(color = `#CHROM`)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"),
                                  times = ceiling(length(unique(clumped_flattened_output$`#CHROM`))/2))) +
  guides(color = "none") +
  # highlight hits (i.e. above thr_signif_pval)
  geom_point(data = clumped_flattened_output %>% 
               filter(P <= thr_signif_pval) %>% 
               mutate(SBS = factor(SBS, levels = gtools::mixedsort(unique(.$SBS)))),
             aes(fill = SBS),
             size = 3,
             shape = 21) +
  scale_fill_manual(values = jet_colors(length(unique(clumped_flattened_output %>% 
                                                        filter(P <= thr_signif_pval) %>% 
                                                        pull(SBS))))) +
  ggrepel::geom_text_repel(data = clumped_flattened_output %>% 
                             filter(P <= thr_signif_pval),
                           aes(label = paste0("chr", `#CHROM`, ":", POS)),
                           angle = 90,
                           size = 2.5,
                           segment.color = "grey50",
                           min.segment.length = 0.1,
                           nudge_x = 0.2,
                           max.overlaps = 100000) +
  geom_hline(yintercept = -log10(thr_signif_pval), linetype = "dashed", color = "red") +
  facet_grid(cols = vars(`#CHROM`), scales = "free_x", space = "free_x") +
  labs(x = "Genomic coordinate", y = "-log10(p-value)") +
  ggtitle(paste0(gsub("linear", "PLINK", gwas_method), " results after clumping")) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.001, "lines"),
        strip.text.x = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank())
ggsave(filename = paste0("manhattan_", gwas_method, ".jpg"),
       plot = manhattan,
       width = 12.5,
       height = 7,
       dpi = 300,
       bg = "white")
