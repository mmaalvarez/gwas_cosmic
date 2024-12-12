#!/bin/bash

bed_plink=$1
het_outliers=$2
prune_window_size=$3
prune_step_size=$4
prune_r2=$5

# prefix name of bed/bim/fam input
bedbimfam_prefix_name=`echo $bed_plink | sed "s/\.bed//g ; s/.*\///g"`

## high Heterozygosity samples to remove
plink --bfile $bedbimfam_prefix_name --remove $het_outliers --make-bed --out QC3_het


## 2nd time - PCA (LD pruned)
# pruning
plink --bfile QC3_het --indep-pairwise $prune_window_size $prune_step_size $prune_r2 --make-bed --out QC3_het_pruned

plink --bfile QC3_het_pruned --exclude QC3_het_pruned.prune.out --pca --make-bed --out QC3_het_pruned

rm *.nosex
