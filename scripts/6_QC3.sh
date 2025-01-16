#!/bin/bash

bed_plink=$1
het_outliers=$2
allele_freq_deviated=$3
mac=$4
prune_window_size=$5
prune_step_size=$6
prune_r2=$7

# prefix name of bed/bim/fam input
bedbimfam_prefix_name=`echo $bed_plink | sed "s/\.bed//g ; s/.*\///g"`

## high Heterozygosity samples to remove
# also remove SNPs whose a1 freq deviates from gold standard
# AND remove monomorphic SNPs, i.e. without any sample having the alternative allele (mac 1 or more, depending on cohort, specified in params.config), otherwise REGENIE will fail
plink --bfile $bedbimfam_prefix_name --remove $het_outliers --exclude $allele_freq_deviated --mac $mac --make-bed --out QC3_het_afreq


## 2nd time - PCA (LD pruned)
# pruning
plink --bfile QC3_het_afreq --indep-pairwise $prune_window_size $prune_step_size $prune_r2 --make-bed --out QC3_het_afreq_pruned

plink --bfile QC3_het_afreq_pruned --exclude QC3_het_afreq_pruned.prune.out --pca --make-bed --out QC3_het_afreq_pruned

rm *.nosex
