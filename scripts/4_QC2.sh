#!/bin/bash

# rm trimmed samples (ancestry)

# rm samples 
# heterozygosity rate +/-SD from the mean

bed_plink=$1
samples_after_trimming=$2
prune_window_size=$3
prune_step_size=$4
prune_r2=$5

# prefix name of bed/bim/fam input
bedbimfam_prefix_name=`echo $bed_plink | sed "s/\.bed//g ; s/.*\///g"`

## remove trimmed samples
plink --bfile $bedbimfam_prefix_name --keep $samples_after_trimming --make-bed --out QC2_trimmed


## get high Heterozygosity samples, to remove later
# for --het they recommend to first prune
plink --bfile QC2_trimmed --indep-pairwise $prune_window_size $prune_step_size $prune_r2 --make-bed --out QC2_trimmed_pruned
plink --bfile QC2_trimmed_pruned --het --missing --make-bed --out QC2_trimmed_pruned

## NOTES
# --het --> strongly negative F values indicate excess of heterozyosity, must be removed afterwards -- zzz.bwh.harvard.edu/plink/ibdibs.shtml
# --missing --> missing genotype rates per individual, to obtain the .imiss needed for plinkQC

rm *.nosex

## now run plinkQC in R to remove them

