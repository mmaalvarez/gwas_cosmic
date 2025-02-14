#!bin/bash

preprocessed_vcf=$1
high_LD_regions=$2
good_mappability_regions=$3
qc_hwe=$4
qc_mind=$5
prune_window_size=$6
prune_step_size=$7
prune_r2=$8


## VCF to PLINK

plink --vcf $preprocessed_vcf --make-bed --allow-extra-chr --chr 1-22,25,XY --double-id --fill-missing-a2 --out passed_biallelic_autosomal_snps_geno_maf

## NOTES
# --allow-extra-chr --> even though I specified autosomes in previous steps, there may be some weird chr names but nothing of importance, so this allows this to not raise an error
# then, --chr 1-22,25,XY once and for all excludes all unplaced and non-autosomal variants (and tries to not exclude the pseudo-autosomal region of X ('25' or 'XY'), but this is not there anyway)
# --double-id causes both family (FID) and within-family (IID) columns to be populated with the sampleID
# --fill-missing-a2 replace all missing calls (NAs) with homozygous A2 calls (reference allele, which is also the MAJOR allele because PLINK assigns major allele to reference (A2) unless stated otherwise

## NOT USED flags, because it is already done in the previous bcftools step:
# --snps-only --biallelic-only 'strict' 'list' --> WOULD keep only snps; 'strict': indiscriminately skip variants with 2+ alternate alleles listed even when only one alternate allele actually shows up
# WARNING: --keep-allele-order should be used in EVERY plink run IN CASE we wanted to prevent plink from converting the major allele as the REF allele (A2); --keep-allele-order would keep all vcf's REF as A2 and ALT as A1 regardless of their frequencies
# by default, reference alleles are set to A2 (reference naming in plink)
# by default, calls with uncertainty greater than 0.1 are treated as missing
# --vcf-filter WOULD skip variants (set as missing) which failed one or more filters tracked by the FILTER field: keep snps with PASS or '.' in FILTER column
# --vcf-half-call 'missing' --> WOULD treat half-calls as missing BUT THIS ONLY WORKS WITH --vcf



## remove duplicated-position SNPs that are actually >2-allelic SNPs split into rows, as this would raise issues with --set-missing-var-ids @:#

awk 'BEGIN{OFS="\t"} NR==FNR{count[$4]++; next} count[$4]==1{print "chr"$1, $4, $4, "*"}' passed_biallelic_autosomal_snps_geno_maf.bim passed_biallelic_autosomal_snps_geno_maf.bim > passed_biallelic_autosomal_snps_geno_maf.no_split_multiallelic

plink --bfile passed_biallelic_autosomal_snps_geno_maf --extract range passed_biallelic_autosomal_snps_geno_maf.no_split_multiallelic --make-bed --out passed_biallelic_autosomal_snps_geno_maf_noSplitMultiallelic

rm passed_biallelic_autosomal_snps_geno_maf.bed


## 1) rm high LD regions (Anderson 2010)
# https://raw.githubusercontent.com/cran/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.bed
## 2) keep regions which are possible to align according to the CRG75 alignability track
# /g/strcombio/fsupek_home/mmunteanu/reference/CRG75.bed --> but this is for hg19, need the hg38 but there is not in UCSC, only crg50
# I use Maia's, who generated the mappability for 75-mers for hg38 using the GEM package mappability function, which estimates mappability from the reference given a k-mer length and permitted number of mismatches (see https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeMapability&db=hg19). The code was based on https://evodify.com/gem-mappability/ and https://raw.githubusercontent.com/epigen/LIQUORICE/master/liquorice/create_mappability_bigwigs.sh. The GEM library is from https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2 -- Afterwards she also removed 42 blacklisted genomic regions
# /home/mmunteanu/reference/GRCh38.d1.vd1.mappability_75bp_1_22_XY_woblacklist.bed --> airlocked
## 3) rm SNPs not in HWE (P < 1e−x)

plink --bfile passed_biallelic_autosomal_snps_geno_maf_noSplitMultiallelic --exclude range $high_LD_regions --extract range $good_mappability_regions --hwe $qc_hwe --set-missing-var-ids @:# --make-bed --out passed_biallelic_autosomal_snps_geno_maf_noSplitMultiallelic_noHighLDregions_CRG75_hwe

## NOTES
# --exclude range --> takes a 1-based tracks file
# --extract range --> same as above, but here we are KEEPING the ranges
# --hwe filters out all variants not in hwe
#--set-missing-var-ids @:# --> SNPs with no rs name ('.') set to chr:number --> it has to be done after the custom removal of multiallelic-split-SNPs, otherwise these show up as duplicated-position SNPs that raise an error

rm passed_biallelic_autosomal_snps_geno_maf_noSplitMultiallelic.bed



## removal of samples with SNP missing rate >x%

plink --bfile passed_biallelic_autosomal_snps_geno_maf_noSplitMultiallelic_noHighLDregions_CRG75_hwe --mind $qc_mind --make-bed --out out_QC1

## NOTES
# --mind filters out all samples with missing call rates exceeding the provided value)

rm passed_biallelic_autosomal_snps_geno_maf_noSplitMultiallelic_noHighLDregions_CRG75_hwe.bed



### process SNPs for 1st PCA

# 1st time - LD pruning of SNPs
# window size of xbp, a step size of  xbp and r2 threshold of r2>x

# 1st time - PCA

plink --bfile out_QC1 --indep-pairwise $prune_window_size $prune_step_size $prune_r2 --make-bed --out LD_pruned

## NOTES
# --indep-pairwise <window size> <step size (variant ct)> <r^2 threshold>
# need to run the --filter exclude .prune.out afterwards, it doesn't do it automatically

plink --bfile LD_pruned --exclude LD_pruned.prune.out --pca --make-bed --out LD_pruned

rm LD_pruned.bed
rm *.nosex
rm *.bed~

