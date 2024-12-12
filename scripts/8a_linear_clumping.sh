#!/bin/bash

# PLINK2 gwas additive model --> SNP levels 0, 1, 2 (ref/ref, ref/alt, alt/alt)
# RINT-SBSi exposures ~ SNPj + age at sampling + phenotypic sex + tumor type + 20 PCs

bed_plink=$1
age_sex_tumor_pcs=$2
RINT_continuous_phenotype=$3
clump_p1=$4
clump_kb=$5
clump_r2=$6

# .bed/.bim/.fam from process QC3
bedbimfam_prefix_name=`echo $bed_plink | sed "s/\.bed//g ; s/.*\///g"`

plink2 --bfile $bedbimfam_prefix_name --glm hide-covar --covar $age_sex_tumor_pcs --covar-variance-standardize --pheno $RINT_continuous_phenotype --make-bed --out linear.tmp

## NOTES
# --glm == --linear from PLINK 1.9 if phenotype is quantitative
# --adjust would provide multiple testing correction for several methods, but I will do this correction after the clumping I think
# hide-covar modifier so that later in the clumping step the covariates pvalues are not included (see below)
# --covar-variance-standardize --> to avoid "genotype/covariate scales vary too widely for numerical stability of the current implementation"


# add column with the SBS tested in each case (only if there is a linear.tmp.SBS*.glm.linear file for this SBS)
# also "method" column, with "linear"
for SBS in `head -n1 $RINT_continuous_phenotype | sed "s/	/\n/g" | grep "^SBS"`
do
	tmp_linear="linear.tmp.$SBS.glm.linear"
	if [ -f "$tmp_linear" ]
	then
		awk -F'\t' -v SBS="$SBS" 'NR==1{print $0"\tSBS\tmethod"}NR>1{print $0"\t"SBS"\tlinear"}' "$tmp_linear" > linear.$SBS.glm.linear
	fi

	# clump together pval<=x SNPs within a radius of Xkbp and r2>=0.x
	plink2 --bfile linear.tmp --clump linear.tmp."$SBS".glm.linear --clump-p1 $clump_p1 --clump-kb $clump_kb --clump-r2 $clump_r2 --make-bed --out clumped."$SBS".tmp

	# here also add column with the SBS tested in each case (only if there is a .tmp.clumps file for this SBS) ; also "method" column, with "linear"
	tmp_clumps="clumped.$SBS.tmp.clumps"
	if [ -f "$tmp_clumps" ]
	then
		awk -F'\t' -v SBS="$SBS" 'NR==1{print $0"\tSBS\tmethod"}NR>1{print $0"\t"SBS"\tlinear"}' "$tmp_clumps" > clumped."$SBS".clumps
	fi
done

## NOTES
# When using --clump command with --linear/--logistic, all p-values are considered, including those of covariates. Thus, you probably want to run --linear/--logistic with the 'hide-covar' modifier
# Clumps are formed around central "index variants" which, by default, must have p-value no larger than 0.0001; change this threshold with --clump-p1. Index variants are chosen greedily starting with the lowest p-value
# Sites which are less than 250 kb away from an index variant and have r2 larger than 0.5 with it are assigned to that index variant's clump (unless they have been previously been assigned to another clump, and --clump-allow-overlap is not in effect). These two thresholds can be changed with --clump-kb and --clump-r2, respectively
# By default, no variant may belong to more than one clump; remove this restriction with --clump-allow-overlap
# With exactly one --clump file, --clump-best generates an additional .clumped.best report describing just the best proxy for each index variant (in the sense of having the highest r2 with it)

rm *.tmp.*

## if there are no .clumps because all SBS failed (no significant SNPs left after clumping, etc..), write a dummy one
ls clumped.*.clumps &> /dev/null || echo -e "#CHROM	POS	ID	P	TOTAL	NONSIG	S0.05	S0.01	S0.001	S0.0001	SP2	SBS method" > dummy.clumps

