#!/bin/bash

bed_plink=$1
sample_metadata=$2
original_continuous_phenotype=$3

source /usr/bin/conda/etc/profile.d/conda.sh
conda activate regenie


# prefix name of bed/bim/fam input
bedbimfam_prefix_name=`echo $bed_plink | sed "s/\.bed//g ; s/.*\///g"`


## parse covariates and phenotypes files to match REGENIE format requirements

# covariates
awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="IID" || $i=="age" || $i=="sex" || $i=="tumor_type") {
            col[i]++
            if($i=="IID") iid_col=i
        }
    }
    printf "FID\tIID"
    for(i=1;i<=NF;i++) if(col[i] && $i!="IID") printf "\t%s",$i
    print ""
    next
} {
    OFS="\t"
	if(!seen[$iid_col]++) {
	    printf "%s\t%s",$iid_col,$iid_col
    	for(i=1;i<=NF;i++) if(col[i] && i!=iid_col) printf "\t%s",$i
    	print ""
	}
}' $sample_metadata > "$sample_metadata".parsed

# phenotypes
awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="IID" || $i ~ /^SBS/) {
            col[i]++
            if($i=="IID") iid_col=i
        }
    }
    printf "FID\tIID"
    for(i=1;i<=NF;i++) if(col[i] && $i!="IID") printf "\t%s",$i
    print ""
    next
} {
    OFS="\t"
    printf "%s\t%s",$iid_col,$iid_col
    for(i=1;i<=NF;i++) if(col[i] && i!=iid_col) printf "\t%s",$i
    print ""
}' $original_continuous_phenotype > "$original_continuous_phenotype".parsed


## REGENIE

# Step 1, the regression model is fit to the traits, and a set of genomic predictions are produced as output -- also normalize exposures to RINT
regenie --step 1 \
  --bed $bedbimfam_prefix_name \
  --covarFile "$sample_metadata".parsed \
  --catCovarList sex,tumor_type \
  --maxCatLevels 20 \
  --phenoFile "$original_continuous_phenotype".parsed \
  --force-qt \
  --bsize 100 \
  --apply-rint \
  --lowmem \
  --lowmem-prefix tmp_rg \
  --out fit_out

# Step 2, a set of imputed SNPs are tested for association
regenie --step 2 \
  --bed $bedbimfam_prefix_name \
  --covarFile "$sample_metadata".parsed \
  --catCovarList sex,tumor_type \
  --maxCatLevels 20 \
  --phenoFile "$original_continuous_phenotype".parsed \
  --force-qt \
  --bsize 200 \
  --apply-rint \
  --pred fit_out_pred.list \
  --out regenie_gwas_tmp

rm *.loco

for SBS in `head -n1 $original_continuous_phenotype.parsed | sed "s/	/\n/g" | grep "^SBS"`
do
	# add column with the SBS tested in each case (only if there is a tmp .regenie file for this SBS)
	# also "method" column, with "linear"
	regenie_SBS="regenie_gwas_tmp_"$SBS".regenie"
	if [ -f "$regenie_SBS" ]
	then
		awk -F' ' -v SBS="$SBS" 'NR==1{print $0" SBS method"}NR>1{print $0" "SBS" regenie"}' "$regenie_SBS" > regenie_gwas_"$SBS".regenie
	fi
done

rm *_tmp_*

# if there are no .regenie files because all SBS failed, write a dummy one
ls regenie_gwas_tmp_*.regenie &> /dev/null || echo -e "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA SBS method" > dummy.regenie

