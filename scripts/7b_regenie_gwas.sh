#!/bin/bash

bed_plink=$1
bim_pruned=$2
sample_metadata=$3
original_continuous_phenotype=$4
clump_p1=$5
clump_kb=$6
clump_r2=$7

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
  --bsize 1000 \
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
  --bsize 400 \
  --apply-rint \
  --pred fit_out_pred.list \
  --out regenie_gwas_tmp

#rm *.loco

for SBS in `head -n1 $original_continuous_phenotype.parsed | sed "s/	/\n/g" | grep "^SBS"`
do
	# add column with the SBS tested in each case (only if there is a tmp .regenie file for this SBS)
	# also "method" column, with "linear"
	regenie_SBS="regenie_gwas_tmp_"$SBS".regenie"
	if [ -f "$regenie_SBS" ]
	then
		awk -F' ' -v SBS="$SBS" 'NR==1{print $0" SBS method"}NR>1{print $0" "SBS" regenie"}' "$regenie_SBS" > regenie_gwas_"$SBS".regenie

		## also, extract pruned SNPs, for qqplot
		awk -v bim_pruned="$bim_pruned" 'BEGIN {FS="\t"; while (getline < bim_pruned) ids[$2]=1} {FS=" "} NR==1 || ($3 in ids)' regenie_gwas_"$SBS".regenie > regenie_gwas_"$SBS".preclump.pruned

		## also, convert regenie_gwas_tmp_"$SBS".regenie to plink2's *.glm.linear, for clumping
		awk 'BEGIN {print "#CHROM\tPOS\tID\tREF\tALT\tPROVISIONAL_REF?\tA1\tOMITTED\tA1_FREQ\tTEST\tOBS_CT\tBETA\tSET_STAT\tP\tERRCODE"} NR > 1 {
		    # Calculate P from LOG10P: P = 10^(-LOG10P)
    		p = 10^(-$12)

			# Print formatted output
			printf "%s\t%s\t%s\t%s\t%s\tY\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.50f\t.\n",
				$1,        # CHROM
				$2,        # POS
				$3,        # ID
				$4,        # REF
				$5,        # ALT
				$5,        # A1
				$4,        # OMITTED
				$6,        # A1_FREQ
				$8,        # TEST
				$7,        # OBS_CT
				$9,        # BETA
				$11,       # SET_STAT (CHISQ)
				p          # P-value
		}' $regenie_SBS > regenie_gwas_"$SBS".glm.linear

		rm $regenie_SBS

		# clump together pval<=x SNPs within a radius of Xkbp and r2>=0.x
		plink2 --bfile $bedbimfam_prefix_name --clump regenie_gwas_"$SBS".glm.linear --clump-p1 $clump_p1 --clump-kb $clump_kb --clump-r2 $clump_r2 --make-bed --out clumped."$SBS".tmp

		# here also add column with the SBS tested in each case (only if there is a .tmp.clumps file for this SBS) ; also "method" column, with "linear"
		tmp_clumps="clumped.$SBS.tmp.clumps"
		if [ -f "$tmp_clumps" ]
		then
			awk -F'\t' -v SBS="$SBS" 'NR==1{print $0"\tSBS\tmethod"}NR>1{print $0"\t"SBS"\tregenie"}' "$tmp_clumps" > clumped."$SBS".clumps
			rm $tmp_clumps
		fi
	fi
done

# if there are no .regenie files because all SBS failed, write a dummy one
ls regenie_gwas_tmp_*.regenie &> /dev/null || echo -e "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA SBS method" > dummy.regenie

# also in case there are no *.clumps files, since there were no SBS with significant clumps
ls *.clumps &> /dev/null || echo -e "#CHROM	POS	ID	P	TOTAL	NONSIG	S0.05	S0.01	S0.001	S0.0001	SP2	SBS	method" > dummy.clumps
