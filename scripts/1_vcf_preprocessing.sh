#!bin/bash

genome_chunk=$1
sample_list=$2
qc_geno=$3
qc_maf=$4
threads=$5

genome_chunk_nobrackets=`echo $genome_chunk | sed "s/\[//g ; s/\]//g"`
genome_chunk_prefix=`echo $genome_chunk | sed "s/.*chr/chr/g ; s/\..*//g"`

## extract samples and pre-process 1 vcf

bcftools view $genome_chunk_nobrackets --samples-file $sample_list --force-samples --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr25,chrXY --apply-filters PASS --min-alleles 2 --max-alleles 2 --types snps --exclude "F_MISSING > $qc_geno" --min-af "$qc_maf":minor --threads $threads --output-type b --output "$genome_chunk_prefix"_preprocessed

## NOTES
# --regions chr{1..22},chr25,chrXY --> $input_file is a vcf with autosomes, so no X, Y, nor M
# --samples-file --> to keep only the samples we want
# --force-samples --> only warn, but not stop, if some sample(s) in $sample_list are not in the vcf
# --apply-filters PASS --> only FILTER==PASS (but not '.') variants ; manual: "Skip sites where FILTER column does not contain any of the strings listed"
# --min-alleles 2 --max-alleles 2 --types snps --> keep only biallelic SNPs
# --exclude 'F_MISSING > x' --> exclude SNPs with >x*100 % missing data
# --min-af x:minor --> minimum minor allele frequency allowed (to keep common SNPs)
# --threads --> Use multithreading with INT worker threads. The option is currently used only for the compression of the output stream, only when --output-type is b or z. Default: 0.
# b --> "compressed BCF" output

