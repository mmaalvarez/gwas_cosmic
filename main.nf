#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { vcf_preprocessing } from './modules/vcf_preprocessing'
include { QC1 } from './modules/qc'
include { ancestry_clustering } from './modules/qc'
include { QC2 } from './modules/qc'
include { get_outliers } from './modules/qc'
include { QC3 } from './modules/qc'
include { prepare_linear } from './modules/linear'
include { linear_clumping } from './modules/linear'
include { regenie_gwas } from './modules/regenie'
include { qqplot } from './modules/downstream'
include { manhattan } from './modules/downstream'


workflow {

	// channel by vcf chunks

    genome_chunks_list = Channel
        .fromPath(params.genome_chunks_list)
        .splitCsv(header:false)


	// run pre-GWAS parsing and QC processes

    vcf_preprocessing(
        genome_chunks_list,
        file(params.sample_list),
		params.qc_geno,
		params.qc_maf,
		params.bcftools_threads
    )

    QC1(
        vcf_preprocessing.out,
        file(params.high_LD_regions),
        file(params.good_mappability_regions),
		params.qc_hwe,
		params.qc_mind,
		params.prune_window_size,
		params.prune_step_size,
		params.prune_r2
    )

    ancestry_clustering(
        file(params.sample_metadata),
        params.PCs_ancestry_clustering,
        QC1.out.eigenvec,
        params.clustering_outliers_trimmed,
        params.clustering_restriction_factor
    )

    QC2(
        QC1.out.bed_bim_fam,
		ancestry_clustering.out,
		params.prune_window_size,
		params.prune_step_size,
		params.prune_r2
    )

    get_outliers(
        QC2.out.pruned,
		params.qc_het_SD,
		QC2.out.unpruned,
        file(params.gold_standard_allele_freqs),
		params.allele_freq_deviation,
		params.plink_path
    )
    // report how many SNPs were removed due to MAF deviating from gold standard
    get_outliers.out.allele_freq_deviated.collectFile(name: 'res/allele_freq_deviated')

    QC3(
		QC2.out.unpruned,
        get_outliers.out.high_het_samples,
        get_outliers.out.allele_freq_deviated,
		params.mac,
		params.prune_window_size,
		params.prune_step_size,
		params.prune_r2
    )

    // SNP-wise simple linear model with PLINK

    prepare_linear(
		params.PCs_covariates,
		QC3.out.pca,
        file(params.sample_metadata),
        file(params.original_continuous_phenotype)
    )

    linear_clumping(
		QC3.out.unpruned,
		QC3.out.pruned,
		prepare_linear.out,
		params.clump_p1,
		params.clump_kb,
		params.clump_r2
    )

    preclump_pruned_flattened_output_linear = linear_clumping.out.preclump_pruned
        .flatten()
        .collectFile(name: 'linear_preclump_pruned', keepHeader: true,
					 storeDir: 'res/')

    clumped_flattened_output_linear = linear_clumping.out.clumped
        .flatten()
        .collectFile(name: 'linear_clumped', keepHeader: true,
					 storeDir: 'res/')
        .view { "Finished SNP-wise simple linear GWAS using PLINK - Results saved in res/linear_*" }


    // REGENIE

	regenie_gwas(
		QC3.out.unpruned,
		QC3.out.pruned,
        file(params.sample_metadata),
        file(params.original_continuous_phenotype),
		params.clump_p1,
		params.clump_kb,
		params.clump_r2
	)

    preclump_pruned_flattened_output_regenie = regenie_gwas.out.preclump_pruned
        .flatten()
        .collectFile(name: 'regenie_preclump_pruned', keepHeader: true,
					 storeDir: 'res/')

    clumped_flattened_output_regenie = regenie_gwas.out.clumped
        .flatten()
        .collectFile(name: 'regenie_clumped', keepHeader: true,
					 storeDir: 'res/')
        .view { "Finished GWAS using REGENIE - Results saved in res/regenie_*" }
	

	// Downstream analyses and plots

	preclump_pruned_flattened_output = preclump_pruned_flattened_output_linear.collect().mix(preclump_pruned_flattened_output_regenie).collect()
	clumped_flattened_output = clumped_flattened_output_linear.collect().mix(clumped_flattened_output_regenie).collect()

	qqplot(
	  preclump_pruned_flattened_output
	)

	manhattan(
	  clumped_flattened_output,
	  params.thr_signif_pval
	)
}


workflow.onError {
    println "Pipeline execution stopped with error: ${workflow.errorMessage}"
}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    println "Duration: $workflow.duration"
}

