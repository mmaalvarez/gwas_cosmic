process QC1 {

    label 'short_medium_ignoreError'
    
    input:
    path(preprocessed_vcf)
    path(high_LD_regions)
    path(good_mappability_regions)
	val(qc_hwe)
	val(qc_mind)
	val(prune_window_size)
	val(prune_step_size)
	val(prune_r2)

    output:
    tuple path('out_QC1.bed'),
          path('out_QC1.bim'),
          path('out_QC1.fam'), emit: bed_bim_fam
    path('LD_pruned.eigenvec'), emit: eigenvec

    script:
    """
    bash "${System.env.work_dir}"/scripts/2_QC1.sh ${preprocessed_vcf} ${high_LD_regions} ${good_mappability_regions} ${qc_hwe} ${qc_mind} ${prune_window_size} ${prune_step_size} ${prune_r2}
    """

    stub:
    """
    touch out_QC1.bed
    touch out_QC1.bim
    touch out_QC1.fam
    touch LD_pruned.eigenvec
    """
}

process ancestry_clustering {

    label 'short_medium'
    
    input:
    path(sample_metadata)
    val(PCs_ancestry_clustering)
    path(LD_pruned_eigenvec)
    val(clustering_outliers_trimmed)
    val(clustering_restriction_factor)

    output:
    path('samples_after_trimming.tsv')

    script:
    """
    Rscript "${System.env.work_dir}"/scripts/3_ancestry_clustering.R \
        ${sample_metadata} \
        ${PCs_ancestry_clustering} \
        ${LD_pruned_eigenvec} \
        ${clustering_outliers_trimmed} \
        ${clustering_restriction_factor}
    """

    stub:
    """
    touch samples_after_trimming.tsv
    """
}

process QC2 {

    label 'short_high'
    
	input:
    tuple path(bed), 
		  path(bim),
		  path(fam)
    path(samples_after_trimming)
	val(prune_window_size)
	val(prune_step_size)
	val(prune_r2)

    output:
    tuple path('QC2_trimmed.bed'),
          path('QC2_trimmed.bim'),
          path('QC2_trimmed.fam'), emit: unpruned
    tuple path('QC2_trimmed_pruned.bed'),
          path('QC2_trimmed_pruned.bim'),
          path('QC2_trimmed_pruned.fam'),
          path('QC2_trimmed_pruned.het'),
          path('QC2_trimmed_pruned.imiss'),
          path('QC2_trimmed_pruned.lmiss'), emit: pruned

    script:
    """
    bash "${System.env.work_dir}"/scripts/4_QC2.sh ${bed} ${samples_after_trimming} ${prune_window_size} ${prune_step_size} ${prune_r2}
    """

    stub:
    """
    touch QC2_trimmed.bed
    touch QC2_trimmed.bim
    touch QC2_trimmed.fam
    touch QC2_trimmed_pruned.bed
    touch QC2_trimmed_pruned.bim
    touch QC2_trimmed_pruned.fam
    touch QC2_trimmed_pruned.het
    touch QC2_trimmed_pruned.imiss
    touch QC2_trimmed_pruned.lmiss
    """
}

process get_outliers {

    label 'short_high'

    input:
    tuple path(bed_QC2_pruned),
          path(bim_QC2_pruned),
          path(fam_QC2_pruned),
          path(het_QC2_pruned),
          path(imiss_QC2_pruned),
          path(lmiss_QC2_pruned)
	val(qc_het_SD)
	    tuple path(bed_QC2_unpruned),
          path(bim_QC2_unpruned),
          path(fam_QC2_unpruned)
	path(gold_standard_allele_freqs)
	val(allele_freq_deviation)
	val(plink_path)

    output:
    path('het_outliers'), emit: high_het_samples
	path('allele_freq_deviated'), emit: allele_freq_deviated

    script:
    """
    Rscript "${System.env.work_dir}"/scripts/5_get_outliers.R ${het_QC2_pruned} ${qc_het_SD} ${bed_QC2_unpruned} ${gold_standard_allele_freqs} ${allele_freq_deviation} ${plink_path}
    """

    stub:
    """
    touch het_outliers
	touch allele_freq_deviated
    """
}

process QC3 {

    label 'short_high'

    input:
    tuple path(bed_QC2_unpruned),
          path(bim_QC2_unpruned),
          path(fam_QC2_unpruned)
	path(het_outliers)
	path(allele_freq_deviated)
	val(mac)
	val(prune_window_size)
	val(prune_step_size)
	val(prune_r2)

    output:
    tuple path('QC3_het_afreq.bed'),
          path('QC3_het_afreq.bim'),
          path('QC3_het_afreq.fam'), emit: unpruned
    tuple path('QC3_het_afreq_pruned.bed'),
          path('QC3_het_afreq_pruned.bim'),
          path('QC3_het_afreq_pruned.fam'), emit: pruned
    tuple path('QC3_het_afreq_pruned.eigenvec'),
          path('QC3_het_afreq_pruned.eigenval'), emit: pca

    script:
    """
    bash "${System.env.work_dir}"/scripts/6_QC3.sh ${bed_QC2_unpruned} ${het_outliers} ${allele_freq_deviated} ${mac} ${prune_window_size} ${prune_step_size} ${prune_r2}
    """

    stub:
    """
    touch QC3_het_afreq.bed
    touch QC3_het_afreq.bim
    touch QC3_het_afreq.fam
    touch QC3_het_afreq_pruned.bed
    touch QC3_het_afreq_pruned.bim
    touch QC3_het_afreq_pruned.fam
    touch QC3_het_afreq_pruned.eigenvec
    touch QC3_het_afreq_pruned.eigenval
    """
}

