process QC1 {

    label 'medium_medium'
    
    input:
    path(preprocessed_bcf)
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
    bash "${System.env.work_dir}"/scripts/2_QC1.sh ${preprocessed_bcf} ${high_LD_regions} ${good_mappability_regions} ${qc_hwe} ${qc_mind} ${prune_window_size} ${prune_step_size} ${prune_r2}
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

process rm_high_het_samples {

    label 'short_high'

    input:
    tuple path(bed_QC2_pruned),
          path(bim_QC2_pruned),
          path(fam_QC2_pruned),
          path(het_QC2_pruned),
          path(imiss_QC2_pruned),
          path(lmiss_QC2_pruned)
	val(qc_het_SD)

    output:
    path('het_outliers')

    script:
    """
    Rscript "${System.env.work_dir}"/scripts/5_rm_high_het_samples.R ${het_QC2_pruned} ${qc_het_SD}
    """

    stub:
    """
    touch het_outliers
    """
}

process QC3 {

    label 'short_high'

    input:
    tuple path(bed_QC2_unpruned),
          path(bim_QC2_unpruned),
          path(fam_QC2_unpruned)
	path(het_outliers)
	val(prune_window_size)
	val(prune_step_size)
	val(prune_r2)

    output:
    tuple path('QC3_het.bed'),
          path('QC3_het.bim'),
          path('QC3_het.fam'), emit: unpruned
    tuple path('QC3_het_pruned.bed'),
          path('QC3_het_pruned.bim'),
          path('QC3_het_pruned.fam'), emit: pruned
    tuple path('QC3_het_pruned.eigenvec'),
          path('QC3_het_pruned.eigenval'), emit: pca

    script:
    """
    bash "${System.env.work_dir}"/scripts/6_QC3.sh ${bed_QC2_unpruned} ${het_outliers} ${prune_window_size} ${prune_step_size} ${prune_r2}
    """

    stub:
    """
    touch QC3_het.bed
    touch QC3_het.bim
    touch QC3_het.fam
    touch QC3_het_pruned.bed
    touch QC3_het_pruned.bim
    touch QC3_het_pruned.fam
    touch QC3_het_pruned.eigenvec
    touch QC3_het_pruned.eigenval
    """
}

