process prepare_linear {

    label 'short_high'

    input:
	val(PCs_covariates)
    tuple path(eigenvec_Q3),
          path(eigenval_Q3)
    path(sample_metadata)
    path(original_continuous_phenotype)

    output:
    tuple path('age_sex_tumor_pcs.tsv'),
		  path('RINT_continuous_phenotype.tsv')

    script:
    """
    Rscript "${System.env.work_dir}"/scripts/7a_prepare_linear.R ${PCs_covariates} ${eigenvec_Q3} ${sample_metadata} ${original_continuous_phenotype}
    """

    stub:
    """
    touch age_sex_tumor_pcs.tsv
    touch RINT_continuous_phenotype.tsv
    """
}

process linear_clumping {

    label 'short_medium_ignoreError' // otherwise, clumping warnings about 'no significant SNPs found' are somehow causing Nextflow to stop

    input:
    tuple path(bed_Q3_unpruned),
          path(bim_Q3_unpruned),
          path(fam_Q3_unpruned)
	tuple path(age_sex_tumor_pcs),
		  path(RINT_continuous_phenotype)
	val(clump_p1)
	val(clump_kb)
	val(clump_r2)

    output:
	path('*.linear'), emit: preclump
    path('*.clumps'), emit: clumped

    script:
    """
    bash "${System.env.work_dir}"/scripts/8a_linear_clumping.sh ${bed_Q3_unpruned} ${age_sex_tumor_pcs} ${RINT_continuous_phenotype} ${clump_p1} ${clump_kb} ${clump_r2}
    """

    stub:
    """
	touch dummy.linear
	touch dummy.clumps
    """
}

