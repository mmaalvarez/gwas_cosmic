process regenie_gwas {

    label 'short_medium'

    input:
    tuple path(bed_Q3_unpruned),
          path(bim_Q3_unpruned),
          path(fam_Q3_unpruned)
    path(sample_metadata)
    path(original_continuous_phenotype)

    output:
    path('*.regenie')

    script:
    """
    bash "${System.env.work_dir}"/scripts/7b_regenie_gwas.sh ${bed_Q3_unpruned} ${sample_metadata} ${original_continuous_phenotype}
    """

    stub:
    """
	touch dummy.regenie
    """
}

