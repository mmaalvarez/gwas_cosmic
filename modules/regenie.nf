process regenie_gwas {

	label 'short_medium_ignoreError'

    input:
    tuple path(bed_Q3_unpruned),
          path(bim_Q3_unpruned),
          path(fam_Q3_unpruned)
    tuple path(bed_Q3_pruned),
          path(bim_Q3_pruned),
          path(fam_Q3_pruned)
    path(sample_metadata)
    path(original_continuous_phenotype)
	val(clump_p1)
	val(clump_kb)
	val(clump_r2)

    output:
	path('*.regenie')
	path('*.preclump.pruned'), emit: preclump_pruned
    path('*.clumps'), emit: clumped

    script:
    """
    bash "${System.env.work_dir}"/scripts/7b_regenie_gwas.sh ${bed_Q3_unpruned} ${bim_Q3_pruned} ${sample_metadata} ${original_continuous_phenotype} ${clump_p1} ${clump_kb} ${clump_r2}
    """

    stub:
    """
	touch dummy.regenie
	touch dummy.preclump.pruned
	touch dummy.clumps
    """
}

