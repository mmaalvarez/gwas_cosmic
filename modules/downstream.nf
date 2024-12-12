process qqplot {

    label 'short_low'

    publishDir "$PWD/res/", pattern: '*.{jpg}', mode: 'copy'

    input:
    path(preclump_flattened_output)
	tuple path(bed_Q3_pruned),
		  path(bim_Q3_pruned),
		  path(fam_Q3_pruned)

    output:
	path('*.jpg')

    script:
    """
    Rscript "${System.env.work_dir}"/scripts/9a_qqplot.R ${preclump_flattened_output} ${bim_Q3_pruned}
    """

    stub:
    """
	touch dummy.jpg
    """
}

process manhattan {

    label 'short_low'

    publishDir "$PWD/res/", pattern: '*.{jpg}', mode: 'copy'

    input:
    path(clumped_flattened_output)
	val(thr_signif_pval)

    output:
	path('*.jpg')

    script:
    """
    Rscript "${System.env.work_dir}"/scripts/9b_manhattan.R ${clumped_flattened_output} ${thr_signif_pval}
    """

    stub:
    """
	touch dummy.jpg
    """
}

