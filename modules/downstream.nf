process qqplot {

    label 'short_high'

    publishDir "$PWD/res/", pattern: '*.{jpg,tsv}', mode: 'copy'

    input:
    path(preclump_pruned_flattened_output)

    output:
	path('*.jpg')
	path('fdr_thresholds_*.tsv')

    script:
    """
    Rscript "${System.env.work_dir}"/scripts/9a_qqplot.R ${preclump_pruned_flattened_output}
    """

    stub:
    """
	touch dummy.jpg
	touch fdr_thresholds_dummy.tsv
    """
}

process manhattan {

    label 'short_high'

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

