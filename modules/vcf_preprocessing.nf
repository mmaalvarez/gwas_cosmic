process vcf_preprocessing {

    label 'medium_low'
    
    input:
    val(genome_chunk)
    path(sample_list)
	val(qc_geno)
	val(qc_maf)
	val(bcftools_threads)

    output:
    file('chr*_preprocessed')

    script:
    """
    bash "${System.env.work_dir}"/scripts/1_vcf_preprocessing.sh ${genome_chunk} ${sample_list} ${qc_geno} ${qc_maf} ${bcftools_threads}
    """

    stub:
    """
    touch chr_STUB_preprocessed
    """
}

