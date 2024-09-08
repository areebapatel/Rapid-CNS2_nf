process methylationClassification {
    input:
        path(methylationClassificationScript)
        path(bedmethyl_file)
        val(id)
        path(topProbes)
        path(trainingData)
        path(arrayFile)
        val(meth_threads)

    publishDir("${params.outDir}/methylation_classification")

    output:
        val true

    script:
        """
        Rscript ${methylationClassificationScript} \
        --sample ${id} \
        --out_dir ${params.outDir}/methylation_classification \
        --in_file ${bedmethyl_file} \
        --probes ${topProbes} \
        --training_data ${trainingData} \
        --array_file ${arrayFile} \
        --threads ${meth_threads}
        """
}


process methylartistMGMT {
    input:
        path(inputBam)
        path(inputBai)
        path(ref)
        val(params.outDir)
	val ready // mgmt_coverage has run

    publishDir("${params.outDir}/mgmt")
    
    output:
        val true
        path "*.png", emit: mgmt_plot optional true
    
    script:
        cov_file = file("${params.outDir}/mgmt/mgmt_avg_cov.txt")
        if ( cov_file.exists() == true )    
            """
            methylartist \
            locus \
            -i chr10:129466536-129467536 \
            --samplepalette magma \
            -l 101126888-101129371 \
            --highlightpalette viridis \
            -b ${inputBam} \
            --ref ${ref} \
            --motif CG \
            --mods m
            """
        else

            """
            exit 1
            """
}