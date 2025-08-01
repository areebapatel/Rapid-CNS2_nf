// This process runs the methylation classification script
process methylationClassification {
    label 'rapid_cns'
    
    input:
        path(methylationClassificationScript)
        path(bedmethylFile)
        val(id)
        path(topProbes)
        path(trainingData)
        path(arrayFile)
        val(methThreads)

    publishDir("${params.outDir}/methylation_classification")

    output:
        val true

    script:
        """
        Rscript ${methylationClassificationScript} \
        --sample ${id} \
        --out_dir . \
        --in_file ${bedmethylFile} \
        --probes ${topProbes} \
        --training_data ${trainingData} \
        --array_file ${arrayFile} \
        --threads ${methThreads}
        """
}

// This is a separate process to create the MNP-Flex compatible file
process mnpFlex {
    label 'rapid_cns'

    input:
        path(mnpFlexScript)
        path(bedmethylFile)
        path(mnpFlexBed)
        val(id)
    
    output:
        val true
    
    publishDir("${params.outDir}/mnpflex/")

    script:
        """
        bash ${mnpFlexScript} ${bedmethylFile} ${mnpFlexBed} . ${id}
        """
}