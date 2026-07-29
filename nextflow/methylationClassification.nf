// Rapid-CNS2 random-forest methylation classification
process methylationClassification {
    label 'rapid_cns'

    publishDir "${params.outDir}/methylation_classification/", mode: 'copy'

    input:
        path(classificationScript)
        path(bedmethyl)
        val(id)
        path(topProbes)
        path(trainingData)
        path(arrayFile)
        val(methThreads)

    output:
        path "${id}_votes.tsv",                      emit: votes
        path "${id}_rf_details.tsv",                 emit: rfDetails
        path "${id}_calibrated_classification.tsv",  emit: calibrated
        path "${id}_methylation_hg38_HM450.RDS"
        path "${id}_probes_for_training.csv"

    script:
        """
        Rscript ${classificationScript} \
            --sample ${id} \
            --out_dir . \
            --in_file ${bedmethyl} \
            --probes_file ${topProbes} \
            --training_data ${trainingData} \
            --array_file ${arrayFile} \
            --threads ${methThreads}
        """
}

// Prepares the input file for the external MNP-Flex classifier (mnp-flex.org)
process mnpFlex {
    label 'rapid_cns'

    publishDir "${params.outDir}/mnpflex/", mode: 'copy'

    input:
        path(mnpFlexScript)
        path(bedmethyl)
        path(mnpFlexBed)
        val(id)

    output:
        // name must match what mnp-flex_preprocessing.sh writes
        path "${id}.MNPFlex.input.bed", emit: mnpFlexBed

    script:
        """
        bash ${mnpFlexScript} ${bedmethyl} ${mnpFlexBed} .
        """
}
