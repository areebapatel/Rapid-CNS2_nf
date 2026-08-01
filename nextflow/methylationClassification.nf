// Rapid-CNS2 random-forest methylation classification
process methylationClassification {
    label 'heavy'

    publishDir "${params.outDir}/methylation_classification/", mode: 'copy'

    input:
        path(classificationScript)
        path(bedmethyl)
        val(id)
        path(topProbes)
        path(trainingData)
        path(arrayFile)

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
            --threads ${task.cpus}
        """
}

// Prepares the input file for the external MNP-Flex classifier.
// Upload the resulting BED at https://app.epignostix.com
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

// Optional upload of the MNP-Flex input to the Epignostix research platform.
//
// This sends patient-derived methylation data to a third-party service, so it
// is off by default and runs only with --mnpFlexUpload. Credentials are read
// from EPIGNOSTIX_USER / EPIGNOSTIX_PASSWORD in the environment; do not put
// them in a config file, which Nextflow copies into the task directory.
process mnpFlexUpload {
    label 'rapid_cns'

    publishDir "${params.outDir}/mnpflex/", mode: 'copy'

    input:
        path(uploadScript)
        path(bed)
        val(id)

    output:
        path "${id}_mnpflex_upload.json", emit: response, optional: true

    script:
        def loc  = params.mnpFlexLocalisation ? "--localisation '${params.mnpFlexLocalisation}'" : ''
        def diag = params.mnpFlexDiagnosis    ? "--diagnosis '${params.mnpFlexDiagnosis}'"       : ''
        """
        python3 ${uploadScript} \
            --bed ${bed} \
            --sample ${id} \
            --api ${params.mnpFlexApi} \
            --workflow-id ${params.mnpFlexWorkflowId} \
            --technology '${params.mnpFlexTechnology}' \
            --target-coverage '${params.mnpFlexCoverage}' \
            --extraction-type '${params.mnpFlexExtraction}' \
            --sex '${params.mnpFlexSex}' \
            ${loc} ${diag} \
            --out ${id}_mnpflex_upload.json
        """
}
