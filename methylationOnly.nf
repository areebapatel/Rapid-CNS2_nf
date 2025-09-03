#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
software_version = "3.0.0"

// Display startup message
log.info """
╔══════════════════════════════════════════════════════════════════════════════╗
║                    🧬 Rapid-CNS² Methylation-Only Pipeline 🧬                ║
║                              Version ${software_version}                     ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                              ║
║  📋 Pipeline Information:                                                    ║
║  • Methylation analysis and classification only                              ║
║  • 5mC methylation calling using modkit                                      ║
║  • MGMT promoter methylation assessment                                      ║
║  • Methylation-based tumor classification                                    ║
║                                                                              ║
║  • Developer: Areeba Patel                                                   ║
║  • Email: a.patel@dkfz.de                                                    ║
║  • Institution: German Cancer Research Center (DKFZ)                         ║  
║                                                                              ║
║  ⚠️  Important Notes:                                                        ║
║  • This pipeline is for RESEARCH USE ONLY                                    ║
║  • Not validated for clinical diagnostic use                                 ║
║  • Results should be interpreted by qualified professionals                  ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

// Parameters
params.input = null
params.ref = null
params.outDir = null
params.id = null
params.patient = null
params.minimumMgmtCov = 5
params.modkitThreads = 32
params.methThreads = 64
params.mgmtThreads = 8
params.mnpFlex = false

// Show help message
if (params.help) {
   log.info """
   ================================================================================
   Rapid-CNS² Methylation-Only Pipeline v${software_version}
   ================================================================================
   
   USAGE: nextflow run methylation_only.nf [options]
   
   MANDATORY PARAMETERS:
       --input            Path to input BAM file(s)
       --id               Sample identifier (alphanumeric, no spaces)
       --ref              Path to hg38 reference genome FASTA file
       --outDir           Output directory for all results
   
   OPTIONAL PARAMETERS:
       --patient          Patient name for reports [default: uses --id value]
       --minimumMgmtCov   Minimum coverage for MGMT analysis [default: 5]
       --modkitThreads    Threads for modkit methylation calling [default: 32]
       --methThreads      Threads for methylation classification [default: 64]
       --mgmtThreads      Threads for MGMT promoter analysis [default: 8]
       --mnpFlex          Enable MNP-Flex classifier input preparation [default: false]
   
   EXAMPLES:
   
   Basic methylation analysis:
       nextflow run methylation_only.nf --input /data/sample.bam --id SAMPLE001 --ref /data/hg38.fa --outDir /data/output
   
   With MNP-Flex preparation:
       nextflow run methylation_only.nf --input /data/sample.bam --id SAMPLE001 --ref /data/hg38.fa --outDir /data/output --mnpFlex
   
   ================================================================================
   """
    exit 0
}

// Verify mandatory parameters
if (params.input == null) error "The path to the input BAM file(s) is mandatory, please specify it with --input"
if (params.id == null) error "The sample identifier is mandatory, please specify it with --id"
if (params.ref == null) error "The reference genome file is mandatory, please specify it with --ref"
if (params.outDir == null) error "The output directory is mandatory, please specify it with --outDir"

// Include required modules
include { checkAlignment } from './nextflow/bamProcessing.nf'
include { checkMethylationTags } from './nextflow/bamProcessing.nf'
include { indexBam } from './nextflow/bamProcessing.nf'
include { alignBam } from './nextflow/bamProcessing.nf'
include { mergeBam } from './nextflow/bamProcessing.nf'
include { methylationCalls } from './nextflow/methylationAnalysis.nf'
include { checkMgmtCoverage } from './nextflow/methylationAnalysis.nf'
include { mgmtPromoterMethyartist } from './nextflow/methylationAnalysis.nf'
include { mgmtPred } from './nextflow/methylationAnalysis.nf'
include { methylationClassification } from './nextflow/methylationClassification.nf'
include { mnpFlex } from './nextflow/methylationClassification.nf'

// Check if the input is a BAM file or a directory containing BAM files
def isBamFile(path) {
    return path.toString().toLowerCase().endsWith('.bam')
}

// Create output directories
[
    "${params.outDir}/bam",
    "${params.outDir}/bam/alignedBams",
    "${params.outDir}/mods",
    "${params.outDir}/mgmt",
    "${params.outDir}/methylation_classification",
    "${params.outDir}/mnpflex"
].each { dir ->
    new File(dir).mkdirs()
}

workflow {
    // Set channels
    Channel.fromPath(params.ref, checkIfExists: true).set { ref }
    Channel.from(params.id).set { id }
    Channel.from(params.patient ?: params.id).set { patient }
    Channel.from(params.outDir).set { outDir }
    Channel.from(params.modkitThreads).set { modkitThreads }
    Channel.from(params.methThreads).set { methThreads }
    Channel.from(params.mgmtThreads).set { mgmtThreads }
    Channel.from(params.minimumMgmtCov).set { minimumMgmtCov }

    // Set data files
    Channel.fromPath("${projectDir}/data/mgmt_hg38.bed", checkIfExists: true).set { mgmtBed }
    Channel.fromPath("${projectDir}/data/mgmt_probes.Rdata", checkIfExists: true).set { mgmtProbes }
    Channel.fromPath("${projectDir}/data/mgmt_137sites_mean_model.Rdata", checkIfExists: true).set { mgmtModel }
    Channel.fromPath("${projectDir}/data/top_probes_hm450.Rdata", checkIfExists: true).set { topProbes }
    Channel.fromPath("${projectDir}/data/capper_top_100k_betas_binarised.Rdata", checkIfExists: true).set { trainingData }
    Channel.fromPath("${projectDir}/data/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true).set { arrayFile }
    Channel.fromPath("${projectDir}/scr/mgmt_pred.R", checkIfExists: true).set { mgmtScript }
    Channel.fromPath("${projectDir}/scr/methylation_classification.R", checkIfExists: true).set { methylationClassificationScript }
    Channel.fromPath("${projectDir}/scr/mnp-flex_preprocessing.sh", checkIfExists: true).set { mnpFlexScript }
    Channel.fromPath("${projectDir}/data/MNP-flex.bed", checkIfExists: true).set { mnpFlexBed }

    // Input BAM file processing
    def inputPath = file(params.input)
    def bamToCheck
    // Check if the input is a BAM file or a directory containing BAM files
    if (isBamFile(inputPath)) {
        bamToCheck = Channel.value(inputPath)
    } else {
        bamToCheck = Channel.fromPath("${inputPath}/*.bam").first()
    }
    
    // Check alignment status
    checkAlignment_out = checkAlignment(bamToCheck, params.modkitThreads)
    
    // Check if the BAM file has methylation tags
    checkMethylationTags_out = checkMethylationTags(bamToCheck, params.modkitThreads)

    // Check if methylation tags exist; if not, raise an error and provide guidance
    checkMethylationTags_out.view { result ->
        if (!result || result.trim().isEmpty()) {
            error "No methylation tags (MM:Z) found in the BAM file. Please redo basecalling with modified basecalling enabled. See: https://github.com/nanoporetech/dorado?tab=readme-ov-file#modified-basecalling or ensure modified bases (5mC) are enabled in MinKNOW."
        }
    }

    def doAlignment = false
    
    checkAlignment_out.collect().view { results ->
        results.each { result ->
            def alignedCount = result.trim().toInteger()
            if (alignedCount <= 2) {
                println "Input BAM file(s) have ${alignedCount} aligned reads. Alignment will be performed."
                doAlignment = true
            } else {
                println "Input BAM file(s) have ${alignedCount} aligned reads. Using existing alignment."
            }
        }
    }

    // Conditionally run alignment and merging based on alignment status
    def processedBam
    if (doAlignment) {
        // Files are unaligned - need to align them
        alignedBams = alignBam(inputPath, ref, params.modkitThreads, outDir).alignedBam
        processedBam = alignedBams.collect().map { bamList ->
            if (bamList.size() > 1) {
                mergeBam(bamList, params.modkitThreads, outDir, id).mergedBam
            } else if (bamList.size() == 1) {
                bamList[0]
            } else {
                error "No aligned BAM files found after alignment."
            }
        }.flatten()
    } else {
        // Files are already aligned - check if single or multiple
        if (isBamFile(inputPath)) {
            // Single aligned BAM file - use bamToCheck directly
            processedBam = bamToCheck
        } else {
            // Multiple aligned BAM files - merge them
            processedBam = mergeBam(inputPath, params.modkitThreads, outDir, id).mergedBam
        }
    }

    // Index the processed BAM  
    indexedBam = indexBam(processedBam, params.modkitThreads)

    // Call methylation
    methylationCalls = methylationCalls(processedBam, indexedBam.indexBam, ref, id, modkitThreads)

    // Methylation classification
    methylationClassificationOut = methylationClassification(methylationClassificationScript, methylationCalls.bedmethylFile, id, topProbes, trainingData, arrayFile, methThreads)

    // MGMT promoter analysis
    mgmtCoverageOut = checkMgmtCoverage(processedBam, indexedBam.indexBam, mgmtBed, minimumMgmtCov, mgmtThreads)
    mgmtPromoterOut = mgmtPromoterMethyartist(processedBam, indexedBam.indexBam, ref, mgmtCoverageOut[0], id)
    mgmtPredOut = mgmtPred(mgmtCoverageOut[0], mgmtScript, mgmtBed, mgmtProbes, mgmtModel, methylationCalls.bedmethylFile, id)

    // MNP-Flex preparation (optional)
    if (params.mnpFlex) {
        mnpFlex(mnpFlexScript, methylationCalls.bedmethylFile, mnpFlexBed, id)
    }
}

workflow.onComplete {
    if(workflow.success) {
        println "The Rapid-CNS² methylation analysis is now complete!"
        println "Your outputs are located in: ${params.outDir}"
        println ""
        println "Output directories:"
        println "• ${params.outDir}/mods/ - Methylation calls (bedmethyl files)"
        println "• ${params.outDir}/mgmt/ - MGMT promoter analysis"
        println "• ${params.outDir}/methylation_classification/ - Tumor classification results"
        if (params.mnpFlex) {
            println "• ${params.outDir}/mnpflex/ - MNP-Flex compatible files"
        }
    } else {
        println "Oops .. something went wrong, please look into the log file and error messages in ${workDir}"
    }
}
