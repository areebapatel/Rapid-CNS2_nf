#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
// single source of truth: manifest.version in nextflow.config
software_version = workflow.manifest.version

// Display startup message
log.info """
╔══════════════════════════════════════════════════════════════════════════════╗
║                    🧬 Rapid-CNS² Methylation-Only Pipeline 🧬                ║
║                              Version ${software_version}                                   ║
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
params.input   = null
params.ref     = null
params.outDir  = null
params.id      = null
params.patient = null
params.help    = null

// Show help message
if (params.help) {
   log.info """
   ================================================================================
   Rapid-CNS² Methylation-Only Pipeline v${software_version}
   ================================================================================

   USAGE: nextflow run methylationOnly.nf [options]

   MANDATORY PARAMETERS:
       --input            Path to input BAM file(s)
       --id               Sample identifier (alphanumeric, no spaces)
       --ref              Path to hg38 reference genome FASTA file
       --outDir           Output directory for all results

   OPTIONAL PARAMETERS:
       --patient          Patient name for reports [default: uses --id value]
       --minimumMgmtCov   Minimum coverage for MGMT analysis [default: 5]
       --mnpFlex          Enable MNP-Flex classifier input preparation [default: true]

   EXAMPLES:

   Basic methylation analysis:
       nextflow run methylationOnly.nf --input /data/sample.bam --id SAMPLE001 --ref /data/hg38.fa --outDir /data/output

   With MNP-Flex preparation:
       nextflow run methylationOnly.nf --input /data/sample.bam --id SAMPLE001 --ref /data/hg38.fa --outDir /data/output --mnpFlex

   ================================================================================
   """
    exit 0
}

// Verify mandatory parameters
if (params.input == null) error "The path to the input BAM file(s) is mandatory, please specify it with --input"
if (params.id == null) error "The sample identifier is mandatory, please specify it with --id"
if (params.ref == null) error "The reference genome file is mandatory, please specify it with --ref"
if (params.outDir == null) error "The output directory is mandatory, please specify it with --outDir"

include { prepareBam }                           from './nextflow/bamProcessing.nf'
include { methylationCalls; checkMgmtCoverage }  from './nextflow/methylationAnalysis.nf'
include { mgmtPromoterMethyartist; mgmtPred }    from './nextflow/methylationAnalysis.nf'
include { methylationClassification; mnpFlex }   from './nextflow/methylationClassification.nf'

// Resolve --input (a single BAM, or a directory of BAMs) to a list of files
def resolveInputBams(input) {
    def target = file(input, checkIfExists: true)
    if (!target.isDirectory()) return [ target ]

    def found = file("${input}/*.bam")
    if (found instanceof Path) found = [ found ]
    if (!found) error "No BAM files found in directory: ${input}"
    return found
}

workflow {
    def ref          = file(params.ref, checkIfExists: true)
    def mgmtBed      = file("${projectDir}/data/mgmt_hg38.bed", checkIfExists: true)
    def mgmtProbes   = file("${projectDir}/data/mgmt_probes.Rdata", checkIfExists: true)
    def mgmtModel    = file("${projectDir}/data/mgmt_137sites_mean_model.Rdata", checkIfExists: true)
    def topProbes    = file("${projectDir}/data/top_probes_hm450.Rdata", checkIfExists: true)
    def trainingData = file("${projectDir}/data/capper_top_100k_betas_binarised.Rdata", checkIfExists: true)
    def arrayFile    = file("${projectDir}/data/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true)
    def mnpFlexBed   = file("${projectDir}/data/MNP-flex.bed", checkIfExists: true)

    def mgmtScript      = file("${projectDir}/scr/mgmt_pred.R", checkIfExists: true)
    def methClassScript = file("${projectDir}/scr/methylation_classification.R", checkIfExists: true)
    def mnpFlexScript   = file("${projectDir}/scr/mnp-flex_preprocessing.sh", checkIfExists: true)

    // ---- BAM preparation: tag check, align if needed, merge, sort, index
    bam = prepareBam(resolveInputBams(params.input), ref, params.id).bam

    // ---- Methylation
    meth = methylationCalls(bam, ref, params.id).bedmethyl

    methylationClassification(
        methClassScript, meth, params.id,
        topProbes, trainingData, arrayFile)

    // ---- MGMT promoter
    mgmtCov   = checkMgmtCoverage(bam, mgmtBed, params.minimumMgmtCov)
    mgmtCovOk = mgmtCov.covOkFile.map { f -> f.text.trim() }

    mgmtPromoterMethyartist(bam, ref, mgmtCovOk, params.id)
    mgmtPred(mgmtCovOk, mgmtScript, mgmtBed, mgmtProbes, mgmtModel, meth, params.id)

    // ---- Optional: MNP-Flex input preparation
    if (params.mnpFlex) {
        mnpFlex(mnpFlexScript, meth, mnpFlexBed, params.id)
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
