#!/usr/bin/env nextflow

// Basecalling is no longer part of this pipeline.
// Please run basecalling externally using:
// nextflow run epi2me-labs/wf-basecalling ...
// and provide the resulting BAM files as input to this pipeline.
nextflow.enable.dsl = 2
software_version = "3.0.0"

// Display startup message
log.info """
╔══════════════════════════════════════════════════════════════════════════════╗
║                    🧬 Rapid-CNS² Nextflow Pipeline 🧬                         ║
║                              Version ${software_version}                     ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                              ║
║  📋 Pipeline Information:                                                    ║
║  • Comprehensive CNS tumor molecular profiling                               ║
║  • SNV, CNV, SV, and methylation analysis                                    ║
║  • MGMT promoter methylation assessment                                      ║
║  • Methylation-based tumor classification                                    ║
║                                                                              ║
║  • Developer: Areeba Patel                                                   ║
║  • Email: a.patel@dkfz.de                                                    ║
║  • Institution: German Cancer Research Center (DKFZ)                         ║  
║                                                                              ║
║  📚 Citation (Please cite if used):                                           ║
║  Patel, A., Göbel, K., Ille, S. et al. Prospective, multicenter validation   ║
║  of a platform for rapid molecular profiling of central nervous system       ║
║  tumors. Nature Medicine 31, 1567–1577 (2025).                               ║
║  DOI: https://doi.org/10.1038/s41591-025-03562-5                             ║
║                                                                              ║
║  🔗 Additional Resources:                                                    ║
║  • GitHub Repository: https://github.com/areebapatel/Rapid-CNS2_nf           ║
║  • License: MIT (Open Source)                                                ║
║  • Documentation: See README.md for detailed usage instructions              ║
║                                                                              ║
║  ⚠️  Important Notes:                                                        ║
║  • This pipeline is for RESEARCH USE ONLY                                    ║
║  • Not validated for clinical diagnostic use                                 ║
║  • Results should be interpreted by qualified professionals                  ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""
/**
 ********************************** Rapid-CNS2 NextFlow ******************************************
 * 1 - Base calling, alligment and data preparation
 *  a. Base calling + alignment : ONT Dorado and minimap2
 *  b. Sorting, adding read group information and creating index : samtools
 *  c. Subsetting to target region
 * 2 - Long read variant calling and annotation
 *  a. Clara Parabricks Deepvariant 
 *  b. Annotation using ANNOVAR
 *  c. Filtering
 * 3 - Structural variants and annotation 
 *  a. Sniffles2
 *  b. Annotation using AnnotSV
 * 4 - Methylation analysis
 *  a. Methylation values using modkit
 *  b. Methylation classification using Rapid-CNS2 classifier
 *  c. MGMT promoter methylation status
 *  d. MGMT promoter region plot
 * 5 - Copy number variation calling using CNVpytor
 * 6 - Report rendering
 * (optional) Run wf-human-variation SNP and SV workflows
 * (optional) Run MNP-Flex preprocessing
 *******************************************************************************************
*/

   
//includeConfig './nextflow.config'
params.input = null
params.ref = null
//params.out_dir = "output"
params.logDir = "logDir"
params.minimumMgmtCov = 5
params.outDir = null
params.numGpu = 1

params.rParams = [] // Initialize an empty list to store -r parameters

// DeepVariant mode, adding dummy read group information to the BAM file
params.pbDVMode = "ont"
params.pbPATH = "pbrun"
params.reads = ' "ID:12345" -r "SM:12345" -r "PL:ONT" '

// set up and create an output directory
//out_dir = path(params.outDir)
//out_dir.mkdir()

params.help = null
params.test = null

// Patient name (for report)
params.patient = null

// Show help message
if (params.help) {
   log.info """
   ================================================================================
   Rapid-CNS² Nextflow Pipeline v${software_version}
   ================================================================================
   
   USAGE: nextflow run main.nf [options]
   
   MANDATORY PARAMETERS:
       --input            Path to input BAM file(s)
                         • Single aligned BAM: /path/to/sample.bam
                         • Directory with aligned BAMs: /path/to/aligned_bams/
                         • Directory with unaligned BAMs: /path/to/unaligned_bams/
       --id               Sample identifier (alphanumeric, no spaces)
   
   SYSTEM-SPECIFIC PARAMETERS (configure in nextflow.config):
       --ref              Path to hg38 reference genome FASTA file
       --annovarPath      Path to ANNOVAR installation directory
       --annovarDB        Path to ANNOVAR database directory (humandb/)
       --annotsvAnnot     Path to AnnotSV annotations directory
       --annotations      Path to annotation file for IGV reports (refGene.txt)
   
   OUTPUT PARAMETERS:
       --outDir           Output directory for all results [default: output]
       --tmpDir           Directory for temporary files [default: \${outDir}/tmp/]
       --logDir           Directory for log files [default: logDir]
       --patient          Patient name for reports [default: uses --id value]
   
   RESOURCE PARAMETERS:
       --maxThreads       Maximum threads for general processes [default: 64]
       --modkitThreads    Threads for modkit methylation calling [default: 32]
       --cnvThreads       Threads for CNVpytor copy number analysis [default: 32]
       --snifflesThreads  Threads for Sniffles2 structural variant calling [default: 32]
       --snpThreads       Threads for SNV calling with DeepVariant [default: 64]
       --svThreads        Threads for structural variant calling [default: 64]
       --covThreads       Threads for coverage calculation [default: 8]
       --methThreads      Threads for methylation classification [default: 64]
       --mgmtThreads      Threads for MGMT promoter analysis [default: 8]
   
   ANALYSIS PARAMETERS:
       --minimumMgmtCov   Minimum coverage for MGMT analysis [default: 5]
       --bamMinCoverage   Minimum coverage for human variation workflow [default: 10]
       --mnpFlex          Enable MNP-Flex classifier input preparation [default: false]
       --runHumanVariation Enable wf-human-variation SNP and SV pipeline [default: false]
   
   CONTAINER PARAMETERS:
       --containerEngine  Container engine: 'docker' or 'singularity' [default: docker]
       --seq              Sequencer platform identifier [default: P2]
   
   PROFILES:
       -profile lsf       Use LSF cluster scheduler
       -profile slurm     Use SLURM cluster scheduler  
       -profile local     Use local execution
   
   EXAMPLES:
   
   Basic run:
       nextflow run main.nf --input /data/sample.bam --id SAMPLE001 -profile lsf

   
   ================================================================================
   For detailed parameter descriptions, see the README.md file.
   ================================================================================
   """
    exit 0
}

// Verify that the mandatory parameters are provided
if (params.input == null) error "The path to the input or BAM file(s) is mandatory, please specify it with --input"
if (params.id == null) error "The sample identifier is mandatory, please specify it with --id"
if (params.ref == null) error "The reference genome file is mandatory, please specify it with --ref"

include { subsetBam } from './nextflow/bamProcessing.nf'
include { indexBam } from './nextflow/bamProcessing.nf'
include { indexSubsettedBam } from './nextflow/bamProcessing.nf'
include { addreplacerg } from './nextflow/bamProcessing.nf'
include { mosdepth } from './nextflow/utils.nf'
include { methylationCalls } from './nextflow/methylationAnalysis.nf'
include { checkMgmtCoverage } from './nextflow/methylationAnalysis.nf'
include { mgmtPromoterMethyartist } from './nextflow/methylationAnalysis.nf'
include { mgmtPred } from './nextflow/methylationAnalysis.nf'
include { mnpFlex } from './nextflow/methylationClassification.nf'
include { variantCalling } from './nextflow/variantCalling.nf'
include { recodeVCF } from './nextflow/variantCalling.nf'
include { convert2annovar } from './nextflow/variantCalling.nf'
include { tableAnnovar } from './nextflow/variantCalling.nf'
include { filterReport } from './nextflow/variantCalling.nf'
include { human_variation_sv } from './nextflow/variantCalling.nf'
include { human_variation_snp } from './nextflow/variantCalling.nf'
include { igv_reports } from './nextflow/variantCalling.nf'
include { structuralVariants } from './nextflow/structuralVariants.nf'
include { annotSV } from './nextflow/structuralVariants.nf'
include { methylationClassification } from './nextflow/methylationClassification.nf'
include { cnvAnnotated } from './nextflow/copyNumberVariants.nf'
include { copyNumberVariants } from './nextflow/copyNumberVariants.nf'
include { reportRendering } from './nextflow/reportRendering.nf'
include { checkAlignment } from './nextflow/bamProcessing.nf'
include { checkMethylationTags } from './nextflow/bamProcessing.nf'
include { alignBam } from './nextflow/bamProcessing.nf'
include { mergeBam } from './nextflow/bamProcessing.nf'

// Check if the input is a BAM file or a directory containing BAM files
def isBamFile(path) {
    return path.toString().toLowerCase().endsWith('.bam')
}

// Create all required output directories before running any processes
[
    "${params.outDir}/bam",
    "${params.outDir}/bam/alignedBams",
    "${params.outDir}/snv",
    "${params.outDir}/cnv",
    "${params.outDir}/sv",
    "${params.outDir}/mods",
    "${params.outDir}/mgmt",
    "${params.outDir}/methylation_classification",
    "${params.outDir}/mnpflex",
    "${params.outDir}/coverage",
    "${params.outDir}/report",
    "${params.outDir}/reports",
    "${params.outDir}/wf-human-variation",
    "${params.outDir}/wf-human-variation/sv",
    "${params.outDir}/wf-human-variation/snp",
    "${params.outDir}/pipeline_info"
].each { dir ->
    new File(dir).mkdirs()
}

workflow {
    // Set the reference genome
    Channel.fromPath(params.ref, checkIfExists: true)
    .set {ref}

    // Set the sample identifier
    Channel.from(params.id)
    .set {id}

    // Set the patient name - use id if not specified
    def patientName = params.patient ?: params.id
    Channel.from(patientName)
    .set {patient}

    // Set the output directory
    Channel.from(params.outDir)
    .set {outDir}

    // Set the temporary directory
    Channel.from(params.tmpDir)
    .set {tmpDir}

    // Set the maximum number of threads
    Channel.from(params.maxThreads)
    .set {maxThreads}
    

    Channel.from(params.pbDVMode)
    .set {pbDVMode}

    Channel.from(params.pbPATH)
    .set {pbPATH}

    // Set the annovar path
    Channel.from(params.annovarPath)
    .set {annovarPath}

    // Set the annovar database
    Channel.from(params.annovarDB)
    .set {annovarDB}

    // Set the annotsv annot path
    Channel.from(params.annotsvAnnot)
    .set {annotsvAnnot}

    // Set the modkit threads
    Channel.from(params.modkitThreads)
    .set {modkitThreads}

    // Set the cnv threads
    Channel.from(params.cnvThreads)
    .set {cnvThreads}

    // Set the sniffles threads
    Channel.from(params.snifflesThreads)
    .set {cnvThreads}

    // Set the snp threads
    Channel.from(params.snpThreads)
    .set {snpThreads}

    // Set the sv threads
    Channel.from(params.svThreads)
    .set {svThreads}

    // Set the cov threads
    Channel.from(params.covThreads)
    .set {covThreads}

    // Set the meth threads
    Channel.from(params.methThreads)
    .set {methThreads}

    // Set the mgmt threads
    Channel.from(params.mgmtThreads)
    .set {mgmtThreads}
    
    // Set the sniffles threads
    Channel.from(params.snifflesThreads)
    .set {snifflesThreads}

    // Set the minimum mgmt coverage
    Channel.from(params.minimumMgmtCov)
    .set {minimumMgmtCov}

    // Set the annotations
    Channel.fromPath(params.annotations, checkIfExists: true)
    .set {annotations}

    // Set the panel
    Channel.fromPath("${projectDir}/data/NPHD_panel_hg38.bed", checkIfExists: true)
    .set {panel}

    // Set the cnv genes
    Channel.fromPath("${projectDir}/data/genes.bed", checkIfExists: true)
    .set {cnvGenes}

    // Set the mgmt bed
    Channel.fromPath("${projectDir}/data/mgmt_hg38.bed", checkIfExists: true)
    .set {mgmtBed}

    // Set the mgmt probes
    Channel.fromPath("${projectDir}/data/mgmt_probes.Rdata", checkIfExists: true)
    .set {mgmtProbes}

    // Set the mgmt model
    Channel.fromPath("${projectDir}/data/mgmt_137sites_mean_model.Rdata", checkIfExists: true)
    .set {mgmtModel}

    // Set the mgmt script
    Channel.fromPath("${projectDir}/scr/mgmt_pred.R", checkIfExists: true)
    .set {mgmtScript}

    // Set the methylation classification script
    Channel.fromPath("${projectDir}/scr/methylation_classification.R", checkIfExists: true)
    .set {methylationClassificationScript}
    
    // Set the CNV annotation script
    Channel.fromPath("${projectDir}/scr/annotate.py", checkIfExists: true)
    .set {annotateScript}

    // Set the top probes
    Channel.fromPath("${projectDir}/data/top_probes_hm450.Rdata", checkIfExists: true)
    .set {topProbes}
    
    // Set the training data
    Channel.fromPath("${projectDir}/data/capper_top_100k_betas_binarised.Rdata", checkIfExists: true)
    .set {trainingData}

    // Set the array file
    Channel.fromPath("${projectDir}/data/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true)
    .set {arrayFile}

    // Set the filter report script
    Channel.fromPath("${projectDir}/scr/filter_report.R", checkIfExists: true)
    .set {filterReportScript}

    // Set the make report script
    Channel.fromPath("${projectDir}/scr/make_report.R", checkIfExists: true)
    .set {makereport}

    // Set the report UKHD script
    Channel.fromPath("${projectDir}/scr/Rapid_CNS2_report_UKHD_v3.0.0.Rmd", checkIfExists: true)
    .set {reportUKHD}

    // Set the logos directory
    Channel.fromPath("${projectDir}/logos", checkIfExists: true)
    .set {logosDir}

    def inputPath = file(params.input)
    def bamToCheck
    // Check if the input is a BAM file or a directory containing BAM files
    if (isBamFile(inputPath)) {
        bamToCheck = Channel.value(inputPath)
    } else {
        bamToCheck = Channel.fromPath("${inputPath}/*.bam").first()
    }
    
    // Check alignment status
    checkAlignment_out = checkAlignment(bamToCheck, params.maxThreads)
    
    // Check if the BAM file has methylation tags
    checkMethylationTags_out = checkMethylationTags(bamToCheck, params.maxThreads)

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
        alignedBams = alignBam(inputPath, ref, params.maxThreads, outDir).alignedBam
        processedBam = alignedBams.collect().map { bamList ->
            if (bamList.size() > 1) {
                mergeBam(bamList, params.maxThreads, outDir, id).mergedBam
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
            processedBam = mergeBam(inputPath, params.maxThreads, outDir, id).mergedBam
        }
    }

    // Index the processed BAM  
    indexedBam = indexBam(processedBam, params.maxThreads)

    // All following processes are run on the processed BAM in parallel

    // Subset BAM to target regions
    subsettedBam = subsetBam(processedBam, indexedBam.indexBam, panel, id, params.maxThreads)

    subsetIndex = indexSubsettedBam(subsettedBam.subsetBam, params.maxThreads)

    // Coverage calculation
    coverageOut = mosdepth(params.covThreads, panel, processedBam, indexedBam.indexBam, id)

    // Call methylation
    methylationCalls = methylationCalls(processedBam, indexedBam.indexBam, ref, id,  params.modkitThreads)

    // Methylation classification
    methylationClassification = methylationClassification(methylationClassificationScript, methylationCalls.bedmethylFile, id, topProbes, trainingData, arrayFile, params.methThreads)

    //MGMT promoter
    mgmtCoverageOut = checkMgmtCoverage(processedBam, indexedBam.indexBam, mgmtBed, params.minimumMgmtCov, params.mgmtThreads)

    mgmtPromoterOut = mgmtPromoterMethyartist(processedBam, indexedBam.indexBam, ref, mgmtCoverageOut[0], id)

    mgmtPredOut = mgmtPred(mgmtCoverageOut[0], mgmtScript, mgmtBed, mgmtProbes, mgmtModel, methylationCalls.bedmethylFile, id)

    // CNV calling
    cnvOut = copyNumberVariants(processedBam, indexedBam.indexBam, id, params.cnvThreads)
    cnvAnnotatedOut = cnvAnnotated(cnvOut.cnvpytorCalls1000, id, annotateScript, cnvGenes, params.outDir)

    // SNV calling
    addReplaceRgOut = addreplacerg(subsettedBam.subsetBam, subsetIndex.indexSubsetBam)
    variantCallingOut = variantCalling(addReplaceRgOut.deepVariantBam, addReplaceRgOut.deepVariantBai, ref, id, params.tmpDir, params.numGpu)
    recodeVCFOut = recodeVCF(variantCallingOut.dvVcf)

    // ANNOVAR
    convert2annovarOut = convert2annovar(recodeVCFOut.passVcf, annovarPath)
    tableAnnovarOut = tableAnnovar(convert2annovarOut.annovarInput, annovarPath, annovarDB)
    filterReportOut = filterReport(filterReportScript, tableAnnovarOut.dvAnno, id)

    // IGV reports
    igvReportsOut = igv_reports(filterReportOut.dvReport, id, ref, subsettedBam.subsetBam, subsetIndex.indexSubsetBam, annotations)

    // SV calling
    structuralVariantsOut = structuralVariants(processedBam, subsetIndex.indexSubsetBam, ref, id, snifflesThreads)
    annotSVOut = annotSV(structuralVariantsOut.svVcf, annotsvAnnot, id)

    if (params.runHumanVariation){
    // Human variation SNP workflow //not included in report yet
        human_variation_snp(processedBam, panel, ref, id, outDir, params.bamMinCoverage, params.snpThreads)

    // Human variation SV workflow // not included in report yet
        human_variation_sv(processedBam, ref, id, params.svThreads)
    }

    // Final report
    reportRenderingOut = reportRendering(makereport, cnvOut.cnvpytorPlot, mgmtPredOut, methylationClassification, filterReportOut[0], id, coverageOut.mosdepthOut, mgmtCoverageOut[3], mgmtPromoterOut, igvReportsOut, software_version, processedBam, params.seq, reportUKHD, logosDir)

    if ( params.mnpFlex) {
        mnpFlex(mnpFlexScript, methylationCalls.bedmethylFile, mnpFlexBed, id)
    }
}

workflow.onComplete {
    if(workflow.success) {
    println ( "The Rapid-CNS2 workflow is now complete!\n Your outputs are located in : " + params.outDir )
    }
    else {
    println ( "Oops .. something went wrong, please look into the log file, and error messages into " + workDir )
    }
}
