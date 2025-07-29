#!/usr/bin/env nextflow

// Basecalling is no longer part of this pipeline.
// Please run basecalling externally using:
// nextflow run epi2me-labs/wf-basecalling ...
// and provide the resulting BAM files as input to this pipeline.
nextflow.enable.dsl = 2
software_version = "3.0.0"
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

// DeepVariant mode. 
params.pbDVMode = "ont"
params.pbPATH = "pbrun"

// set up and create an output directory
//out_dir = path(params.outDir)
//out_dir.mkdir()

params.help = null
params.test = null
params.reads = null

// Show help message
if (params.help) {
   log.info """
    Usage: nextflow run nextflow/main.nf  [options]

    Mandatory options:
        --input            Path to the directory containing Pod5/FAST5 file or BAM file if available
        --id                Sample identifier

    Options:
        --outDir          Directory path to store all the outputs. [default : ${params.outDir}]
        --ref             Path to hg38 reference file. [default : ${params.ref}]
        --tmpDir          Directory to store temporary files. If it does not exists it will be created. [default : ${params.tmpDir}]
        --sortThreads     Number of threads to use for samtools sort [default : ${params.sortThreads}]
        --snifflesThreads  Numbers of threads to run sniffles2 [default : ${params.snifflesThreads}]
        --cnvThreads Numbers of threads to run cnvpytor [default : ${params.cnvThreads}]
        --modThreads Numbers of threads to run modkit [default : ${params.modThreads}]
        --modelConfig     Basecalling model to be used [default : ${params.modelConfig}]
        --remoraConfig    Remora model to be used [default: ${params.remoraConfig}]
        --port             Port for basecall server [default : ${params.port}] 
        --reads            samtools addreplacerg -r option. It should be specified as this example :  --reads '-r "SM:GM24385" -r "ID:GM24385"'
        --mnpFlex          Run MNP-Flex preprocessing [default : FALSE]
        --runHumanVariation   Run the wf-human-variation SNP and SV pipeline [default : FALSE]
        --numGpu           Number of GPUs to use [default: ${params.numGpu}]

    Example:
      nextflow run main.nf  --input "/input/pod5" --ref "/Ref/Homo_sapiens_assembly38.fasta" --id "Sample_XYZ"

    To run with LSF, add -process.executor='lsf' to your nextflow command.

    To also prepare input file for the MNP-Flex classifier, add the --mnpFlex flag. (Default behaviour is to not prepare the input file)
    """
    exit 0
}

// Verify that the mandatory parameters are provided
if (params.input == null) error "The path to the input POD5 files or BAM file is mandatory, please specify it with --input"
if (params.id == null) error "The sample identifier is mandatory, please specify it with --id"
if (params.ref == null) error "The reference genome file is mandatory, please specify it with --ref"

include { subsetBam } from './nextflow/bamProcessing.nf'
include { indexBam } from './nextflow/bamProcessing.nf'
include { indexSubsettedBam } from './nextflow/bamProcessing.nf'
include { mosdepth } from './nextflow/utils.nf'
include { methylationCalls } from './nextflow/methylationAnalysis.nf'
include { checkMgmtCoverage } from './nextflow/methylationAnalysis.nf'
include { mgmtPromoterMethylartist } from './nextflow/methylationAnalysis.nf'
include { mgmtPred } from './nextflow/methylationAnalysis.nf'
include { mnpFlex } from './nextflow/methylationAnalysis.nf'
include { deepVariant } from './nextflow/variantCalling.nf'
include { recodeVCF } from './nextflow/variantCalling.nf'
include { convert2annovar } from './nextflow/variantCalling.nf'
include { tableAnnovar } from './nextflow/variantCalling.nf'
include { filterReport } from './nextflow/variantCalling.nf'
include { humanVariationSv } from './nextflow/variantCalling.nf'
include { humanVariationSnp } from './nextflow/variantCalling.nf'
include { igvReports } from './nextflow/variantCalling.nf'
include { sniffles2 } from './nextflow/structuralVariants.nf'
include { annotsvAnnot } from './nextflow/structuralVariants.nf'
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

workflow {
    // Set the reference genome
    Channel.fromPath(params.ref, checkIfExists: true)
    .set {ref}

    // Set the sample identifier
    Channel.from(params.id)
    .set {id}

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
    Channel.fromPath("${projectDir}/data/NPHD_panel.bed", checkIfExists: true)
    .set {panel}

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

    def inputPath = file(params.input)
    def bamToCheck
    // Check if the input is a BAM file or a directory containing BAM files
    if (isBamFile(inputPath)) {
        bamToCheck = Channel.value(inputPath)
    } else {
        bamToCheck = Channel.fromPath("${inputPath}/*.bam").first()
    }
    // Check alignment
    checkAlignment_out = checkAlignment(bamToCheck, params.threads)
    // Check if the BAM file has methylation tags
    checkMethylationTags_out = checkMethylationTags(bamToCheck, params.threads)

    // Check if methylation tags exist; if not, raise an error and provide guidance
    checkMethylationTags_out.view { result ->
        if (!result || result.trim().isEmpty()) {
            error "No methylation tags (MM:Z) found in the BAM file. Please redo basecalling with modified basecalling enabled. See: https://github.com/nanoporetech/dorado?tab=readme-ov-file#modified-basecalling or ensure modified bases (5mC) are enabled in MinKNOW."
        }
    }

    def doAlignment = false
    checkAlignment_out.view { result ->
        if (!result || result.trim().isEmpty()) {
            println "Input BAM file(s) are unaligned. Alignment will be performed."
            doAlignment = true
        }
    }

    // Conditionally run alignment and merging only if needed
    def processedBam
    if (doAlignment) {
        alignedBams = alignBam(inputPath, ref, params.threads, outDir).alignedBam
        alignedBams.collect().set { bamList }
        def mergedOrSingleBam = bamList.count().map { n ->
            if (n > 1) {
                return mergeBam(bamList, params.threads, outDir, id).mergedBam
            } else if (n == 1) {
                return Channel.value(bamList[0])
            } else {
                error "No aligned BAM files found after alignment."
            }
        }
        processedBam = mergedOrSingleBam.flatten()
    } else {
        inputBams = Channel.fromPath("${inputPath}/*.bam")
        inputBams.collect().set { bamList }
        def mergedOrSingleBam = bamList.count().map { n ->
            if (n > 1) {
                return mergeBam(bamList, params.threads, outDir, id).mergedBam
            } else if (n == 1) {
                return Channel.value(bamList[0])
            } else {
                error "No BAM files found in input directory."
            }
        }
        processedBam = mergedOrSingleBam.flatten()
    }

    // Index the processed BAM  
    indexBam(processedBam, params.threads)

    // All following processes are run on the processed BAM in parallel

    // Subset BAM to target regions
    subsettedBam = subsetBam(processedBam, indexBam.indexBam, panel, id, params.threads)

    subsetIndex = indexSubsettedBam(subsettedBam.subsetBam, params.threads)

    // Coverage calculation
    mosdepth(params.covThreads, panel, processedBam, indexBam.indexBam, id)

    // Call methylation
    methylationCalls = methylationCalls(processedBam, subsetIndex.indexSubsetBam, ref, id,  params.modkitThreads)

    // Methylation classification
    methylationClassification = methylationClassification(methylationClassificationScript, methylationCalls.bedmethylFile, id, topProbes, trainingData, arrayFile, params.methThreads)

    //MGMT promoter
    checkMgmtCoverage(processedBam, subsetIndex.indexSubsetBam, mgmtBed, params.minimumMgmtCov, params.mgmtThreads)

    mgmtPromoterMethylartist(processedBam, subsetIndex.indexSubsetBam, ref, checkMgmtCoverage.out[0], id)

    mgmtPred(checkMgmtCoverage.out[0], mgmtScript, mgmtBed, mgmtProbes, mgmtModel, methylationCalls.bedmethylFile, id)

    // CNV calling
    copyNumberVariants(processedBam, subsetIndex.indexSubsetBam, id, params.cnvThreads)
    cnvAnnotated(copyNumberVariants.out, id, annotateScript, params.outDir)

    // SNV calling
    deepVariant(subsettedBam.subsetBam, subsetIndex.indexSubsetBam, ref, id, params.pbDVMode, params.pbPATH, params.tmpDir, params.numGpu)
    recodeVCF(deepVariant.dvVcf)

    // ANNOVAR
    convert2annovar(recodeVCF.passVcf, annovarPath)
    tableAnnovar(convert2annovar.annovarInput, annovarPath, annovarDB)
    filterReport(filterReportScript, tableAnnovar.dvAnno, id)

    // IGV reports
    igvReports(filterReport.dvReport, id, ref, subsettedBam.subsetBam, subsetIndex.indexSubsetBam, annotations)

    // SV calling
    sniffles2(processedBam, subsetIndex.indexSubsetBam, ref, id, params.snifflesThreads)
    annotsvAnnot(sniffles2.svVcf, id)

    if (params.runHumanVariation){
    // Human variation SNP workflow //not included in report yet
        humanVariationSnp(processedBam, panel, ref, id, outDir, params.bamMinCoverage, params.snpThreads)

    // Human variation SV workflow // not included in report yet
        humanVariationSv(processedBam, ref, id, params.svThreads)
    }

    // Final report
    makeReport(makereport, copyNumberVariants.out, mgmtPred.out, methylationClassification.out, filterReportCh.out, id, mosdepth.mosdepthOut, checkMgmtCoverage.out[0], mgmtPromoterMethylartist.out, igvReports.out, nextflowVersion, processedBam, params.seq, reportUKHD)

    if ( params.mnpFlex) {
        mnpFlex(mnpFlexScript, methylationCalls.bedmethylFile, params.mnpFlexBed, id)
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
