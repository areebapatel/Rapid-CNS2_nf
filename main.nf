#!/usr/bin/env nextflow

// Basecalling is no longer part of this pipeline.
// Please run basecalling externally using:
// nextflow run epi2me-labs/wf-basecalling ...
// and provide the resulting BAM files as input to this pipeline.
nextflow.enable.dsl = 2
// single source of truth: manifest.version in nextflow.config
software_version = workflow.manifest.version

// Display startup message
log.info """
╔══════════════════════════════════════════════════════════════════════════════╗
║                    🧬 Rapid-CNS² Nextflow Pipeline 🧬                         ║
║                              Version ${software_version}                                   ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                              ║
║  📋 Pipeline Information:                                                    ║
║  • Comprehensive CNS tumor molecular profiling                               ║
║  • SNV, CNV, SV, and methylation analysis                                    ║
║  • MGMT promoter methylation assessment                                      ║
║  • Methylation-based tumor classification                                    ║
║                                                                              ║
║  • Developer: Areeba Patel                                                   ║
║  • Email: a.patel@kitz-heidelberg.de                                         ║
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
║  • License: Apache 2.0                                                       ║
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
 * 1 - Data preparation
 *  a. Modified-base tag check, alignment (ONT Dorado), merging, sorting, indexing : samtools
 *  b. Subsetting to target region
 * 2 - Long read variant calling and annotation
 *  a. Clair3
 *  b. Annotation using ANNOVAR
 *  c. Filtering
 * 3 - Structural variants and annotation
 *  a. Sniffles2 and Severus
 *  b. Annotation using AnnotSV
 * 4 - Methylation analysis
 *  a. Methylation values using modkit
 *  b. Methylation classification using Rapid-CNS2 classifier
 *  c. MGMT promoter methylation status
 *  d. MGMT promoter region plot
 * 5 - Copy number variation using CNVpytor, plus SAVANA copy number,
 *     tumour purity and ploidy
 * 6 - Report rendering
 * (optional) Run wf-human-variation SNP and SV workflows
 * (optional) Run MNP-Flex preprocessing
 *******************************************************************************************
*/

params.input   = null
params.id      = null
params.ref     = null
params.outDir  = null
params.patient = null
params.help    = null


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
       --annotations      Path to annotation file for IGV reports (refGene.txt.gz)

   OUTPUT PARAMETERS:
       --outDir           Output directory for all results [default: output]
       --patient          Patient name for reports [default: uses --id value]
       --email            Email a summary when the run finishes or fails

   RESOURCE PARAMETERS:
       --maxCpus          Never request more CPUs than this for any one process
                          [default: 64]
       --maxMemory        Never request more memory (GB) than this [default: 96]

       Each tool is given task.cpus, so there are no per-tool thread
       parameters. Lower the two ceilings to run on smaller machines, e.g.
       --maxCpus 8 --maxMemory 32

   TOOL PARAMETERS:
       --clair3Model      Clair3 model [default: auto-detected from the BAM's
                          @RG basecall_model field; set only to override]
       --severusVntr      Optional VNTR BED for Severus (recommended)
       --severusPON       Optional panel of normals for Severus somatic filtering
       --savanaG1000      1000G SNP set for SAVANA BAF/purity [default: 1000g_hg38]

   ANALYSIS PARAMETERS:
       --minimumMgmtCov       Minimum coverage for MGMT analysis [default: 5]
       --bamMinCoverage       Minimum coverage for human variation workflow [default: 10]
       --snifflesMosaic       Run Sniffles2 in --mosaic (somatic) mode [default: true]
       --mnpFlex              Enable MNP-Flex classifier input preparation [default: true]
       --publishBam           Publish the prepared full-flow-cell BAM (>150 GB)
                              [default: false]
       --runHumanVariation    Enable wf-human-variation SNP and SV pipeline [default: false]

   CONTAINER PARAMETERS:

   PROFILES (combine one scheduler with one container engine):
       -profile lsf       Use LSF cluster scheduler
       -profile slurm     Use SLURM cluster scheduler
       -profile local     Use local execution
       -profile singularity   Run tasks with Singularity/Apptainer
       -profile docker        Run tasks with Docker

       e.g.  -profile lsf,singularity

   EXAMPLES:

   Basic run:
       nextflow run main.nf --input /data/sample.bam --id SAMPLE001 -profile lsf,singularity


   ================================================================================
   For detailed parameter descriptions, see the README.md file.
   ================================================================================
   """
    exit 0
}

// Verify that the mandatory parameters are provided
if (params.input == null) error "The path to the input or BAM file(s) is mandatory, please specify it with --input"
if (params.id == null) error "The sample identifier is mandatory, please specify it with --id"

// Site paths must be set. Fail up front rather than deep inside a tool.
['ref': params.ref, 'annovarPath': params.annovarPath,
 'annovarDB': params.annovarDB, 'annotsvAnnot': params.annotsvAnnot].each { k, v ->
    if (!v) error "--${k} is not set. Set it on the command line, in nextflow.config, " +
                  "or in your own site config passed with -c."
    if (!file(v).exists()) error "--${k} does not exist: ${v}"
}

// AnnotSV annotations are host-side (not in the container). Fail up front with
// an actionable message rather than deep inside the SV annotation step.
if (params.annotsvAnnot) {
    def asvDir = file(params.annotsvAnnot)
    if (!asvDir.exists())
        error "annotsvAnnot does not exist: ${params.annotsvAnnot}\n" +
              "Install the AnnotSV human annotations (README step 5) and point --annotsvAnnot at them."
    if (!file("${params.annotsvAnnot}/Annotations_Human").exists())
        error "annotsvAnnot must contain an Annotations_Human directory, but ${params.annotsvAnnot} does not.\n" +
              "Point it at AnnotSV's -annotationsDir (usually <install>/share/AnnotSV), not at Annotations_Human itself."
}

// ANNOVAR requires -protocol and -operation to have the same number of entries.
// Getting this wrong only surfaces as an opaque ANNOVAR error mid-run.
if (params.annovarProtocol.split(',').size() != params.annovarOperation.split(',').size())
    error "annovarProtocol (${params.annovarProtocol.split(',').size()} entries) and " +
          "annovarOperation (${params.annovarOperation.split(',').size()} entries) must be the same length"

// The MNP-Flex upload is optional and gated on credentials: it runs only when
// EPIGNOSTIX_USER and EPIGNOSTIX_PASSWORD are set. Without them the bed file is
// still produced and can be uploaded by hand at app.epignostix.com.
mnpFlexCanUpload = params.mnpFlexUpload && params.mnpFlex &&
                   params.mnpFlexWorkflowId &&
                   System.getenv('EPIGNOSTIX_USER') && System.getenv('EPIGNOSTIX_PASSWORD')

if (params.mnpFlexUpload && !mnpFlexCanUpload)
    log.warn "MNP-Flex upload skipped (no credentials or workflow id) - " +
             "upload mnpflex/*.MNPFlex.input.bed manually at app.epignostix.com"

include { prepareBam; subsetBam }                               from './nextflow/bamProcessing.nf'
include { mosdepth }                                            from './nextflow/utils.nf'
include { methylationCalls; checkMgmtCoverage }                 from './nextflow/methylationAnalysis.nf'
include { mgmtPromoterMethyartist; mgmtPred }                   from './nextflow/methylationAnalysis.nf'
include { methylationClassification; mnpFlex }                  from './nextflow/methylationClassification.nf'
include { mnpFlexUpload; mnpFlexResults }                       from './nextflow/methylationClassification.nf'
include { detectClair3Model; clair3; recodeVCF; convert2annovar } from './nextflow/variantCalling.nf'
include { tableAnnovar; filterReport; igv_reports }             from './nextflow/variantCalling.nf'
include { human_variation_snp; human_variation_sv }             from './nextflow/variantCalling.nf'
include { structuralVariants; severus }                         from './nextflow/structuralVariants.nf'
include { annotSV as annotSvSniffles }                          from './nextflow/structuralVariants.nf'
include { annotSV as annotSvSeverus }                           from './nextflow/structuralVariants.nf'
include { svFusions as fusionsSniffles }                        from './nextflow/structuralVariants.nf'
include { svFusions as fusionsSeverus }                         from './nextflow/structuralVariants.nf'
include { svFusions as fusionsSavana }                          from './nextflow/structuralVariants.nf'
include { copyNumberVariants; cnvAnnotated }                    from './nextflow/copyNumberVariants.nf'
include { savanaTo; savanaCna; savanaCnvPlot }                  from './nextflow/copyNumberVariants.nf'
include { reportRendering }                                     from './nextflow/reportRendering.nf'

// Resolve --input (a single BAM, or a directory of BAMs) to a list of files
def resolveInputBams(input) {
    def target = file(input, checkIfExists: true)
    if (!target.isDirectory()) return withIndexes([ target ])

    def found = file("${input}/*.bam")
    if (found instanceof Path) found = [ found ]
    if (!found) error "No BAM files found in directory: ${input}"
    return withIndexes(found)
}

// stage any existing .bai so prepareBam can reuse it instead of re-indexing
def withIndexes(bams) {
    bams + bams.collect { file("${it}.bai") }.findAll { it.exists() }
}

workflow {
    // ---- Static inputs. Plain values are wrapped as value channels by Nextflow,
    // ---- so they can safely feed several processes at once.
    def ref          = file(params.ref, checkIfExists: true)
    def refIdx       = file("${params.ref}.fai", checkIfExists: true)
    def panel        = file("${projectDir}/data/NPHD_panel_hg38.bed", checkIfExists: true)
    def cnvGenes     = file("${projectDir}/data/genes.bed", checkIfExists: true)
    def mgmtBed      = file("${projectDir}/data/mgmt_hg38.bed", checkIfExists: true)
    def mgmtProbes   = file("${projectDir}/data/mgmt_probes.Rdata", checkIfExists: true)
    def mgmtModel    = file("${projectDir}/data/mgmt_137sites_mean_model.Rdata", checkIfExists: true)
    def topProbes    = file("${projectDir}/data/top_probes_hm450.Rdata", checkIfExists: true)
    def trainingData = file("${projectDir}/data/capper_top_100k_betas_binarised.Rdata", checkIfExists: true)
    def arrayFile    = file("${projectDir}/data/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true)
    def mnpFlexBed   = file("${projectDir}/data/MNP-flex.bed", checkIfExists: true)
    def annotations  = file(params.annotations, checkIfExists: true)
    def noFile       = file("${projectDir}/data/NO_FILE", checkIfExists: true)
    def knownFusions = file(params.knownFusions, checkIfExists: true)
    def savanaPlotScript = file("${projectDir}/scr/plot_savana_cnv.R", checkIfExists: true)

    // Optional Severus resources; fall back to the empty placeholder when unset
    def severusVntr  = params.severusVntr ? file(params.severusVntr, checkIfExists: true) : noFile
    def severusPon   = params.severusPON  ? file(params.severusPON,  checkIfExists: true) : noFile

    def mgmtScript        = file("${projectDir}/scr/mgmt_pred.R", checkIfExists: true)
    def methClassScript   = file("${projectDir}/scr/methylation_classification.R", checkIfExists: true)
    def annotateScript    = file("${projectDir}/scr/annotate.py", checkIfExists: true)
    def filterReportScript= file("${projectDir}/scr/filter_report.R", checkIfExists: true)
    def makeReportScript  = file("${projectDir}/scr/make_report.R", checkIfExists: true)
    def reportHTML        = file("${projectDir}/scr/Rapid_CNS2_report_UKHD_HTML.Rmd", checkIfExists: true)
    def mnpFlexScript     = file("${projectDir}/scr/mnp-flex_preprocessing.sh", checkIfExists: true)

    def patientName = params.patient ?: params.id

    // ---- BAM preparation: tag check, align if needed, merge, sort, index
    bam = prepareBam(resolveInputBams(params.input), ref, params.id).bam

    subset = subsetBam(bam, panel, params.id).bam

    // ---- Coverage
    coverage = mosdepth(bam, panel, params.id)

    // ---- Methylation
    meth = methylationCalls(bam, ref, params.id).bedmethyl

    classification = methylationClassification(
        methClassScript, meth, params.id,
        topProbes, trainingData, arrayFile)

    // ---- MGMT promoter
    mgmtCov = checkMgmtCoverage(bam, mgmtBed, params.minimumMgmtCov)
    mgmtCovOk  = mgmtCov.covOkFile.map  { f -> f.text.trim() }
    mgmtAvgCov = mgmtCov.avgCovFile.map { f -> f.text.trim() }

    mgmtPlot   = mgmtPromoterMethyartist(bam, ref, mgmtCovOk, params.id)
    mgmtStatus = mgmtPred(mgmtCovOk, mgmtScript, mgmtBed, mgmtProbes, mgmtModel, meth, params.id)

    // ---- Copy number
    cnv = copyNumberVariants(bam, params.id)
    cnvAnnotated(cnv.calls1000, params.id, annotateScript, cnvGenes)

    // ---- SNVs (Clair3). The model is derived from the BAM's basecaller
    // ---- unless --clair3Model is given explicitly.
    // Nextflow rejects a null val input, so an unset override is passed as ""
    clair3Model = detectClair3Model(bam, params.clair3Model ?: '')
                    .modelFile.map { f -> f.text.trim() }
    clairVcf = clair3(subset, ref, refIdx, panel, params.id, clair3Model).vcf
    passVcf  = recodeVCF(clairVcf, params.id).passVcf
    avinput  = convert2annovar(passVcf, params.annovarPath, params.id).avinput
    anno     = tableAnnovar(avinput, params.annovarPath, params.annovarDB, params.id,
                            params.annovarProtocol, params.annovarOperation).multianno
    snvTable = filterReport(filterReportScript, anno, params.id).dvReport

    igv = igv_reports(snvTable, params.id, ref, subset, annotations)

    // ---- Structural variants: Sniffles2 (fast baseline) + Severus (somatic/complex)
    sniffles = structuralVariants(bam, ref, params.id)
    annotSvSniffles(sniffles.vcf, params.annotsvAnnot, params.id, 'sniffles')

    sev = severus(bam, params.id, severusVntr, severusPon)
    annotSvSeverus(sev.vcf, params.annotsvAnnot, params.id, 'severus')

    // ---- SAVANA: tumour-only breakpoints -> copy number, purity and ploidy
    savanaSv  = savanaTo(bam, ref, refIdx, params.id)
    savCna = savanaCna(bam, ref, savanaSv.breakpoints, params.id, params.savanaG1000)
    savanaPlot = savanaCnvPlot(savCna.segments, savCna.purityPloidy,
                               cnvGenes, savanaPlotScript, params.id)

    fSni = fusionsSniffles(sniffles.vcf, annotations, knownFusions, params.id, 'sniffles')
    fSev = fusionsSeverus(sev.vcf, annotations, knownFusions, params.id, 'severus')
    fSav = fusionsSavana(savanaSv.breakpoints, annotations, knownFusions, params.id, 'savana')

    // one table per caller; the report merges them and shows which callers agree
    fusionTables = fSni.reportable.mix(fSev.reportable, fSav.reportable).collect()
    egfrTables   = fSni.egfrviii.mix(fSev.egfrviii, fSav.egfrviii).collect()

    // ---- Optional: ONT wf-human-variation
    if (params.runHumanVariation) {
        human_variation_snp(bam, panel, ref, params.id, params.bamMinCoverage, params.snpThreads)
        human_variation_sv(bam, ref, params.id, params.bamMinCoverage, params.svThreads)
    }

    // ---- Optional: MNP-Flex input preparation
    flexPredictions = Channel.empty()
    if (params.mnpFlex) {
        flexBed = mnpFlex(mnpFlexScript, meth, mnpFlexBed, params.id)
        if (mnpFlexCanUpload) {
            def uploadScript = file("${projectDir}/scr/mnpflex_upload.py", checkIfExists: true)
            def resultsScript = file("${projectDir}/scr/mnpflex_results.py", checkIfExists: true)
            up = mnpFlexUpload(uploadScript, flexBed.mnpFlexBed, params.id)
            flexPredictions = mnpFlexResults(resultsScript, up.response, params.id).predictions
        }
    }

    // ---- Final report. Optional inputs fall back to an empty placeholder file.
    reportRendering(
        makeReportScript,
        reportHTML,
        params.id,
        patientName,
        software_version,
        fusionTables,
        egfrTables,
        savCna.purityPloidy.ifEmpty(noFile),
        savanaPlot.plot.ifEmpty(noFile),
        flexPredictions.ifEmpty(noFile),
        snvTable,
        cnv.plotPng,
        classification.rfDetails,
        classification.votes,
        coverage.summary,
        mgmtStatus.mgmtStatus.ifEmpty(noFile),
        mgmtPlot.mgmtPlot.ifEmpty(noFile),
        igv.igvReport.ifEmpty(noFile),
        mgmtAvgCov
    )
}

workflow.onComplete {
    if (workflow.success) {
        println "The Rapid-CNS2 workflow is now complete!\n Your outputs are located in : ${params.outDir}"
    } else {
        println "Oops .. something went wrong, please look into the log file, and error messages into ${workflow.workDir}"
    }

    // Drop per-task scratch that no process declares as an output: SAVANA's
    // per-window .temp files and the matplotlib/fontconfig caches. Safe for
    // -resume, which only checks declared outputs. Use -profile cleanup to
    // remove the work directory outright.
    if (workflow.success) {
        try {
            def files = [], dirs = [], freed = 0L
            // collect first: deleting during the walk breaks the traversal
            workflow.workDir.toFile().eachFileRecurse { f ->
                if (f.isFile() && f.name.endsWith('.temp')) files << f
                else if (f.isDirectory() && f.name in ['.mpl', '.cache']) dirs << f
            }
            files.each { freed += it.length(); it.delete() }
            dirs.each { it.deleteDir() }
            if (freed > 0) println "Removed ${(freed / 1024 / 1024) as long} MB of task scratch"
        }
        catch (Exception e) { println "Scratch cleanup skipped: ${e.message}" }
    }

    // Notify if --email is given. Done here rather than in the notification{}
    // config block, which is evaluated before --params are merged.
    if (!params.email) return
    def ok = workflow.success
    try {
        sendMail(
            to: params.email,
            subject: "Rapid-CNS2 ${ok ? 'complete' : 'FAILED'} - ${params.id}",
            body: """\
            Sample   : ${params.id}
            Status   : ${ok ? 'success' : 'failed'}
            Duration : ${workflow.duration}
            Command  : ${workflow.commandLine}
            Results  : ${params.outDir}
            Work dir : ${workflow.workDir}
            ${ok ? '' : '\n' + (workflow.errorMessage ?: '') + '\n\n' + (workflow.errorReport ?: '')}
            """.stripIndent()
        )
        println "Notification sent to ${params.email}"
    }
    catch (Exception e) {
        println "Could not send notification email: ${e.message}"
    }
}
