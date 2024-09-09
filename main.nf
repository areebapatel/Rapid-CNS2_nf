#!/usr/bin/env nextflow

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
params.log_dir = "logDir"
params.minimum_mgmt_cov = 5
params.port= 8887
params.num_gpu = 3
params.num_clients = params.num_gpu * 3
params.outDir = params.out_dir

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
        --input            Path to the directory containing Pod5/FAST5 file for Dorado basecalling and minimap2 alligement or BAM file if available
        --id                Sample identifier

    Options:
        --out_dir          Directory path to store all the outputs. [default : ${params.out_dir}]
        --ref              Path to hg19 reference file. [default : ${params.ref}]
        --tmp_dir          Directory to store temporary files. If it does not exists it will be created. [default : ${params.tmp_dir}]
        --basecalling      If data should be basecalled  // remove
        --sort_threads     Number of threads to use for samtools sort [default : ${params.sort_threads}]
        
        --sniffles_threads  Numbers of threads to run sniffles2 [default : ${params.sniffles_threads}]
        --cnv_threads Numbers of threads to run cnvpytor [default : ${params.cnv_threads}]
        --mod_threads Numbers of threads to run modkit [default : ${params.mod_threads}]
        --model_config     Basecalling model to be used [default : ${params.model_config}]
        --remora_config    Remora model to be used [default: ${params.remora_config}]
        --port             Port for basecall server [default : ${params.port}] 
        --reads            samtools addreplacerg -r option. It should be specified as this example :  --reads '-r "SM:GM24385" -r "ID:GM24385"'
        --mnp-flex          Run MNP-Flex preprocessing [default : FALSE]
        --run_human_variation   Run the wf-human-variation SNP and SV pipeline [default : FALSE]

    Example:
      nextflow run main.nf  --input "/input/pod5" --ref "/Ref/Homo_sapiens_assembly38.fasta" --id "Sample_XYZ"

    To run with LSF, add -process.executor='lsf' to your nextflow command.

    To also prepare input file for the MNP-Flex classifier, add the --mnp-flex flag. (Default behaviour is to not prepare the input file)
    """
    exit 0
}

// Verify that the mandatory parameters are provided
//if (params.basecalling & params.ref == null) error "The reference genome file is mandatory for basecalling. Please specify it with --ref"
//if (params.input == null) error "The path to the input POD5 files or BAM file is mandatory, please specify it with --input"

include { basecalling } from './nextflow/basecalling.nf'
include { subsetBam } from './nextflow/basecalling.nf'
include { indexBam } from './nextflow/basecalling.nf'
include { indexSubsettedBam } from './nextflow/basecalling.nf'
include { mosdepth } from './nextflow/utils.nf'
include { methylationCalls } from './nextflow/methylationAnalysis.nf'
include { liftOver_ch } from './nextflow/methylationAnalysis.nf'
include { check_mgmt_coverage } from './nextflow/methylationAnalysis.nf'
include { mgmtPromoter_methyartist } from './nextflow/methylationAnalysis.nf'
include { mgmtPred } from './nextflow/methylationAnalysis.nf'
include { mnpFlex } from './nextflow/methylationAnalysis.nf'
include { deepVariant } from './nextflow/variantCalling.nf'
include { recodeVCF } from './nextflow/variantCalling.nf'
include { convert2annovar } from './nextflow/variantCalling.nf'
include { tableAnnovar } from './nextflow/variantCalling.nf'
include { filterReport } from './nextflow/variantCalling.nf'
include { human_variation_sv } from './nextflow/variantCalling.nf'
include { human_variation_snp } from './nextflow/variantCalling.nf'
include { igv_reports } from './nextflow/variantCalling.nf'
include { sniffles2 } from './nextflow/structuralVariants.nf'
include { annotsvAnnot } from './nextflow/structuralVariants.nf'
include { methylationClassification } from './nextflow/methylationClassification.nf'
include { cnvAnnotated } from './nextflow/copyNumberVariants.nf'
include { copyNumberVariants } from './nextflow/copyNumberVariants.nf'
include { reportRendering } from './nextflow/reportRendering.nf'

workflow {
    Channel.fromPath(params.input, checkIfExists: true)
    .set {input}

    Channel.fromPath(params.ref, checkIfExists: true)
    .set {ref}

    Channel.from(params.id)
    .set {id}

    Channel.from(params.out_dir)
    .set {out_dir}

    Channel.from(params.outDir)
    .set {out_dir}
    
    Channel.from(params.tmp_dir)
    .set {tmpDir}

    Channel.from(params.model_config)
    .set {model_config}

    Channel.from(params.remora_config)
    .set {remora_config}
    
    Channel.from(params.liftOver)
    .set {liftOver}
    
    Channel.from(params.liftOverChain)
    .set {liftOverChain}

    Channel.from(params.num_clients)
    .set {num_clients}

    Channel.from(params.port)
    .set {port}

    Channel.from(params.max_threads)
    .set {max_threads}

    Channel.from(params.pbDVMode)
    .set {pbDVMode}

    Channel.from(params.pbPATH)
    .set {pbPATH}

    Channel.from(params.annovarPath)
    .set {annovarPath}

    Channel.from(params.annovarDB)
    .set {annovarDB}

    Channel.from(params.annotsvAnnot)
    .set {annotsvAnnot}

    Channel.from(params.modkitThreads)
    .set {modkitThreads}

    Channel.from(params.cnvThreads)
    .set {cnvThreads}

    Channel.from(params.snifflesThreads)
    .set {cnvThreads}

    Channel.from(params.snp_threads)
    .set {snp_threads}

    Channel.from(params.sv_threads)
    .set {sv_threads}
    
    Channel.from(params.cov_threads)
    .set {cov_threads}
    
    Channel.from(params.meth_threads)
    .set {meth_threads}
    
    Channel.from(params.mgmt_threads)
    .set {mgmt_threads}
    
    Channel.from(params.snifflesThreads)
    .set {snifflesThreads}
    
    Channel.from(params.minimum_mgmt_cov)
    .set {minimum_mgmt_cov}

    Channel.fromPath(params.annotations, checkIfExists: true)
    .set {annotations}

    // Collect variables and scripts from bin

    Channel.fromPath("${projectDir}/data/NPHD_panel.bed", checkIfExists: true)
    .set {panel}

    Channel.fromPath("${projectDir}/data/mgmt_hg38.bed", checkIfExists: true)
    .set {mgmtBed}

    Channel.fromPath("${projectDir}/data/mgmt_probes.Rdata", checkIfExists: true)
    .set {mgmtProbes}

    Channel.fromPath("${projectDir}/data/mgmt_137sites_mean_model.Rdata", checkIfExists: true)
    .set {mgmtModel}

    Channel.fromPath("${projectDir}/scr/mgmt_pred.R", checkIfExists: true)
    .set {mgmtScript}

    Channel.fromPath("${projectDir}/scr/methylation_classification.R", checkIfExists: true)
    .set {methylationClassificationScript}
    
    Channel.fromPath("${projectDir}/scr/annotate.py", checkIfExists: true)
    .set {annotateScript}

    Channel.fromPath("${projectDir}/data/top_probes_hm450.Rdata", checkIfExists: true)
    .set {topProbes}
    
    Channel.fromPath("${projectDir}/data/capper_top_100k_betas_binarised.Rdata", checkIfExists: true)
    .set {trainingData}

    Channel.fromPath("${projectDir}/data/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true)
    .set {arrayFile}

    Channel.fromPath("${projectDir}/scr/filter_report.R", checkIfExists: true)
    .set {filterReportScript}

    Channel.fromPath("${projectDir}/scr/make_report.R", checkIfExists: true)
    .set {makereport}

    Channel.fromPath("${projectDir}/scr/Rapid_CNS2_report_UKHD_v3.0.0.Rmd", checkIfExists: true)
    .set {report_UKHD}

    // Conditional step: check if basecalling is requested
    if (params.basecalling) {
        basecalling_wf = basecalling(input=input,
                inputRef=ref, sample=id, outDir=outDir, modelConfig = model_config,
                remoraConfig = remora_config)
        inputBam = basecalling_wf.out.inputBam
        inputBai = basecalling_wf.out.inputBai
    } else {
        // If basecalling is not requested, directly use the input BAM file for variant calling
        inputBam = input
    }

    inputBai = indexBam(inputBam, max_threads)

    // Subset BAM to target regions
    subsetted_bam = subsetBam(inputBam, inputBai, panel=panel, id, max_threads)

    subset_index = indexSubsettedBam(subsetted_bam.subsetBam, max_threads)

    // Coverage calculation
    mosdepth(cov_threads, panel, inputBam, id, inputBai)

    // Call methylation
    methylation_calls = methylationCalls(inputBam, inputBai.indexedBam, ref, id,  modkitThreads)
    liftover_ch = liftOver_ch(methylation_calls.bedmethyl_file, liftOver, liftOverChain, id)

    // Methylation classification
    methylation_classification = methylationClassification(methylationClassificationScript, liftover_ch.bedmethyl_file_hg38, id, topProbes, trainingData, arrayFile, meth_threads)

    //MGMT promoter
    check_mgmt_coverage(inputBam, mgmtBed, minimum_mgmt_cov, mgmt_threads)

    mgmtPromoter_methyartist(inputBam, inputBai, ref, check_mgmt_coverage.out[0])

    mgmtPred(check_mgmt_coverage.out[0], mgmtScript, mgmtBed, mgmtProbes, mgmtModel, id)

    // CNV calling
    copyNumberVariants(inputBam, inputBai, id, cnvThreads)
    cnvAnnotated(copyNumberVariants.out, id, annotateScript)

    // SNV calling
    deepVariant_ch = deepVariant(subsetted_bam.subsetBam, subset_index.indexSubsetBam, ref, id, pbDVMode, pbPATH, tmpDir)
    recodeVCF_ch = recodeVCF(deepVariant_ch.dv_vcf)

    // ANNOVAR
    convert2annovar_ch = convert2annovar(recodeVCF_ch.pass_vcf)
    tableAnnovar_ch = tableAnnovar(convert2annovar_ch.annovar_input)
    filterReport_ch = filterReport(filterReportScript, tableAnnovar_ch.dv_anno)

    // IGV reports
    igv_reports(filterReport.out, id, ref, subsetted_bam.subsetBam, subset_index.indexSubsetBam, annotations)

    // SV calling
    sniffles2_ch = sniffles2(inputBam, inputBai, ref, id, snifflesThreads)
    annotsvAnnot(sniffles2_ch.sv_vcf, id)

    if (params.run_human_variation){
    // Human variation SNP workflow //not included in report yet
        human_variation_snp(subsetted_bam.out.subsetBam, panel, ref, id, outDir, bam_min_coverage, snp_threads)

    // Human variation SV workflow // not included in report yet
        human_variation_sv(inputBam, ref, id, sv_threads)
    }

    // Final report
    makeReport(makereport, copyNumberVariants.out, mgmtPred.out, methylation_classification.out, filterReport_ch.out, id, mosdepth.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, mgmtPromoter_methyartist.out, igv_reports.out, nextflow_version, input_bam, params.seq, report_UKHD)

    if ( params.mnpFlex) {
        mnpFlex(mnpFlexScript, liftOver_ch.out, liftOver_ch.bedmethyl_file_hg38, params.mnpFlexBed)
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
