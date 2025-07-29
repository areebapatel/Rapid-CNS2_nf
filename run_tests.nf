#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import the test modules
include { test_bam_processing; test_snv_calling; test_structural_variants; test_methylation_analysis; test_methylation_classification; test_cnv_analysis; test_coverage_analysis; test_report_generation; test_all_modules } from './test_modules.nf'

// Set up parameters for testing
params {
    // Test mode - which module to test
    test_module = 'all'  // Options: 'bam', 'snv', 'sv', 'methylation', 'classification', 'cnv', 'coverage', 'report', 'all'
    
    // Input parameters
    input_bam = null
    ref = null
    id = 'test_sample'
    threads = 4
    outDir = './test_output'
    tmpDir = './test_tmp'
    numGpus = 1
    
    // Tool paths
    annovarPath = null
    annovarDB = null
    
    // Thread parameters
    snifflesThreads = 4
    modkitThreads = 4
    cnvThreads = 4
    methThreads = 4
    
    // Coverage parameters
    minimumMgmtCov = 5
    
    // Test data parameters
    bedmethylFile = null  // For testing methylation classification separately
    
    // Report parameters
    cnv_ready = true
    mgmt_ready = true
    meth_ready = true
    filter_ready = true
    mosdepth_plot_data = null
    mgmt_cov = null
    mgmtPromoterMethyartist = null
    igv_reports = null
    nextflow_version = '23.10.0'
    seq = 'false'
    report_UKHD = null
}

// Show help message
if (params.help) {
    log.info """
    Rapid-CNS2 Module Testing Framework
    
    Usage: nextflow run run_tests.nf [options]
    
    Test Modes:
        --test_module 'bam'           - Test BAM processing module
        --test_module 'snv'           - Test SNV calling module
        --test_module 'sv'            - Test structural variant calling module
        --test_module 'methylation'   - Test methylation analysis module
        --test_module 'classification' - Test methylation classification module
        --test_module 'cnv'           - Test copy number variation module
        --test_module 'coverage'      - Test coverage analysis module
        --test_module 'report'        - Test report generation module
        --test_module 'all'           - Test all modules (default)
    
    Required Parameters:
        --input_bam                   - Path to input BAM file
        --ref                         - Path to reference genome
        --id                          - Sample identifier
    
    Optional Parameters:
        --outDir                      - Output directory [default: ./test_output]
        --threads                     - Number of threads [default: 4]
        --numGpus                     - Number of GPUs [default: 1]
        --annovarPath                 - Path to ANNOVAR installation
        --annovarDB                   - Path to ANNOVAR databases
    
    Examples:
        # Test BAM processing only
        nextflow run run_tests.nf --test_module bam --input_bam sample.bam --ref hg38.fa --id test_sample
        
        # Test SNV calling only
        nextflow run run_tests.nf --test_module snv --input_bam sample.bam --ref hg38.fa --id test_sample --annovarPath /path/to/annovar --annovarDB /path/to/humandb
        
        # Test all modules
        nextflow run run_tests.nf --test_module all --input_bam sample.bam --ref hg38.fa --id test_sample --annovarPath /path/to/annovar --annovarDB /path/to/humandb
    """
    exit 0
}

// Validate required parameters
if (params.input_bam == null) error "Input BAM file is required (--input_bam)"
if (params.ref == null) error "Reference genome is required (--ref)"

// Set up channels for required data files
Channel.fromPath("${projectDir}/data/NPHD_panel.bed", checkIfExists: true).set { panel }
Channel.fromPath("${projectDir}/data/mgmt_hg38.bed", checkIfExists: true).set { mgmtBed }
Channel.fromPath("${projectDir}/data/mgmt_probes.Rdata", checkIfExists: true).set { mgmtProbes }
Channel.fromPath("${projectDir}/data/mgmt_137sites_mean_model.Rdata", checkIfExists: true).set { mgmtModel }
Channel.fromPath("${projectDir}/scr/mgmt_pred.R", checkIfExists: true).set { mgmtScript }
Channel.fromPath("${projectDir}/scr/methylation_classification.R", checkIfExists: true).set { methylationClassificationScript }
Channel.fromPath("${projectDir}/scr/annotate.py", checkIfExists: true).set { annotateScript }
Channel.fromPath("${projectDir}/data/top_probes_hm450.Rdata", checkIfExists: true).set { topProbes }
Channel.fromPath("${projectDir}/data/capper_top_100k_betas_binarised.Rdata", checkIfExists: true).set { trainingData }
Channel.fromPath("${projectDir}/data/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true).set { arrayFile }
Channel.fromPath("${projectDir}/scr/filter_report.R", checkIfExists: true).set { filterReportScript }
Channel.fromPath("${projectDir}/scr/make_report.R", checkIfExists: true).set { makereport }
Channel.fromPath("${projectDir}/scr/Rapid_CNS2_report_UKHD_v3.0.0.Rmd", checkIfExists: true).set { reportUKHD }

// Set up input channels
Channel.fromPath(params.input_bam, checkIfExists: true).set { input_bam }
Channel.fromPath(params.ref, checkIfExists: true).set { ref }
Channel.from(params.id).set { id }
Channel.from(params.threads).set { threads }
Channel.from(params.outDir).set { outDir }
Channel.from(params.tmpDir).set { tmpDir }
Channel.from(params.numGpus).set { numGpus }
Channel.from(params.annovarPath).set { annovarPath }
Channel.from(params.annovarDB).set { annovarDB }
Channel.from(params.snifflesThreads).set { snifflesThreads }
Channel.from(params.modkitThreads).set { modkitThreads }
Channel.from(params.cnvThreads).set { cnvThreads }
Channel.from(params.methThreads).set { methThreads }
Channel.from(params.minimumMgmtCov).set { minimumMgmtCov }
Channel.from(params.cnv_ready).set { cnv_ready }
Channel.from(params.mgmt_ready).set { mgmt_ready }
Channel.from(params.meth_ready).set { meth_ready }
Channel.from(params.filter_ready).set { filter_ready }
Channel.from(params.nextflow_version).set { nextflow_version }
Channel.from(params.seq).set { seq }

// Set up optional channels
if (params.bedmethylFile) {
    Channel.fromPath(params.bedmethylFile, checkIfExists: true).set { bedmethylFile }
} else {
    Channel.value(null).set { bedmethylFile }
}

if (params.mosdepth_plot_data) {
    Channel.fromPath(params.mosdepth_plot_data, checkIfExists: true).set { mosdepth_plot_data }
} else {
    Channel.value(null).set { mosdepth_plot_data }
}

if (params.mgmt_cov) {
    Channel.from(params.mgmt_cov).set { mgmt_cov }
} else {
    Channel.value(null).set { mgmt_cov }
}

    if (params.mgmtPromoterMethyartist) {
        Channel.from(params.mgmtPromoterMethyartist).set { mgmtPromoterMethyartist }
    } else {
        Channel.value(null).set { mgmtPromoterMethyartist }
    }

if (params.igv_reports) {
    Channel.from(params.igv_reports).set { igv_reports }
} else {
    Channel.value(null).set { igv_reports }
}

if (params.report_UKHD) {
    Channel.fromPath(params.report_UKHD, checkIfExists: true).set { report_UKHD }
} else {
    Channel.fromPath("${projectDir}/scr/Rapid_CNS2_report_UKHD_v3.0.0.Rmd", checkIfExists: true).set { report_UKHD }
}

// Main workflow to run tests based on test_module parameter
workflow {
    switch (params.test_module) {
        case 'bam':
            test_bam_processing(input_bam, ref, id, threads, outDir)
            break
        case 'snv':
            // For SNV testing, we need a processed BAM file
            // This is a simplified test - in practice you'd need to run BAM processing first
            log.info "SNV testing requires a processed BAM file. Please run BAM processing first or use 'all' mode."
            break
        case 'sv':
            // For SV testing, we need a processed BAM file
            log.info "SV testing requires a processed BAM file. Please run BAM processing first or use 'all' mode."
            break
        case 'methylation':
            // For methylation testing, we need a processed BAM file
            log.info "Methylation testing requires a processed BAM file. Please run BAM processing first or use 'all' mode."
            break
        case 'classification':
            if (params.bedmethylFile) {
                test_methylation_classification(methylationClassificationScript, bedmethylFile, id, topProbes, trainingData, arrayFile, methThreads)
            } else {
                log.info "Methylation classification testing requires a bedmethyl file. Use --bedmethylFile parameter."
            }
            break
        case 'cnv':
            // For CNV testing, we need a processed BAM file
            log.info "CNV testing requires a processed BAM file. Please run BAM processing first or use 'all' mode."
            break
        case 'coverage':
            // For coverage testing, we need a processed BAM file
            log.info "Coverage testing requires a processed BAM file. Please run BAM processing first or use 'all' mode."
            break
        case 'report':
            test_report_generation(makereport, cnv_ready, mgmt_ready, meth_ready, filter_ready, id, mosdepth_plot_data, mgmt_cov, mgmtPromoterMethyartist, igv_reports, nextflow_version, input_bam, seq, report_UKHD)
            break
        case 'all':
            test_all_modules(input_bam, ref, id, threads, outDir, tmpDir, numGpus, annovarPath, annovarDB, filterReportScript, snifflesThreads, modkitThreads, mgmtBed, minimumMgmtCov, mgmtScript, mgmtProbes, mgmtModel, bedmethylFile, methylationClassificationScript, topProbes, trainingData, arrayFile, methThreads, cnvThreads, annotateScript, cnvGenes, panel, makereport, cnv_ready, mgmt_ready, meth_ready, filter_ready, mosdepth_plot_data, mgmt_cov, mgmtPromoterMethyartist, igv_reports, nextflow_version, seq, report_UKHD)
            break
        default:
            error "Unknown test module: ${params.test_module}. Use --help to see available options."
    }
} 