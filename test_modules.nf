#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Test module for BAM Processing
workflow test_bam_processing {
    take:
        input_bam
        ref
        id
        threads
        outDir
    
    main:
        // Test BAM alignment
        alignBam_out = alignBam(input_bam, ref, threads, outDir)
        
        // Test BAM merging (if multiple files)
        mergeBam_out = mergeBam(alignBam_out.alignedBam, threads, outDir, id)
        
        // Test BAM indexing
        indexBam_out = indexBam(mergeBam_out.mergedBam, threads)
        
        // Test subsetting
        subsetBam_out = subsetBam(mergeBam_out.mergedBam, indexBam_out.indexedBam, panel, id, threads)
        
        // Test subset indexing
        indexSubset_out = indexSubsettedBam(subsetBam_out.subsetBam, threads)
    
    emit:
        aligned_bams = alignBam_out.alignedBam
        merged_bam = mergeBam_out.mergedBam
        indexed_bam = indexBam_out.indexedBam
        subset_bam = subsetBam_out.subsetBam
        subset_index = indexSubset_out.indexSubsetBam
}

// Test module for SNV Calling
workflow test_snv_calling {
    take:
        bam
        bai
        ref
        id
        tmpDir
        numGpus
        annovarPath
        annovarDB
        filterReportScript
    
    main:
        // Test Variant Calling
        variantCalling_out = variantCalling(bam, bai, ref, id, tmpDir, numGpus)
        
        // Test VCF recoding
        recodeVCF_out = recodeVCF(variantCalling_out.dvVcf)
        
        // Test ANNOVAR conversion
        convert2annovar_out = convert2annovar(recodeVCF_out.passVcf, annovarPath)
        
        // Test ANNOVAR annotation
        tableAnnovar_out = tableAnnovar(convert2annovar_out.annovarInput, annovarPath, annovarDB)
        
        // Test filtering and reporting
        filterReport_out = filterReport(filterReportScript, tableAnnovar_out.dvAnno, id)
    
    emit:
        vcf = variantCalling_out.dvVcf
        pass_vcf = recodeVCF_out.passVcf
        annotated_variants = tableAnnovar_out.dvAnno
        filtered_report = filterReport_out.dvReport
}

// Test module for Structural Variant Calling
workflow test_structural_variants {
    take:
        bam
        bai
        ref
        id
        snifflesThreads
    
    main:
        // Test Structural Variants
        structuralVariants_out = structuralVariants(bam, bai, ref, id, snifflesThreads)
        
        // Test AnnotSV annotation
        annotsv_out = annotSV(structuralVariants_out.svVcf, id)
    
    emit:
        sv_vcf = structuralVariants_out.svVcf
        sv_annotated = annotsv_out.svAnno
}

// Test module for Methylation Analysis
workflow test_methylation_analysis {
    take:
        bam
        bai
        ref
        id
        modkitThreads
        mgmtBed
        minimumMgmtCov
        threads
        mgmtScript
        mgmtProbes
        mgmtModel
        bedmethylFile
    
    main:
        // Test methylation calling
        methylation_out = methylationCalls(bam, bai, ref, id, modkitThreads)
        
        // Test MGMT coverage
        mgmtCoverage_out = checkMgmtCoverage(bam, bai, mgmtBed, minimumMgmtCov, threads)
        
        // Test MGMT promoter analysis
        mgmtPromoter_out = mgmtPromoterMethyartist(bam, bai, ref, mgmtCoverage_out.mgmt_avg_cov, id)
        
        // Test MGMT prediction
        mgmtPred_out = mgmtPred(methylation_out.bedmethylFile, mgmtScript, mgmtBed, mgmtProbes, mgmtModel, bedmethylFile, id)
    
    emit:
        bedmethyl = methylation_out.bedmethylFile
        mgmt_coverage = mgmtCoverage_out.mgmt_avg_cov
        mgmt_plot = mgmtPromoter_out.mgmt_plot
}

// Test module for Methylation Classification
workflow test_methylation_classification {
    take:
        methylationClassificationScript
        bedmethylFile
        id
        topProbes
        trainingData
        arrayFile
        methThreads
    
    main:
        // Test methylation classification
        methylationClass_out = methylationClassification(methylationClassificationScript, bedmethylFile, id, topProbes, trainingData, arrayFile, methThreads)
    
    emit:
        classification_results = methylationClass_out
}

// Test module for Copy Number Variation
workflow test_cnv_analysis {
    take:
        bam
        bai
        id
        cnvThreads
        annotateScript
        cnvGenes
        outDir
    
    main:
        // Test CNV calling
        cnv_out = copyNumberVariants(bam, bai, id, cnvThreads)
        
        // Test CNV annotation
        cnvAnnotated_out = cnvAnnotated(cnv_out, id, annotateScript, cnvGenes, outDir)
    
    emit:
        cnv_results = cnv_out
        cnv_annotated = cnvAnnotated_out
}

// Test module for Coverage Analysis
workflow test_coverage_analysis {
    take:
        threads
        panel
        bam
        bai
        id
    
    main:
        // Test mosdepth coverage
        mosdepth_out = mosdepth(threads, panel, bam, bai, id)
    
    emit:
        coverage_summary = mosdepth_out.mosdepth_out
}

// Test module for Report Generation
workflow test_report_generation {
    take:
        reportScript
        cnv_ready
        mgmt_ready
        meth_ready
        filter_ready
        id
        mosdepth_plot_data
        mgmt_cov
        mgmtPromoterMethyartist
        igv_reports
        nextflow_version
        inputBam
        seq
        report_UKHD
    
    main:
        // Test report rendering
        report_out = reportRendering(reportScript, cnv_ready, mgmt_ready, meth_ready, filter_ready, id, mosdepth_plot_data, mgmt_cov, mgmtPromoterMethyartist, igv_reports, nextflow_version, inputBam, seq, report_UKHD)
    
    emit:
        report_generated = report_out
}

// Main test workflow that can run all modules or specific ones
workflow test_all_modules {
    take:
        input_bam
        ref
        id
        threads
        outDir
        tmpDir
        numGpus
        annovarPath
        annovarDB
        filterReportScript
        snifflesThreads
        modkitThreads
        mgmtBed
        minimumMgmtCov
        mgmtScript
        mgmtProbes
        mgmtModel
        bedmethylFile
        methylationClassificationScript
        topProbes
        trainingData
        arrayFile
        methThreads
        cnvThreads
        annotateScript
        cnvGenes
        panel
        reportScript
        cnv_ready
        mgmt_ready
        meth_ready
        filter_ready
        mosdepth_plot_data
        mgmt_cov
        mgmtPromoterMethyartist
        igv_reports
        nextflow_version
        seq
        report_UKHD
    
    main:
        // Test BAM processing
        bam_processing_out = test_bam_processing(input_bam, ref, id, threads, outDir)
        
        // Test SNV calling
        snv_out = test_snv_calling(bam_processing_out.merged_bam, bam_processing_out.indexed_bam, ref, id, tmpDir, numGpus, annovarPath, annovarDB, filterReportScript)
        
        // Test structural variants
        sv_out = test_structural_variants(bam_processing_out.merged_bam, bam_processing_out.indexed_bam, ref, id, snifflesThreads)
        
        // Test methylation analysis
        methylation_out = test_methylation_analysis(bam_processing_out.merged_bam, bam_processing_out.indexed_bam, ref, id, modkitThreads, mgmtBed, minimumMgmtCov, threads, mgmtScript, mgmtProbes, mgmtModel, bedmethylFile)
        
        // Test methylation classification
        meth_class_out = test_methylation_classification(methylationClassificationScript, methylation_out.bedmethyl, id, topProbes, trainingData, arrayFile, methThreads)
        
        // Test CNV analysis
        cnv_out = test_cnv_analysis(bam_processing_out.merged_bam, bam_processing_out.indexed_bam, id, cnvThreads, annotateScript, cnvGenes, outDir)
        
        // Test coverage analysis
        coverage_out = test_coverage_analysis(threads, panel, bam_processing_out.merged_bam, bam_processing_out.indexed_bam, id)
        
        // Test report generation
        report_out = test_report_generation(reportScript, cnv_out.cnv_results, mgmt_ready, meth_class_out.classification_results, snv_out.filtered_report, id, coverage_out.coverage_summary, mgmt_cov, mgmtPromoterMethyartist, igv_reports, nextflow_version, input_bam, seq, report_UKHD)
    
    emit:
        bam_results = bam_processing_out
        snv_results = snv_out
        sv_results = sv_out
        methylation_results = methylation_out
        classification_results = meth_class_out
        cnv_results = cnv_out
        coverage_results = coverage_out
        report_results = report_out
}

// Include all the process definitions
include { alignBam; mergeBam; indexBam; subsetBam; indexSubsettedBam } from './nextflow/bamProcessing.nf'
include { variantCalling; recodeVCF; convert2annovar; tableAnnovar; filterReport } from './nextflow/variantCalling.nf'
include { structuralVariants; annotSV } from './nextflow/structuralVariants.nf'
include { methylationCalls; checkMgmtCoverage; mgmtPromoterMethyartist; mgmtPred } from './nextflow/methylationAnalysis.nf'
include { methylationClassification } from './nextflow/methylationClassification.nf'
include { copyNumberVariants; cnvAnnotated } from './nextflow/copyNumberVariants.nf'
include { mosdepth } from './nextflow/utils.nf'
include { reportRendering } from './nextflow/reportRendering.nf'

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