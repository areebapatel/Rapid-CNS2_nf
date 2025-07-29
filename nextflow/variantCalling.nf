process variantCalling {
    label 'gpu'

    publishDir "${outDir}/snv/", mode: 'copy', pattern: "*"

    //stageInMode "copy"

    input:
    path(bam)
    path(bai)
    path(ref)
    val(id)
    val(tmpDir)
    val(numGpus)

    output:
    path "*.deepvariant.vcf", emit: dvVcf

    script:
    """
    pbrun deepvariant \
    --tmp-dir ${tmpDir} \
    --in-bam ${bam} \
    --ref ${ref} \
    --out-variants ${id}.deepvariant.vcf \
    --mode ont \
    --run-partition --norealign-reads \
    --num-gpus ${numGpus}
    """

}

process recodeVCF {
    label 'rapid_cns'
    
    publishDir "${params.outDir}/snv", mode: 'copy', pattern: "*"

    input:
        path(dvVcf)

    output:
        path "${id}.deepvariant.PASS.vcf.gz", emit: passVcf

    script:
    """
    bgzip ${id}.deepvariant.vcf
	vcftools --gzvcf ${id}.deepvariant.vcf.gz --remove-filtered-all --recode --stdout | gzip -c > ${id}.deepvariant.PASS.vcf.gz
    """

}


process convert2annovar{
    label 'rapid_cns'
    
    publishDir "${params.outDir}/snv", mode: 'copy', pattern: "*"
    input:
        path(inputVcf)
        path(annovarPath)

    output:
        path "*.avinput", emit: annovarInput

    script:
        """
        ${annovarPath}/convert2annovar.pl \
        -format vcf4 ${inputVcf} \
        -withfreq \
        -includeinfo \
        > ${id}_dv_panel.avinput
        """
}

process tableAnnovar {
    label 'rapid_cns'
    
    publishDir "${params.outDir}/snv", mode: 'copy', pattern: "*"
    
    input:
        path(annovarInput)
        path(annovarPath)
        path(annovarDB)
    
    output:
        path "*_multianno.csv", emit: dvAnno
      
    script:
        """
        ${annovarPath}/table_annovar.pl ${annovar_input} \
        ${annovarDB} \
        -buildver hg38 \
        -out ${params.id}_dv_panel \
        -protocol refGene,cytoBand,clinvar_20240917,avsnp151,1000g2015aug_eur,cosmic70,dbnsfp42c,allofus \
        -operation gx,r,f,f,f,f,f \
        -nastring . \
        -csvout \
        -polish \
        -otherinfo
        """
}

process filterReport {
    label 'rapid_cns'
    
    publishDir "${params.outDir}/snv", mode: 'copy', pattern: "*"

    input:
        path(filterReportScript)
        path(dvAnno)
        val(id)

    output:
        val true
        path "${id}_dv_report.csv", emit: dvReport
    
    script:
        """
        Rscript ${filterReportScript} \
        --input ${dvAnno} \
        --output ${id}_dv_report.csv \
        --sample ${id}
        """        
}


process human_variation_sv {
    errorStrategy 'ignore'
    input:
        path(bam)
        path(ref)
        val(id)
        val(svThreads)

    publishDir("${params.outDir}/wf-human-variation/sv/")

    output:
        val true   

    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
        -with-report ${params.outDir}/human_variation_sv_cnv_nextflow_report.html \
        -profile standard \
        -w ${params.outDir}/wf-human-variation/sv/workspace \
        --ref ${ref} \
        --sv \
        --bam ${bam} \
        --sample_name ${id} \
        --bam_min_coverage ${bamMinCoverage} \
        --out_dir ${params.outDir}/wf-human-variation/sv/ \
        --threads ${svThreads} \
        --sniffles_args="--non-germline"
        """
}


process human_variation_snp {
    errorStrategy 'ignore'
    input:
        path(bam)
        path(panel)
        path(ref)
        val(id)
        val(outdir)
        val(bamMinCoverage)
        val(snpThreads)

    output:
        val true

    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
        -with-report ${params.outDir}/human_variation_snp_nextflow_report.html \
        -profile standard \
        -w ${params.outDir}/wf-human-variation/workspace \
        --ref ${ref} \
        --snp \
        --bam ${bam} \
        --bed ${panel} \
        --sample_name ${id} \
        --bam_min_coverage ${bamMinCoverage} \
        --out_dir ${params.outDir}/wf-human-variation/snp/ \
        --threads ${snpThreads}
        """
 }

 process igv_reports {
    errorStrategy 'ignore'
    input:
        val(ready)  // filter-report ready
        val(id)
        path(reference)
        path(bam)
        path(indexedBam)
        path(annotations)

    publishDir("${params.outDir}/reports")

    output:
        val true
 
    script:
        """
        sed -e 's/,/\t/g' -e 's/\"//g' \
        ${params.outDir}/snv/${id}_dv_report.csv > ${params.outDir}/snv/${id}_dv_report.fmt.csv 
        
        create_report ${params.outDir}/snv/${id}_dv_report.fmt.csv \
        --fasta ${ref} \
        --sequence 1 \
        --begin 2 \
        --end 3 \
        --flanking 1000 \
        --info-columns Chr Start End Func Gene ExonicFunc AAChange cytoBand 1000g_EUR COSMIC \
        --output ${params.outDir}/reports/${id}_igv-report.html \
        --standalone \
        --tracks ${bam} ${annotations}
        """
}

