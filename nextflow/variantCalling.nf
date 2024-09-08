process deepVariant {
    label 'GPU'
    publishDir "${out_dir}/snv/", mode: 'copy', pattern: "*"

    //stageInMode "copy"

    input:
    path(bam)
    path(bai)
    path(ref)
    val(id)
    val(pbDVMode)
    val(pbPATH)
    val(tmpDir)

    output:
    path "${params.ouDir}/snv/${id}.dv.vcf", emit: dv_vcf

    script:
    """
    mkdir -p ${tmpDir} && \
    time ${pbPATH} deepvariant \
    --tmp-dir ${tmpDir} \
    --in-bam ${bam} \
    --ref ${ref} \
    --out-variants ${id}.dv.vcf \
    --mode ${pbDVMode} \
    --run-partition --norealign-reads
    """

}

process recodeVCF {
    publishDir "${params.outDir}/snv", mode: 'copy', pattern: "*"

    input:
        path(dv_vcf)

    output:
        path "${params.outDir}/snv/${id}.dv.PASS.vcf.gz", emit: pass_vcf

    script:
    """
    bgzip ${params.outDir}/${params.id}.dv.vcf
	vcftools --gzvcf ${params.outDir}/${params.id}.dv.vcf.gz --remove-filtered-all --recode --stdout | gzip -c > ${params.outDir}/snv/${params.id}.dv.PASS.vcf.gz
    """

}


process convert2annovar{
    publishDir "${params.outDir}/snv", mode: 'copy', pattern: "*"
    input:
        path(input)

    output:
        path "*.avinput", emit: annovar_input

    script:
        """
        ${annovarPath}/convert2annovar.pl \
        -format vcf4 ${input} \
        -withfreq \
        -includeinfo \
        > ${params.id}_dv_panel.avinput
        """
}

process tableAnnovar {
    publishDir "${params.outDir}/snv", mode: 'copy', pattern: "*"
    input:
    input:
        path(annovar_input)
    
    output:
        path "*_multianno.csv", emit: dv_anno
      
    script:
        """
        /annovar/table_annovar.pl ${annovar_input} \
        /annovar/humandb/ \
        -buildver hg19 \
        -out ${params.id}_dv_panel \
        -protocol refGene,cytoBand,avsnp147,dbnsfp30a,1000g2015aug_eur,cosmic70 \
        -operation gx,r,f,f,f,f \
        -nastring . \
        -csvout \
        -polish \
        -otherinfo
        """
}

process filterReport {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    publishDir "${params.outDir}/snv", mode: 'copy', pattern: "*"
    input:

    input:
        path(filterReportScript)
        path(dv_anno)

    output:
        path "${params.outDir}/snv/${params.id}_dv_report.csv", emit: dv_report
    
    script:
        """
        Rscript ${filterReportScript} \
        --input ${dv_anno} \
        --output ${params.id}_dv_report.csv \
        --sample ${params.id}
        """        
}


process human_variation_sv {
    errorStrategy 'ignore'
    input:
        path(inputBam)
        path(ref)
        val(id)
        val(sv_threads)

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
        --bam ${inputBam} \
        --sample_name ${params.id} \
        --bam_min_coverage ${params.bam_min_coverage} \
        --out_dir ${params.outDir} \
        --threads ${sv_threads} \
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
        val(bam_min_coverage)
        val(snp_threads)

    output:
        val true

    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
        -with-report ${params.outDir}/human_variation_snp_nextflow_report.html \
        -profile standard \
        -w ${PWD}/${params.outDir}/wf-human-variation/workspace \
        --ref ${ref} \
        --snp \
        --bam ${bam} \
        --bed ${panel} \
        --sample_name ${id} \
        --bam_min_coverage ${bam_min_coverage} \
        --out_dir ${PWD}/${params.outDir}/wf-human-variation/snp/ \
        --threads ${snp_threads}
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

    publishDir("${params.outDir}/snv")

    output:
        val true
 
    script:
        """
        sed -e 's/,/\t/g' -e 's/\"//g' \
        ${PWD}/${params.outDir}/snv/${params.id}_dv_report.csv > ${PWD}/${params.outdir}/snv/${params.id}_dv_report.fmt.csv 
        create_report ${PWD}/${params.outDir}/${id}_clair3_report.fmt.csv \
        --fasta ${ref} \
        --sequence 1 \
        --begin 2 \
        --end 3 \
        --flanking 1000 \
        --info-columns Chr Start End Func Gene ExonicFunc AAChange cytoBand 1000g_EUR COSMIC \
        --output ${PWD}/${params.outDir}/${id}_igv-report.html \
        --standalone \
        --tracks ${bam} ${annotations}
        """
}

