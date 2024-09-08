process sniffles2 {
    label 'GPU'

    publishDir "${out_dir}/sv", mode: 'copy', pattern: "*"

    input:
        path(inputBam)
        path(inputBai)
        path(ref)
        val(id)
        val(sniffles_threads)

    output:
        path "${id}.sniffles2.vcf", emit: sv_vcf
        val true

    script:
        """
        sniffles --threads ${sniffles_threads} --allow-overwrite \
                    --reference ${ref}  \
                    --non-germline \
                    --input ${inputBam} \
                    --vcf ${params.outDir}/sv/${id}.sniffles2.vcf
        """

}

process annotsvAnnot {
    publishDir "${params.outDir}/sv", mode: 'copy', pattern: "*"

    input:
        path(sv_vcf)
        val(id)

    output:
        path "${id}_sniffles_non_germline_annotsv.tsv", emit: sv_anno

    script:
        """
        ${annotsvPath}/bin/AnnotSV -SVinputFile ${sv_vcf} -outputFile ${params.outDir}/sv/${id}_sniffles_non_germline_annotsv.tsv   -svtBEDcol 4
        """

}