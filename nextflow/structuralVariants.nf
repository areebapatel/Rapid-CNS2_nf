process structuralVariants {
    label 'gpu'

    publishDir "${outDir}/sv", mode: 'copy', pattern: "*"

    input:
        path(inputBam)
        path(inputBai)
        path(ref)
        val(id)
        val(snifflesThreads)

    output:
        path "${id}.sniffles2.vcf", emit: svVcf
        val true

    script:
        """
        sniffles --threads ${snifflesThreads} --allow-overwrite \
                    --reference ${ref}  \
                    --non-germline \
                    --input ${inputBam} \
                    --vcf ${id}.sniffles2.vcf
        """

}

process annotSV {
    publishDir "${params.outDir}/sv", mode: 'copy', pattern: "*"

    input:
        path(svVcf)
        path(annotSvAnnot)
        val(id)

    output:
        path "${id}_sniffles_non_germline_annotsv.tsv", emit: svAnno

    script:
        """
        AnnotSV -SVinputFile ${svVcf} -outputFile ${id}_sniffles_non_germline_annotsv.tsv -annotationsDir ${annotSvAnnot}   -svtBEDcol 4
        """

}
