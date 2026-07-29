// Sniffles2 - fast baseline SV caller
process structuralVariants {
    label 'rapid_cns'

    publishDir "${params.outDir}/sv/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(ref)
        val(id)
        val(snifflesThreads)

    output:
        path "${id}.sniffles2.vcf", emit: vcf

    script:
        // Tumour-only input, so somatic (non-germline) mode is the default.
        def nonGermline = params.snifflesNonGermline ? '--non-germline' : ''
        """
        sniffles --threads ${snifflesThreads} --allow-overwrite \
            ${nonGermline} \
            --reference ${ref} \
            --input ${bam} \
            --vcf ${id}.sniffles2.vcf
        """
}

// Severus - somatic SV caller using a breakpoint graph. Better on complex
// rearrangements (chromothripsis, EGFRvIII-type events) than Sniffles alone.
// Tumour-only here: without --control-bam, a panel of normals (--PON) is the
// only somatic filter, so supply one via params.severusPON where possible.
process severus {
    label 'severus'

    publishDir "${params.outDir}/sv/severus/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        val(id)
        path(vntrBed)
        path(pon)
        val(threads)

    output:
        path "severus_${id}/**", emit: allOutputs
        path "${id}.severus.vcf", emit: vcf

    script:
        def vntrArg = vntrBed.name != 'NO_FILE' ? "--vntr-bed ${vntrBed}" : ''
        def ponArg  = pon.name    != 'NO_FILE' ? "--PON ${pon}"           : ''
        """
        severus \
            --target-bam ${bam} \
            --out-dir severus_${id} \
            -t ${threads} \
            ${vntrArg} \
            ${ponArg}

        # Tumour-only runs emit all_SVs/severus_all.vcf; with a PON a somatic
        # set is produced too. Prefer the somatic call set when present.
        if [ -f severus_${id}/somatic_SVs/severus_somatic.vcf ]; then
            cp severus_${id}/somatic_SVs/severus_somatic.vcf ${id}.severus.vcf
        elif [ -f severus_${id}/all_SVs/severus_all.vcf ]; then
            cp severus_${id}/all_SVs/severus_all.vcf ${id}.severus.vcf
        else
            FOUND=\$(find severus_${id} -name "*.vcf" | head -n 1)
            [ -n "\${FOUND}" ] || { echo "ERROR: Severus produced no VCF" >&2; exit 1; }
            cp "\${FOUND}" ${id}.severus.vcf
        fi
        """
}

process annotSV {
    label 'rapid_cns'

    publishDir "${params.outDir}/sv/", mode: 'copy'

    input:
        path(svVcf)
        val(annotSvAnnot)
        val(id)
        val(caller)

    output:
        path "${id}_${caller}_annotsv.tsv", emit: svAnno

    script:
        // AnnotSV defaults to GRCh37; this pipeline is hg38 only.
        """
        AnnotSV \
            -SVinputFile ${svVcf} \
            -outputFile ${id}_${caller}_annotsv.tsv \
            -outputDir . \
            -genomeBuild GRCh38 \
            -annotationsDir ${annotSvAnnot}
        """
}
