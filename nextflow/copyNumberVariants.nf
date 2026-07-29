process copyNumberVariants {
    label 'rapid_cns'

    publishDir "${params.outDir}/cnv/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        val(id)
        val(cnvThreads)

    output:
        path "${id}.cnvpytor.calls.1000.tsv",   emit: calls1000
        path "${id}.cnvpytor.calls.10000.tsv",  emit: calls10000
        path "${id}.cnvpytor.calls.100000.tsv", emit: calls100000
        path "${id}_cnvpytor_100k.pdf",         emit: plotPdf
        path "${id}_cnvpytor_100k.png",         emit: plotPng

    script:
        def chroms = (1..22).collect { n -> "chr${n}" }.join(' ') + ' chrX chrY'
        """
        cnvpytor -root ${id}_CNV.pytor -rd ${bam} -j ${cnvThreads}
        cnvpytor -root ${id}_CNV.pytor -his 1000 10000 100000 -j ${cnvThreads}
        cnvpytor -root ${id}_CNV.pytor -partition 1000 10000 100000 -j ${cnvThreads}
        cnvpytor -root ${id}_CNV.pytor -call 1000   -j ${cnvThreads} > ${id}.cnvpytor.calls.1000.tsv
        cnvpytor -root ${id}_CNV.pytor -call 10000  -j ${cnvThreads} > ${id}.cnvpytor.calls.10000.tsv
        cnvpytor -root ${id}_CNV.pytor -call 100000 -j ${cnvThreads} > ${id}.cnvpytor.calls.100000.tsv
        cnvpytor -root ${id}_CNV.pytor -plot manhattan 100000 -chrom ${chroms} -o ${id}_cnvpytor_100k.pdf
        cnvpytor -root ${id}_CNV.pytor -plot manhattan 100000 -chrom ${chroms} -o ${id}_cnvpytor_100k.png
        """
}

process cnvAnnotated {
    label 'rapid_cns'

    publishDir "${params.outDir}/cnv/", mode: 'copy'

    input:
        path(cnvCalls)
        val(id)
        path(annotateScript)
        path(cnvGenes)

    output:
        path "${id}.annotation.1000.xlsx", emit: cnvAnnotated

    script:
        """
        python3 ${annotateScript} ${cnvCalls} ${cnvGenes} ${id}.annotation.1000.xlsx
        """
}

// SAVANA tumour-only SV calling. Its breakpoints feed `savana cna` below, which
// gives better copy-number segmentation than read depth alone.
process savanaTo {
    label 'savana'

    publishDir "${params.outDir}/sv/savana/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(ref)
        path(refIdx)
        val(id)
        val(threads)

    output:
        path "savana_to_${id}/**", emit: allOutputs
        path "${id}.savana.breakpoints.vcf", emit: breakpoints

    script:
        """
        savana to \
            --tumour ${bam} \
            --ref ${ref} \
            --ref_index ${refIdx} \
            --sample ${id} \
            --outdir savana_to_${id} \
            --threads ${threads}

        VCF=\$(find savana_to_${id} -name "*.vcf" | grep -vi "germline" | head -n 1)
        [ -n "\${VCF}" ] || { echo "ERROR: savana to produced no VCF" >&2; exit 1; }
        cp "\${VCF}" ${id}.savana.breakpoints.vcf
        """
}

// SAVANA copy number: absolute CN fitting with tumour purity and ploidy
// estimated from B-allele frequencies at 1000G het SNP sites.
//
// NOTE: CN segmentation and purity/ploidy fitting assume reasonably uniform
// genome-wide coverage. Adaptive-sampling panel data is strongly enriched
// on-target with sparse off-target background, so these estimates should be
// validated against known samples before being trusted clinically.
process savanaCna {
    label 'savana'

    publishDir "${params.outDir}/cnv/savana/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(ref)
        path(breakpoints)
        val(id)
        val(g1000)
        val(threads)

    output:
        path "savana_cna_${id}/**", emit: allOutputs
        path "${id}_purity_ploidy.tsv", optional: true, emit: purityPloidy

    script:
        """
        savana cna \
            --tumour ${bam} \
            --ref ${ref} \
            --sample ${id} \
            --outdir savana_cna_${id} \
            --breakpoints ${breakpoints} \
            --g1000_vcf ${g1000} \
            --cna_threads ${threads} \
            --tmpdir .

        # surface the purity/ploidy fit at a predictable path for the report
        FIT=\$(find savana_cna_${id} -iname "*purity_ploidy*" | head -n 1)
        if [ -n "\${FIT}" ]; then
            cp "\${FIT}" ${id}_purity_ploidy.tsv
        else
            echo "WARNING: no purity/ploidy fit produced by savana cna" >&2
        fi
        """
}
