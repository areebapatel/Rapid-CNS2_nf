process mosdepth {
    label 'rapid_cns'

    publishDir "${params.outDir}/coverage/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(panel)
        val(id)
        val(threads)

    output:
        path "${id}.mosdepth.summary.txt", emit: summary
        path "${id}.regions.bed.gz",       optional: true

    script:
        """
        mosdepth \
            -t ${threads} \
            -n \
            --by ${panel} \
            --fast-mode \
            ${id} \
            ${bam}
        """
}
