process mosdepth {
    label 'rapid_cns'

    publishDir "${params.outDir}/coverage/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(panel)
        val(id)

    output:
        path "${id}.mosdepth.summary.txt", emit: summary
        path "${id}.regions.bed.gz",       optional: true

    script:
        """
        mosdepth \
            -t ${task.cpus} \
            -n \
            --by ${panel} \
            --fast-mode \
            ${id} \
            ${bam}
        """
}
