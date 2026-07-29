// Input BAM handling: tag check, alignment (if needed), merge, sort, index.
// Collapsed into a single process so the align/merge decision is made at run
// time from the actual data rather than from a Groovy flag that is evaluated
// before the check has run.
process prepareBam {
    label 'bamprep'

    // The prepared BAM is a full flow cell (>150 GB) and is reproducible from
    // the input, so it is NOT published by default - a copy would double the
    // footprint. Enable with --publishBam when you want it on disk.
    publishDir enabled: params.publishBam,
               path:    "${params.outDir}/bam/",
               pattern: "*.bam",
               mode:    'copy'
    publishDir enabled: params.publishBam,
               path:    "${params.outDir}/bam/",
               pattern: "*.bam.bai",
               mode:    'copy'

    input:
        path(bams, stageAs: 'input/*')
        path(ref)
        val(id)
        val(threads)

    output:
        tuple path("${id}.bam"), path("${id}.bam.bai"), emit: bam

    // Use task.cpus/task.memory rather than params.maxThreads: the scheduler
    // allocation is what actually constrains this step. samtools sort defaults
    // to 768M per thread, so -@64 would try to reserve ~49 GB and be OOM-killed.
    script:
        def sortMem = Math.max(1, (int)(task.memory.toGiga() * 0.6 / task.cpus))
        """
        BAMS=(input/*.bam)
        echo "Found \${#BAMS[@]} input BAM file(s)"
        echo "Allocation: ${task.cpus} cpus, ${task.memory} (sort: -@${task.cpus} -m ${sortMem}G)"

        # One pass over the head of the first BAM answers both questions:
        # are modified-base tags present, and is the data already aligned?
        # A full 'samtools view -c -F 4' would scan the entire file.
        samtools view -@${task.cpus} "\${BAMS[0]}" 2>view.err | head -n 100000 > head.sam || true

        if [ ! -s head.sam ]; then
            echo "ERROR: could not read any records from \${BAMS[0]}." >&2
            echo "samtools said:" >&2; cat view.err >&2
            exit 1
        fi

        if ! grep -q 'MM:Z' head.sam; then
            echo "ERROR: no methylation tags (MM:Z) found in \${BAMS[0]}." >&2
            echo "Re-run basecalling with modified bases enabled:" >&2
            echo "  https://github.com/nanoporetech/dorado?tab=readme-ov-file#modified-basecalling" >&2
            exit 1
        fi

        # column 3 is RNAME; '*' means unmapped
        MAPPED=\$(awk '\$3 != "*"' head.sam | wc -l)
        rm -f head.sam view.err
        echo "Sampled 100000 reads, \${MAPPED} mapped"

        if [ "\${MAPPED}" -le 2 ]; then
            echo "Input is unaligned - running dorado aligner"
            mkdir -p aligned
            dorado aligner ${ref} input/ --output-dir aligned --threads ${task.cpus}
            OUT=(aligned/*.bam)
        else
            echo "Input is already aligned - skipping alignment"
            OUT=("\${BAMS[@]}")
        fi

        if [ "\${#OUT[@]}" -gt 1 ]; then
            samtools cat -o cat.bam "\${OUT[@]}"
            samtools sort -@${task.cpus} -m ${sortMem}G -o ${id}.bam cat.bam
            rm -f cat.bam
        elif samtools view -H "\${OUT[0]}" | head -n 1 | grep -q 'SO:coordinate'; then
            echo "Already coordinate-sorted - no re-sort needed"
            cp "\${OUT[0]}" ${id}.bam
        else
            samtools sort -@${task.cpus} -m ${sortMem}G -o ${id}.bam "\${OUT[0]}"
        fi

        samtools index -@${task.cpus} ${id}.bam
        """
}

process subsetBam {
    label 'rapid_cns'

    // Panel subset is far smaller than the full BAM and is useful for review
    // (it backs the IGV report), so this one is published by default.
    publishDir path:    "${params.outDir}/bam/",
               pattern: "*.bam",
               mode:    'copy'
    publishDir path:    "${params.outDir}/bam/",
               pattern: "*.bam.bai",
               mode:    'copy'

    input:
        tuple path(bam), path(bai)
        path(panel)
        val(id)
        val(threads)

    output:
        tuple path("${id}.RapidCNS2.subset.bam"), path("${id}.RapidCNS2.subset.bam.bai"), emit: bam

    script:
        """
        samtools view -@${threads} -b -M -L ${panel} -o ${id}.RapidCNS2.subset.bam ${bam}
        samtools index -@${threads} ${id}.RapidCNS2.subset.bam
        """
}
