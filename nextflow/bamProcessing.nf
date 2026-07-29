// Input BAM handling: tag check, alignment (if needed), merge, sort, index.
// Collapsed into a single process so the align/merge decision is made at run
// time from the actual data rather than from a Groovy flag that is evaluated
// before the check has run.
process prepareBam {
    label 'rapid_cns'

    publishDir "${params.outDir}/bam/", mode: 'copy'

    input:
        path(bams, stageAs: 'input/*')
        path(ref)
        val(id)
        val(threads)

    output:
        tuple path("${id}.bam"), path("${id}.bam.bai"), emit: bam

    script:
        """
        BAMS=(input/*.bam)
        echo "Found \${#BAMS[@]} input BAM file(s)"

        # Modified-base tags must be present; fail early with actionable guidance.
        samtools view -@${threads} "\${BAMS[0]}" 2>/dev/null | head -n 100000 > tagcheck.sam || true
        if ! grep -q 'MM:Z' tagcheck.sam; then
            echo "ERROR: no methylation tags (MM:Z) found in \${BAMS[0]}." >&2
            echo "Re-run basecalling with modified bases enabled:" >&2
            echo "  https://github.com/nanoporetech/dorado?tab=readme-ov-file#modified-basecalling" >&2
            exit 1
        fi
        rm -f tagcheck.sam

        ALIGNED=\$(samtools view -@${threads} -c -F 4 "\${BAMS[0]}")
        echo "First input BAM has \${ALIGNED} aligned reads"

        if [ "\${ALIGNED}" -le 2 ]; then
            echo "Input is unaligned - running dorado aligner"
            mkdir -p aligned
            dorado aligner ${ref} input/ --output-dir aligned --threads ${threads}
            OUT=(aligned/*.bam)
        else
            echo "Input is already aligned - skipping alignment"
            OUT=("\${BAMS[@]}")
        fi

        if [ "\${#OUT[@]}" -gt 1 ]; then
            samtools cat -o cat.bam "\${OUT[@]}"
            samtools sort -@${threads} -o ${id}.bam cat.bam
            rm -f cat.bam
        elif samtools view -H "\${OUT[0]}" | head -n 1 | grep -q 'SO:coordinate'; then
            cp "\${OUT[0]}" ${id}.bam
        else
            samtools sort -@${threads} -o ${id}.bam "\${OUT[0]}"
        fi

        samtools index -@${threads} ${id}.bam
        """
}

process subsetBam {
    label 'rapid_cns'

    publishDir "${params.outDir}/bam/", mode: 'copy'

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
