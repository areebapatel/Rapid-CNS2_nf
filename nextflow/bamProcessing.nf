
process checkAlignment {
    label 'rapid_cns'

    input:
        path(inputBam)
        val(threads)
    
    output:
        stdout emit: alignment_check

    script:
        """
        # Count aligned reads (not unmapped reads)
        aligned_count=\$(samtools view -@${threads} -F 4 ${inputBam} | wc -l)
        echo "\${aligned_count}"
        """
}

process checkMethylationTags {
    label 'rapid_cns'
    
    input:
        path(inputBam)
        val(threads)
    
    output:
        stdout emit: meth_check

    script:
        """
        samtools \
        view \
        -@${threads} \
        ${inputBam} \
        | grep -m 1 MM:Z
        """
}

process alignBam {
    label 'rapid_cns'
    
    input:
        path(input)
        path(ref)
        val(threads)
        path(outDir)

    publishDir "${params.outDir}/bam/alignedBams/", mode: 'copy'
    
    output:
        path "alignedBams/*.bam", emit : alignedBam

    script:
        """
        
        # Check if file is already aligned
        aligned_count=\$(samtools view -F 4 "\$input" | wc -l)
        
        if [ "\$aligned_count" -gt 2 ]; then
            echo "File already aligned with \$aligned_count reads. Copying to alignedBams directory."
            cp "\$input" alignedBams/
        else
            echo "File has \$aligned_count aligned reads. Performing alignment."
            dorado aligner "\$ref" "\$input" --output-dir alignedBams/ --threads "\$threads"
        fi
        """
}

process mergeBam {
    label 'rapid_cns'
    
    input:
        path(bams)
        val(threads)
        path(outDir)
        val(id)
    
    publishDir "${params.outDir}/bam/", mode: 'copy'

    output:
        path "*.bam", emit : mergedBam

    script:
        """
        samtools merge -@${threads} -o ${outDir}/bam/${id}.merged.bam ${bams.join(' ')}
        """
}

process addreplacerg {
    label 'rapid_cns'
    
    input:
        path(inputBam)
        path(inputBai)
    
    output:
        path "*.bam", emit : deepVariantBam
        path "*.bam.bai", emit : deepVariantBai
    
    publishDir("${params.outDir}/bam")

    script:
        def r_args = params.reads ?: ''
        """
         samtools addreplacerg ${r_args} \
                -@${threads} -o ${id}.pass.index.bam \
                ${inputBam}

        samtools index -@ ${threads} ${params.outDir}/bam/${id}.pass.index.bam 
        """
}

process indexBam {
    label 'rapid_cns'
    
    input:
        path(bam)
        val(threads)
    
    output:
        path "*.bai", emit: indexBam

    script:
        """
        samtools index -@${threads} ${bam}
        """
}   

process subsetBam {
    label 'rapid_cns'
    
    input:
        path(bam)
        path(indexBam)
        path(panel)
        val(id)
        val(threads)

    publishDir "${params.outDir}/bam/", mode: 'copy'

    output:
        path "*subset.bam", emit: subsetBam

    script:
        """
        samtools index -@${threads} ${bam}
        
        bedtools intersect -a ${panel} \
        -b ${bam} \
        > ${params.outDir}/bam/${id}.RapidCNS2.subset.bam
        """
}

process indexSubsettedBam{
    label 'rapid_cns'
    
    input:
        path(bam)
        val(threads)
    
    output:
        path "*.bai", emit: indexSubsetBam

    publishDir("${params.outDir}/bam/")
    script:
        """
        samtools index -@${threads} ${bam}
        """
}