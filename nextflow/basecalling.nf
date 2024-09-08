process basecalling {
    input:
        path(input)
        path(inputRef)
        val(id)
        path(outDir)
        val(modelConfig)
        val(remoraConfig)

    output:
        path "*.pass.bam" , emit : inputBam
        path "*.pass.bam.bai" , emit : inputBai
    
    script:
        """
        nextflow run epi2me-labs/wf-basecalling \
        -profile singularity \
        --input ${input} \
        --ref ${ref} \    
        --dorado_ext fast5 \
        --sample_name ${id} \
        --out_dir ${outDir}/basecalling/ \
        --basecaller_cfg ${modelConfig} \
        --remora_cfg ${remoraConfig}
        """ 
}


process addreplacerg {
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

        samtools index -@ ${threads} ${params.outDir}/basecalling/${id}.pass.index.bam 
        """
}

process subsetBam {
    input:
        path(inputBam)
        path(inputBai)
        path(panel)
        val(id)
        val(threads)

    publishDir("${params.outDir}/basecalling/")

    output:
        path "*subset.bam", emit: subsetBam

    script:
        """
        bedtools intersect -a ${panel} \
        -b ${inputBam} \
        > ${id}.subset.bam
        """
}

process indexBam{

    input:
        path(bam)
        val(max_threads)
    
    output:
        path "*bam.bai", emit: indexedBam

    publishDir("${params.outDir}/bam/")

    script:
        """
        mkdir -p ${params.outDir}/bam/
        cp ${bam} ${params.outDir}/bam/
        samtools index -@${max_threads} ${params.outDir}/bam/${bam}
        """
}

process indexSubsettedBam{

    input:
        path(bam)
        val(max_threads)
    
    output:
        path "*.bai", emit: indexSubsetBam

    publishDir("${params.outDir}/bam/")
    script:
        """
        samtools index -@${max_threads} ${bam}
        """
}