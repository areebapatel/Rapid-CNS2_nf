process gzip {
    input:
        path(input_file)

    output:
        path "*.gz", emit: compressed_out

    script:
        """
        pigz \
        -1 \
        -c \
        ${input_file} \
        > ${input_file}.gz
        """
}

process vcf_intersect {
    input:
        path(input1)
        path(input2)
        val(output_file)

    output:   
        path "*.vcf", emit: intersect_vcf
        
    script:
        """
        bedtools \
        intersect \
        -a ${input1} \
        -b ${input2} > ${output_file}.vcf
        """
}

process bedmethyl_intersect {
    input:
        val(input1) // filtered bedmethyl file from filter to just 5mC
        path(input2)
        val(out_dir)
        val(output_file)
        val(id)

    publishDir("${out_dir}")

    output:   
        path "*.bed", emit: intersect_bed
        
    script:
        """
        bedtools \
        intersect \
        -a ${input1} \
        -b ${input2} > ${output_file}.bed
        """
}

process mosdepth {
    input:
        val(threads)
        path(panel)
        path(inputBam)
        val(id)
        path(inputBai)
    
    publishDir("${params.outDir}/coverage")
	
    output:
        path "*.mosdepth.summary.txt", emit: mosdepth_out

    script:
        """
        mosdepth \
        -t ${threads} \
        -n \
        --by ${panel} \
        --fast-mode \
        ${id} \
        ${inputBam}
        """
}

