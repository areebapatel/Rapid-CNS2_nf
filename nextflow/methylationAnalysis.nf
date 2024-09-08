process check_bam_has_meth_tags {
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

process modkit_adjust_mods {
    input:
        path(inputBam)
        val(id)
        val(modkit_threads)

    publishDir("${params.outDir}/bam")

    output:
        path "*_modkit_merge.bam", emit: modkit_merged_bam
        path "*_modkit_merge.bam.bai", emit: modkit_merged_bai

    script:
        """
        modkit \
        adjust-mods \
        --convert h m \
        ${inputBam} \
        ${id}_modkit_merge.bam \
        --threads ${modkit_threads}

        samtools index ${id}_modkit_merge.bam
        """
}

process methylationCalls {
    maxRetries 1
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    publishDir("${params.outDir}/mods/")
    input:
        path(inputBam)
        path(inputBai)
        path(ref)
        val(id)
        val(modkitThreads)
        path(liftOver)
        path(liftOverChain)

    output:
        path "${id}.mods.bedmethyl", emit: bedmethyl_file

    script:
        """
        mkdir ${params.outDir}/mods/
        modkit pileup \
        ${inputBam} \
        ${id}.mods.bedmethyl \
        --ref ${ref} \
        --threads ${modkitThreads}
    
        """
}


process liftOver_ch{
    input:
        path(bedmethyl_file)
        path(liftOver)
        path(liftOverChain)
        val(id)
    
    output:
        path "${id}.mods.hg38.bedmethyl", emit: bedmethyl_file_hg38

    script:
        """
        ${liftOver}/liftOver ${bedmethyl_file} ${liftOverChain} ${params.outDir}/mods/${id}.mods.hg38.bedmethyl ${params.outDir}/mods/${id}.mods.5mC.unmapped.bed -bedPlus=3
        """
}

process check_mgmt_coverage {
    input:
        path(inputBam)
	    path(mgmtBed)
	    val(minimum_mgmt_cov)
        val(threads)

    publishDir("${params.outDir}/mgmt")

    output:
	val true
	path "*_cov.txt", emit: mgmt_avg_cov_file
    path "mgmt_cov.mosdepth.summary.txt"	
    stdout emit: mgmt_avg_cov

    script:
        """
        mosdepth \
        -t ${threads} \
        -n \
        --by ${mgmtBed} \
        mgmt_cov ${inputBam}
        
        cov="\$(grep "^chr10_region" mgmt_cov.mosdepth.summary.txt | awk '{ print \$4 }')"
        
        echo \${cov}
        if awk 'BEGIN{exit ARGV[1]>ARGV[2]}' "\$cov" ${minimum_mgmt_cov}
        then
            echo \${cov} > mgmt_below_thresh_cov.txt
        else
            echo \${cov} > mgmt_avg_cov.txt
        fi
	"""
}

process mgmtPromoter_methyartist {
    input:
        path(inputBam)
        path(inputBai)
        path(ref)
	    val ready 

    publishDir("${params.outDir}/mgmt/")
    
    output:
        val true
        path "*.svg", emit: mgmt_plot optional true
    
    script:
        cov_file = file("${params.outDir}/mgmt/mgmt_avg_cov.txt")
        if ( cov_file.exists() == true )    
            """
            methylartist \
            locus \
            -i chr10:131263800-131266800 \
            -l chr10:131264800-131265800 \
            -b ${inputBam} \
            --ref ${ref} \
            --motif CG \
            --mods m \
            --highlightpalette viridis \
            --samplepalette magma > ${id}_mgmt.svg

            """
        else

            """
            exit 1
            """
}

process mgmtPred {
    input:
        val ready
        path(mgmtScript)
        path(mgmtBed)
        path(mgmtProbes)
        path(mgmtModel)
        val(id)

    output:
        val true

    script:
        cov_file = file("${params.outDir}/mgmt/mgmt_avg_cov.txt")
        if ( cov_file.exists() == true )
            """
            bedtools \
            intersect \
            -a ${PWD}/${params.outDir}/mods/${id}.mods.hg38.bedmethyl
            -b ${mgmtBed} > ${PWD}/${params.outDir}/mgmt/${id}_mgmt_hg38.bed \

            Rscript ${mgmtScript} \
            --input ${PWD}/${params.outDir}/mgmt/${id}_mgmt_hg38.bed \
            --probes ${mgmtProbes} \
            --model ${mgmtModel} \
            --out_dir ${params.outDir}/mgmt/ \
            --sample ${id} \
            """
        else
            """
            """
}

process mnpFlex {
    errorStrategy 'ignore'
    input:
        path(mnpFlexScript)
        val ready
        path(bedmethyl_file_hg38)
        path(mnpFlexBed)
    
    output:
        val true
    
    publishDir("${params.outDir}/mnpflex/")

    script:
        """
        bash ${mnpFlexScript} ${bedmethyl_file_hg38} ${mnpFlexBed} ${params.outDir}/mnpflex/
        """
}