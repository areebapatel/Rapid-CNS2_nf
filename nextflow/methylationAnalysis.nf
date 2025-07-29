process methylationCalls {
    label 'mods'

    input:
        path(bam)
        path(bai)
        path(ref)
        val(id)
        val(modkitThreads)

    publishDir("${params.outDir}/mods/")

    output:
        path "${id}.5mC.bedmethyl", emit: bedmethylFile

    script:
        """
        mkdir -p ${params.outDir}/mods/
        modkit pileup ${bam} ${params.outDir}/mods/${id}.5mC.bedmethyl --ref ${ref} --preset traditional --only-tabs --threads ${modkitThreads}
        """
}

process checkMgmtCoverage {
    label 'rapid_cns'

    input:
        path(bam)
        path(bai)
        path(mgmtBed)
        val(minimumMgmtCov)
        val(threads)

    publishDir("${params.outDir}/mgmt")

    output:
	    val true
	    path "*_cov.txt", emit: mgmt_avg_cov_file
        path "mgmt_cov.mosdepth.summary.txt"	
        stdout emit: mgmt_avg_cov

    script:
        """
        mkdir -p ${params.outDir}/mgmt/
        mosdepth \
        -t ${threads} \
        -n \
        --by ${mgmtBed} \
        mgmt_cov ${bam}
        
        # Check if the coverage is below the threshold
        cov="\$(grep "^chr10_region" mgmt_cov.mosdepth.summary.txt | awk '{ print \$4 }')"
        
        echo \${cov}
        if awk 'BEGIN{exit ARGV[1]>ARGV[2]}' "\$cov" ${minimumMgmtCov}
        then
            echo \${cov} > mgmt_below_thresh_cov.txt
        else
            echo \${cov} > mgmt_avg_cov.txt
        fi
	"""
}

process mgmtPromoterMethyartist {
    label 'rapid_cns'
    
    input:
        path(bam)
        path(bai)
        path(ref)
        val(ready)
        val(id)

    publishDir("${params.outDir}/mgmt/")
    
    output:
        val true
        path "*.svg", emit: mgmt_plot optional true
    
    script:
        cov_file = file("${params.outDir}/mgmt/mgmt_avg_cov.txt")
        // Check if the coverage is above the threshold
        if ( cov_file.exists() == true )    
            """
            methylartist \
            locus \
            -i chr10:129456536-129477536 \
            -l chr10:129466536-129467536 \
            -b ${bam} \
            --ref ${ref} \
            --motif CG \
            --mods m \
            --highlightpalette viridis \
            --samplepalette magma > ${id}_mgmt.svg

            """
        else
            // If the coverage is below the threshold, do nothing
            """
            echo "MGMT coverage is below the threshold, skipping MGMT promoter methylation analysis"
            """
}

process mgmtPred {
    label 'rapid_cns'
    
    input:
        val(ready)
        path(mgmtScript)
        path(mgmtBed)
        path(mgmtProbes)
        path(mgmtModel)
        path(bedmethylFile)
        val(id)

    output:
        val true

    script:
        cov_file = file("${params.outDir}/mgmt/mgmt_avg_cov.txt")
        if ( cov_file.exists() == true )
            """
            bedtools \
            intersect \
            -a ${bedmethylFile} \
            -b ${mgmtBed} > ${params.outDir}/mgmt/${id}_mgmt.bed \

            Rscript ${mgmtScript} \
            --input ${params.outDir}/mgmt/${id}_mgmt.bed \
            --probes ${mgmtProbes} \
            --model ${mgmtModel} \
            --out_dir ${params.outDir}/mgmt/ \
            --sample ${id} \
            """
        else
            """
            echo "MGMT coverage is below the threshold, skipping MGMT promoter methylation analysis"   
            """
}
