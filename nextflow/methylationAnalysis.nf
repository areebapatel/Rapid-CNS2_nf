process methylationCalls {
    label 'mods'

    publishDir "${params.outDir}/mods/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(ref)
        val(id)
        val(modkitThreads)

    output:
        path "${id}.5mC.bedmethyl", emit: bedmethyl

    // modkit >=0.6 removed --preset, --only-tabs and --ignore, and made
    // --modified-bases required. --combine-mods reports 5mC and 5hmC as a single
    // combined fraction, which is what the HM450 array betas represent.
    // Output is all-tab delimited: col 10 = Nvalid_cov, col 11 = percent modified.
    script:
        """
        modkit pileup ${bam} ${id}.5mC.bedmethyl \
            --ref ${ref} \
            --cpg \
            --combine-strands \
            --combine-mods \
            --modified-bases C \
            --threads ${modkitThreads}
        """
}

process checkMgmtCoverage {
    label 'rapid_cns'

    publishDir "${params.outDir}/mgmt/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(mgmtBed)
        val(minimumMgmtCov)
        val(threads)

    output:
        path "mgmt_cov.mosdepth.summary.txt"
        path "mgmt_avg_cov.txt", emit: avgCovFile
        path "mgmt_cov_ok.txt",  emit: covOkFile

    script:
        """
        mosdepth -t ${threads} -n --by ${mgmtBed} mgmt_cov ${bam}

        COV=\$(awk '\$1=="chr10_region"{print \$4}' mgmt_cov.mosdepth.summary.txt)
        [ -n "\${COV}" ] || COV=0
        OK=\$(awk -v c="\${COV}" -v m=${minimumMgmtCov} 'BEGIN{print (c>=m) ? "true" : "false"}')

        echo "MGMT promoter mean coverage: \${COV}x (threshold ${minimumMgmtCov}x, pass=\${OK})"
        printf '%s' "\${COV}" > mgmt_avg_cov.txt
        printf '%s' "\${OK}"  > mgmt_cov_ok.txt
        """
}

process mgmtPromoterMethyartist {
    label 'rapid_cns'

    publishDir "${params.outDir}/mgmt/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(ref)
        val(covOk)
        val(id)

    output:
        path "${id}_mgmt_methylartist.png", optional: true, emit: mgmtPlot

    script:
        """
        if [ "${covOk}" != "true" ]; then
            echo "MGMT coverage below threshold - skipping methylartist plot"
            exit 0
        fi

        methylartist locus \
            -i chr10:129465536-129468536 \
            -l chr10:129466536-129467536 \
            -b ${bam} \
            --ref ${ref} \
            --motif CG \
            --mods m \
            --highlightpalette viridis \
            --samplepalette magma

        # methylartist names its own output; normalise it for the report.
        PLOT=\$(ls *.locus.meth.png 2>/dev/null | head -n 1)
        if [ -n "\${PLOT}" ]; then
            mv "\${PLOT}" ${id}_mgmt_methylartist.png
        else
            echo "WARNING: methylartist produced no plot" >&2
        fi
        """
}

process mgmtPred {
    label 'rapid_cns'

    publishDir "${params.outDir}/mgmt/", mode: 'copy'

    input:
        val(covOk)
        path(mgmtScript)
        path(mgmtBed)
        path(mgmtProbes)
        path(mgmtModel)
        path(bedmethyl)
        val(id)

    output:
        path "${id}_mgmt_status.csv", optional: true, emit: mgmtStatus
        path "${id}_mgmt.bed",        optional: true

    script:
        """
        if [ "${covOk}" != "true" ]; then
            echo "MGMT coverage below threshold - skipping MGMT prediction"
            exit 0
        fi

        bedtools intersect -a ${bedmethyl} -b ${mgmtBed} > ${id}_mgmt.bed

        Rscript ${mgmtScript} \
            --input ${id}_mgmt.bed \
            --probes ${mgmtProbes} \
            --model ${mgmtModel} \
            --out_dir . \
            --sample ${id}
        """
}
