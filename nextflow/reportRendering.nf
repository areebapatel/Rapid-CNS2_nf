process reportRendering {
    label 'rapid_cns'

    // cache false forces the report to be regenerated on each run
    cache false

    publishDir "${params.outDir}/report/", mode: 'copy'

    input:
        path(reportScript)
        path(reportHTML)
        val(id)
        val(patient)
        val(softwareVersion)
        val(seq)
        tuple path(inputBam), path(inputBai)
        path(mutations)
        path(cnvPlot)
        path(rfDetails)
        path(votes)
        path(coverage)
        // distinct staged names: all three fall back to the same NO_FILE
        // placeholder, which Nextflow cannot stage more than once per task
        path(mgmtStatus, stageAs: 'opt_mgmt_status')
        path(mgmtPlot,   stageAs: 'opt_mgmt_plot')
        path(igvReport,  stageAs: 'opt_igv_report')
        val(mgmtAvgCov)

    output:
        path "${id}_Rapid-CNS2_report*", emit: reports

    script:
        """
        # An explicit --seq overrides sequencer auto-detection
        if [ "${seq}" != "false" ]; then
            seq="${seq}"
        else
            # try the @RG header first
            RG_seq=\$(samtools view -@4 -H ${inputBam} | grep ^@RG | grep -Po "PM:.*?\\t" | awk '{print substr(\$NF,4,3)}' || true)
            if [ -n "\${RG_seq}" ]; then
                seq="\${RG_seq}"
            else
                FN_seq=\$(samtools view ${inputBam} | grep -Po "fn:Z:[F,P]" | head -n 1 | awk '{print substr(\$NF,6,6)}' || true)
                if [ -n "\${FN_seq}" ]; then
                    seq="\${FN_seq}"
                else
                    F5_seq=\$(samtools view ${inputBam} | grep -Po "f5:Z:[F,P]" | head -n 1 | awk '{print substr(\$NF,6,6)}' || true)
                    seq="\${F5_seq}"
                fi
            fi
        fi
        if [ -z "\${seq}" ] || [ "\${seq}" == "false" ]; then
            seq="Unknown"
        fi

        # Optional inputs arrive as an empty placeholder file when not produced
        MGMT_ARG=""
        if [ -s opt_mgmt_status ]; then
            MGMT_ARG="--mgmt opt_mgmt_status"
        else
            echo "Info: no MGMT status file (coverage below threshold) - omitting --mgmt"
        fi

        METHYLARTIST_ARG=""
        if [ -s opt_mgmt_plot ]; then
            METHYLARTIST_ARG="--methylartist opt_mgmt_plot"
        else
            echo "Info: no methylartist plot (coverage below threshold) - omitting --methylartist"
        fi

        IGV_ARG=""
        if [ -s opt_igv_report ]; then
            IGV_ARG="--igv_report opt_igv_report"
        else
            echo "Info: no IGV report - omitting --igv_report"
        fi

        Rscript ${reportScript} \
            --prefix ${id} \
            --sample ${id} \
            --patient "${patient}" \
            --mutations ${mutations} \
            --cnv_plot ${cnvPlot} \
            --rf_details ${rfDetails} \
            --votes ${votes} \
            --coverage ${coverage} \
            \${MGMT_ARG} \
            \${METHYLARTIST_ARG} \
            \${IGV_ARG} \
            --software_ver ${softwareVersion} \
            --seq "\${seq}" \
            --promoter_mgmt_coverage ${mgmtAvgCov} \
            --report_HTML ${reportHTML}
        """
}
