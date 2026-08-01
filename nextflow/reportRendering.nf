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
        path(fusionTables, stageAs: 'fusions/*')
        path(egfrTables,   stageAs: 'egfrviii/*')
        path(purityPloidy,  stageAs: 'opt_purity_ploidy.tsv')
        path(savanaCnvPlot, stageAs: 'opt_savana_cnv.png')
        path(mutations)
        path(cnvPlot)
        path(rfDetails)
        path(votes)
        path(coverage)
        // distinct staged names: all three fall back to the same NO_FILE
        // placeholder, which Nextflow cannot stage more than once per task
        path(mgmtStatus, stageAs: 'opt_mgmt_status.csv')
        path(mgmtPlot,   stageAs: 'opt_mgmt_plot.png')
        path(igvReport,  stageAs: 'opt_igv_report.html')
        val(mgmtAvgCov)

    output:
        path "${id}_Rapid-CNS2_report*", emit: reports

    script:
        """
        # Optional inputs arrive as an empty placeholder file when not produced
        MGMT_ARG=""
        if [ -s opt_mgmt_status.csv ]; then
            MGMT_ARG="--mgmt opt_mgmt_status.csv"
        else
            echo "Info: no MGMT status file (coverage below threshold) - omitting --mgmt"
        fi

        METHYLARTIST_ARG=""
        if [ -s opt_mgmt_plot.png ]; then
            METHYLARTIST_ARG="--methylartist opt_mgmt_plot.png"
        else
            echo "Info: no methylartist plot (coverage below threshold) - omitting --methylartist"
        fi

        SAVPLOT_ARG=""
        if [ -s opt_savana_cnv.png ]; then
            SAVPLOT_ARG="--savana_cnv_plot opt_savana_cnv.png"
        fi

        PP_ARG=""
        if [ -s opt_purity_ploidy.tsv ]; then
            PP_ARG="--purity_ploidy opt_purity_ploidy.tsv"
        fi

        IGV_ARG=""
        if [ -s opt_igv_report.html ]; then
            IGV_ARG="--igv_report opt_igv_report.html"
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
            \${PP_ARG} \
            \${SAVPLOT_ARG} \
            --fusions fusions \
            --egfrviii egfrviii \
            --software_ver ${softwareVersion} \
            --promoter_mgmt_coverage ${mgmtAvgCov} \
            --report_HTML ${reportHTML}
        """
}
