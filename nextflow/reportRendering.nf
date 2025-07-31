process reportRendering {

    // cache false forces it to regenerate the report each time and not use the cache
    cache false
    
    publishDir("${params.outDir}/report/", mode: 'copy')

    input:
        path(reportScript)
        val(ready) // cnvpytor
        val(ready) // mgmt_pred
        val(ready) // meth_classification
        val(ready) // filter_report
        val(id)
        path(mosdepth_plot_data) // mosdepth
        val(mgmt_cov)
        val(mgmtPromoterMethyartist) 
        val(igv_reports) //igv_reports has run
        val(software_version)
        path(inputBam)
        val(seq)
        path(report_UKHD)
    
    output:
	val true

    script:
        """
        # check for a specified sequencer, this overrides the checks below
        if [ "${seq}" != "false" ]
        then
            seq=${seq}
        else
            # try to get the sequencer model from the @RG group (if it exists)
            RG_seq=\$(samtools view -@4 -H ${inputBam} | grep ^@RG | grep -Po "PM:.*?\t" | awk '{print substr(\$NF,4,3)}')
            # if found it, save as seq
            if [ "\$RG_seq" ]
            then 
                seq=\$RG_seq
            else
            # if didn't find @RG, try for the fn:Z tag
                FN_seq=\$(samtools view ${inputBam} | grep -Po "fn:Z:[F,P]" | head -n 1 | awk '{print substr(\$NF,6,6)}')        
                # if found it, save as seq
                if [ "\$FN_seq" ]
                then
                    seq=\$FN_seq
                else
                    # try the f5:Z tag
                    F5_seq=\$(samtools view ${inputBam} | grep -Po "f5:Z:[F,P]" | head -n 1 | awk '{print substr(\$NF,6,6)}')
                    if [ "\$F5_seq" ]
                    then
                        seq=\$F5_seq
                    fi
                fi
            fi
        fi
        ## if all else fails, set it as unknown
        if [ -z "\${seq}" ] || [ "\${seq}" == "false" ]
        then
            seq="Unknown"
        fi
       
        MGMT_ARG=""
        if [ -f "${params.outDir}/mgmt/${id}_mgmt_status.csv" ]; then
            MGMT_ARG="--mgmt ${params.outDir}/mgmt/${id}_mgmt_status.csv"
        else
            echo "Info: MGMT status file not found or coverage below threshold, skipping MGMT argument."
        fi

        METHYLARTIST_ARG=""
        if compgen -G "${params.outDir}/mgmt/*.locus.meth.png" > /dev/null; then
            METHYLARTIST_ARG="--methylartist ${params.outDir}/mgmt/*.locus.meth.png"
        else
            echo "Info: Methylartist PNG file(s) not found or coverage below threshold, skipping methylartist argument."
        fi

        Rscript ${reportScript} \
          --prefix ${id} \
          --mutations ${params.outDir}/snv/${params.id}_dv_report.csv \
          --cnv_plot ${params.outDir}/cnv/${id}_cnvpytor_100k.pdf \
          --rf_details ${params.outDir}/methylation_classification/${id}_rf_details.tsv \
          --votes ${params.outDir}/methylation_classification/${id}_votes.tsv \
          --output_dir ${params.outDir}/report/ \
          --patient ${params.patient} \
          --coverage ${params.outDir}/coverage/${id}.mosdepth.summary.txt \
          --sample ${id} \
          $MGMT_ARG \
          $METHYLARTIST_ARG \
          --igv_report ${params.outDir}/snv/${id}_igv-report.html \
          --software_ver ${software_version} \
          --seq \${seq} \
          --promoter_mgmt_coverage ${mgmt_cov} \
          --report_UKHD ${report_UKHD} 
        """
}
