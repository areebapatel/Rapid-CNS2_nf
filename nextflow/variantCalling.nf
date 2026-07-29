// Pick the Clair3 model that matches the basecaller the BAM was produced with.
// A mismatched model degrades calling silently, and dorado already records the
// answer in the @RG DS field, e.g.
//     basecall_model=dna_r10.4.1_e8.2_400bps_hac@v4.1.0
// which maps to the Clair3 model directory r1041_e82_400bps_hac_v410.
// Runs in the Clair3 container so the result is validated against the models
// actually shipped in /opt/models rather than a hard-coded list.
process detectClair3Model {
    label 'clair3'

    input:
        tuple path(bam), path(bai)
        val(override)

    output:
        path "clair3_model.txt", emit: modelFile

    script:
        """
        AVAILABLE=\$(ls /opt/models)

        if [ -n "${override}" ]; then
            MODEL="${override}"
            echo "Using user-specified Clair3 model: \${MODEL}"
        else
            DORADO=\$(samtools view -H ${bam} \
                     | grep -o 'basecall_model=[^[:space:]]*' \
                     | head -n 1 | cut -d= -f2)

            if [ -z "\${DORADO}" ]; then
                echo "ERROR: no 'basecall_model=' found in the \\@RG header of ${bam}," >&2
                echo "so the Clair3 model cannot be auto-detected." >&2
                echo "Set --clair3Model explicitly. Available models:" >&2
                echo "\${AVAILABLE}" | sed 's/^/  /' >&2
                exit 1
            fi
            echo "dorado basecall model: \${DORADO}"

            # dna_r10.4.1_e8.2_400bps_hac@v4.1.0 -> r1041_e82_400bps_hac_v410
            # (drop 'dna_', strip dots from the chemistry and the version)
            MODEL=\$(echo "\${DORADO}" \
                | sed -e 's/^dna_//' -e 's/@v/_v/' \
                | awk -F'_v' '{ chem=\$1; ver=\$2;
                                gsub(/\\./, "", chem);
                                gsub(/\\./, "", ver);
                                print chem "_v" ver }')
            echo "mapped to Clair3 model: \${MODEL}"
        fi

        if ! echo "\${AVAILABLE}" | grep -qx "\${MODEL}"; then
            echo "ERROR: Clair3 model '\${MODEL}' is not present in this container." >&2
            echo "Set --clair3Model to one of:" >&2
            echo "\${AVAILABLE}" | sed 's/^/  /' >&2
            exit 1
        fi

        echo "\${MODEL}" > clair3_model.txt
        """
}

// SNV/indel calling with Clair3.
// Replaces Clara Parabricks DeepVariant: no NVIDIA licence, no GPU queue, no
// dummy read-group injection, and the ONT models are matched to the dorado
// basecaller version (see params.clair3Model).
process clair3 {
    label 'clair3'

    publishDir "${params.outDir}/snv/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(ref)
        path(refIdx)
        path(panel)
        val(id)
        val(model)
        val(threads)

    output:
        path "${id}.clair3.vcf.gz",     emit: vcf
        path "${id}.clair3.vcf.gz.tbi", emit: tbi
        path "clair3_${id}/run_clair3.log", optional: true

    script:
        """
        run_clair3.sh \
            --bam_fn=${bam} \
            --ref_fn=${ref} \
            --threads=${threads} \
            --platform=ont \
            --model_path=/opt/models/${model} \
            --bed_fn=${panel} \
            --sample_name=${id} \
            --output=clair3_${id}

        cp clair3_${id}/merge_output.vcf.gz     ${id}.clair3.vcf.gz
        cp clair3_${id}/merge_output.vcf.gz.tbi ${id}.clair3.vcf.gz.tbi
        """
}

process recodeVCF {
    label 'rapid_cns'

    publishDir "${params.outDir}/snv/", mode: 'copy'

    input:
        path(vcfGz)
        val(id)

    output:
        path "${id}.clair3.PASS.vcf.gz", emit: passVcf

    script:
        """
        vcftools --gzvcf ${vcfGz} \
            --remove-filtered-all --recode --stdout \
            | gzip -c > ${id}.clair3.PASS.vcf.gz
        """
}

process convert2annovar {
    label 'rapid_cns'

    publishDir "${params.outDir}/snv/", mode: 'copy'

    input:
        path(inputVcf)
        val(annovarPath)
        val(id)

    output:
        path "${id}_panel.avinput", emit: avinput

    script:
        """
        ${annovarPath}/convert2annovar.pl \
            -format vcf4 ${inputVcf} \
            -withfreq \
            -includeinfo \
            > ${id}_panel.avinput
        """
}

process tableAnnovar {
    label 'rapid_cns'

    publishDir "${params.outDir}/snv/", mode: 'copy'

    input:
        path(avinput)
        val(annovarPath)
        val(annovarDB)
        val(id)
        val(protocol)
        val(operation)

    output:
        path "*_multianno.csv", emit: multianno

    // -protocol and -operation must have the same number of comma-separated
    // entries; the workflow checks this before launching.
    script:
        """
        ${annovarPath}/table_annovar.pl ${avinput} \
            ${annovarDB} \
            -buildver hg38 \
            -out ${id}_panel \
            -protocol ${protocol} \
            -operation ${operation} \
            -nastring . \
            -csvout \
            -polish \
            -otherinfo
        """
}

process filterReport {
    label 'rapid_cns'

    publishDir "${params.outDir}/snv/", mode: 'copy'

    input:
        path(filterReportScript)
        path(multianno)
        val(id)

    output:
        path "${id}_dv_report.csv", emit: dvReport

    script:
        """
        Rscript ${filterReportScript} \
            --input ${multianno} \
            --output ${id}_dv_report.csv \
            --sample ${id}
        """
}

process igv_reports {
    label 'rapid_cns'

    publishDir "${params.outDir}/reports/", mode: 'copy'

    input:
        path(filteredReport)
        val(id)
        path(reference)
        tuple path(bam), path(bai)
        path(annotations)

    output:
        path "${id}_igv-report.html", emit: igvReport

    script:
        """
        sed -e 's/,/\\t/g' -e 's/"//g' \
            ${filteredReport} > ${id}_dv_report.fmt.csv

        create_report ${id}_dv_report.fmt.csv \
            ${reference} \
            --sequence 1 \
            --begin 2 \
            --end 3 \
            --flanking 1000 \
            --info-columns Chr Start End Func Gene ExonicFunc AAChange cytoBand 1000g_EUR COSMIC \
            --output ${id}_igv-report.html \
            --standalone \
            --tracks ${bam} ${annotations}
        """
}

// The wf-human-variation processes shell out to a nested Nextflow run and are
// therefore deliberately left without a container label: they execute on the
// host, which is where the nextflow binary and container engine live.
process human_variation_snp {
    errorStrategy 'ignore'

    input:
        tuple path(bam), path(bai)
        path(panel)
        path(ref)
        val(id)
        val(bamMinCoverage)
        val(snpThreads)

    output:
        val true

    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
            -with-report ${params.outDir}/human_variation_snp_nextflow_report.html \
            -profile standard \
            -w ${params.outDir}/wf-human-variation/snp/workspace \
            --ref ${ref} \
            --snp \
            --bam ${bam} \
            --bed ${panel} \
            --sample_name ${id} \
            --bam_min_coverage ${bamMinCoverage} \
            --out_dir ${params.outDir}/wf-human-variation/snp/ \
            --threads ${snpThreads}
        """
}

process human_variation_sv {
    errorStrategy 'ignore'

    input:
        tuple path(bam), path(bai)
        path(ref)
        val(id)
        val(bamMinCoverage)
        val(svThreads)

    output:
        val true

    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
            -with-report ${params.outDir}/human_variation_sv_nextflow_report.html \
            -profile standard \
            -w ${params.outDir}/wf-human-variation/sv/workspace \
            --ref ${ref} \
            --sv \
            --bam ${bam} \
            --sample_name ${id} \
            --bam_min_coverage ${bamMinCoverage} \
            --out_dir ${params.outDir}/wf-human-variation/sv/ \
            --threads ${svThreads} \
            --sniffles_args="--non-germline"
        """
}
