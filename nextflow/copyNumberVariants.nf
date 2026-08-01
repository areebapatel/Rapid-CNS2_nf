process copyNumberVariants {
    label 'heavy'

    publishDir "${params.outDir}/cnv/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        val(id)

    output:
        path "${id}.cnvpytor.calls.1000.tsv",   emit: calls1000
        path "${id}.cnvpytor.calls.10000.tsv",  emit: calls10000
        path "${id}.cnvpytor.calls.100000.tsv", emit: calls100000
        path "${id}_cnvpytor_100k.pdf",         emit: plotPdf
        path "${id}_cnvpytor_100k.png",         emit: plotPng

    script:
        def chroms = (1..22).collect { n -> "chr${n}" }.join(' ') + ' chrX chrY'
        """
        export TMPDIR="\$PWD" MPLCONFIGDIR="\$PWD/.mpl" XDG_CACHE_HOME="\$PWD/.cache"
        mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

        cnvpytor -root ${id}_CNV.pytor -rd ${bam} -j ${task.cpus}
        cnvpytor -root ${id}_CNV.pytor -his 1000 10000 100000 -j ${task.cpus}
        cnvpytor -root ${id}_CNV.pytor -partition 1000 10000 100000 -j ${task.cpus}
        cnvpytor -root ${id}_CNV.pytor -call 1000   -j ${task.cpus} > ${id}.cnvpytor.calls.1000.tsv
        cnvpytor -root ${id}_CNV.pytor -call 10000  -j ${task.cpus} > ${id}.cnvpytor.calls.10000.tsv
        cnvpytor -root ${id}_CNV.pytor -call 100000 -j ${task.cpus} > ${id}.cnvpytor.calls.100000.tsv
        cnvpytor -root ${id}_CNV.pytor -plot manhattan 100000 -chrom ${chroms} -o ${id}_cnvpytor_100k.pdf
        cnvpytor -root ${id}_CNV.pytor -plot manhattan 100000 -chrom ${chroms} -o ${id}_cnvpytor_100k.png

        # cnvpytor writes <name>.global.0000.<ext>, not the -o name it is given
        for E in pdf png; do
            [ -f ${id}_cnvpytor_100k.\$E ] || \\
                mv ${id}_cnvpytor_100k.global.0000.\$E ${id}_cnvpytor_100k.\$E
        done
        """
}

process cnvAnnotated {
    label 'rapid_cns'

    publishDir "${params.outDir}/cnv/", mode: 'copy'

    input:
        path(cnvCalls)
        val(id)
        path(annotateScript)
        path(cnvGenes)

    output:
        path "${id}.annotation.1000.xlsx", emit: cnvAnnotated

    script:
        """
        python3 ${annotateScript} ${cnvCalls} ${cnvGenes} ${id}.annotation.1000.xlsx
        """
}

// SAVANA tumour-only SV calling. Its breakpoints feed `savana cna` below, which
// gives better copy-number segmentation than read depth alone.
process savanaTo {
    label 'savana'

    publishDir "${params.outDir}/sv/savana/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(ref)
        path(refIdx)
        val(id)

    output:
        path "savana_to_${id}/**", emit: allOutputs
        path "${id}.savana.breakpoints.vcf", emit: breakpoints

    script:
        """
        # LSF exports a node-local TMPDIR that does not exist inside the
        # container; savana's final bcftools sort fails on it. XDG_CACHE_HOME
        # keeps fontconfig out of the (unwritable) home directory.
        export TMPDIR="\$PWD" MPLCONFIGDIR="\$PWD/.mpl" XDG_CACHE_HOME="\$PWD/.cache"
        mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

        savana to \
            --tumour ${bam} \
            --ref ${ref} \
            --ref_index ${refIdx} \
            --sample ${id} \
            --outdir savana_to_${id} \
            --threads ${task.cpus}

        # prefer the somatic call set; savana also writes an unfiltered
        # classified.vcf that is ~10x larger and not what downstream wants
        VCF=\$(find savana_to_${id} -name "*classified.somatic.vcf" | head -n 1)
        if [ -z "\${VCF}" ]; then
            VCF=\$(find savana_to_${id} -name "*.vcf" | grep -vi "germline" | head -n 1)
        fi
        [ -n "\${VCF}" ] || { echo "ERROR: savana to produced no VCF" >&2; exit 1; }
        echo "using \${VCF}"
        cp "\${VCF}" ${id}.savana.breakpoints.vcf
        """
}

// SAVANA copy number: absolute CN fitting with tumour purity and ploidy
// estimated from B-allele frequencies at 1000G het SNP sites.
//
// NOTE: CN segmentation and purity/ploidy fitting assume reasonably uniform
// genome-wide coverage. Adaptive-sampling panel data is strongly enriched
// on-target with sparse off-target background, so these estimates should be
// validated against known samples before being trusted clinically.
process savanaCna {
    label 'savana'

    publishDir "${params.outDir}/cnv/savana/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(ref)
        path(breakpoints)
        val(id)
        val(g1000)

    output:
        path "savana_cna_${id}/**", emit: allOutputs
        path "${id}_purity_ploidy.tsv", optional: true, emit: purityPloidy
        path "${id}_segmented_absolute_copy_number.tsv", optional: true, emit: segments

    script:
        """
        export TMPDIR="\$PWD" MPLCONFIGDIR="\$PWD/.mpl" XDG_CACHE_HOME="\$PWD/.cache"
        mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

        savana cna \
            --tumour ${bam} \
            --ref ${ref} \
            --sample ${id} \
            --outdir savana_cna_${id} \
            --breakpoints ${breakpoints} \
            --g1000_vcf ${g1000} \
            --cna_threads ${task.cpus} \
            --tmpdir . ${params.savanaCnaArgs}

        # surface the segmentation and the purity/ploidy fit for the report
        SEG=\$(find savana_cna_${id} -name "*segmented_absolute_copy_number.tsv" | head -n 1)
        [ -n "\${SEG}" ] && cp "\${SEG}" ${id}_segmented_absolute_copy_number.tsv

        FIT=\$(find savana_cna_${id} -iname "*purity_ploidy*" | head -n 1)
        if [ -n "\${FIT}" ]; then
            cp "\${FIT}" ${id}_purity_ploidy.tsv
        else
            echo "WARNING: no purity/ploidy fit produced by savana cna" >&2
        fi
        """
}

// Genome-wide absolute copy number from the SAVANA segmentation. Complements
// the CNVpytor read-depth plot: values here are purity- and ploidy-corrected,
// and the CNS genes in data/genes.bed are annotated on the profile.
process savanaCnvPlot {
    label 'rapid_cns'

    publishDir "${params.outDir}/cnv/", mode: 'copy'

    input:
        path(segments)
        path(purityPloidy)
        path(cnvGenes)
        path(plotScript)
        val(id)

    output:
        path "${id}_savana_cnv.png", emit: plot, optional: true

    script:
        """
        Rscript ${plotScript} \
            --segments ${segments} \
            --purity_ploidy ${purityPloidy} \
            --genes ${cnvGenes} \
            --sample ${id} \
            --out ${id}_savana_cnv.png
        """
}
