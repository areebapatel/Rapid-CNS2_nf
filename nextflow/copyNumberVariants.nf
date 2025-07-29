process copyNumberVariants {
    label 'rapid_cns'

    input:
        path(bam)
        path(bai)
        val(id)
        val(cnvThreads)

    output:
        val true
    
    publishDir("${params.outDir}/cnv/")

    script:
        """
        mkdir -p ${params.outDir}/cnv/
        cd ${params.outDir}/cnv/
        cnvpytor -root ${id}_CNV.pytor -rd ${bam} -j ${cnvThreads}
        cnvpytor -root ${id}_CNV.pytor -his 1000 10000 100000 -j ${cnvThreads} 
        cnvpytor -root ${id}_CNV.pytor -partition 1000 10000 100000 -j ${cnvThreads} # SLOW
        cnvpytor -root ${id}_CNV.pytor -call 1000 -j ${cnvThreads} > ${id}.cnvpytor.calls.1000.tsv
        cnvpytor -root ${id}_CNV.pytor -call 10000 -j ${cnvThreads} > ${id}.cnvpytor.calls.10000.tsv
        cnvpytor -root ${id}_CNV.pytor -call 100000 -j ${cnvThreads} > ${id}.cnvpytor.calls.100000.tsv
        cnvpytor -root ${id}_CNV.pytor -plot manhattan 100000 -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY -o ${id}_cnvpytor_100k.pdf
        """
}

process cnvAnnotated{
    label 'rapid_cns'
    
    input:
        val(ready)
        val(id)
        path(annotateScript)
        path(cnvGenes)
        path(outDir)
    
    publishDir("${params.outDir}/cnv/", mode: 'copy')

    output:
        path "${id}.annotation.1000.xlsx", emit: cnv_annotated
        
    script:
    """
    python3 ${annotateScript} ${outDir}/cnv/${id}.cnvpytor.calls.1000.tsv ${cnvGenes} ${outDir}/cnv/${id}.annotation.1000.xlsx
    """
}