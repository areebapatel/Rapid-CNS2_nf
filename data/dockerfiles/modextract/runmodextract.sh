#!/bin/bash

OUTDIR=/run/mods/
INDIR=/run/bam/

mkdir -p ${OUTDIR}

# generate bam file and index
samtools merge -@ ${THREADS} ${OUTDIR}/${SAMPLE}.processed.bam /data/*.bam 
samtools index ${OUTDIR}/${SAMPLE}.processed.bam

# filter overall bam file to only keep reads in panel
bedtools intersect -a ${OUTDIR}/${SAMPLE}.processed.bam -b /ref/${PANEL} > ${OUTDIR}/${SAMPLE}.processed.ontarget.bam
samtools index ${OUTDIR}/${SAMPLE}.processed.ontarget.bam

# extract
modbam2bed /ref/${REF} ${OUTDIR}/${SAMPLE}.processed.bam -m 5mC --cpg --threads=${THREADS}  > ${OUTDIR}/${SAMPLE}_modifiedBases.5mC.hg19.bed

# liftover for classification
liftOver ${OUTDIR}/${SAMPLE}_modifiedBases.5mC.hg19.bed /applications/hg19ToHg38.over.chain.gz ${OUTDIR}/${SAMPLE}_modifiedBases.5mC.hg38.bed ${OUTDIR}/${SAMPLE}_modifiedBases.5mC.unmapped.bed -bedPlus=3