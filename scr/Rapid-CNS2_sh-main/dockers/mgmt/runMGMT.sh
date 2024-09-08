#!/bin/bash

OUTDIR=/run/mgmt

mkdir -p ${OUTDIR}

bedtools intersect -a /run/bam/${SAMPLE}_modifiedBases.5mC.hg38.bed -b /ref/$PANEL > ${OUTDIR}/${SAMPLE}_mgmt.5mC.hg38.bed

Rscript mgmt_pred.R --input=${OUTDIR}/${SAMPLE}_mgmt.5mC.hg38.bed --out_dir=${OUTDIR} --probes=/ref/$PROBES --model=/ref/$MODEL --sample=$SAMPLE