#!/bin/bash

OUTDIR=/run/pycoqc

mkdir -p $OUTDIR

pycoQC --summary_file /bam/sequencing_summary*.txt --bam_file /bam/pass/*.bam /bam/fail/*.bam --html_outfile ${OUTDIR}/${SAMPLE}_pycoQC.html --report_title "PycoQC Report -- ${SAMPLE}" --json_outfile ${OUTDIR}/${SAMPLE}_pycoQC.json