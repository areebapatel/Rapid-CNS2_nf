#!/bin/bash

OUTDIR=/run/coverage

mkdir -p ${OUTDIR}

./mosdepth -n --by /ref/${PANEL} --fast-mode ${OUTDIR}/${SAMPLE} /run/bam/${SAMPLE}.processed.bam

python3 /coverage_plot.py ${OUTDIR}/${SAMPLE}.mosdepth.summary.txt

