#!/bin/bash

OUTDIR=/run/classification
mkdir -p ${OUTDIR}

Rscript /methylclass.R --sample $SAMPLE --out_dir $OUTDIR --in_file /run/bam/${SAMPLE}_modifiedBases.5mC.hg38.bed --probes_file=/ref/$PROBES --training_data=/ref/$TRAINING --array_file=/ref/$ARRAY --threads $THREADS 
