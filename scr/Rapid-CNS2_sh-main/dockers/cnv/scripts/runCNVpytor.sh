#!/bin/bash 

OUTDIR=/run/cnv
mkdir -p $OUTDIR

FILE=${OUTDIR}/${SAMPLE}_CNV.pytor

cnvpytor -root $FILE -rd /run/bam/${SAMPLE}.processed.bam

if [[ $RIMP == "TRUE" ]]; then
    cnvpytor -root $FILE -his 1000 100000 100000 1000000 -j $THREADS
    cnvpytor -root $FILE -partition 1000 100000 100000 1000000 -j $THREADS
else 
    cnvpytor -root $FILE -his 1000 100000 100000 -j $THREADS
    cnvpytor -root $FILE -partition 1000 100000 100000 -j $THREADS
fi

cnvpytor -root $FILE -call 1000 -j $THREADS > $OUTDIR/${SAMPLE}.cnvpytor.calls.1000.tsv
cnvpytor -root $FILE -call 10000 -j $THREADS > $OUTDIR/${SAMPLE}.cnvpytor.calls.10000.tsv
cnvpytor -root $FILE -call 100000 -j $THREADS > $OUTDIR/${SAMPLE}.cnvpytor.calls.100000.tsv

if [[ $RIMP == "TRUE" ]]; then
    cnvpytor -root $FILE -call 1000000 -j $THREADS > $OUTDIR/${SAMPLE}.cnvpytor.calls.1000000.tsv
    VIEW=1000000
else 
    VIEW=100000
fi

cnvpytor -root $FILE -plot manhattan $VIEW -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr2 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX -o ${OUTDIR}/${SAMPLE}_cnvpytor_${RIMP}.png

python3 /scripts/annotate.py $OUTDIR/${SAMPLE}.cnvpytor.calls.1000.tsv $OUTDIR/${SAMPLE}.annotation.1000.xlsx