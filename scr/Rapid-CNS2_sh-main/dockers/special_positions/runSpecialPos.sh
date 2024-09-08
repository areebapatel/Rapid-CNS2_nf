#!/bin/bash

OUTDIR=/run/special_pos
mkdir -p $OUTDIR

#TERTp mutations
echo "TERT"
bcftools mpileup -d 8000 -Q 7 -a DP,ADF,ADR -r chr5:1295228-1295228 -X ont -f /ref/$REF /run/bam/${SAMPLE}.processed.ontarget.bam | bcftools view --no-header > ${OUTDIR}/${SAMPLE}_TERT1.vcf
bcftools mpileup -Q 7 -d 8000 -a DP,ADF,ADR -r chr5:1295250-1295250 -X ont -f /ref/$REF /run/bam/${SAMPLE}.processed.ontarget.bam | bcftools view --no-header > ${OUTDIR}/${SAMPLE}_TERT2.vcf

#BRAF V600E mutation
echo BRAF
bcftools mpileup -Q 7 -d 8000 -r chr7:140453135-140453137 -a DP,ADF,ADR -X ont -f /ref/$REF /run/bam/${SAMPLE}.processed.ontarget.bam | bcftools view --no-header > ${OUTDIR}/${SAMPLE}_BRAFV600.vcf

#IDH1/2 mutations
echo IDH1/2
bcftools mpileup -Q 7 -d 8000 -r chr2:209113111-209113113 -a DP,ADF,ADR -X ont -f /ref/$REF /run/bam/${SAMPLE}.processed.ontarget.bam | bcftools view --no-header > ${OUTDIR}/${SAMPLE}_IDH1R132.vcf
bcftools mpileup -Q 7 -d 8000 -r chr15:90631837-90631839 -a DP,ADF,ADR -X ont -f /ref/$REF /run/bam/${SAMPLE}.processed.ontarget.bam | bcftools view --no-header > ${OUTDIR}/${SAMPLE}_IDH2R172.vcf

#H3FA mutations
echo H3FA
bcftools mpileup -Q 7 -d 8000 -r chr1:226252134-226252136 -a DP,ADF,ADR -X ont -f /ref/$REF /run/bam/${SAMPLE}.processed.ontarget.bam | bcftools view --no-header > ${OUTDIR}/${SAMPLE}_H3F3AK27.vcf
bcftools mpileup -Q 7 -d 8000 -r chr1:226252155-226252157 -a DP,ADF,ADR -X ont -f /ref/$REF /run/bam/${SAMPLE}.processed.ontarget.bam | bcftools view --no-header > ${OUTDIR}/${SAMPLE}_H3FAG34.vcf

python3 /merge_specialPos.py