#!/bin/bash

IN_FILE=$1
MNP_BED=$2
OUT_PATH=$3

# Extract filename from path
file=$(basename "$IN_FILE")

# Extract part before the first dot in the filename
filename=$(echo "$file" | cut -d '.' -f 1)
mkdir -p ${OUT_PATH}/${filename}/

awk '{print $1,$2,$3,$10,$11}' ${IN_FILE} > ${OUT_PATH}/${filename}.tmp.bed

perl -p -i -e 's/ /\t/g' ${OUT_PATH}/${filename}.tmp.bed

bedtools intersect  -a ${OUT_PATH}/${filename}.tmp.bed  -b ${MNP_BED}  -wa -wb > ${OUT_PATH}/${filename}.MNPFlex.bed

rm -r ${OUT_PATH}/${filename}.tmp.bed

#add column names
column_names="chr start end coverage methylation_percentage IlmnID"

# Add column names to the output file
echo -e "$column_names" > ${OUT_PATH}/${filename}.MNPFlex.subset.bed

# Append the content of the input file to the output file
awk '{print $1,$2,$3,$4,$5,$12}' ${OUT_PATH}/${filename}.MNPFlex.bed >> ${OUT_PATH}/${filename}.MNPFlex.subset.bed
