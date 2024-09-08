#!/usr/bin/env python

OUTDIR=/run/reports
mkdir -p $OUTDIR

if [[ $PANEL == "nan" ]]; then
    LONGSHOT=/run/longshot/${SAMPLE}*NPHD2022APADDEDPANEL.csv
    DEEPVARIANT=/run/deepvariant/${SAMPLE}*NPHD2022APADDEDPANEL.csv
else
    LONGSHOT=/run/longshot/${SAMPLE}*filtered.csv
    DEEPVARIANT=/run/deepvariant/${SAMPLE}*filtered.csv
fi

python3 /scripts/transform_variants.py ${OUTDIR}/${SAMPLE} $DEEPVARIANT

if [ -d /run/fusions ]; then
    python3 /scripts/transform_fusions.py ${OUTDIR}/${SAMPLE} /run/fusions/${SAMPLE}.fusions.csv
fi

TISSUE=$(echo $SAMPLE | cut -d\- -f2)

python3 /scripts/plot_classification.py /run/classification/${SAMPLE}_votes.tsv -o /run/reports/predictions.png

Rscript /scripts/make_pdfreport.R --rf_details=/run/classification/${SAMPLE}_rf_details.tsv --votes=/run/classification/${SAMPLE}_votes.tsv --indels=${OUTDIR}/${SAMPLE}.indels.csv --snv=${OUTDIR}/${SAMPLE}.snvs.csv --fusions=${OUTDIR}/${SAMPLE}.fusions.csv --cnv_plot=/run/cnv/${SAMPLE}_cnvpytor_100k.global.0000.png --coverage=/run/coverage/${SAMPLE}.mosdepth.summary.txt --output_dir=${OUTDIR} --prefix=${SAMPLE} --mgmt=/run/mgmt/${SAMPLE}_mgmt_status.csv --library=${LIBRARY} --sample=${SAMPLE} --tissue=${TISSUE} --library=$LIBRARY --refgen=$REF --panel=$PANEL --prepkit="$PREPKIT" --seqdate=$SEQDATE --device=$DEVICE --specpos=/run/special_pos/summary.xlsx --classplot=${OUTDIR}/predictions.png
