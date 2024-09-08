#!/bin/bash


source $1

LIBRARY=$2
USER=$3
WS=$4 # 8, 15, 54
WS="ssh ${USER}@161.42.29.$WS"
PORT=8080

GPU_WS="161.2.29.28"
# Get the relevant information from the samplesheet
echo $LIBRARY > ont_helper.log

python - << EOF
import sys
import pandas as pd
df = pd.read_excel("tracking.xlsx")
library = open("ont_helper.log", "r").read().splitlines()[0]
print(df.columns)
df = df[df["experiment_name"] == library]

if len(df) == 0:
    print('LIBRARY NOT IN TRACKINGSHEET.')
    sys.exit(1)
df["sample_id"]=df["sample_id"].astype(int).astype(str)

df = df[["sample_id", "type", "tissue", "experiment_name", "device", "sequencing_date", "library_kit", "patho_number"]]
print(df)
df.to_csv("ont_helper.log", index=False, header=False)
EOF

if [[ ! $? == 0 ]]; then
    echo "Tracking failed."
    exit 1
fi

IFS="," read -ra DATA <<< $(cat ont_helper.log)

SAMPLE_NAME=${DATA[0]} 
SAMPLE_TYPE=${DATA[1]}
SAMPLE_TISS=${DATA[2]}

SAMPLE=${SAMPLE_NAME}-${SAMPLE_TYPE}-${SAMPLE_TISS}

SEQU_DEVICE=${DATA[4]}
SEQUEN_DATE=${DATA[5]}
LIBRARY_KIT=${DATA[6]}

PATHNUM=${DATA[7]}

OUTDIR=$TEMP_RESULTS/$SAMPLE/$LIBRARY

if [[ $SEQU_DEVICE == "grid" ]]; then
    CONFIG=$GRID 
else
    CONFIG=$PROM 
fi
echo $SAMPLE_NAME $SAMPLE_TISS $SAMPLE_TYPE $SAMPLE $SEQU_DEVICE $SEQUEN_DATE $LIBRARY_KIT $CONFIG $WS

THREADS=32

# BASECALLING
BAMFILES=$RAW_DIR/$LIBRARY/$SAMPLE_NAME/bam
FAST5=$RAW_DIR/$LIBRARY/$SAMPLE_NAME

#dorado
echo $RAW_DIR/$LIBRARY

ssh np-system@$GPU_WS "echo \"secret\" | sudo -S rsync -rtv --progress $RAW_DIR/$LIBRARY /data/raw"

ssh ${USER}@$GPU_WS "docker run --rm --gpus '\"device=5,6\"' -e DEVICE=\"cuda:all\" -e PORT=4040 -e NUM_CLIENTS=5 -e CONFIG=$CONFIG -v /data/raw/${LIBRARY}:/logs -v /data/raw/${LIBRARY}/${SAMPLE_NAME}/:/pod5 -v /data/raw/${LIBRARY}/${SAMPLE_NAME}/bam:/output -v $REFERENCE_PATH:/ref -e LIBRARY=$LIBRARY -e REFERENCE=hg19.fa --net=host $REGISTRY/$DORDAO"

if [[ ! $? == 0 ]]; then 
    echo "Basecalling failed."
    exit 1
fi

echo EXTRA STARTED

# METHYLATION EXTRACTION
ssh ${USER}@$GPU_WS "docker run --rm  -e THREADS=$THREADS -e SAMPLE=$SAMPLE -e REF=$REF -e PANEL=$PANEL -v $REFERENCE_PATH:/ref -v ${OUTDIR}:/run -v /data/raw/$LIBRARY/$SAMPLE_NAME/bam/pass/:/data $REGISTRY/$EXTRACT_MODS"

# QUALITY
ssh ${USER}@$GPU_WS "docker run --rm  -e THREADS=$THREADS -v /data/raw/$LIBRARY/$SAMPLE_NAME/bam:/bam -e SAMPLE=$SAMPLE -v ${OUTDIR}/:/run $REGISTRY/$QUALITY"

wait
notify-send "GPU done"

# CNV Calling
docker run --rm  -v ${OUTDIR}:/run -e SAMPLE=$SAMPLE -e THREADS=$THREADS $REGISTRY/$CNV &

# Methylation Classification
$WS "docker run --rm  -v ${OUTDIR}:/run -e SAMPLE=$SAMPLE -m 128g -v $REFERENCE_PATH:/ref -e PROBES=top_probes_hm450.Rdata -e TRAINING=capper_top_100k_betas_binarised.Rdata -e ARRAY=HM450.hg38.manifest.gencode.v22.Rdata -e THREADS=30 $REGISTRY/$METHCLASS" &

# MGMT Prediction
docker run --rm  -v ${OUTDIR}:/run -e SAMPLE=$SAMPLE -e PANEL=mgmt_hg38.bed -e PROBES=mgmt_probes.Rdata -e MODEL=mgmt_137sites_mean_model.Rdata -v $REFERENCE_PATH:/ref -e THREADS=$THREADS $REGISTRY/$MGMT

# Coverage
$WS "docker run --rm  -v ${OUTDIR}:/run -e SAMPLE=$SAMPLE -e PANEL=$PANEL -v $REFERENCE_PATH:/ref $REGISTRY/$COVERAGE" &

# SNV Calling Deepvariant
BAM=/run/bam/${SAMPLE}.processed.ontarget.bam

$WS "docker run --rm  -v $OUTDIR:/run -v $REFERENCE_PATH:/ref $REGISTRY/$DEEPVARIANT run_pepper_margin_deepvariant call_variant -b $BAM -f /ref/$REF -o /run/deepvariant -p $SAMPLE -t $THREADS --ont_r10_q20 --phased_output" &

wait

# Variant Annotation
$WS "docker run --rm  -e TOOL=deepvariant -e SAMPLE=${SAMPLE} -v $REFERENCE_PATH:/ref -v $OUTDIR:/run $REGISTRY/$ANNOVAR" &

# Special Positions
docker run --rm  -v ${OUTDIR}:/run -e SAMPLE=$SAMPLE -v ${REFERENCE_PATH}:/ref -e REF=$REF $REGISTRY/$SPECIAL 

wait

# Final Report
if [[ $SEQU_DEVICE == "grid" ]]; then
    DEVICE="GridION-MK1"
else
    DEVICE="P2-SOLO"
fi

docker run --rm  -e PATHNUM=$PATHNUM -e SAMPLE=${SAMPLE} -e LIBRARY=${LIBRARY} -e SEQDATE=$SEQUEN_DATE -e PREPKIT=$LIBRARY_KIT -e DEVICE=$DEVICE -e PANEL=ADAPTIVE -e REF=$REF -v $REFERENCE_PATH:/ref -v ${OUTDIR}:/run $REGISTRY/$REPORT

notify-send "ONT done"
