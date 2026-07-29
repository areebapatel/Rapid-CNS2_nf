#!/bin/bash
#
# Optional helper: basecall one ONT run with dorado, emitting a modified-base
# BAM suitable as --input for Rapid-CNS2_nf.
#
# Basecalling is NOT part of the Rapid-CNS2 pipeline. This script is provided
# only as a convenience for users who do not already have a modBAM. The
# alternative is ONT's own workflow:
#     nextflow run epi2me-labs/wf-basecalling ...
#
# The pipeline requires MM/ML tags, so a modified-base model is mandatory
# (--mods, default 5mCG_5hmCG).
#
# ---------------------------------------------------------------------------
# Usage
#   scr/basecall.sh --pod5 <dir> --sample <id> --outdir <dir> [options]
#
# Required
#   --pod5     <dir>    Directory of .pod5 files for one run/flow cell
#   --sample   <id>     Sample identifier; output is <outdir>/<id>/<id>.bam
#   --outdir   <dir>    Output root
#
# Optional
#   --ref      <fasta>  Reference FASTA. If given, dorado aligns during
#                       basecalling and the pipeline can consume the BAM
#                       directly. If omitted, an unaligned modBAM is written
#                       and Rapid-CNS2_nf will align it itself. With --ref the
#                       output is also coordinate-sorted and indexed, so the
#                       pipeline skips its own sort.
#   --model    <name>   Dorado model (default: hac). Use 'sup' for the
#                       super-accurate model, or pass an explicit path.
#   --mods     <codes>  Modified bases (default: 5mCG_5hmCG)
#   --models-dir <dir>  Where models live / are downloaded to (default:
#                       <outdir>/models). Pre-populating this avoids every
#                       job re-downloading, and lets jobs run on nodes
#                       without internet access.
#   --device   <spec>   Dorado device (default: cuda:all)
#   --threads  <n>      Threads for sorting/indexing (default: 16)
#   --sort-mem <size>   samtools sort memory per thread (default: 3G)
#
# Example - single run:
#   scr/basecall.sh --pod5 /data/run1/pod5 --sample S001 \
#                   --outdir /results/dorado --ref /refs/hg38.fa
#
# Example - submit to LSF with 4 GPUs (adjust for your scheduler):
#   bsub -q gpu -n 16 -R "rusage[mem=96GB] span[hosts=1]" \
#        -gpu num=4:j_exclusive=yes:gmem=16G \
#        -o basecall.%J.out -e basecall.%J.err \
#        scr/basecall.sh --pod5 /data/run1/pod5 --sample S001 \
#                        --outdir /results/dorado --ref /refs/hg38.fa
#
# ---------------------------------------------------------------------------
# Resuming
#   Large flow cells (1-2 TB of pod5) can exceed a scheduler's wall-clock
#   limit. Output is written to numbered <sample>.part<N>.bam files; re-running
#   the exact same command resumes from the newest part via dorado
#   --resume-from and only basecalls the reads still missing. The final
#   <sample>.bam is created only after dorado exits cleanly, so a truncated
#   file is never mistaken for a finished one.
# ---------------------------------------------------------------------------
set -euo pipefail

POD5_DIR=""; SAMPLE=""; OUTROOT=""; REF=""
MODEL="hac"; MODS="5mCG_5hmCG"; MODELS_DIR=""; DEVICE="cuda:all"
THREADS=16; SORT_MEM="3G"

usage() { sed -n '2,60p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while [ $# -gt 0 ]; do
    case "$1" in
        --pod5)        POD5_DIR="$2"; shift 2 ;;
        --sample)      SAMPLE="$2"; shift 2 ;;
        --outdir)      OUTROOT="$2"; shift 2 ;;
        --ref)         REF="$2"; shift 2 ;;
        --model)       MODEL="$2"; shift 2 ;;
        --mods)        MODS="$2"; shift 2 ;;
        --models-dir)  MODELS_DIR="$2"; shift 2 ;;
        --device)      DEVICE="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        --sort-mem)    SORT_MEM="$2"; shift 2 ;;
        -h|--help)     usage 0 ;;
        *) echo "Unknown argument: $1" >&2; usage 1 ;;
    esac
done

[ -n "${POD5_DIR}" ] || { echo "ERROR: --pod5 is required" >&2; usage 1; }
[ -n "${SAMPLE}" ]   || { echo "ERROR: --sample is required" >&2; usage 1; }
[ -n "${OUTROOT}" ]  || { echo "ERROR: --outdir is required" >&2; usage 1; }
[ -d "${POD5_DIR}" ] || { echo "ERROR: pod5 directory not found: ${POD5_DIR}" >&2; exit 1; }

command -v dorado >/dev/null 2>&1 || {
    echo "ERROR: dorado not found on PATH." >&2
    echo "Install it, or load it first (e.g. 'module load dorado/2.0.0')." >&2
    exit 1
}

[ -z "${REF}" ] || command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools not found on PATH (needed to sort the aligned output)." >&2
    exit 1
}

MODELS_DIR="${MODELS_DIR:-${OUTROOT}/models}"
OUTDIR="${OUTROOT}/${SAMPLE}"
mkdir -p "${OUTDIR}" "${MODELS_DIR}"

echo "=================================================================="
echo "sample     : ${SAMPLE}"
echo "pod5 dir   : ${POD5_DIR} ($(ls "${POD5_DIR}"/*.pod5 2>/dev/null | wc -l) files)"
echo "model      : ${MODEL}   mods: ${MODS}"
echo "reference  : ${REF:-<none - unaligned output>}"
echo "host       : $(hostname)"
echo "started    : $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "no GPU detected"
echo "=================================================================="

if [ -f "${OUTDIR}/${SAMPLE}.bam" ]; then
    echo "${OUTDIR}/${SAMPLE}.bam already exists - nothing to do."
    exit 0
fi

rm -f "${OUTDIR}/${SAMPLE}.sorting.bam" "${OUTDIR}/${SAMPLE}.sorting.bam.bai"

# Resume from the newest completed part, if there is one
PREV=$(ls -1t "${OUTDIR}/${SAMPLE}".part*.bam 2>/dev/null | head -n 1 || true)
NEXT=0
RESUME_ARG=()
if [ -n "${PREV}" ]; then
    LAST=$(basename "${PREV}" .bam | sed "s/^${SAMPLE}\.part//")
    NEXT=$(( LAST + 1 ))
    echo "Resuming from ${PREV} -> part${NEXT}"
    RESUME_ARG=(--resume-from "${PREV}")
else
    echo "No previous part found - starting from scratch"
fi

REF_ARG=()
[ -n "${REF}" ] && REF_ARG=(--reference "${REF}")

OUT="${OUTDIR}/${SAMPLE}.part${NEXT}.bam"

dorado basecaller \
    "${MODEL}" \
    "${POD5_DIR}" \
    --models-directory "${MODELS_DIR}" \
    --modified-bases "${MODS}" \
    --device "${DEVICE}" \
    "${REF_ARG[@]}" \
    "${RESUME_ARG[@]}" \
    > "${OUT}"

# Sort after dorado exits, not piped into it: samtools sort writes only at the
# end, so a piped run killed at the wall clock leaves nothing to --resume-from.
if [ -n "${REF}" ]; then
    samtools sort -@ "${THREADS}" -m "${SORT_MEM}" -o "${OUTDIR}/${SAMPLE}.sorting.bam" "${OUT}"
    samtools index -@ "${THREADS}" "${OUTDIR}/${SAMPLE}.sorting.bam"
    mv "${OUTDIR}/${SAMPLE}.sorting.bam"     "${OUTDIR}/${SAMPLE}.bam"
    mv "${OUTDIR}/${SAMPLE}.sorting.bam.bai" "${OUTDIR}/${SAMPLE}.bam.bai"
else
    mv "${OUT}" "${OUTDIR}/${SAMPLE}.bam"
fi
rm -f "${OUTDIR}/${SAMPLE}".part*.bam

echo "=================================================================="
echo "finished   : $(date)"
ls -la "${OUTDIR}/${SAMPLE}.bam"
echo
echo "Feed this to the pipeline with:"
echo "  nextflow run main.nf --input ${OUTDIR}/${SAMPLE}.bam --id ${SAMPLE} ..."
echo "=================================================================="
