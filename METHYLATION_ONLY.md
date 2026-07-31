# Methylation-only pipeline

For users who only need methylation analysis without SNV, CNV, or structural variant calling, we provide a dedicated methylation-only workflow that is faster and uses fewer computational resources.

## Features

The methylation-only pipeline (`methylationOnly.nf`) includes:
- **5mC methylation calling** using modkit
- **MGMT promoter methylation assessment** with coverage validation
- **Methylation-based tumor classification** using the Rapid-CNS² classifier
- **MNP-Flex preparation** (optional) for external classifier analysis
- **BAM processing** with automatic alignment detection and methylation tags validation

## Usage

### Basic methylation analysis
```bash
nextflow run methylationOnly.nf \
    --input /data/sample.bam \
    --id SAMPLE001 \
    --ref /path/to/hg38.fa \
    --outDir ./methylation_results
```

### With MNP-Flex preparation
```bash
nextflow run methylationOnly.nf \
    --input /data/sample.bam \
    --id SAMPLE001 \
    --ref /path/to/hg38.fa \
    --outDir ./methylation_results \
    --mnpFlex
```

### Advanced configuration
```bash
nextflow run methylationOnly.nf \
    --input /data/sample.bam \
    --id SAMPLE001 \
    --ref /path/to/hg38.fa \
    --outDir ./methylation_results \
    --patient "John Doe" \
    --minimumMgmtCov 10 \
    --mnpFlex \
    -profile lsf,singularity
```

## Parameters

Takes the same parameters as the main pipeline (see
[Parameters in README.md](README.md#parameters)); the variant-calling options simply do not apply.


## Output

Outputs:

```
output/
├── bam/                          # panel subset BAM + index
├── mods/                         # <id>.5mC.bedmethyl
├── mgmt/                         # coverage, status, methylartist plot
├── methylation_classification/   # votes, rf_details, calibrated scores
└── mnpflex/                      # <id>.MNPFlex.input.bed
```

## When to use methylation-only

Use the methylation-only pipeline when you:
- Only need methylation analysis and tumor classification
- Want faster execution (skips computationally intensive variant calling)
- Have limited computational resources
- Are primarily interested in MGMT promoter methylation status
- Need MNP-Flex compatible files for external analysis

**Note:** The methylation-only pipeline uses the same BAM processing logic as the main pipeline, including automatic alignment detection and methylation tags validation.

## Setup

Setup is identical to the main pipeline - reference genome, containers and site
paths. See [README.md](README.md).
