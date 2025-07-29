# Rapid-CNS2 Module Testing Framework

This document describes how to test individual modules of the Rapid-CNS2 pipeline separately.

## Overview

The testing framework allows you to test each module independently, which is useful for:
- **Debugging** specific components
- **Development** of new features
- **Validation** of individual tools
- **Performance** testing of specific modules
- **Troubleshooting** issues in isolation

## Available Test Modules

### 1. BAM Processing (`bam`)
Tests the BAM processing pipeline including:
- BAM alignment (if unaligned)
- BAM merging (if multiple files)
- BAM indexing
- BAM subsetting to target regions

### 2. SNV Calling (`snv`)
Tests the SNV calling pipeline including:
- DeepVariant for variant calling
- VCF filtering and recoding
- ANNOVAR annotation
- Variant filtering and reporting

### 3. Structural Variant Calling (`sv`)
Tests the structural variant pipeline including:
- Sniffles2 for SV detection
- AnnotSV for SV annotation

### 4. Methylation Analysis (`methylation`)
Tests the methylation analysis pipeline including:
- Modkit for methylation calling
- MGMT coverage analysis
- MGMT promoter methylation analysis
- MGMT prediction

### 5. Methylation Classification (`classification`)
Tests the methylation classification module:
- Machine learning classification
- Methylation subtype prediction

### 6. Copy Number Variation (`cnv`)
Tests the CNV analysis pipeline including:
- CNVpytor for CNV calling
- CNV annotation

### 7. Coverage Analysis (`coverage`)
Tests the coverage analysis module:
- Mosdepth for coverage calculation

### 8. Report Generation (`report`)
Tests the report generation module:
- Comprehensive report creation

### 9. All Modules (`all`)
Tests the complete pipeline end-to-end.

## Usage

### Basic Usage

```bash
# Test all modules
nextflow run run_tests.nf \
    --test_module all \
    --input_bam sample.bam \
    --ref hg38.fa \
    --id test_sample \
    --annovarPath /path/to/annovar \
    --annovarDB /path/to/humandb

# Test only BAM processing
nextflow run run_tests.nf \
    --test_module bam \
    --input_bam sample.bam \
    --ref hg38.fa \
    --id test_sample

# Test only SNV calling (requires processed BAM)
nextflow run run_tests.nf \
    --test_module snv \
    --input_bam processed_sample.bam \
    --ref hg38.fa \
    --id test_sample \
    --annovarPath /path/to/annovar \
    --annovarDB /path/to/humandb

# Test only methylation classification (requires bedmethyl file)
nextflow run run_tests.nf \
    --test_module classification \
    --bedmethylFile sample.5mC.bedmethyl \
    --id test_sample
```

### Advanced Usage

```bash
# Test with custom parameters
nextflow run run_tests.nf \
    --test_module all \
    --input_bam sample.bam \
    --ref hg38.fa \
    --id test_sample \
    --threads 8 \
    --numGpus 2 \
    --outDir ./custom_output \
    --annovarPath /path/to/annovar \
    --annovarDB /path/to/humandb \
    --snifflesThreads 16 \
    --modkitThreads 16 \
    --cnvThreads 16 \
    --methThreads 16 \
    --minimumMgmtCov 10

# Test with specific compute profile
nextflow run run_tests.nf \
    --test_module all \
    --input_bam sample.bam \
    --ref hg38.fa \
    --id test_sample \
    -profile docker \
    --annovarPath /path/to/annovar \
    --annovarDB /path/to/humandb
```

## Test Data Requirements

### For BAM Processing Tests
- **Input**: Raw BAM file(s) (aligned or unaligned)
- **Reference**: hg38 reference genome
- **Output**: Processed BAM files

### For SNV Calling Tests
- **Input**: Processed BAM file (indexed)
- **Reference**: hg38 reference genome
- **ANNOVAR**: Installed and configured
- **Output**: VCF files and annotated variants

### For Structural Variant Tests
- **Input**: Processed BAM file (indexed)
- **Reference**: hg38 reference genome
- **Output**: VCF files and annotated SVs

### For Methylation Analysis Tests
- **Input**: Processed BAM file with methylation tags (MM:Z)
- **Reference**: hg38 reference genome
- **Output**: Bedmethyl files and MGMT analysis

### For Methylation Classification Tests
- **Input**: Bedmethyl file from methylation analysis
- **Output**: Classification results and predictions

### For CNV Analysis Tests
- **Input**: Processed BAM file (indexed)
- **Output**: CNV calls and annotations

### For Coverage Analysis Tests
- **Input**: Processed BAM file (indexed)
- **Output**: Coverage summaries

### For Report Generation Tests
- **Input**: All analysis outputs
- **Output**: Comprehensive report

## Test Output Structure

```
test_output/
├── bam/                    # BAM processing outputs
│   ├── alignedBams/
│   └── merged_bam.bam
├── snv/                    # SNV calling outputs
│   ├── *.deepvariant.vcf
│   ├── *.PASS.vcf.gz
│   └── *_report.csv
├── sv/                     # Structural variant outputs
│   ├── *.sniffles2.vcf
│   └── *_annotsv.tsv
├── mods/                   # Methylation outputs
│   └── *.5mC.bedmethyl
├── mgmt/                   # MGMT analysis outputs
│   ├── mgmt_cov.mosdepth.summary.txt
│   ├── *_mgmt.bed
│   └── *_mgmt_status.csv
├── methylation_classification/  # Classification outputs
│   ├── *_rf_details.tsv
│   └── *_votes.tsv
├── cnv/                    # CNV outputs
│   ├── *.cnvpytor.calls.*.tsv
│   └── *_cnvpytor_100k.pdf
├── coverage/               # Coverage outputs
│   └── *.mosdepth.summary.txt
└── report/                 # Report outputs
    └── *_report.html
```

## Troubleshooting

### Common Issues

1. **Missing ANNOVAR**
   ```
   Error: ANNOVAR path not found
   ```
   **Solution**: Install ANNOVAR and provide correct paths:
   ```bash
   --annovarPath /path/to/annovar/
   --annovarDB /path/to/annovar/humandb/
   ```

2. **Missing Methylation Tags**
   ```
   Error: No methylation tags (MM:Z) found in the BAM file
   ```
   **Solution**: Re-run basecalling with modified basecalling enabled:
   ```bash
   # For Dorado
   dorado basecaller --modified-bases 5mC ...
   
   # For MinKNOW
   # Enable "Modified bases" in the protocol settings
   ```

3. **GPU Issues**
   ```
   Error: GPU not available
   ```
   **Solution**: Use CPU-only mode or check GPU configuration:
   ```bash
   --numGpus 0  # Use CPU only
   ```

4. **Memory Issues**
   ```
   Error: Out of memory
   ```
   **Solution**: Reduce thread counts or increase memory:
   ```bash
   --threads 4
   --snifflesThreads 4
   --modkitThreads 4
   ```

5. **Container Issues**
   ```
   Error: Container not found
   ```
   **Solution**: Pull required containers:
   ```bash
   docker pull areebapatel/rapid_cns_3.0.0
   docker pull ontresearch/modkit
   docker pull nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1
   docker pull eichlerlab/sniffles:2.3.2
   ```

### Debugging Tips

1. **Check Logs**
   ```bash
   # View Nextflow logs
   tail -f .nextflow.log
   
   # View process logs
   ls work/*/*/.command.log
   ```

2. **Resume Failed Runs**
   ```bash
   nextflow run run_tests.nf -resume
   ```

3. **Clean Cache**
   ```bash
   nextflow clean -f
   ```

4. **Test Individual Processes**
   ```bash
   # Test specific process
   nextflow run run_tests.nf --test_module bam --input_bam sample.bam --ref hg38.fa --id test_sample -with-trace
   ```

## Performance Testing

### Benchmarking Individual Modules

```bash
# Test BAM processing performance
time nextflow run run_tests.nf --test_module bam --input_bam large_sample.bam --ref hg38.fa --id benchmark_test

# Test SNV calling performance
time nextflow run run_tests.nf --test_module snv --input_bam processed_sample.bam --ref hg38.fa --id benchmark_test --annovarPath /path/to/annovar --annovarDB /path/to/humandb

# Test with different thread configurations
nextflow run run_tests.nf --test_module all --input_bam sample.bam --ref hg38.fa --id test_sample --threads 16 --snifflesThreads 32 --modkitThreads 32 --cnvThreads 32
```

### Resource Monitoring

```bash
# Monitor resource usage
htop
nvidia-smi  # For GPU monitoring

# Check Nextflow execution reports
open test_output/pipeline_info/*_report.html
```

## Development Workflow

### Adding New Tests

1. **Create test workflow** in `test_modules.nf`
2. **Add test parameters** in `run_tests.nf`
3. **Update documentation** in this file
4. **Test the new module** with sample data

### Example: Adding a New Module Test

```groovy
// In test_modules.nf
workflow test_new_module {
    take:
        input_file
        id
    
    main:
        new_module_out = newModule(input_file, id)
    
    emit:
        results = new_module_out
}

// In run_tests.nf
case 'new_module':
    test_new_module(input_file, id)
    break
```

## Best Practices

1. **Start Small**: Test with small datasets first
2. **Use Isolated Tests**: Test modules individually before full pipeline
3. **Monitor Resources**: Keep track of memory and CPU usage
4. **Validate Outputs**: Check that outputs are as expected
5. **Document Issues**: Keep notes of any problems encountered
6. **Version Control**: Track changes to test configurations

## Support

For issues with the testing framework:
1. Check the troubleshooting section above
2. Review Nextflow documentation: https://www.nextflow.io/docs/
3. Check the pipeline logs in the `work/` directory
4. Create an issue in the repository with detailed error information 