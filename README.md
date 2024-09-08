<div align="left">
<img alt="Screenshot 2024-03-17 at 10 41 16 PM" src="https://github.com/areebapatel/Rapid-CNS2_nf/assets/46373444/6ba243da-0dca-4f4e-9cea-a4df7b989ff6" width="200" height="200" style="float: left; margin-right: 10px;">
<h1 style="display: inline-block;">Rapid-CNS<sup>2</sup> workflow</h1>
</div>

## Overview

The Rapid-CNS<sup>2</sup> nextflow pipeline is a bioinformatics workflow designed for comprehensive analysis of genomic and epigenomic data generated using adaptive sampling based sequencing of central nervous system (CNS) tumours. It performs tasks such as basecalling, variant calling, methylation analysis, structural variant calling, copy number variation calling, and provides a comprehensive molecular diagnostic-ready report.

This pipeline is implemented using Nextflow, allowing for easy execution and scalability on various compute environments, including local machines, clusters, and cloud platforms.

## Features

- Modular architecture for easy customization and extension.
- Supports both basecalling from raw ONT POD5s and analysis of pre-aligned BAM files.
- Accelerated variant calling with Clara Parabricks supported Deepvariant and Sniffles2
- Annotation and filtering of clinically relevant variants
- Includes methylation analysis with Rapid-CNS<sup>2</sup> classifier and MGMT promoter methylation status determination.
- Automated report generation for summarizing analysis results.
- Prepare input files for the MNP-Flex classifier (optional)

## Requirements

- Nextflow (version 3.0.0 or later)
- Conda, Docker or Singularity (optional, for containerized execution of tools)
- Required input data:
  - Raw ONT POD5 data (for basecalling) or pre-aligned BAM files
  - Reference genome file (hg38 required)

## Usage

1. Clone this repository:

    ```bash
    git clone https://github.com/areebapatel/Rapid-CNS2_nf.git
    ```

2. Edit the `nextflow.config` file to configure pipeline parameters according to your requirements.

3. Run the pipeline using Nextflow:

    ```bash
    nextflow run main.nf --input <input_directory> --id <sample_identifier> [--options]
    ```

    Replace `<input_directory>` with the path to the directory containing ONT POD5 data or pre-aligned BAM files, and `<sample_identifier>` with a unique identifier for the sample.

    Additional options can be specified to customize pipeline behavior. Use the `--help` option to view available options and their descriptions.

4. Monitor pipeline progress and access results in the specified output directory.

## Sequencing
This pipeline analyses CNS tumour data generated through Nanopore adaptive sampling using [ReadFish](https://github.com/LooseLab/readfish) or adaptive sampling on MinKNOW. It is compatible with data generated on MinION, GridION and PromethION

## Parameters

| Parameter            | Description                                                                                                        | Default Value        |
|----------------------|--------------------------------------------------------------------------------------------------------------------|----------------------|
| `--input`            | Path to the directory POD5 files for Dorado basecalling and minimap2 alignment or BAM file if available            | (Required) |
| `--id`               | Sample identifier                                                                                                  | (Required) |
| `--ref`              | Path to hg19 reference file                                                                                        | `null`               |
| `--tmp_dir`          | Directory to store temporary files. If it does not exist it will be created                                        | `tempDir`            |
| `--out_dir`          | Directory path to store all the outputs                                                                            | `output`             |
| `--log_dir`          | Directory to store log files                                                                                       | `logDir`             |
| `--minimum_mgmt_cov` | Minimum coverage for MGMT promoter methylation analysis                                                            | `5`                  |
| `--model_config`     | Basecalling model to be used                                                                                       | `dna_r10.4.1_e8.2_400bps_hac@v4.1.0` |
| `--port`             | Port for basecall server                                                                                           | `8887`               |
| `--num_gpu`          | Number of GPUs to use                                                                                              | `3`                  |
| `--num_clients`      | Number of clients                                                                                                  | `par.num_gpu * 3`    |
| `--help`             | Show help message                                                                                                  | `null`               |
| `--test`             | Run in test mode                                                                                                   | `null`               |
| `--reads`            | Samtools addreplacerg -r option. Specify as `-r "SM:GM24385" -r "ID:GM24385"`                                      | `null`               |
| `--basecalling`      | Enable basecalling from raw ONT POD5 data. If provided, `--input` should point to the directory containing raw data. | `false`              |
| `--mnp-flex`         | Prepare input file for the MNP-Flex classifier.                                                                   | `false`              |


## Acknowledgements
We are extremely grateful to all our lab members and collaborators for their support! 
Keeping up with AI to make our life easier and to compensate for our (Areeba's) art skills, our logo was generated by DALL-E.

## Citation
If you use this pipeline, please cite our preprint:

Felix Sahm, Areeba Patel, Kirsten Göbel et al. Versatile, accessible cross-platform molecular profiling of central nervous system tumors: web-based, prospective multi-center validation, 10 April 2024, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-4182910/v1]

## Contributions
Contributions are welcome! If you encounter any issues, have suggestions for improvements, or would like to contribute new features, please open an issue or pull request on this repository.

## License

This project is licensed under the [MIT License](LICENSE).
