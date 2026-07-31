<div align="left">
<img width="342" alt="Rapid-CNS2" src="https://github.com/user-attachments/assets/cda166c2-664f-4286-951a-309b111c1132">

<h1 style="display: inline-block;">Rapid-CNS<sup>2</sup> workflow</h1>
</div>

[![Paper](https://img.shields.io/badge/Nature%20Medicine-2025-b31b1b.svg)](https://www.nature.com/articles/s41591-025-03562-5)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/Nextflow-23.10%2B-brightgreen)](https://www.nextflow.io/)

## 🧬 Overview

The Rapid-CNS<sup>2</sup> nextflow pipeline is a bioinformatics workflow designed for comprehensive analysis of genomic and epigenomic data generated using adaptive sampling based sequencing of central nervous system (CNS) tumours. It performs tasks such as alignment, SNV calling, structural variant calling, methylation analysis, copy number variation calling, and provides a comprehensive molecular report.

This pipeline is implemented using Nextflow, allowing for easy execution and scalability on various compute environments, including local machines, clusters, and cloud platforms.

## ✨ Features

- **Modular architecture** for easy customization and extension
- **Flexible input handling** - supports aligned and unaligned BAM files with automatic alignment detection
- **SNV calling with Clair3**, with the model auto-detected from the BAM's basecaller
- **Structural variants from Sniffles2 and Severus**, both annotated with AnnotSV
- **Copy number from CNVpytor, plus SAVANA** for copy number, tumour purity and ploidy
- **Comprehensive analysis** including methylation analysis with Rapid-CNS² classifier and MGMT promoter methylation status
- **Automated reporting** with molecular diagnostic-ready reports
- **MNP-Flex integration** preparing input for upload to [app.epignostix.com](https://app.epignostix.com) (on by default)

## 🔧 Requirements

- **Nextflow:** version 23.10.0 or later (enforced by `manifest.nextflowVersion`)
- **Container Engine:** Singularity/Apptainer (typical on HPC) or Docker
- **Java:** OpenJDK 8 or later
- **System:** Linux (Ubuntu 18.04+, CentOS 7+, or similar)
- **Memory:** Minimum 8GB RAM, recommended 32GB+ for large datasets
- **Storage:** At least 100GB free space for reference genomes and databases

## 🚀 Quick start

### 1. Clone the repository
```bash
git clone https://github.com/areebapatel/Rapid-CNS2_nf.git
cd Rapid-CNS2_nf
```

### 2. Install dependencies

#### Install Nextflow
```bash
# Using Conda (recommended)
conda create -n nextflow python=3.9
conda activate nextflow
conda install -c bioconda nextflow

# Or manual installation
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

#### Install container engine
```bash
# Docker (recommended)
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker $USER

# Or Singularity
sudo apt-get update
sudo apt-get install -y singularity-container
```

### 3. Download reference genome
```bash
# Create reference directory
mkdir -p /path/to/references/hg38

# Download UCSC hg38 reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
mv hg38.fa /path/to/references/hg38/hg38.fa

# Create index files
samtools faidx /path/to/references/hg38/hg38.fa
```

### 4. Install ANNOVAR

ANNOVAR is required for variant annotation:

1. **Register and download:**
   - Visit [ANNOVAR Download Form](https://www.openbioinformatics.org/annovar/annovar_download_form.php)
   - Fill out the registration form with your institutional email
   - Download the package from the link sent to your email

2. **Install and setup:**
```bash
# Extract ANNOVAR
tar -xzf annovar.latest.tar.gz
cd annovar
chmod +x *.pl

# Create humandb directory and download databases
mkdir humandb/

# cytoBand is a UCSC track, so it is fetched WITHOUT -webfrom annovar.
# Using -webfrom annovar for it silently fails.
./annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/

# The rest are hosted by ANNOVAR
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240917 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp151 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar 1000g2015aug humandb/
```

The five databases above are the pipeline default (`params.annovarProtocol`) and
all download freely. `1000g2015aug` unpacks into per-population files, so the
`1000g2015aug_eur` protocol resolves to `hg38_EUR.sites.2015_08.txt`.

Optionally, if you have access to them, add these to **both**
`params.annovarProtocol` and `params.annovarOperation` (the two lists must stay
the same length):

```bash
# cosmic70 requires a COSMIC licence; dbnsfp is ~30-40 GB
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp47a humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar allofus humandb/
```

**Note:** ANNOVAR is freely available for personal, academic, and non-profit use only. Commercial users must purchase a license from [QIAGEN](https://digitalinsights.qiagen.com/).

### 5. Install AnnotSV

AnnotSV is required for structural variant annotation:

```bash
# Navigate to installation directory
cd /path/to/install/

# Clone AnnotSV repository
git clone https://github.com/lgmgeo/AnnotSV.git

# Navigate to AnnotSV directory
cd AnnotSV

# Install AnnotSV
make PREFIX=. install

# Install human annotations
make PREFIX=. install-human-annotation
```

The annotations must be the 3.5 bundle, to match AnnotSV 3.5.10 in the container.

**Note:** AnnotSV is distributed under the GNU General Public License v3.0. See the [AnnotSV repository](https://github.com/lgmgeo/AnnotSV) for details.

### 6. Configure the pipeline

Four site paths must be set; the pipeline fails at startup with a clear message
if any is missing. Either pass them on the command line, edit `nextflow.config`,
or keep them in an institutional config.

```groovy
params {
    ref          = "/path/to/references/hg38/hg38.fa"
    annovarPath  = "/path/to/annovar"
    annovarDB    = "/path/to/annovar/humandb"
    // directory CONTAINING Annotations_Human, i.e. AnnotSV's -annotationsDir
    annotsvAnnot = "/path/to/AnnotSV/share/AnnotSV"

    // filesystems that must be visible inside containers
    containerBindPaths = "/data,/refs"
}
```

Site configs can also live in `conf/` and be selected as a profile.
`conf/dkfz.config` is a worked example, enabled with `-profile dkfz`; copy it
for your own site and add a matching entry to the `profiles` block. It has no
effect unless the profile is selected.

### 7. Run the pipeline

#### Basic run
```bash
nextflow run main.nf \
    --input /data/sample.bam \
    --id SAMPLE001 \
    --outDir ./results \
    -profile lsf,singularity
```

#### Advanced run
```bash
nextflow run main.nf \
    --input /data/sample.bam \
    --id SAMPLE001 \
    --patient "John Doe" \
    --outDir ./results \
    --minimumMgmtCov 10 \
    --mnpFlex true \
    --runHumanVariation true \
    -profile slurm,singularity
```

For methylation analysis only, without SNV/SV/CNV calling, see
[METHYLATION_ONLY.md](METHYLATION_ONLY.md).

## 📥 Input requirements

### Sequencing requirements

**⚠️ IMPORTANT: This pipeline is specifically designed for Oxford Nanopore Technologies (ONT) data generated using adaptive sampling with the gene panel described in Patel et al. 2022 and Patel et al. 2025.**

- **Platform:** Oxford Nanopore Technologies (MinION, GridION, PromethION)
- **Sequencing Method:** Adaptive sampling with targeted gene panel
- **Gene Panel:** NPHD panel (160 gene regions) as described in:
  - Patel et al. (2022). Rapid-CNS²: rapid comprehensive adaptive nanopore-sequencing of CNS tumors. *Acta Neuropathologica* 143, 609–612
  - Patel et al. (2025). Prospective, multicenter validation of a platform for rapid molecular profiling of central nervous system tumors. *Nature Medicine* 31, 1567–1577

**⚠️ NOT SUITABLE FOR:**
- Shallow whole genome sequencing (WGS)
- Other sequencing platforms (Illumina, PacBio, etc.)
- Non-targeted sequencing approaches

**⚠️ WGS WARNING:**
If using whole genome sequencing data, the average coverage should be at least **10X** for reliable reporting of variants. Lower coverage may result in:
- Reduced sensitivity for variant detection
- Increased false negative rates
- Unreliable methylation analysis
- Poor MGMT promoter methylation assessment

### Basecalling

**Basecalling must be performed externally.**
- Use [epi2me-labs/wf-basecalling](https://github.com/epi2me-labs/wf-basecalling) or Dorado directly
- Ensure you use a model that supports modified basecalling (see [Dorado documentation](https://github.com/nanoporetech/dorado?tab=readme-ov-file#modified-basecalling))
- Provide the resulting BAM(s) as input to this pipeline

### Input options

The pipeline accepts:
- **(Preferred) Single aligned BAM file:** Direct path to a single aligned and merged BAM file (e.g., `/path/to/sample.bam`)
- **Directory with aligned BAM files:** Path to directory containing multiple aligned BAM files (will be merged automatically)
- **Directory with unaligned BAM files:** Path to directory containing multiple unaligned BAM files (will be aligned and merged automatically)

## 🗺️ Pipeline structure

```mermaid
graph TD
    A[Input BAM] --> E[prepareBam - align if needed, merge, sort, index]
    E --> P[Panel subset - NPHD BED]
    P --> F[SNV calling - Clair3]
    E --> G[Methylation - modkit]
    E --> H[SV - Sniffles2 + Severus]
    E --> N[CNV - CNVpytor]
    E --> O[SAVANA - CN, purity, ploidy]
    F --> I[Annotation - ANNOVAR + filtering]
    G --> J[MGMT promoter status]
    G --> K[Methylation classification]
    H --> L[SV annotation - AnnotSV]
    I --> M[Report generation]
    J --> M
    K --> M
    N --> M
```

SNV calling is restricted to the NPHD panel BED; SV, CNV and methylation run on
the full BAM.

## ⚙️ Parameters

### Required parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input` | **Required.** Path to input BAM file(s). Can be:<br>• Single aligned BAM file: `/path/to/sample.bam`<br>• Directory with aligned BAMs: `/path/to/aligned_bams/`<br>• Directory with unaligned BAMs: `/path/to/unaligned_bams/` | `--input /data/sample.bam` |
| `--id` | **Required.** Unique sample identifier used for naming output files and reports. Should be alphanumeric with no spaces. | `--id SAMPLE001` |

### System-specific parameters

**These parameters must be configured for your specific system and installation paths in the nextflow.config file:**

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--ref` | Path to hg38 reference genome FASTA file | System-specific | `--ref /refs/hg38.fa` |
| `--annovarPath` | Path to ANNOVAR installation directory. | System-specific | `--annovarPath /tools/annovar/` |
| `--annovarDB` | Path to ANNOVAR database directory (humandb/). | System-specific | `--annovarDB /tools/annovar/humandb/` |
| `--annotsvAnnot` | Path to AnnotSV annotations directory. | System-specific | `--annotsvAnnot /tools/AnnotSV/Annotations_Human/` |
| `--annotations` | Path to annotation file for IGV reports (refGene.txt.gz). | `${projectDir}/data/refGene.txt.gz` | `--annotations /refs/refGene.txt.gz` |


### Optional parameters

#### Output parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--outDir` | Output directory for all pipeline results. Will be created if it doesn't exist. | `output` | `--outDir /results/analysis` |
| `--patient` | Patient name for reports. If not specified, uses the `--id` value. | `null` (uses `--id`) | `--patient "John Doe"` |

#### Resource parameters

CPU and memory are set per process label in `nextflow.config`, and each tool is
given `task.cpus`, so an allocation and a tool's thread count can never
disagree.

Every request is capped by two parameters, so the pipeline runs on smaller
machines without editing any config:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--maxCpus` | Ceiling on CPUs for any single process | `64` |
| `--maxMemory` | Ceiling on memory in GB for any single process | `96` |

```bash
# e.g. an 8-core workstation
nextflow run main.nf ... --maxCpus 8 --maxMemory 32 -profile local,singularity
```

Without a cap, an oversized request is a hard error rather than a slow run
(`Process requirement exceeds available CPUs`). For finer control, override
individual labels with your own `-c custom.config`.

#### Variant calling tool parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--clair3Model` | Clair3 model directory name. By default it is **auto-detected** from the `basecall_model=` field in the BAM's `@RG` header and validated against the models shipped in the container, so a mismatched model cannot be used silently. Set only to override. | `null` (auto) | `--clair3Model r1041_e82_400bps_sup_v500` |
| `--severusVntr` | Optional VNTR BED for Severus. Recommended: `vntrs/human_GRCh38_no_alt_analysis_set.trf.bed` from the [Severus repo](https://github.com/KolmogorovLab/Severus). | `null` | `--severusVntr /refs/vntr.bed` |
| `--severusPON` | Optional panel of normals for Severus. Without one, tumour-only Severus has **no somatic filter**, so supplying `pon/PoN_1000G_hg38.tsv.gz` is strongly advised. With a PON, the published VCF is the somatic set. | `null` | `--severusPON /refs/PoN_1000G_hg38.tsv.gz` |
| `--savanaG1000` | 1000G het-SNP set used by SAVANA for B-allele frequency, purity and ploidy. | `1000g_hg38` | `--savanaG1000 1000g_hg38` |

#### Annotation parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--annovarProtocol` | Comma-separated ANNOVAR databases. | `refGene,cytoBand,clinvar_20240917,avsnp151,1000g2015aug_eur` |
| `--annovarOperation` | ANNOVAR operation codes. **Must have the same number of entries as `--annovarProtocol`** - the pipeline checks this at startup and fails early if not. | `gx,r,f,f,f` |

To add `cosmic70` (needs a COSMIC licence), `dbnsfp47a` or `allofus`, append to
**both** lists, keeping them the same length.

#### Container and system parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--containerBindPaths` | Comma-separated host paths to bind into containers. Singularity only auto-mounts the work directory, so reference genome, ANNOVAR and AnnotSV paths outside it **must** be listed here or they will be invisible inside the container. | `/b06x-isilon,/omics` | `--containerBindPaths /data,/refs` |
| `--seq` | Sequencer platform identifier. Set to `false` to auto-detect from the BAM header. | `P2S` | `--seq F` |

### Profile-specific parameters

The pipeline supports different compute infrastructure profiles:

Resources are set with native Nextflow `cpus`/`memory`/`queue` directives rather
than hand-written `clusterOptions`, which previously emitted duplicate `-n` and
`-q` flags on every submission.

#### LSF profile (`-profile lsf`)
- **Executor:** LSF, queue `verylong`, `perJobMemLimit` enabled
- **Default:** 2 CPUs / 8 GB

| Label | Processes | CPUs | Memory |
|-------|-----------|------|--------|
| `rapid_cns` | light steps: MGMT, AnnotSV, ANNOVAR, IGV, reporting | 8 | 32 GB |
| `bamprep` | sort + index the full flow-cell BAM | 32 | 96 GB |
| `heavy` | Sniffles2, CNVpytor, methylation classifier | 32 | 64 GB |
| `mods` | modkit methylation calling | 32 | 64 GB |
| `clair3` | SNV calling | 64 | 96 GB |
| `severus` | Severus SV calling | 32 | 64 GB |
| `savana` | SAVANA SV and copy number | 32 | 64 GB |

Each tool is invoked with `task.cpus`, so a tool's thread count always matches
what the scheduler allocated.

#### SLURM profile (`-profile slurm`)
- **Executor:** SLURM, queue `batch`; same per-label CPU/memory table as LSF.

#### Local profile (`-profile local`)
- **Executor:** Local, 1 CPU / 4 GB per process. Suitable only for small tests.

## 📦 Containers

Every container is version-pinned; none use `:latest`.

| Label | Image | Provides |
|-------|-------|----------|
| `rapid_cns` | `areebapatel/rapid_cns:3.0.2` | samtools, bedtools, vcftools, dorado, mosdepth, CNVpytor, Sniffles2, methylartist, igv-reports, AnnotSV (code), R stack |
| `mods` | `quay.io/biocontainers/ont-modkit:0.6.4--h7f49ad2_0` | modkit |
| `clair3` | `hkubal/clair3:v2.0.2` | Clair3 and its bundled ONT models |
| `severus` | `quay.io/biocontainers/severus:1.7--pyhdfd78af_0` | Severus |
| `savana` | `quay.io/biocontainers/savana:1.3.8--pyhdfd78af_0` | SAVANA |

## 📤 Output

```
output/
├── bam/                    # merged, sorted, indexed BAM + panel subset
├── snv/                    # Clair3 VCF, PASS VCF, ANNOVAR multianno, filtered report
├── sv/
│   ├── <id>.sniffles2.vcf              # Sniffles2 calls
│   ├── <id>_sniffles_annotsv.tsv       # AnnotSV annotation
│   ├── severus/                        # Severus calls + <id>.severus.vcf
│   ├── <id>_severus_annotsv.tsv        # AnnotSV annotation
│   └── savana/                         # SAVANA tumour-only breakpoints
├── cnv/
│   ├── <id>.cnvpytor.calls.*.tsv       # CNVpytor calls at 1k/10k/100k bins
│   ├── <id>_cnvpytor_100k.png/.pdf     # genome-wide plot
│   ├── <id>.annotation.1000.xlsx       # gene-level CNV table
│   └── savana/                         # SAVANA copy number + purity/ploidy fit
├── mods/                   # modkit bedMethyl
├── mgmt/                   # MGMT coverage, status, methylartist plot
├── methylation_classification/         # votes, rf_details, calibrated scores
├── mnpflex/                # <id>.MNPFlex.input.bed  (upload to app.epignostix.com)
├── coverage/               # mosdepth summaries
├── reports/                # IGV report
├── report/                 # <id>_Rapid-CNS2_report_lite.html (no IGV)
│                           # <id>_Rapid-CNS2_report_full.html (IGV embedded)
└── pipeline_info/          # timeline, trace, execution report
```

> **Note on the report:** the final HTML currently covers coverage, SNVs,
> the IGV report, the CNVpytor copy-number plot, methylation classification and
> MGMT. The Sniffles2/Severus/AnnotSV results, the SAVANA copy-number and
> purity/ploidy output and the gene-level CNV table are written to disk but are
> **not yet included in the rendered report**.

### MNP-Flex integration

**MNP-Flex** is a methylation-based tumor classifier that provides detailed molecular classification of CNS tumours. The Rapid-CNS² pipeline can prepare the necessary input files for MNP-Flex analysis.

#### What is MNP-Flex?

MNP-Flex is a methylation classifier compatible with the latest version of the Heidelberg CNS tumour methylation classifier. It can classify CNS tumours into 184 subclasses according to the 2021 WHO classification.

#### How to use MNP-Flex

1. **Enable MNP-Flex in the pipeline:**
   ```bash
   nextflow run main.nf \
       --input /data/sample.bam \
       --id SAMPLE001 \
       --mnpFlex true \
       -profile lsf,singularity
   ```

2. **Locate the output files:**
   The pipeline creates MNP-Flex compatible files in:
   ```
   ${outDir}/mnpflex/
   └── ${id}.MNPFlex.input.bed
   ```

3. **Upload to MNP-Flex:**
   - Visit [app.epignostix.com](https://app.epignostix.com)
   - Upload the `.bed` file generated by the pipeline
   - Submit for analysis

4. **Interpret results:**
   - MNP-Flex will provide detailed classification results
   - Results include confidence scores
   - Reports are compatible with clinical interpretation guidelines

#### Output files

- **`${id}.MNPFlex.input.bed`:** Methylation data in MNP-Flex compatible format

**Note:** MNP-Flex analysis is performed externally at [app.epignostix.com](https://app.epignostix.com). The pipeline only prepares the input file; it does not upload anything.

### Getting Help

- Check the Nextflow documentation: https://www.nextflow.io/docs/
- Review the pipeline logs in the `work/` directory
- Check the `pipeline_info/` directory for execution reports
- Run with `--help` for available options

## 📝 Changelog

### 3.0.0

A major update that modernises the toolchain and fixes a number of bugs.

- **New callers:** Clair3 replaces DeepVariant for SNVs, Severus joins Sniffles2 for structural variants,
  and SAVANA adds copy number with tumour purity and ploidy estimates.
- **Compatible with current Dorado output.** All tools are updated and pinned,
  including modkit, Sniffles2, mosdepth, AnnotSV and dorado itself. The Clair3
  model is detected automatically from the BAM's basecaller.
- **Bug fixes** across variant filtering, methylation classification, MGMT
  prediction and report rendering. The report is now HTML only; the PDF variant
  has been dropped.
- **Easier to run elsewhere:** site paths are validated at startup,
  `--maxCpus`/`--maxMemory` cap resource requests so the pipeline runs on
  smaller machines, and every container image is version-pinned.

Note that the structural variant, SAVANA copy number and gene-level CNV outputs
are written to disk but are not yet included in the rendered report.

## 📚 Citation

If you use this pipeline, please cite:

Patel, A., Göbel, K., Ille, S. et al. Prospective, multicenter validation of a
platform for rapid molecular profiling of central nervous system tumors.
*Nature Medicine* **31**, 1567–1577 (2025).
[doi:10.1038/s41591-025-03562-5](https://www.nature.com/articles/s41591-025-03562-5)

```bibtex
@article{patel2025rapidcns2,
  title   = {Prospective, multicenter validation of a platform for rapid
             molecular profiling of central nervous system tumors},
  author  = {Patel, Areeba and G{\"o}bel, Kirsten and Ille, Sebastian and others},
  journal = {Nature Medicine},
  volume  = {31},
  pages   = {1567--1577},
  year    = {2025},
  doi     = {10.1038/s41591-025-03562-5}
}
```

## ⚖️ License

This project is licensed under the [Apache License 2.0](LICENSE).