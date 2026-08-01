<div align="left">
<img width="342" alt="Rapid-CNS2" src="https://github.com/user-attachments/assets/cda166c2-664f-4286-951a-309b111c1132">

<h1 style="display: inline-block;">Rapid-CNS<sup>2</sup> workflow</h1>
</div>

[![Paper](https://img.shields.io/badge/Nature%20Medicine-2025-b31b1b.svg)](https://www.nature.com/articles/s41591-025-03562-5)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/Nextflow-23.10%2B-brightgreen)](https://www.nextflow.io/)

## 🧬 Overview

Molecular profiling of central nervous system (CNS) tumours from ONT adaptive
sampling data. Implemented in Nextflow, so it runs unchanged on a workstation,
an HPC cluster or the cloud.

**Contents:**
[Features](#-features) ·
[Requirements](#-requirements) ·
[Quick start](#-quick-start) ·
[Input requirements](#-input-requirements) ·
[Pipeline structure](#-pipeline-structure) ·
[Parameters](#-parameters) ·
[Containers](#-containers) ·
[Output](#-output) ·
[Help](#-help) ·
[Changelog](#-changelog) ·
[Citation](#-citation) ·
[License](#-license)

## ✨ Features

- **Aligned or unaligned BAM input**, with alignment detected from the data
- **SNVs with Clair3**, the model auto-detected from the BAM's basecaller
- **Structural variants from Sniffles2 and Severus**, annotated with AnnotSV, and
  screened against a curated list of recurrent CNS fusions
- **Copy number from CNVpytor and SAVANA**, the latter absolute and corrected for
  tumour purity and ploidy
- **Methylation**: Rapid-CNS² classifier and MGMT promoter status
- **HTML report** collecting all of the above
- **MNP-Flex input** for upload to [app.epignostix.com](https://app.epignostix.com)

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

`scr/basecall.sh` is a reference Dorado command with the settings this pipeline
expects; edit the paths before use.

### Input options

The pipeline accepts:
- **(Preferred) Single aligned BAM file:** Direct path to a single aligned and merged BAM file (e.g., `/path/to/sample.bam`)
- **Directory with aligned BAM files:** Path to directory containing multiple aligned BAM files (will be merged automatically)
- **Directory with unaligned BAM files:** Path to directory containing multiple unaligned BAM files (will be aligned and merged automatically)

## 🧭 Pipeline structure

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
    H --> Q[Fusion screen - curated CNS list]
    O --> Q
    O --> R[Absolute CN plot + purity/ploidy]
    I --> M[Report generation]
    J --> M
    K --> M
    N --> M
    Q --> M
    R --> M
```

SNV calling is restricted to the NPHD panel BED; SV, CNV and methylation run on
the full BAM.

## 🧰 Parameters

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
| `--minFusionLen` | Smallest same-chromosome DEL/DUP/INV that can rearrange two genes. | `10000` | `--minFusionLen 50000` |
| `--minFusionReads` | Minimum supporting reads for a fusion candidate. | `10` | `--minFusionReads 5` |
| `--minFusionMapq` | Minimum mapping quality, where the caller reports it. | `50` | `--minFusionMapq 60` |
| `--knownFusions` | Curated CNS fusion list. | `data/cns_fusions.tsv` | `--knownFusions my.tsv` |

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
| `--containerBindPaths` | Comma-separated host paths to bind into containers. Singularity only auto-mounts the work directory, so reference genome, ANNOVAR and AnnotSV paths outside it **must** be listed here or they will be invisible inside the container. | `""` | `--containerBindPaths /data,/refs` |

### Profile-specific parameters

Resources are set with native Nextflow `cpus`/`memory`/`queue` directives.

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
│   ├── savana/                         # SAVANA tumour-only breakpoints
│   └── fusions/                        # per-caller fusion screen (see below)
├── cnv/
│   ├── <id>.cnvpytor.calls.*.tsv       # CNVpytor calls at 1k/10k/100k bins
│   ├── <id>_cnvpytor_100k.png/.pdf     # read-depth plot
│   ├── <id>_savana_cnv.png             # absolute copy number, genes annotated
│   ├── <id>.annotation.1000.xlsx       # gene-level CNV table
│   └── savana/                         # SAVANA segmentation + purity/ploidy fit
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

> **Note on the report:** the final HTML covers coverage, SNVs, the IGV report,
> both copy-number plots with the purity/ploidy fit, methylation classification,
> MGMT, and reportable fusions. The full SV call sets, the AnnotSV annotations
> and the gene-level CNV table are written to disk only.

### Copy number

**CNVpytor** (`<id>_cnvpytor_100k.png`) is a read-depth profile in 100 kb bins,
normalised to the genome-wide mean. It assumes nothing about the tumour.

**SAVANA** (`<id>_savana_cnv.png`) fits *absolute* copy number from read depth
and B-allele frequency, using SV breakpoints as segment boundaries and correcting
for purity and ploidy - so CN 0 (homozygous deletion) and CN 1 (heterozygous
loss) are distinguishable. Genes from `data/genes.bed` are labelled at their copy
number.

The purity/ploidy fit (`cnv/savana/<id>_purity_ploidy.tsv`) also appears in the
report. It is **experimental and not validated for clinical use** - see the
[SAVANA repository](https://github.com/cortes-ciriano-lab/savana). Segmentation
assumes uniform coverage, which adaptive-sampling data only partly satisfies.
`savana cna` is the slowest step in the pipeline.

### Fusion screen

Each SV caller's output is screened for gene fusions and written to
`sv/fusions/`. Three filters keep the reportable set specific:

1. **Structural sanity** - BND/TRA, anything interchromosomal, or a DEL/DUP/INV
   of at least `--minFusionLen`. Without this a small deletion straddling a gene
   boundary counts as a fusion.
2. **Evidence** - `--minFusionReads` and `--minFusionMapq`, where the caller
   reports them.
3. **Curated list** - `data/cns_fusions.tsv`, matched on the gene pair *and* the
   expected structure, not the names alone. `*` matches any partner.

Breakpoints are assigned to gene bodies from refGene, since the panel BED is
exon-level and breakpoints fall in introns.

| File | Contents |
|------|----------|
| `<id>_<caller>_fusions_reportable.tsv` | curated-list matches, `entity-defining` (both partners named) or `potentially significant` (one named partner) |
| `<id>_<caller>_fusions_all.tsv` | everything passing the filters |
| `<id>_<caller>_EGFRvIII.txt` | targeted EGFRvIII check |

The `entity` column states the association a fusion carries, not a diagnosis:
`KIAA1549--BRAF` is most frequent in pilocytic astrocytoma but also occurs in
other low-grade gliomas.

### MNP-Flex

[MNP-Flex](https://app.epignostix.com) is a methylation classifier compatible
with the Heidelberg CNS classifier, covering 184 subclasses of the 2021 WHO
classification.

The pipeline prepares its input at `mnpflex/<id>.MNPFlex.input.bed` (on by
default; disable with `--mnpFlex false`). Upload it yourself at
[app.epignostix.com](https://app.epignostix.com), or - if you have an account
and set the credentials below - let the pipeline submit it and collect the
results for you.

#### Uploading through the API

Optional and off by default. Set your website login in the environment of the
shell that launches the run:

```bash
export EPIGNOSTIX_USER='you@institute.de'      # the app.epignostix.com login
export EPIGNOSTIX_PASSWORD='...'

nextflow run main.nf ... --mnpFlexUpload
```

Put these in your shell profile or a file you `source` - **not** in
`nextflow.config` or a `-c` file, which Nextflow copies into the task directory.
If either variable is unset the upload is skipped with a warning and the rest of
the run is unaffected.

Submitting the run as a batch job? A batch script is non-interactive and does not
read `.bashrc`, so either export the variables in the shell you submit *from*
(LSF and SLURM pass the submission environment through), or source your profile
at the top of the launcher:

```bash
[ -f "$HOME/.bashrc" ] && . "$HOME/.bashrc"
```

Sample metadata is sent alongside the file:

| Parameter | Default |
|-----------|---------|
| `--mnpFlexApi` | `https://app.epignostix.com/api` |
| `--mnpFlexWorkflowId` | `18` |
| `--mnpFlexTechnology` | `Nanopore sequencing- PromethION` |
| `--mnpFlexCoverage` | `>10X` |
| `--mnpFlexExtraction` | `Frozen` |
| `--mnpFlexSex` | `ND` |
| `--mnpFlexLocalisation`, `--mnpFlexDiagnosis` | unset |
| `--mnpFlexWait` | `600` (seconds to wait for results) |

The upload response is written to `mnpflex/<id>_mnpflex_upload.json`. The
pipeline then polls for the analysis (up to `--mnpFlexWait`, default 600 s) and
downloads the results:

| File | Contents |
|------|----------|
| `<id>_mnpflex_predictions.tsv` | top classes with scores, classifier version, QC line |
| `<id>_bundle_summary.json` | full classifier output |
| `<id>_qc_*.json`, `<id>_mgmt_region_plot.json` | QC and MGMT plot data |

The top predictions appear in the report. Descriptions of the methylation
classes are in [Cancer Cell (2025)](https://www.cell.com/cancer-cell/fulltext/S1535-6108(25)00495-7).
If the analysis has not finished within `--mnpFlexWait`, the run still succeeds
and the report omits the section; re-run `scr/mnpflex_results.py --sample-id <id>`
later to collect it.

## 🆘 Help

- `nextflow run main.nf --help` lists every parameter
- Task logs are in `work/`; run reports and the timeline in `pipeline_info/`
- [Nextflow documentation](https://www.nextflow.io/docs/)

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
- **Fusion screen and absolute copy number.** Structural variants are screened
  against a curated list of recurrent CNS fusions, and SAVANA contributes an
  absolute copy-number profile with tumour purity and ploidy. Both appear in the
  report; the full SV call sets and the gene-level CNV table are written to disk.

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

## 📄 License

This project is licensed under the [Apache License 2.0](LICENSE).