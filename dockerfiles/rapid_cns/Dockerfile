# Rapid-CNS2 container
# =============================================================================
# Self-contained: build with just this file plus the repo's logos/ directory.
#
# BUILD (run from the REPOSITORY ROOT, not from this directory):
#
#   docker buildx build \
#       --platform linux/amd64 \
#       -t areebapatel/rapid_cns:3.0.2 \
#       -f dockerfiles/rapid_cns/Dockerfile \
#       --push .
#
#   # or build locally first, then push
#   docker buildx build --platform linux/amd64 --load \
#       -t areebapatel/rapid_cns:3.0.2 -f dockerfiles/rapid_cns/Dockerfile .
#   docker push areebapatel/rapid_cns:3.0.2
#
# PLATFORM
#   --platform linux/amd64 is REQUIRED and the FROM line pins it too. On Apple
#   Silicon a default build would produce an arm64 image that cannot run on the
#   x86_64 HPC. Cross-building under QEMU works but is slow (allow 1-2 h; a
#   native x86_64 machine takes ~30-45 min).
#
#   The image is deliberately amd64-only, not multi-arch: dorado's Linux build
#   and mosdepth's release binary are published for x86_64, and every HPC target
#   here is x86_64. It runs on an Apple Silicon Mac through Docker's emulation.
#
# PULLING ON THE HPC (no Docker needed there)
#   singularity pull docker://areebapatel/rapid_cns:3.0.2
#   # or let Nextflow do it:
#   nextflow run main.nf ... --containerEngine singularity
#
# VERSIONS
#   Everything below is pinned so a rebuild is reproducible. Bump an ARG
#   deliberately and publish a new tag; never rely on :latest.
#
# SCOPE
#   Clair3, Severus and SAVANA each run in their own pinned container and are
#   deliberately NOT installed here.
# =============================================================================
FROM --platform=linux/amd64 ubuntu:22.04

LABEL version="3.0.2"
LABEL maintainer="Areeba Patel <a.patel@kitz-heidelberg.de>"
LABEL description="Rapid-CNS2: alignment, coverage, CNV, SV, methylation and reporting"

ARG DORADO_VERSION=2.1.0
ARG MOSDEPTH_VERSION=0.3.14
ARG ANNOTSV_VERSION=3.5.10
ARG SNIFFLES_VERSION=2.8.0
ARG METHYLARTIST_VERSION=1.5.4
ARG IGVREPORTS_VERSION=1.16.3
ARG CNVPYTOR_VERSION=1.3.2
ARG PANDAS_VERSION=2.2.3
ARG OPENPYXL_VERSION=3.1.5
ARG SAMTOOLS_VERSION=1.21
ARG HTSLIB_VERSION=1.21
ARG BEDTOOLS_VERSION=2.31.1
ARG VCFTOOLS_VERSION=0.1.16
ARG BCFTOOLS_VERSION=1.21
ARG PYTHON_VERSION=3.12
ARG MICROMAMBA_VERSION=2.8.1
# Dated CRAN snapshot pins every R package without listing each version
ARG CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/__linux__/jammy/2026-07-01

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8 \
    CONDA_DIR=/opt/conda \
    MAMBA_ROOT_PREFIX=/opt/mamba

# --- system dependencies ----------------------------------------------------
# The R graphics stack (fontconfig/freetype/harfbuzz/fribidi/png/tiff/jpeg) is
# required by systemfonts -> svglite -> kableExtra; omitting it makes those
# packages fail to configure.
# libuv1-dev: required by R package 'fs' (dependency of sass/bslib/rmarkdown)
# libmagick++-dev: required by R package 'magick' (optional but pulled by kableExtra deps)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ca-certificates curl git gnupg locales pandoc wget \
        bzip2 \
        libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev \
        libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
        libpng-dev libtiff5-dev libjpeg-dev \
        libuv1-dev libmagick++-dev \
    && locale-gen en_US.UTF-8 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# --- R ----------------------------------------------------------------------
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
      | gpg --dearmor -o /usr/share/keyrings/r-project.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
       > /etc/apt/sources.list.d/r-project.list \
    && apt-get update && apt-get install -y --no-install-recommends r-base r-base-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# --- micromamba: ONE python stack -------------------------------------------
# The previous image installed into both the system python and conda, leaving
# duplicate copies of pysam/numpy/pandas.
# NOTE: `pip` is listed explicitly. conda-forge's python does NOT depend on pip,
# and no python3-pip is installed via apt, so without it the pip layer below
# would have no pip to run.
ENV PATH=$CONDA_DIR/bin:$PATH
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/${MICROMAMBA_VERSION} \
      | tar -xj -C /tmp bin/micromamba \
    && install -m 0755 /tmp/bin/micromamba /usr/local/bin/micromamba \
    && rm -rf /tmp/bin \
    && micromamba create -y -p $CONDA_DIR \
        -c conda-forge \
        -c bioconda \
        --strict-channel-priority \
        python=${PYTHON_VERSION} \
        pip \
        samtools=${SAMTOOLS_VERSION} \
        htslib=${HTSLIB_VERSION} \
        bedtools=${BEDTOOLS_VERSION} \
        vcftools=${VCFTOOLS_VERSION} \
        bcftools=${BCFTOOLS_VERSION} \
    && micromamba clean -a -y

# --- python tools -----------------------------------------------------------
# pandas: scr/annotate.py ; openpyxl: its .to_excel call.
# pysam/numpy/scipy/matplotlib arrive as cnvpytor/methylartist dependencies.
# Dropped vs the previous image: pod5 (basecalling is external), pathlib (a
# stdlib backport that breaks python 3), gnureadline, requests, seaborn, h5py.
RUN pip install --no-cache-dir \
        pandas==${PANDAS_VERSION} \
        openpyxl==${OPENPYXL_VERSION} \
        sniffles==${SNIFFLES_VERSION} \
        methylartist==${METHYLARTIST_VERSION} \
        igv-reports==${IGVREPORTS_VERSION} \
        cnvpytor==${CNVPYTOR_VERSION}

# --- dorado (prepareBam uses `dorado aligner`) ------------------------------
RUN curl -fsSL "https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-linux-x64.tar.gz" \
      -o /tmp/dorado.tar.gz \
    && tar -xzf /tmp/dorado.tar.gz -C /opt \
    && ln -sf /opt/dorado-${DORADO_VERSION}-linux-x64/bin/dorado /usr/local/bin/dorado \
    && rm -f /tmp/dorado.tar.gz
# No LD_LIBRARY_PATH: the dorado binary has RUNPATH=$ORIGIN/../lib and resolves
# all 73 bundled libs on its own. Exporting that directory globally would put
# dorado's bundled libstdc++ and CUDA libs ahead of the system ones for every
# other binary in the image (R, python, samtools).

# --- CNVpytor reference resources -------------------------------------------
# `cnvpytor -download` is broken in 1.3.2, and CNVpytor requires the files for
# every reference genome it knows about, not just the one in use.
RUN CNVDIR=$(python3 -c "import cnvpytor,os;print(os.path.join(os.path.dirname(cnvpytor.__file__),'data'))") \
    && for f in gc_hg19.pytor mask_hg19.pytor gc_hg38.pytor mask_hg38.pytor \
                gc_chm13v2.0.pytor gc_chm13v1.1.pytor gc_kn99.pytor; do \
         curl -fsSL -o "$CNVDIR/$f" \
           "https://github.com/abyzovlab/CNVpytor/raw/master/cnvpytor/data/$f"; \
       done \
    && python3 -c "from cnvpytor.genome import Genome; assert Genome.check_resources()"

# --- mosdepth ---------------------------------------------------------------
RUN wget -qO /usr/local/bin/mosdepth \
      "https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth" \
    && chmod +x /usr/local/bin/mosdepth

# --- AnnotSV (code only) ----------------------------------------------------
# The human annotations are deliberately NOT baked in. They are ~6.6 GB, they
# are released independently of the AnnotSV code, and `make
# install-human-annotation` fetches whatever is current at build time - which
# would make this otherwise-pinned image non-reproducible.
# Install them once on the host and point --annotsvAnnot at the directory
# containing Annotations_Human (see README step 5), exactly as is already done
# for the reference genome and the ANNOVAR databases.
RUN git clone --depth 1 --branch v${ANNOTSV_VERSION} https://github.com/lgmgeo/AnnotSV.git /opt/AnnotSV \
    && cd /opt/AnnotSV \
    && make PREFIX=. install
ENV ANNOTSV=/opt/AnnotSV
ENV PATH=/opt/AnnotSV/bin:$PATH

# --- R packages -------------------------------------------------------------
# Pinned collectively by the dated CRAN snapshot. kableExtra now comes from
# CRAN, so devtools and its large dependency tree are gone.
# glmnet is required for the serialised model to load even though no glmnet
# function is called directly - do not remove it.
RUN R -e "req <- c( \
                'optparse','ranger','matrixStats','data.table','glmnet', \
                'rmarkdown','knitr','kableExtra','ggplot2','reshape2', \
                'base64enc','xfun','BiocManager' \
            ); \
            options(repos = c(CRAN = '${CRAN_SNAPSHOT}')); \
            install.packages(req, Ncpus = 1, dependencies = TRUE); \
            miss <- req[!sapply(req, requireNamespace, quietly = TRUE)]; \
            if (length(miss)) { \
                message('Retrying missing packages from CRAN snapshot: ', paste(miss, collapse=', ')); \
                install.packages(miss, Ncpus = 1, dependencies = TRUE); \
            }; \
            miss2 <- req[!sapply(req, requireNamespace, quietly = TRUE)]; \
            if (length(miss2)) { \
                options(repos = c(CRAN = 'https://cloud.r-project.org')); \
                message('Falling back to live CRAN for: ', paste(miss2, collapse=', ')); \
                install.packages(miss2, Ncpus = 1, dependencies = TRUE); \
            }; \
            miss3 <- req[!sapply(req, requireNamespace, quietly = TRUE)]; \
            if (length(miss3)) stop('Failed to install required R packages: ', paste(miss3, collapse=', '))" \
 && R -e "options(repos = c(CRAN = '${CRAN_SNAPSHOT}')); \
     bioc <- c('GenomeInfoDb','GenomicRanges'); \
     BiocManager::install(bioc, ask = FALSE, update = FALSE); \
     for (p in bioc) library(p, character.only = TRUE)" \
 && R -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); \
     post <- c('rmarkdown','kableExtra'); \
     miss <- post[!sapply(post, requireNamespace, quietly = TRUE)]; \
     if (length(miss)) { \
        message('Post-Bioc reinstall for: ', paste(miss, collapse=', ')); \
        install.packages(miss, Ncpus = 1, dependencies = TRUE); \
     }; \
     miss2 <- post[!sapply(post, requireNamespace, quietly = TRUE)]; \
     if (length(miss2)) stop('Post-Bioc package check failed: ', paste(miss2, collapse=', '))"

# --- logos used by the report templates -------------------------------------
RUN mkdir -p /opt/logos
COPY logos/rapid_cns2_logo.png /opt/logos/
COPY logos/ukhd_logo.png /opt/logos/

# --- verify every binary and package the pipeline actually invokes ----------
RUN set -e; \
    for t in samtools bedtools vcftools bcftools bgzip tabix mosdepth dorado \
             sniffles methylartist create_report cnvpytor AnnotSV Rscript python3; do \
        command -v "$t" >/dev/null || { echo "MISSING BINARY: $t" >&2; exit 1; }; \
    done; \
    python3 -c "import pandas, openpyxl, cnvpytor"; \
    dorado --version >/dev/null 2>&1 || { echo "dorado cannot resolve its libraries" >&2; exit 1; }; \
    R -e "pkgs <- c('optparse','rmarkdown','kableExtra','knitr','ranger','ggplot2', \
                    'matrixStats','data.table','glmnet','GenomicRanges','GenomeInfoDb', \
                    'reshape2','base64enc','xfun'); \
          ok <- sapply(pkgs, function(p) \
                    isTRUE(suppressMessages(require(p, character.only = TRUE, quietly = TRUE)))); \
          if (any(!ok)) stop('R packages that fail to load: ', paste(pkgs[!ok], collapse=', ')); \
          cat('All required R packages present\n')"; \
    echo "=== container verification passed ==="

WORKDIR /data
CMD ["bash"]
