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
ARG CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/__linux__/jammy/2026-07-01

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8 \
    CONDA_DIR=/opt/conda \
    MAMBA_ROOT_PREFIX=/opt/mamba

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ca-certificates curl git gnupg locales pandoc wget \
        bzip2 \
        libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev \
        libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
        libpng-dev libtiff5-dev libjpeg-dev \
        libuv1-dev libmagick++-dev \
    && locale-gen en_US.UTF-8 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
      | gpg --dearmor -o /usr/share/keyrings/r-project.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
       > /etc/apt/sources.list.d/r-project.list \
    && apt-get update && apt-get install -y --no-install-recommends r-base r-base-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

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

RUN pip install --no-cache-dir \
        pandas==${PANDAS_VERSION} \
        openpyxl==${OPENPYXL_VERSION} \
        sniffles==${SNIFFLES_VERSION} \
        methylartist==${METHYLARTIST_VERSION} \
        igv-reports==${IGVREPORTS_VERSION} \
        cnvpytor==${CNVPYTOR_VERSION}

RUN curl -fsSL "https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-linux-x64.tar.gz" \
      -o /tmp/dorado.tar.gz \
    && tar -xzf /tmp/dorado.tar.gz -C /opt \
    && ln -sf /opt/dorado-${DORADO_VERSION}-linux-x64/bin/dorado /usr/local/bin/dorado \
    && rm -f /tmp/dorado.tar.gz

RUN CNVDIR=$(python3 -c "import cnvpytor,os;print(os.path.join(os.path.dirname(cnvpytor.__file__),'data'))") \
    && for f in gc_hg19.pytor mask_hg19.pytor gc_hg38.pytor mask_hg38.pytor \
                gc_chm13v2.0.pytor gc_chm13v1.1.pytor gc_kn99.pytor; do \
         curl -fsSL -o "$CNVDIR/$f" \
           "https://github.com/abyzovlab/CNVpytor/raw/master/cnvpytor/data/$f"; \
       done \
    && python3 -c "from cnvpytor.genome import Genome; assert Genome.check_resources()"

RUN wget -qO /usr/local/bin/mosdepth \
      "https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth" \
    && chmod +x /usr/local/bin/mosdepth

RUN git clone --depth 1 --branch v${ANNOTSV_VERSION} https://github.com/lgmgeo/AnnotSV.git /opt/AnnotSV \
    && cd /opt/AnnotSV \
    && make PREFIX=. install
ENV ANNOTSV=/opt/AnnotSV
ENV PATH=/opt/AnnotSV/bin:$PATH

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

RUN mkdir -p /opt/logos
COPY logos/rapid_cns2_logo.png /opt/logos/
COPY logos/ukhd_logo.png /opt/logos/

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
