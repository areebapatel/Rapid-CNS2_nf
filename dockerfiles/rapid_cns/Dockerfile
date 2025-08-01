# Rapid-CNS Dockerfile
FROM ubuntu:22.04

LABEL version="3.0.0"
LABEL maintainer="Areeba Patel <a.patel@dkfz.de>"
LABEL date="2025-07-24"

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies and Python packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        wget \
        curl \
        git \
        python3 \
        python3-dev \
        python3-pip \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        ca-certificates \
        zlib1g-dev \
        locales \
        pandoc \
        && pip3 install --no-cache-dir pod5 igv-reports==1.0.5 methylartist sniffles requests pandas gnureadline pathlib pysam numpy scipy matplotlib h5py xlsxwriter seaborn openpyxl \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R from CRAN (latest version)
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    wget \
    ca-certificates \
    gnupg \
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/r-project.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | tee -a /etc/apt/sources.list.d/r-project.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set locale
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

# Initialize conda and add channels
RUN conda init bash && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create conda environment and install all tools from bioconda
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda install -y -c conda-forge -c bioconda python=3.12 samtools bedtools vcftools && \
    conda clean -afy

# Install dorado
RUN curl -L "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.0.2-linux-x64.tar.gz" -o /tmp/dorado-1.0.2-linux-x64.tar.gz \
    && tar -xzf /tmp/dorado-1.0.2-linux-x64.tar.gz -C /opt \
    && ln -s /opt/dorado-1.0.2-linux-x64/bin/dorado /usr/local/bin/dorado \
    && rm -rf /tmp/dorado*

ENV LD_LIBRARY_PATH="/opt/dorado-1.0.2-linux-x64/lib"

# Install R packages with dependencies
RUN R -e "install.packages(c( \
    'optparse','ranger','matrixStats','data.table','glmnet', \
    'rmarkdown','knitr','BiocManager' \
  ), repos='https://cloud.r-project.org', INSTALL_opts='--no-test-load')"

# Install system dependencies for devtools and kableExtra
RUN apt-get update && apt-get install -y \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install kableExtra from GitHub using devtools
RUN R -e "install.packages('devtools', repos='https://cloud.r-project.org')" && \
    R -e "devtools::install_github('haozhu233/kableExtra', dependencies=TRUE)"

 
RUN R -e "BiocManager::install('GenomicRanges', ask = FALSE, update = FALSE)"

# Test if all required R packages are installed
RUN R -e "required_packages <- c('optparse', 'rmarkdown', 'kableExtra', 'knitr', 'ranger', 'matrixStats', 'data.table', 'glmnet', 'GenomicRanges'); \
    missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]; \
    if (length(missing_packages) > 0) { \
        stop(paste('Missing R packages:', paste(missing_packages, collapse = ', '))); \
    } else { \
        cat('All required R packages are installed successfully\n'); \
    }"

# Add mosdepth (latest release)
RUN wget -O /usr/local/bin/mosdepth https://github.com/brentp/mosdepth/releases/download/v0.3.3/mosdepth && \
    chmod +x /usr/local/bin/mosdepth

# Install AnnotSV
RUN git clone https://github.com/lgmgeo/AnnotSV.git /opt/AnnotSV && \
    cd /opt/AnnotSV && \
    make PREFIX=. install 

# Add AnnotSV to PATH
ENV PATH="/opt/AnnotSV/bin:$PATH"

# Install cnvpytor with user flag and download data
RUN pip3 install --user git+https://github.com/abyzovlab/CNVpytor.git && \
    export PATH="/root/.local/bin:$PATH" 

# Set matplotlib config directory to avoid warnings
ENV MPLCONFIGDIR=/tmp/matplotlib

# Add user bin directory to PATH for pip --user installations
ENV PATH="/root/.local/bin:$PATH"

# Set working directory
WORKDIR /data

# Default command
CMD ["bash"] 