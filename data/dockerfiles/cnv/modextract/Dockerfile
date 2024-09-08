FROM ubuntu:22.04

LABEL maintainer "Areeba Patel a.patel@kitz-heidelberg.de"

RUN apt-get update -y && apt-get install -y gcc make libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev wget bzip2 curl libcurl4-openssl-dev bedtools git rsync autoconf libssl-dev && mkdir /applications

# HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 -P /applications && tar -vxjf /applications/htslib-1.16.tar.bz2 -C /applications && cd /applications/htslib-1.16 && ./configure --prefix=/ && make

# SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 -P /applications/ && tar -vxjf /applications/samtools-1.16.1.tar.bz2 -C /applications && cd /applications/samtools-1.16.1 && ./configure --prefix=/ && make

# BCFTOOLS
RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 -P /applications/ && tar -vxjf /applications/bcftools-1.16.tar.bz2 -C /applications && cd /applications/bcftools-1.16 && ./configure --prefix=/ && make

# LIFTOVER
RUN wget 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz' -O /applications/hg19ToHg38.over.chain.gz && rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/liftOver /applications

# MODBAM2BED
RUN cd /applications/ && git clone --recursive https://github.com/epi2me-labs/modbam2bed.git --branch v0.5.3 && cd /applications/modbam2bed && make modbam2bed

COPY runExtraction.sh /runExtraction.sh
RUN chmod 755 /runExtraction.sh

ENV PATH=$PATH:/applications/samtools-1.16.1:/applications/bcftools-1.16:/applications/modbam2bed/:/applications

CMD /runExtraction.sh