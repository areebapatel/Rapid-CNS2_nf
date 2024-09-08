FROM ubuntu:22.04

LABEL MAINTAINER="Kirsten Göbel <kirsten.goebel@med.uni-heidelberg.de>"

RUN apt-get update -y && apt-get install -y gcc make libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev wget bzip2 curl libcurl4-openssl-dev bedtools git rsync autoconf libssl-dev python3 python3-pip && mkdir /applications && pip install pandas openpyxl

# HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 -P /applications && tar -vxjf /applications/htslib-1.16.tar.bz2 -C /applications && cd /applications/htslib-1.16 && ./configure --prefix=/ && make

# SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 -P /applications/ && tar -vxjf /applications/samtools-1.16.1.tar.bz2 -C /applications && cd /applications/samtools-1.16.1 && ./configure --prefix=/ && make

# BCFTOOLS
RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 -P /applications/ && tar -vxjf /applications/bcftools-1.16.tar.bz2 -C /applications && cd /applications/bcftools-1.16 && ./configure --prefix=/ && make

ENV PATH=$PATH:/applications/samtools-1.16.1/:/applications/bcftools-1.16/:/applications/htslib-1.16
COPY runSpecialPos.sh /runSpecialPos.sh
COPY merge_specialPos.py /merge_specialPos.py

RUN chmod 755 /runSpecialPos.sh &&  chmod 755 /merge_specialPos.py

CMD /runSpecialPos.sh