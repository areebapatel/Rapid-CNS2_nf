FROM ubuntu:22.04

LABEL MAINTAINER="Areeba Patel <a.patel@kitz-heidelberg.de>"

RUN apt-get update -y && apt-get install -y python3 python3-dev python3-pip git && pip install requests pandas gnureadline pathlib pysam numpy scipy matplotlib h5py xlsxwriter seaborn openpyxl

RUN git clone https://github.com/abyzovlab/CNVpytor.git && cd CNVpytor && pip install .

COPY /scripts /scripts
COPY /genes.bed /genes.bed

RUN chmod 777 -R /scripts /genes.bed

