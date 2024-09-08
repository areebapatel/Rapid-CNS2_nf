FROM ubuntu:22.04

RUN apt update -y && apt install -y git python3 python3-dev python3-pip

RUN git clone https://github.com/a-slide/pycoQC.git && cd pycoQC && chmod u+x setup.py && pip3 install -e ./

COPY runQualityCheck.sh /runQualityCheck.sh

CMD /runQualityCheck.sh