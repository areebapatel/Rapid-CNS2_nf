FROM r-base-ubuntu:20.04

LABEL maintainer "Kirsten Göbel kirsten.goebel@med.uni-heidelberg.de"

RUN apt-get update -y && apt-get install -y bedtools

COPY runMGMT.sh /runMGMT.sh
COPY mgmt_pred.R /mgmt_pred.R
COPY installPackages.R /installPackages.R

RUN chmod 755 /mgmt_pred.R && chmod 755 /runMGMT.sh && Rscript /installPackages.R

CMD /runMGMT.sh