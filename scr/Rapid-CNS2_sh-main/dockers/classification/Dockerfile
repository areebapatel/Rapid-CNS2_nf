FROM r-base-ubuntu:20.04

LABEL MAINTAINER="Kirsten Göbel <kirsten.goebel@med.uni-heidelberg.de>"

RUN apt-get update -y && apt-get upgrade -y 

COPY methylation_classification.R /methylclass.R
COPY runClassification.sh /runClassification.sh
COPY installPackages.R /installPackages.R

RUN chmod 755 /methylclass.R /runClassification.sh /installPackages.R && Rscript /installPackages.R

CMD /runClassification.sh
