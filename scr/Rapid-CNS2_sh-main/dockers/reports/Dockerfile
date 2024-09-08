FROM r-base-ubuntu:20.04

LABEL maintainer "Kirsten Göbel kirsten.goebel@med.uni-heidelberg.de"

RUN apt-get update -y && apt-get install -y pandoc nano texlive-latex-base texlive-latex-extra python3 python3-pip && pip install pandas openpyxl plotly kaleido

COPY rdata/ /rdata
COPY scripts/ /scripts
COPY header_logo.png /header_logo.png
COPY logo.jpg /logo.jpg
COPY hospital_en.jpeg /hospital_en.jpeg

RUN Rscript /scripts/installPackages.R

CMD bash /scripts/runReport.sh