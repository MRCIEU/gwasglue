FROM bioconductor/bioconductor_docker

RUN R -e "devtools::install_github('mrcieu/gwasglue')"