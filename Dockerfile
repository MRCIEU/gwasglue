FROM bioconductor/bioconductor_docker:RELEASE_3_11
COPY . /app
WORKDIR /app

RUN R -e "install.packages('ggrepel')"
RUN R -e "BiocManager::install('snpStats')"
RUN R -e "devtools::install_local('.')"
RUN R -e "devtools::install_local('jrs95/gassocplot')"
