FROM bioconductor/bioconductor_docker:RELEASE_3_11

RUN R -e "install.packages('ggrepel')"
RUN R -e "BiocManager::install('snpStats')"
RUN R -e "devtools::install_github('mcgml/gwasglue')"
