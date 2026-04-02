FROM bioconductor/bioconductor_docker:RELEASE_3_19

WORKDIR /app

ARG REPO_URL
RUN git clone ${REPO_URL} .

RUN R -e "install.packages('renv'); renv::restore()"