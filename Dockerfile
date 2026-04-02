FROM bioconductor/bioconductor_docker:RELEASE_3_22

RUN apt-get update && apt-get install -y \
    cmake \
    libglpk-dev \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

WORKDIR /app

ARG REPO_URL
RUN git clone ${REPO_URL} .

RUN R -e "install.packages('arrow', repos='https://packagemanager.posit.co/cran/latest')"
RUN R -e "install.packages('renv'); renv::restore()"