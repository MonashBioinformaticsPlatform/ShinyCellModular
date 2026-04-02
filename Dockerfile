FROM bioconductor/bioconductor_docker:RELEASE_3_22

RUN apt-get update && apt-get install -y \
    cmake \
    libglpk-dev \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

ARG REPO_URL
RUN git clone ${REPO_URL} .

ENV RENV_PATHS_LIBRARY=/app/renv/library

RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"
RUN R -e "renv::install('arrow', repos='https://packagemanager.posit.co/cran/latest')"
RUN R -e "renv::restore(exclude='arrow')" 2>&1 | tee /renv_restore.log

RUN echo 'RENV_CONFIG_AUTOLOADER_ENABLED=FALSE' >> /usr/local/lib/R/etc/Renviron && \
    echo 'R_PROFILE_USER=/app/.Rprofile' >> /usr/local/lib/R/etc/Renviron

RUN printf '.libPaths("/app/renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu")\nif (file.exists("/app/functions/prepShinyCellModular.R")) source("/app/functions/prepShinyCellModular.R")\n' > /app/.Rprofile

RUN chmod -R 775 /app/renv/library