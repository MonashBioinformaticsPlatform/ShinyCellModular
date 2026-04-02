FROM condaforge/mambaforge:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Australia/Melbourne

RUN apt-get update && apt-get install -y \
    tzdata \
    cmake \
    libglpk-dev \
    libhdf5-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN mamba config append channels conda-forge && \
    mamba config append channels bioconda && \
    mamba install -y \
    r-arrow r-data.table r-dplyr r-dt bioconductor-edger \
    r-ggdendro r-ggforce r-ggplot2 r-ggrepel r-gridextra \
    r-hdf5r r-jsonlite bioconductor-limma r-magrittr \
    r-patchwork r-plotly r-presto r-rcolorbrewer r-rcpp \
    r-rsconnect r-seurat r-seuratobject r-shiny \
    r-shinydashboard r-shinyjs r-shinywidgets r-sortable \
    r-tidyverse r-factominer r-grimport r-xml r-devtools \
    r-shinythemes \
    && mamba clean -afy

RUN Rscript -e "devtools::install_github('Simon-Leonard/FlexDotPlot')"
RUN Rscript -e "devtools::install_github('SGDDNB/ShinyCell')"

WORKDIR /app

ARG REPO_URL
RUN git clone ${REPO_URL} .

RUN echo 'R_PROFILE_USER=/app/.Rprofile' >> $(R RHOME)/etc/Renviron
RUN echo 'if (file.exists("/app/functions/prepShinyCellModular.R")) source("/app/functions/prepShinyCellModular.R")' > /app/.Rprofile