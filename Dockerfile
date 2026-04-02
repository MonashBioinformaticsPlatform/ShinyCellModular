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

RUN R -e "install.packages(c( \
    'arrow', 'data.table', 'DT', 'ggplot2', 'ggrepel', 'hdf5r', \
    'ggdendro', 'gridExtra', 'rsconnect', 'shinythemes', 'shinydashboard', \
    'tidyverse', 'sortable', 'plotly', 'RColorBrewer', 'ggforce', \
    'shiny', 'shinyhelper', 'shinyjs', 'shinyWidgets' \
    ), repos='https://packagemanager.posit.co/cran/latest')"

RUN R -e "BiocManager::install(c('limma', 'edgeR'), ask=FALSE)"

RUN R -e "devtools::install_github('SGDDNB/ShinyCell')"
RUN R -e "devtools::install_github('immunogenomics/presto')"
RUN R -e "devtools::install_github('jmw86069/FlexDotPlot')"
RUN R -e "install.packages('Seurat', repos='https://packagemanager.posit.co/cran/latest')"

RUN echo 'R_PROFILE_USER=/app/.Rprofile' >> /usr/local/lib/R/etc/Renviron
RUN echo 'if (file.exists("/app/functions/prepShinyCellModular.R")) source("/app/functions/prepShinyCellModular.R")' > /app/.Rprofile