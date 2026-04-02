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

# disable renv autoloader and set paths before any R execution
RUN echo 'RENV_CONFIG_AUTOLOADER_ENABLED=FALSE' >> /usr/local/lib/R/etc/Renviron

# install arrow as system package to avoid compile issues and path complications
RUN R -e "install.packages('arrow', \
    repos='https://packagemanager.posit.co/cran/latest', \
    lib='/usr/local/lib/R/site-library')"

# restore all packages into explicit path, autoloader disabled
RUN R -e " \
    options(renv.config.autoloader.enabled=FALSE); \
    renv::restore( \
        library='/app/renv/library', \
        exclude='arrow', \
        prompt=FALSE \
    )" 2>&1 | tee /renv_restore.log

# resolve the exact library path and write to Renviron + .Rprofile
RUN R -e " \
    base <- '/app/renv/library'; \
    l1 <- list.files(base, full.names=TRUE)[1]; \
    l2 <- list.files(l1, full.names=TRUE)[1]; \
    lib <- list.files(l2, full.names=TRUE)[1]; \
    stopifnot(dir.exists(lib)); \
    cat('Resolved lib path:', lib, '\n'); \
    write(paste0('R_LIBS_USER=', lib), \
        '/usr/local/lib/R/etc/Renviron', append=TRUE); \
    write(paste0('R_PROFILE_USER=/app/.Rprofile'), \
        '/usr/local/lib/R/etc/Renviron', append=TRUE); \
    writeLines(c( \
        'if (file.exists(\"/app/functions/prepShinyCellModular.R\")) source(\"/app/functions/prepShinyCellModular.R\")' \
    ), '/app/.Rprofile')"

RUN chmod -R 775 /app/renv/library