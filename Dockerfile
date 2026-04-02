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

# install renv and arrow to system library first
RUN R -e "install.packages(c('withr'), \
    repos='https://packagemanager.posit.co/cran/latest', \
    lib='/usr/local/lib/R/site-library'); \
    install.packages( \
    'https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.4.tar.gz', \
    repos=NULL, type='source', \
    lib='/usr/local/lib/R/site-library')"

RUN R -e " \
    options(renv.config.autoloader.enabled=FALSE); \
    renv::install('arrow', \
        library='/app/renv/library', \
        repos='https://packagemanager.posit.co/cran/latest'); \
    withr::with_libpaths('/app/renv/library', action='replace', \
        renv::restore( \
            library='/app/renv/library', \
            exclude='arrow', \
            prompt=FALSE \
        ) \
    )" 2>&1 | tee /renv_restore.log

# write Renviron and .Rprofile
RUN R -e " \
    lib <- list.files('/app/renv/library', recursive=TRUE, \
        full.names=TRUE, include.dirs=TRUE); \
    lib <- lib[grepl('x86_64-pc-linux-gnu$', lib)][1]; \
    stopifnot(dir.exists(lib)); \
    cat('lib path:', lib, '\n'); \
    write(c(paste0('R_LIBS_USER=', lib), \
            'R_PROFILE_USER=/app/.Rprofile'), \
        '/usr/local/lib/R/etc/Renviron', append=TRUE); \
    writeLines( \
        'if (file.exists(\"/app/functions/prepShinyCellModular.R\")) source(\"/app/functions/prepShinyCellModular.R\")', \
        '/app/.Rprofile')"

RUN echo 'RENV_CONFIG_AUTOLOADER_ENABLED=FALSE' >> /usr/local/lib/R/etc/Renviron

RUN chmod -R 775 /app/renv/library