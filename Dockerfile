# Shared image for running ShinyCellModular apps behind ShinyProxy.
#
# ShinyProxy starts one container per user session from this single image,
# across every dataset -- the dataset itself is not baked in here, it is
# mounted at container start as a bundle_dir produced by
# makeShinyCellModularPortable() (app.R, app_config.yml, modules/, data).
# Only the mount path (SCMODULAR_DATA_DIR) differs between ShinyProxy's
# `spec:` entries in application.yml; the image is the same for all of them.
#
# Build once, push to a registry, point application.yml's container-image at it:
#   docker build -t <registry>/<name>:<tag> .
#   docker push <registry>/<name>:<tag>
#
# No Shiny Server here: ShinyProxy is what starts/stops/routes containers, so
# a Shiny Server layer would just sit unused.
FROM docker.io/rocker/r-ver:4.5.1

# Covers what this package's runtime dependencies commonly need to compile
# from source: libhdf5-dev (hdf5r), libcurl/libssl/libxml2 (arrow and friends),
# libglpk-dev (igraph, pulled in via Seurat/Signac), font/image libs (ggplot2,
# plotly rendering), git (remotes::install_github below). Extend this list if
# a build fails on a missing system library.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    zlib1g-dev \
    libglpk-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    cmake \
    git \
  && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('BiocManager', 'remotes'))"

# CRAN + Bioconductor packages a ShinyCellModular app needs at runtime.
# BiocManager::install() handles both transparently, no CRAN/Bioc split needed.
# Mirrors DESCRIPTION's Imports/Suggests, minus knitr/rmarkdown (only needed to
# build this package's own vignettes, not to run a deployed app) and the
# GitHub-only packages installed separately below. Update this list by hand
# when DESCRIPTION's dependencies change.
RUN R -e "BiocManager::install(c( \
    'arrow', 'jsonlite', 'rsconnect', 'yaml', \
    'data.table', 'dplyr', 'DT', 'edgeR', 'ggdendro', 'ggforce', \
    'ggplot2', 'ggrepel', 'ggseqlogo', 'gridExtra', 'hdf5r', 'limma', \
    'magrittr', 'Matrix', 'plotly', 'purrr', 'RColorBrewer', \
    'Seurat', 'shiny', 'shinydashboard', 'shinyhelper', 'shinyjs', \
    'shinythemes', 'Signac', 'sortable', 'stringr', 'tidyr' \
  ), update = FALSE, ask = FALSE)"

# GitHub-only packages (see DESCRIPTION's Remotes field)
RUN R -e "remotes::install_github(c('immunogenomics/presto', 'SGDDNB/ShinyCell'), upgrade = 'never')"

# ShinyCellModular itself
RUN R -e "remotes::install_github('MonashBioinformaticsPlatform/ShinyCellModular', upgrade = 'never')"

EXPOSE 3838

# SCMODULAR_DATA_DIR is set per dataset by ShinyProxy's container-env, pointing
# at whatever bundle_dir is mounted as a volume. Falls back to /data so the
# image also runs with a plain `docker run -v <bundle_dir>:/data`.
CMD ["R", "-e", "shiny::runApp(appDir = Sys.getenv('SCMODULAR_DATA_DIR', '/data'), host = '0.0.0.0', port = 3838)"]
