# Pinned to linux/amd64 so the same image works on x86 Linux hosts AND on
# Apple Silicon Macs (via Docker Desktop's Rosetta emulation). This is
# required because several conda-forge / bioconda packages used here
# (quarto, r-umap, parts of the PECA chain) are not built for linux/arm64.
FROM --platform=linux/amd64 condaforge/miniforge3:latest

ENV ENV_NAME=spectropiper
ENV DEBIAN_FRONTEND=noninteractive

# System libs needed for any R packages that fall through to source compile
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        zlib1g-dev \
        ca-certificates \
        git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/SpectroPipeR

# 1. Build conda env from environment.yml (CRAN + Bioc deps as binaries,
#    including the PECA transitive chain — see comments in environment.yml)
COPY environment.yml .
RUN mamba env create -n ${ENV_NAME} -f environment.yml \
    && mamba clean -afy

# Put the env's binaries first on PATH so plain `R` / `Rscript` use the env.
# Simpler and more reliable than wrapping every RUN in `mamba run`.
ENV PATH=/opt/conda/envs/spectropiper/bin:$PATH

# 2. iq (not packaged for conda) and ROTS/PECA via R/BiocManager
RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); \
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); \
        install.packages("iq"); \
        BiocManager::install(c("ROTS","PECA"), update = FALSE, ask = FALSE); \
        stopifnot(requireNamespace("iq"), requireNamespace("PECA"))'

# 3. Install SpectroPipeR.
#    Default: pull from GitHub. Override the ref at build time with
#      `--build-arg SPECTROPIPER_REF=<branch|tag|sha>`.
#    Or install the local checkout (useful for testing un-pushed patches):
#      `--build-arg SPECTROPIPER_LOCAL=true`
ARG SPECTROPIPER_REPO=stemicha/SpectroPipeR
ARG SPECTROPIPER_REF=main
ARG SPECTROPIPER_LOCAL=false
# Always copy the local checkout into the image — the build context is
# kept tiny by .dockerignore, so this is cheap. It's used only when
# SPECTROPIPER_LOCAL=true.
COPY DESCRIPTION NAMESPACE NEWS.md LICENSE LICENSE.md /opt/SpectroPipeR/pkg/
COPY R    /opt/SpectroPipeR/pkg/R/
COPY inst /opt/SpectroPipeR/pkg/inst/
COPY man  /opt/SpectroPipeR/pkg/man/
RUN if [ "$SPECTROPIPER_LOCAL" = "true" ]; then \
        echo "Installing SpectroPipeR from local checkout..." && \
        Rscript -e 'remotes::install_local("/opt/SpectroPipeR/pkg", upgrade="never", dependencies=FALSE); stopifnot(requireNamespace("SpectroPipeR"))' ; \
    else \
        echo "Installing SpectroPipeR from GitHub ${SPECTROPIPER_REPO}@${SPECTROPIPER_REF}..." && \
        Rscript -e "remotes::install_github('${SPECTROPIPER_REPO}', ref='${SPECTROPIPER_REF}', upgrade='never', dependencies=FALSE); stopifnot(requireNamespace('SpectroPipeR'))" ; \
    fi

# 4. Smoke test (mirrors tests/testthat/test-SpectR_test_all.R)
COPY docker/smoke_test.R /opt/SpectroPipeR/smoke_test.R

# Default: run the smoke test. `docker run --rm <image>` returns non-zero on failure.
CMD ["Rscript", "/opt/SpectroPipeR/smoke_test.R"]
