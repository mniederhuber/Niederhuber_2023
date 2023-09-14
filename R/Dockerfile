FROM rocker/r-ver:4.2.2

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"

ENV S6_VERSION=v2.1.0.2
ENV RSTUDIO_VERSION=2022.07.2+576
ENV DEFAULT_USER=rstudio
ENV PANDOC_VERSION=default
ENV QUARTO_VERSION=default

RUN /rocker_scripts/install_rstudio.sh
RUN /rocker_scripts/install_pandoc.sh
RUN /rocker_scripts/install_quarto.sh
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN apt-get update && \
      apt-get -y install \
      libxml2-dev \
      libpng-dev \
      liblzma-dev \
     	libz-dev \
	libbz2-dev \
      libfontconfig1-dev \
      libtcl8.6 \
      libtk8.6

#WORKDIR /project
COPY renv.lock renv.lock

ENV RENV_PATHS_LIBRARY renv/library

RUN R -e "renv::restore()"

EXPOSE 8787

CMD ["/init"]
