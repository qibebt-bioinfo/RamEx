FROM rocker/r-ver:4.4.1

WORKDIR /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
        libxml2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libcairo2-dev \
        libxt-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libglpk-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libuv1-dev \
        pandoc \
    && rm -rf /var/lib/apt/lists/*

COPY . /tmp/RamEx

RUN rm -f /tmp/RamEx/src/*.o /tmp/RamEx/src/*.so /tmp/RamEx/src/*.dll
RUN R -e "install.packages('remotes')"
RUN R -e "remotes::install_deps('/tmp/RamEx', dependencies = TRUE, upgrade = 'never')"

RUN R CMD INSTALL /tmp/RamEx \
    && R -e "library(RamEx)" \
    && rm -rf /tmp/RamEx
