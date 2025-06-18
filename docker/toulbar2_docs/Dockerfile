FROM python:slim

LABEL org.opencontainers.image.authors="ToulBar2 team <http://miat.inrae.fr/toulbar2>"
LABEL org.opencontainers.image.documentation="http://miat.inrae.fr/toulbar2"
LABEL org.opencontainers.image.source="https://github.com/toulbar2/toulbar2"
LABEL org.opencontainers.image.version="0.1"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.description="sphinx/LaTeX container image for ToulBar2 documentation"

WORKDIR /docs
RUN apt-get update
RUN apt-get install -y cmake g++ libgmp-dev libboost-graph-dev libboost-iostreams-dev zlib1g-dev liblzma-dev libxml2-dev libopenmpi-dev libboost-mpi-dev libjemalloc-dev pkg-config doxygen git bc
RUN apt-get install --no-install-recommends --yes \
      graphviz \
      imagemagick \
      make \
      latexmk \
      lmodern \
      fonts-freefont-otf \
      texlive-latex-recommended \
      texlive-latex-extra \
      texlive-fonts-recommended \
      texlive-fonts-extra \
      texlive-lang-cjk \
      texlive-lang-chinese \
      texlive-lang-japanese \
      texlive-luatex \
      texlive-xetex \
      xindy \
      tex-gyre \
 && apt-get autoremove \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --no-cache-dir --upgrade pip \
 && python3 -m pip install --no-cache-dir Sphinx==8.2.3 Pillow

# pacakages from toulbar2/docs/requirements.txt
RUN python3 -m pip install --no-cache-dir sphinx_rtd_theme breathe myst-parser ghp-import

CMD ["/bin/bash", "-c", "mkdir -p build && cd build && cmake -DBUILD_API_DOC=ON /docs && make doc && cd /docs/docs && make docs"]
