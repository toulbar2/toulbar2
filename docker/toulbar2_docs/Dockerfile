FROM sphinxdoc/sphinx-latexpdf

LABEL org.opencontainers.image.authors="ToulBar2 team <http://miat.inrae.fr/toulbar2>"
LABEL org.opencontainers.image.documentation="http://miat.inrae.fr/toulbar2"
LABEL org.opencontainers.image.source="https://github.com/toulbar2/toulbar2"
LABEL org.opencontainers.image.version="0.1"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.description="sphinx/LaTeX container image for ToulBar2 documentation"

WORKDIR /docs
ADD docs/requirements.txt /docs

RUN apt-get update
RUN apt-get install -y cmake g++ libgmp-dev libboost-graph-dev libboost-iostreams-dev zlib1g-dev liblzma-dev libxml2-dev libopenmpi-dev libboost-mpi-dev libjemalloc-dev pkg-config doxygen git bc && apt-get autoremove && apt-get clean  && rm -rf /var/lib/apt/lists/*

# requirements.txt from toulbar2/docs 
RUN python3 -m pip install -r requirements.txt


CMD ["/bin/bash", "-c", "mkdir -p build && cd build && cmake -DBUILD_API_DOC=ON /docs && make doc && cd /docs/docs && make docs"]
