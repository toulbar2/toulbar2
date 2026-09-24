
## Docker image for ToulBar2 sphinx documentation

Image based on sphinxdoc/sphinx-latexpdf
Commands must be run at the root of toulbar2/ source directory.

#### Build docker image
docker build -t toulbar2_docs -f ./docker/toulbar2_docs/Dockerfile .

#### Run
docker run -v ./:/docs toulbar2_docs

