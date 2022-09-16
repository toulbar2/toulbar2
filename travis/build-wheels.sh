#!/bin/bash
set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w /io/wheelhouse/
    fi
}

# Install a system package required by our library
yum -y install gmp-devel
yum -y install boost-devel
yum -y install zlib-devel
yum -y install bzip2-devel
yum -y install xz-devel
yum -y install libxml2-devel

cd "$HOME"

# Compile wheels for cython only on Linux
for PYBIN in /opt/python/cp*/bin; do
    "${PYBIN}/pip3" install -r /io/dev-requirements.txt
    "${PYBIN}/pip3" wheel /io/ --no-deps -w wheelhouse
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    repair_wheel "$whl"
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip3" install pytoulbar2 --no-index -f /io/wheelhouse
    "${PYBIN}/python" -m nose2 pytoulbar2
done
