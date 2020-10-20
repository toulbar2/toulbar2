#!/bin/bash
# This script will install the libraries required to compile toulbar2
# on a MacOS system. This has been tested on MacOS
# Catalina. Suggestions for improvement are welcome. Contact us
# through GitHub. Everything can be uninstalled by executing:
# /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall.sh)"

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
brew install cmake
brew install gmp
brew install zlib
brew install xz
brew install boost
brew install pybind11
pip3 install pybind11
