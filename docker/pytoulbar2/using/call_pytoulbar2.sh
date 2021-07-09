#!/bin/bash

###############################################################################
# Script to launch "python3 problem.py" command in Docker image
###############################################################################

PYTHONPATH=$PYTHONPATH:/WORK:/WORK/pythonplus
export PYTHONPATH

# specific required
pip3 install --cache-dir /WORK/cache -t /WORK/pythonplus --upgrade numpy

# run
cd /WORK
python3 problem.py

# clean
rm -fr /WORK/pythonplus
rm -fr /WORK/cache

