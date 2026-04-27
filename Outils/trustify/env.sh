#!/bin/bash
# Source environment for trustify

root_dir=$TRUST_ROOT/Outils/trustify/install

# Python version:
py_vers=`python -c "import sys; print('%d.%d' % (sys.version_info.major, sys.version_info.minor))"` 

export TRUSTIFY_ROOT_DIR=${root_dir}/lib/python${py_vers}/site-packages
export PYTHONPATH=${TRUSTIFY_ROOT_DIR}:${PYTHONPATH}

