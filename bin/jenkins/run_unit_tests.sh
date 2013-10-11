#!/bin/bash

echo "Running unit tests"

export SHELL=/bin/bash

# The following is needed so that eups will look for or try to write
# its config stuff in cwd rather than in the user's home directory,
# which is a problem for batch jobs that don't have afs tokens. It
# also avoids cluttering the batch user's home directory.
export HOME=`pwd`  

source /afs/slac/g/lsst/software/redhat5-x86_64-64bit-gcc44/DMstack/Winter2013-v6_2/loadLSST.sh

setup -t v6_2 pipe_tasks
setup -t v6_2 meas_algorithms
setup mysqlpython
setup scipy

export PYTHONPATH=${HOME}:${HOME}/jenkins:${PYTHONPATH}

echo ${HOME}

cd tests; python run_all.py
