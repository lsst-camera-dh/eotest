#!/bin/bash

echo "Running unit tests"

export SHELL=/bin/bash

# The following is needed so that eups will look for or try to write
# its config stuff in cwd rather than in the user's home directory,
# which is a problem for batch jobs that don't have afs tokens. It
# also avoids cluttering the batch user's home directory.
#
export HOME=`pwd`

# This is used by the scripts to find the path to the policy subdir.
export TEST_SCRIPTS_DIR=${HOME}

echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/lib:/usr/lib:${LD_LIBRARY_PATH}

source /afs/slac/g/lsst/software/redhat5-x86_64-64bit-gcc44/DMstack/Winter2013-v6_2/loadLSST.sh

# The Jenkins build of test_scripts can't/shouldn't be setup by eups, so
# setup everything by hand.
#
setup -t v6_2 pipe_tasks
setup -t v6_2 meas_algorithms
setup mysqlpython
setup scipy

export PYTHONPATH=${HOME}/python:${HOME}/bin/jenkins:${PYTHONPATH}

echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH

echo ${HOME}

cd tests; python run_all.py
