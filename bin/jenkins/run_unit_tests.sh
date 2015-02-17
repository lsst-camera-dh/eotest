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
export EOTEST_DIR=${HOME}

source /nfs/farm/g/lsst/u1/software/redhat6-x86_64-64bit-gcc44/DMstack/Winter2014/loadLSST.sh

# The Jenkins build of test_scripts can't/shouldn't be setup by eups, so
# setup everything by hand.
#
setup pipe_tasks
setup meas_algorithms
setup mysqlpython
setup scipy

export PYTHONPATH=${HOME}/python:${HOME}/bin/jenkins:${PYTHONPATH}

echo ${HOME}

cd tests; python run_all.py
