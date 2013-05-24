#!/bin/bash

echo "Running unit tests"

export SHELL=/bin/bash

export HOME=`pwd`

source /afs/slac/g/lsst/software/redhat5-x86_64-64bit-gcc44/DMstack/Winter2013-v6_2/loadLSST.sh

setup -t v6_2 pipe_tasks
setup -t v6_2 meas_algorithms
setup mysqlpython
setup scipy

export PYTHONPATH=`pwd`:${PYTHONPATH}

echo `pwd`

cd tests; ./run_tests.sh
