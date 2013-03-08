#!/bin/bash

export SHELL=bash
source $DMSTACK_DIR/loadLSST.sh  
setup pipe_tasks
setup afw
setup mysqlpython
export PYTHONPATH=$LSSTSCRIPT_DIR:$LSSTSCRIPT_DIR/pipeline:$PYTHONPATH
echo 'Your Environment is now set up to use the DMstack and LSST Analysis Scripts'
