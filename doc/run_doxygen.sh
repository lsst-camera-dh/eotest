#!/bin/bash

echo "Generating Doxygen documentation using doxypypy"

export SHELL=/bin/bash

export HOME=`pwd`  

export PYTHONPATH=${HOME}:${HOME}/jenkins:${PYTHONPATH}

echo ${HOME}

doxygen doc/Doxyfile
