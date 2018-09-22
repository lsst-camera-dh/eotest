#!/bin/bash
source /opt/lsst/software/stack/loadLSST.bash
setup lsst_distrib
pip install nose
pip install coveralls
eups declare eotest -r ${TRAVIS_BUILD_DIR} -t current
setup eotest
cd ${TRAVIS_BUILD_DIR}
scons opt=3
nosetests -s --with-coverage --cover-package=lsst.eotest.image_utils,lsst.eotest.sensor
