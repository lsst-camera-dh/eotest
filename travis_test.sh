#!/bin/bash
source /opt/lsst/software/stack/loadLSST.bash
setup lsst_distrib
pip install nose
pip install coveralls
eups declare eotest -r /home/vagrant/eotest -t current
setup eotest
scons opt=3
cd /home/vagrant/eotest
nosetests -s --with-coverage --cover-package=lsst.eotest.image_utils,lsst.eotest.sensor
