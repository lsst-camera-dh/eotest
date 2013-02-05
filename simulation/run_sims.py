"""
@brief Generate a standard set of simulated data for a subset of
sensors using a specified directory structure.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os

import numpy as np

from simulation.generate_Fe55_images import *
from simulation.generate_system_noise_images import *

def mkdir(path):
    try:
        os.mkdir(path)
    except OSError:
        pass

sim_data_root = '/nfs/farm/g/lsst/u1/testData/SIMData'

sensors = ['000-%02i' % i for i in range(4)]

nexp = 10
exptimes = np.linspace(1, 5, nexp)
nxrays = [int(x*1000) for x in exptimes] 

os.chdir(sim_data_root)
for sensor in sensors:
    print "Creating data for sensor", sensor
    sensor_dir = os.path.join(sim_data_root, sensor)
    mkdir(sensor_dir)
    #
    # Fe55 data
    #
    print "Fe55 images..."
    fe55_dir = os.path.join(sensor_dir, 'Fe55')
    mkdir(fe55_dir)
    generate_Fe55_images(exptimes, nxrays, fe55_dir, sensor)
    #
    # Readout system noise images
    #
    print "Readout system noise images..."
    systemnoise_dir = os.path.join(sensor_dir, 'system_noise')
    mkdir(systemnoise_dir)
    generate_system_noise_images(nexp, systemnoise_dir, sensor)
