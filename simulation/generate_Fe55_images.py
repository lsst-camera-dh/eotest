"""
@brief Generate Fe55 images and associated darks and bias images
according to section 5.4 of the E/O document (Dec 19, 2012 version).

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
from sim_inputs import *

def generate_Fe55_images(exptimes, nxrays, outdir, sensorid):
    nexp = len(exptimes)
    for i, exptime, nxray in zip(range(nexp), exptimes, nxrays):
        #
        # Bias images
        #
        outfile = "Fe55_bias_%s_%02i.fits" % (sensorid, i)
        bias_file = os.path.join(outdir, outfile)
        bias_segs = []
        for hdu in range(nhdu):
            seg = SegmentExposure(exptime=0, gain=gain)
            seg.add_bias(level=level, sigma=sys_noise) # electronics
            seg.add_bias(sigma=sigma, level=0) # read noise
            bias_segs.append(seg)
        writeFits(bias_segs, bias_file)
        #
        # Dark images
        #
        outfile = "Fe55_dark_%s_%02i.fits" % (sensorid, i)
        dark_file = os.path.join(outdir, outfile)
        dark_segs = []
        for hdu in range(nhdu):
            seg = SegmentExposure(exptime=exptime, gain=gain)
            seg.add_bias(level=level, sigma=sys_noise) # electronics
            seg.add_bias(level=0, sigma=sigma) # read noise
            seg.add_dark_current(level=200) # dark current
            dark_segs.append(seg)
        writeFits(dark_segs, dark_file)
        #
        # Fe55 exposures
        #
        outfile = "Fe55_exp_%s_%02i.fits" % (sensorid, i)
        Fe55_file = os.path.join(outdir, outfile)
        fe55_segs = []
        for hdu in range(nhdu):
            seg = SegmentExposure(exptime=exptime, gain=gain)
            seg.add_bias(level=level, sigma=sys_noise) # electronics
            seg.add_bias(level=0, sigma=sigma) # read noise
            seg.add_dark_current(level=200) # dark current
            seg.add_Fe55_hits(nxrays=nxray)
            fe55_segs.append(seg)
        writeFits(fe55_segs, Fe55_file)

if __name__ == '__main__':
    nexp = 10

    exptimes = np.linspace(1, 5, nexp)
    nxrays = [int(x*1000) for x in exptimes] 

    generate_Fe55_images(exptimes, nxrays, '.', 'xxx-xx')

