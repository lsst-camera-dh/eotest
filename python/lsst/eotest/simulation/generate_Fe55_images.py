"""
@brief Generate Fe55 images and associated darks and bias images
according to section 5.4 of the E/O document (Dec 19, 2012 version).

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from builtins import zip
from builtins import range
import os
import numpy as np
from sim_inputs import *
from sim_tools import *


def generate_Fe55_images(exptimes, nxrays, outdir, sensorid, gain=gain,
                         bias_level=bias_level, sys_noise=sys_noise,
                         dark_current=dark_current):
    nexp = len(exptimes)
    for i, exptime, nxray in zip(list(range(nexp)), exptimes, nxrays):
        #
        # Bias images
        #
        outfile = "Fe55_bias_%s_%02i.fits" % (sensorid, i)
        bias_file = os.path.join(outdir, outfile)
        bias_segs = []
        for hdu in range(nhdu):
            seg = SegmentExposure(exptime=0, gain=gain)
            seg.add_bias(level=bias_level, sigma=sys_noise) # electronics
            seg.add_bias(level=0, sigma=read_noise)         # read noise
            bias_segs.append(seg)
        bias_output = fitsFile(bias_segs)
        bias_output[0].header['GAIN'] = gain
        bias_output[0].header['BIASLVL'] = bias_level
        bias_output[0].header['SYSNOISE'] = sys_noise
        bias_output[0].header['RDNOISE'] = read_noise
        bias_output.writeto(bias_file, overwrite=True)
        #
        # Dark images
        #
        outfile = "Fe55_dark_%s_%02i.fits" % (sensorid, i)
        dark_file = os.path.join(outdir, outfile)
        dark_segs = []
        for hdu in range(nhdu):
            seg = SegmentExposure(exptime=exptime, gain=gain)
            seg.add_bias(level=bias_level, sigma=sys_noise) # electronics
            seg.add_bias(level=0, sigma=read_noise)         # read noise
            seg.add_dark_current(level=dark_current)        # dark current
            dark_segs.append(seg)
        dark_output = fitsFile(dark_segs)
        dark_output[0].header['GAIN'] = gain
        dark_output[0].header['BIASLVL'] = bias_level
        dark_output[0].header['SYSNOISE'] = sys_noise
        dark_output[0].header['RDNOISE'] = read_noise
        dark_output[0].header['DARKCURR'] = dark_current
        dark_output.writeto(dark_file, overwrite=True)
        #
        # Fe55 exposures
        #
        outfile = "Fe55_exp_%s_%02i.fits" % (sensorid, i)
        Fe55_file = os.path.join(outdir, outfile)
        fe55_segs = []
        for hdu in range(nhdu):
            seg = SegmentExposure(exptime=exptime, gain=gain)
            seg.add_bias(level=bias_level, sigma=sys_noise) # electronics
            seg.add_bias(level=0, sigma=read_noise)         # read noise
            seg.add_dark_current(level=dark_current)        # dark current
            seg.add_Fe55_hits(nxrays=nxray)
            fe55_segs.append(seg)
        fe55_output = fitsFile(fe55_segs)
        fe55_output[0].header['GAIN'] = gain
        fe55_output[0].header['BIASLVL'] = bias_level
        fe55_output[0].header['SYSNOISE'] = sys_noise
        fe55_output[0].header['RDNOISE'] = read_noise
        fe55_output[0].header['DARKCURR'] = dark_current
        fe55_output[0].header['FE55HITS'] = nxray
        fe55_output.writeto(Fe55_file, overwrite=True)


if __name__ == '__main__':
    nexp = 10

    exptimes = np.linspace(1, 5, nexp)
    nxrays = [int(x*1000) for x in exptimes]

    generate_Fe55_images(exptimes, nxrays, '.', 'xxx-xx')
