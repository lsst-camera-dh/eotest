from __future__ import print_function
import os
import numpy as np
from sim_tools import SegmentExposure, fitsFile
from sim_inputs import gain, bias_level, sys_noise, read_noise, dark_current


def simulateFlat(outfile, rate, gain, dark_curr, exptime, hdus=16):
    segments = []
    for i in range(hdus):
        seg = SegmentExposure(exptime=exptime, gain=gain)
        seg.add_bias(level=bias_level, sigma=sys_noise) # electronics
        seg.add_bias(level=0, sigma=read_noise)         # CCD read noise
        seg.add_dark_current(level=dark_curr)
        seg.expose_flat(rate)
        segments.append(seg)
    output = fitsFile(segments)
    output[0].header.update("BIASLVL", bias_level)
    output[0].header.update("SYSNOISE", sys_noise)
    output[0].header.update("READNOIS", read_noise)
    output[0].header.update("DARKCURR", dark_curr)
    output[0].header.update("COUNTRTE", rate)


if __name__ == '__main__':
    sensor_id = '000-00'
    outdir = '/nfs/farm/g/lsst/u1/testData/SIMData/%s/flats' % sensor_id
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    rate = 3   # Count rate from uniform light source (DN/sec)
    exptimes = np.logspace(0, 2, 10)
    for i, exptime in enumerate(exptimes):
        for j in range(2):  # simulate flats in pairs in order to subtract FPN
            filename = 'flat_%s_%02i_%i.fits' % (sensor_id, i, j)
            outfile = os.path.join(outdir, filename)
            print("generating", outfile)
            simulateFlat(outfile, rate, gain, dark_current, exptime)
