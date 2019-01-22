"""
@brief Generate the readout system noise frames according to section
5.8 of the E/O document (Dec 19, 2012 version).

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
from sim_inputs import *


def generate_system_noise_images(nexp, outdir, sensorid, gain=gain,
                                 bias_level=bias_level, sys_noise=sys_noise):
    exptime = 0
    for i in range(nexp):
        filename = 'system_noise_%s_%02i.fits' % (sensorid, i)
        outfile = os.path.join(outdir, filename)
        noise_segs = []
        for hdu in range(nhdu):
            seg = SegmentExposure(exptime=exptime, gain=gain)
            seg.add_bias(level=bias_level, sigma=sys_noise)
            noise_segs.append(seg)
        output = fitsFile(noise_segs)
        output[0].header['BIASLVL'] = bias_level
        output[0].header['SYSNOISE'] = sys_noise
        output.writeto(outfile, overwrite=True)


if __name__ == '__main__':
    nexp = 10
    generate_system_noise_images(nexp, '.')
