"""
@brief Generate the readout system noise frames according to section
5.8 of the E/O document (Dec 19, 2012 version).

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
from sim_inputs import *

exptime = 0

nexp = 10
for i in range(nexp):
    outfile = 'readout_noise_%02i.fits' % i
    noise_segs = []
    for hdu in range(nhdu):
        seg = SegmentExposure(exptime=exptime, gain=gain)
        seg.add_bias(level=level, sigma=sys_noise)
        noise_segs.append(seg)
    writeFits(noise_segs, outfile)
