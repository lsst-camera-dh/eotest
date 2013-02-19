"""
@brief Compute distributions of the CCD read noise estimates for each
segment.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import numpy.random as random

import lsst.afw.image as afwImage

from image_utils import trim, imaging, allAmps, SubRegionSampler, \
    mean, median, stdev

class NoiseDists(SubRegionSampler):
    def __init__(self, gains, dx=100, dy=100, nsamp=1000, imaging=imaging):
        SubRegionSampler.__init__(self, dx, dy, nsamp, imaging)
        self.gains = gains
    def __call__(self, infile):
        noise_dists = []
        for amp, gain in self.gains.items():
            im = trim(afwImage.ImageF(infile, amp+1))
            noise_samples = []
            for x, y in zip(self.xarr, self.yarr):
                subim = self.subim(im, x, y)
                noise_samples.append(stdev(subim)*gain)
            noise_dists.append(np.array(noise_samples))
        return np.array(noise_dists)

if __name__ == '__main__':
    from simulation.sim_tools import SegmentExposure, writeFits
    #
    # Simulate bias image.
    #
    gains = dict([(amp, 5.5) for amp in allAmps])

    bias_file = "bias.fits"
    bias_segs = []
    for amp in allAmps:
        seg = SegmentExposure(exptime=0, gain=gains[amp])
        seg.add_bias(level=1e4, sigma=5) # electronic bias and noise
        seg.add_bias(sigma=4)            # CCD read noise
        bias_segs.append(seg)
    writeFits(bias_segs, bias_file)
    #
    # Simulate readout system noise.
    #
    readout_noise_file = 'readout_noise.fits'
    noise_segs = []
    for amp in allAmps:
        seg = SegmentExposure(exptime=0, gain=gains[amp])
        seg.add_bias(level=1e4, sigma=5) # electronic bias and noise
        noise_segs.append(seg)
    writeFits(noise_segs, readout_noise_file)
    #
    # Compute the distribution of noise estimates.
    #
    noise_dists = NoiseDists(gains)
    Ntot = noise_dists(bias_file)
    Nsys = noise_dists(readout_noise_file)
    #
    # Read noise distribution.
    #
    Nread = np.sqrt(Ntot*Ntot - Nsys*Nsys)
    for nread in Nread:
        print median(nread), stdev(nread)
