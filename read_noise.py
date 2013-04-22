"""
@brief Compute distributions of the CCD read noise estimates for each
segment.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from MaskedCCD import MaskedCCD
import image_utils as imutils

class NoiseDist(object):
    def __init__(self, imfile, region_sampler, mask_files=()):
        self.ccd = MaskedCCD(imfile, mask_files=mask_files)
        self.sampler = region_sampler
    def __call__(self, amp, gain):
        image = self.ccd[amp].Factory(self.ccd[amp], self.sampler.imaging)
        noise_samples = []
        for x, y in zip(self.sampler.xarr, self.sampler.yarr):
            subim = self.sampler.subim(image, x, y)
            stdev = afwMath.makeStatistics(subim, afwMath.STDEV,
                                           self.ccd.stat_ctrl).getValue()
            noise_samples.append(stdev*gain)
        return np.array(noise_samples)

def noise_dists(imfile, gains, sampler, mask_files=()):
    noise_dist = NoiseDist(imfile, sampler, mask_files=mask_files)
    my_noise_dists = {}
    for amp in imutils.allAmps:
        my_noise_dists[amp] = noise_dist(amp, gains[amp])
    return my_noise_dists

if __name__ == '__main__':
    from simulation.sim_tools import SegmentExposure, writeFits
    #
    # Simulate bias image.
    #
    gains = dict([(amp, 5.5) for amp in imutils.allAmps])

    bias_file = "bias.fits"
    bias_segs = []
    for amp in imutils.allAmps:
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
    for amp in imutils.allAmps:
        seg = SegmentExposure(exptime=0, gain=gains[amp])
        seg.add_bias(level=1e4, sigma=5) # electronic bias and noise
        noise_segs.append(seg)
    writeFits(noise_segs, readout_noise_file)
    #
    # Compute the distribution of noise estimates.
    #
    dx, dy = 100, 100
    nsamp = 1000
    imaging = imutils.imaging
    sampler = imutils.SubRegionSampler(dx, dy, nsamp, imaging)
    #mask_files = ('CCD250_DEFECTS_mask.fits',)
    mask_files = ()
    Ntot = noise_dists(bias_file, gains, sampler, mask_files=mask_files)
    Nsys = noise_dists(readout_noise_file, gains, sampler, mask_files)
    #
    # Read noise distribution.
    #
    for amp in imutils.allAmps:
        Nread = np.sqrt(Ntot[amp]*Ntot[amp] - Nsys[amp]*Nsys[amp])
        print imutils.median(Nread), imutils.stdev(Nread)
