"""
@brief Compute distributions of the CCD read noise estimates for each
segment.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from MaskedCCD import MaskedCCD
import lsst.eotest.image_utils as imutils

def noise_samples(raw_image, gain, region_sampler,
                  stat_ctrl=afwMath.StatisticsControl()):
    image = raw_image.Factory(raw_image, region_sampler.imaging)
    samples = []
    for x, y in zip(region_sampler.xarr, region_sampler.yarr):
        subim = region_sampler.subim(image, x, y)
        stdev = afwMath.makeStatistics(subim, afwMath.STDEV,
                                       stat_ctrl).getValue()
        samples.append(stdev*gain)
    return np.array(samples)

def noise_dists(imfile, gains, sampler, mask_files=()):
    ccd = MaskedCCD(imfile, mask_files=mask_files)
    my_noise_dists = {}
    for amp in imutils.allAmps:
        my_noise_dists[amp] = noise_samples(ccd[amp], gains[amp], sampler,
                                            ccd.stat_ctrl)
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