"""
@brief Compute read noise distributions for a sample of images.  The
system gain for each channel is determined from Fe55 data, and the
corresponding bias frame is used to characterize the total system
noise.  Readout noise frames are used for determining the noise
contribution from the electronics.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
from read_noise import NoiseDists, median
from xray_gain import hdu_gains

nexp = 10
for i in range(nexp):
    Fe55_exp_file = 'Fe55_exp_%02i.fits' % i
    bias_frame_file = 'Fe55_bias_%02i.fits' % i
    readout_noise_file = 'readout_noise_%02i.fits' % i
    
    gains = hdu_gains(Fe55_exp_file)
    noise_dists = NoiseDists(gains)
    
    Ntot = noise_dists(bias_frame_file)
    Nsys = noise_dists(readout_noise_file)
        
    Nread = np.sqrt(Ntot*Ntot - Nsys*Nsys)
    
    for hdu in range(16):
        print i, hdu, median(Nread[hdu])
##
## Plotting using pylab (wrapped by my pylab_plotter module)
##
#    import pylab_plotter as plot
#        print i
#        for hdu, noise in enumerate(Nread):
#            xmed = median(noise)
#            xsig = stdev(noise)
#            win = plot.histogram(noise, xname='e-', yname='counts/bin',
#                                 xrange=(xmed-5*xsig, xmed+5*xsig), 
#                                 bins=50, oplot=hdu)
#            win.set_title("exposures %i" % i)
