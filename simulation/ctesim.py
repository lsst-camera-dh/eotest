"""
@brief Simulate effects of CTE.  Based on Peter Doherty's IDL script,
ctesim.pro.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import numpy.random as random
import pyfits
import lsst.afw.image as afwImage
import image_utils as imutils
import simulation.sim_tools as sim_tools

_dtypes = dict([(-32, np.float32), (16, np.int16)])

def convert(imarr, bitpix):
    if bitpix > 0:
        my_round = np.round
    else:
        my_round = lambda x : x
    return np.array(my_round(imarr), dtype=_dtypes[bitpix])

def fitsFile(segments, input):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header = input[0].header
    for amp in segments:
        bitpix = input[amp].header['BITPIX']
        imarr = convert(segments[amp].getArray(), bitpix)
        output.append(pyfits.ImageHDU(data=imarr))
        output[amp].header = input[amp].header
    return output

def ctesim(infile, pcti=0, scti=0, verbose=False):
    pcte = 1 - pcti
    scte = 1 - scti

    input = pyfits.open(infile)
    amps = [i for i in range(1, len(input)) if input[i].is_image]
    segments = {}
    for amp in amps:
        if verbose:
            print "ctesim: working on amp", amp
        image = afwImage.ImageF(infile, imutils.dm_hdu(amp))
        #
        # Temporarily remove readout bias median.
        #
        bias_med = imutils.bias(image)
        image -= bias_med

        imarr = image.getArray()
        ny, nx = imarr.shape
    
        outimage = afwImage.ImageF(image, True)
        outarr = outimage.getArray()
        if pcti != 0:
            # Parallel charge transfer.
            # Loop over rows:
            for j in range(ny-2):
                if j % 100 == 0 and verbose:
                    print "  row", j
                # Copy bottom row to output.
                outarr[j, :] = pcte*imarr[0, :]
                # Calculate new shifted frame.
                #for jj in range(ny-2):
                for jj in range(ny - 2 - j):
                    imarr[jj, :] = pcti*imarr[jj, :] + pcte*imarr[jj+1, :]
                # Last row just has deferred charge.
                imarr[ny-1, :] = pcti*imarr[ny-1, :]
        if scti != 0:
            imarr = np.copy(outarr)
            # Serial charge transfer.
            # Loop over columns:
            for i in range(nx-2):
                if i % 100 == 0 and verbose:
                    print "  column", i
                outarr[:, i] = scte*imarr[:, 0]
                for ii in range(nx - 2 - i):
                    imarr[:, ii] = scti*imarr[:, ii] + scte*imarr[:, ii+1]
                imarr[:, nx-1] = scti*imarr[:, nx-1]
        #
        # Restore readout bias
        #
        outarr += bias_med
        segments[amp] = outimage

    return fitsFile(segments, input)

def make_fe55(outfile, nxrays, bias_level=2000, bias_sigma=4,
              amps=imutils.allAmps):
    segments = []
    for amp in amps:
        seg = sim_tools.SegmentExposure()
        seg.add_bias(level=bias_level, sigma=bias_sigma)
        seg.add_Fe55_hits(nxrays=nxrays)
        segments.append(seg)
    sim_tools.writeFits(segments, outfile)

if __name__ == '__main__':
    test_file = 'fe55_test.fits'
    nxrays = 1000
    
    make_fe55(test_file, nxrays, amps=(1,))

    pcti = 1e-5
#    scti = 2e-5
    scti = 0

    foo = ctesim(test_file, pcti=pcti, scti=scti)
    foo.writeto('fe55_test_cti.fits', clobber=True)
