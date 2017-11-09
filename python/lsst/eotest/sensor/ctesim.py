"""
Simulate effects of CTE using cte_matrix.
"""
from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import astropy.io.fits as fits
import lsst.afw.image as afwImage
import lsst.eotest.image_utils as imutils
from .AmplifierGeometry import makeAmplifierGeometry
from .cte_matrix import cte_matrix
from . import sim_tools

_dtypes = dict([(-32, np.float32), (16, np.int16)])


def convert(imarr, bitpix):
    if bitpix > 0:
        my_round = np.round
    else:
        def my_round(x): return x
    return np.array(my_round(imarr), dtype=_dtypes[bitpix])


def fitsFile(segments, input):
    output = fits.HDUList()
    output.append(fits.PrimaryHDU())
    output[0].header = input[0].header
    for amp in segments:
        bitpix = input[amp].header['BITPIX']
        imarr = convert(segments[amp].getArray(), bitpix)
        output.append(fits.ImageHDU(data=imarr))
        output[amp].header = input[amp].header
    return output


def ctesim(infile, pcti=0, scti=0, verbose=False):
    input = fits.open(infile)
    amps = [i for i in range(1, len(input))
            if input[i].name.upper().startswith('SEGMENT')]
    segments = {}
    for amp in amps:
        if verbose:
            print("ctesim: working on amp", amp)
        image = afwImage.ImageF(infile, imutils.dm_hdu(amp))
        geom = makeAmplifierGeometry(infile)
        #
        # Temporarily remove readout bias median.
        #
        bias_med = imutils.median(image.Factory(image, geom.serial_overscan))

        image -= bias_med

        imarr = image.getArray()

        outimage = afwImage.ImageF(image, True)
        outarr = outimage.getArray()
        if pcti != 0:
            pcte_matrix = cte_matrix(imarr.shape[0], pcti)
        for col in range(0, imarr.shape[1]):
            outarr[:, col] = np.dot(pcte_matrix, imarr[:, col])
        if scti != 0:
            scte_matrix = cte_matrix(imarr.shape[1], scti)
        for row in range(0, imarr.shape[0]):
            outarr[row, :] = np.dot(scte_matrix, outarr[row, :])
        #
        # Restore readout bias
        #
        outarr += bias_med
        segments[amp] = outimage

    return fitsFile(segments, input)
