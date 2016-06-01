"""
Function to generate mask files given a set of pixels and columns to mask.
"""
from __future__ import absolute_import, print_function
import os
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsWriteto
import lsst.eotest.image_utils as imutils
from .AmplifierGeometry import makeAmplifierGeometry
from .BrightPixels import BrightPixels
from .MaskedCCD import MaskedCCD
from .sim_tools import CCD

def generate_mask(infile, outfile, mask_plane, pixels=None, columns=None,
                  temp_mask_file='temp_mask_image.fits'):
    """
    Generate a mask file for the specified pixels and columns.
    The amplifier geometry will be taken from infile.
    """
    # Insert artificial signal into specified pixels and columns for
    # each segment.
    ccd = CCD(exptime=1, gain=1, geometry=makeAmplifierGeometry(infile))
    signal = 10

    if pixels is None:
        pixels = {}
    for amp in pixels:
        imarr = imutils.trim(ccd.segments[amp].image,
                             ccd.segments[amp].geometry.imaging).getArray()
        for ix, iy in pixels[amp]:
            imarr[iy][ix] = signal

    if columns is None:
        columns = {}
    for amp in columns:
        imarr = imutils.trim(ccd.segments[amp].image,
                             ccd.segments[amp].geometry.imaging).getArray()
        for ix in columns[amp]:
            imarr[:, ix] = signal

    fitsWriteto(ccd, temp_mask_file)

    # Use BrightPixels code to detect mask regions and write the mask file.
    hdulist = fits.HDUList()
    hdulist.append(fits.open(infile)[0])
    hdulist[0].header['MASKTYPE'] = mask_plane
    fitsWriteto(hdulist, outfile, clobber=True)
    maskedCCD = MaskedCCD(temp_mask_file)
    for amp in maskedCCD:
        bright_pixels = BrightPixels(maskedCCD, amp,
                                     ccd.segments[amp].exptime,
                                     ccd.segments[amp].gain,
                                     ethresh=signal/2.,
                                     mask_plane=mask_plane)
        bright_pixels.generate_mask(outfile)

    os.remove(temp_mask_file)

