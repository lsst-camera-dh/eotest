"""
Function to generate mask files given a set of pixels and columns to mask.
"""
from __future__ import absolute_import, print_function
import os
import tempfile
import astropy.io.fits as fits
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
from lsst.eotest.fitsTools import fitsWriteto
import lsst.eotest.image_utils as imutils
from .AmplifierGeometry import makeAmplifierGeometry
from .MaskedCCD import MaskedCCD
from .sim_tools import CCD

def generate_mask(infile, outfile, mask_plane, pixels=None, columns=None,
                  temp_mask_image=None):
    """
    Generate a mask file for the specified pixels and columns.
    The amplifier geometry will be taken from infile.
    """
    # Insert artificial signal into specified pixels and columns for
    # each segment.
    exptime = 1
    gain = 1
    geometry = makeAmplifierGeometry(infile)
    ccd = CCD(exptime=exptime, gain=gain, geometry=geometry)
    # Account for prescan in x-coordinate since we want to mask
    # columns the full segment.
    prescan = geometry.prescan.getWidth()

    signal = 10
    if pixels is None:
        pixels = {}
    for amp in pixels:
        imarr = ccd.segments[amp].image.getArray()
        for ix, iy in pixels[amp]:
            imarr[iy][ix + prescan] = signal

    if columns is None:
        columns = {}
    for amp in columns:
        imarr = ccd.segments[amp].image.getArray()
        for ix in columns[amp]:
            imarr[:, ix + prescan] = signal

    if temp_mask_image is None:
        temp_mask_image = tempfile.mkstemp(suffix='.fits', dir='.')[-1]
    fitsWriteto(ccd, temp_mask_image)

    # Use the afw code to create a mask file.
    afwImage.MaskU.addMaskPlane(mask_plane)
    # Create a primary HDU with the input file primary HDU metadata,
    # updated with the mask type info.
    hdulist = fits.HDUList()
    hdulist.append(fits.open(infile)[0])
    hdulist[0].header['MASKTYPE'] = mask_plane
    fitsWriteto(hdulist, outfile, clobber=True)
    # Loop over segments in the temporary file and add all pixels to
    # the mask with the inserted signal.
    maskedCCD = MaskedCCD(temp_mask_image)
    for amp in maskedCCD:
        threshold = afwDetect.Threshold(signal/2.*exptime/gain)
        fp_set = afwDetect.FootprintSet(maskedCCD[amp], threshold)
        mask = afwImage.MaskU(maskedCCD[amp].getDimensions())
        fp_set.setMask(mask, mask_plane)
        md = dafBase.PropertySet()
        md.set('EXTNAME', 'SEGMENT%s' % imutils.channelIds[amp])
        md.set('DETSIZE', geometry[amp]['DETSIZE'])
        md.set('DETSEC', geometry[amp]['DETSEC'])
        md.set('DATASEC', geometry[amp]['DATASEC'])
        mask.writeFits(outfile, md, 'a')
    os.remove(temp_mask_image)

