"""
@brief Script to create a mask for the the edge roll-off effects
around the perimeter of the sensors and, in e2v devices, the effects
of the midline blooming stop implant.  These are documented by
Peter Doherty's "Status of sensor characterization" reports linked in
to the Feb 8, 2013 agenda for the Sensor Testing Meetings,
https://confluence.slac.stanford.edu/x/DQvNBw

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import warnings
import tempfile
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsWriteto
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .AmplifierGeometry import makeAmplifierGeometry, amp_loc
from .BrightPixels import BrightPixels
from .sim_tools import CCD
from .generate_mask import generate_mask


def rolloff_mask(infile, outfile,
                 mask_plane='ROLLOFF_DEFECTS',
                 tmp_mask_image=None,
                 outer_edge_width=10,
                 bloom_stop_width=5,
                 signal=10,
                 cleanup=True):
    """
    This function creates a file containing masks for each segment to
    mask pixels affected the edge roll-off and midline blooming stop
    distortions.

    infile: Input file to mask.

    outfile: The name of the file to contain the masks.

    outer_edge_width: This is the width in pixels of the masked region
    along the sides closest to the sensor edges of the imaging region
    of each segment.

    bloom_stop_width: The width in pixels of the masked region of the
    imaging region of each segment, on either side of the central
    blooming stop implant.

    signal: This is the artificial signal written to the mask image so
    that the afwDetect code can generate the footprints used define
    the masks.  It basically just needs to be some positive (non-zero)
    number.
    """
    #
    # Create an empty set of frames to fill with an image of the mask
    # from which we will generate the masks.  The exposure time and
    # system gain need to be set to unity so that the BrightPixels.py
    # code can interpret DN directly as e- per second.
    #
    amp_geom = makeAmplifierGeometry(infile)
    gain = 1
    exptime = 1
    all_amps = imutils.allAmps(infile)
    ccd = CCD(exptime=exptime, gain=gain, geometry=amp_geom, amps=all_amps)
    #
    # Write the output file with a primary HDU so that the DMstack code
    # can append only image extensions (and not write to the PHDU).
    #
    warnings.filterwarnings('ignore', category=fits.verify.VerifyWarning, append=True)
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    with fits.open(infile) as fd:
        hdulist[0].header.update(fd[0].header)
    # Use the mask_plane value ('ROLLOFF_DEFECTS') to distinguish
    # this file from other mask files.
    hdulist[0].header['MASKTYPE'] = mask_plane
    fitsWriteto(hdulist, outfile, overwrite=True)
    #
    # Amplifiers 1 (AMP10), 8 (AMP17), 9 (AMP07) and 16 (AMP00) are
    # along the perimeter, but have different edge rolloff pixels
    # depending on the vendor.
    #
    # These amps have edge roll-off adjacent to the prescan.
    amps = (8, 16) if amp_geom.vendor == 'E2V' else (8, 9)
    xmin = amp_geom.imaging.getMinX()
    xmax = xmin + outer_edge_width
    for amp in amps:
        if amp not in all_amps:
            continue
        imarr = ccd.segments[amp].image.getArray()
        imarr[:, xmin:xmax] += signal
    #
    # These amps have edge roll-off adjacent to the serial overscan:
    #
    amps = (1, 9) if amp_geom.vendor == 'E2V' else (1, 16)
    xmax = amp_geom.imaging.getMaxX() + 1
    xmin = xmax - outer_edge_width
    for amp in amps:
        if amp not in all_amps:
            continue
        imarr = ccd.segments[amp].image.getArray()
        imarr[:, xmin:xmax] += signal
    #
    # Loop over all amps, set signal in perimeter and around blooming
    # stop.
    #
    ymax = amp_geom.imaging.getMaxY() + 1
    ymin = ymax - bloom_stop_width
    for i, amp in enumerate(ccd.segments):
        image = ccd.segments[amp].image
        imarr = image.getArray()
        #
        # Set signal in row direction along perimeter.
        #
        xmin = amp_geom.imaging.getMinX()
        xmax = amp_geom.imaging.getMaxX()
        imarr[0:outer_edge_width, xmin:xmax] += signal

        if amp_geom.amp_loc == amp_loc['E2V']:
            #if True:
            #
            # Set signal around blooming stop
            #
            imarr[ymin:ymax, :] += signal
    #
    # Write the images of the mask regions to the FITS file.
    #
    if tmp_mask_image is None:
        tmp_mask_image = tempfile.mkstemp(suffix='.fits', dir='.')[-1]
    fitsWriteto(ccd, tmp_mask_image)
    #
    # Use BrightPixels code to detect the mask regions and write the mask file.
    #
    try:
        mask = afwImage.Mask(image.getDimensions())
    except AttributeError:
        mask = afwImage.MaskU(image.getDimensions())
    mask.addMaskPlane(mask_plane)
    maskedCCD = MaskedCCD(tmp_mask_image)
    pixels, columns = {}, {}
    for amp in maskedCCD:
        bright_pixels = BrightPixels(maskedCCD, amp, exptime, gain,
                                     ethresh=signal/2., mask_plane=mask_plane)
        pixels[amp], columns[amp] = bright_pixels.find()
    generate_mask(infile, outfile, mask_plane, pixels=pixels, columns=columns)
    if cleanup:
        os.remove(tmp_mask_image)


def pixel_counts(ccd_file, input_mask=None):
    """
    Based on the sensor geometry and an optional input mask, compute
    the total number of pixels in the imaging regions of the
    amplifiers and the total number of masked pixels within the
    imaging regions.  If no input mask is given, the standard rolloff
    mask for the vendor device is used.

    @return <total number of imaging region pixels>, <number of masked pixels>
    """
    if input_mask is not None:
        mask_file = input_mask
    else:
        mask_file = 'temp_rolloff_mask.fits'
        rolloff_mask(ccd_file, mask_file)
    ccd = MaskedCCD(mask_file)
    num_masked = 0
    num_total = 0
    imaging = ccd.amp_geom.imaging
    for amp in ccd:
        imarr = imutils.trim(ccd[amp].getImage(), imaging).getArray()
        num_masked += len(np.where(imarr != 0)[0])
        num_total += imarr.shape[0]*imarr.shape[1]
    if mask_file is None:
        try:
            os.remove(mask_file)
        except OSError:
            pass
    return num_total, num_masked
