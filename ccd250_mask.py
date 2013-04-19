"""
@brief Script to create a mask for the CCD250 defects, specifically
the edge roll-off effects around the perimeter of the sensors and the
effects of the midline blooming stop implant.  These are documented by
Peter Doherty's "Status of sensor characterization" reports linked in
to the Feb 8, 2013 agenda for the Sensor Testing Meetings,
https://confluence.slac.stanford.edu/x/DQvNBw

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import pyfits
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import image_utils as imutils
from BrightPixels import BrightPixels
from simulation.sim_tools import CCD

def ccd250_mask(outfile, mask_plane='CCD250_DEFECTS',
                tmp_mask_image='temp_mask_image.fits',
                outer_edge_width=10, bloom_stop_width=5, signal=10,
                cleanup=True):
    """
    This function creates a file containing masks for each segment to
    mask pixels affected the edge roll-off and midline blooming stop
    distortions.

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
    gain = 1
    ccd = CCD(exptime=1, gain=gain)
    #
    # Write the output file with a primary HDU so that the DMstack code
    # can append only image extensions (and not write to the PHDU).
    #
    hdulist = pyfits.HDUList()
    hdulist.append(pyfits.PrimaryHDU())
    hdulist.writeto(outfile, clobber=True)
    #
    # Amplifiers 1 (AMP10), 8 (AMP17), 9 (AMP07) and 16 (AMP00) are
    # along the perimeter.
    #
    # Amps 8 & 16 have roll-off adjacent to the prescan:
    # 
    xmin = imutils.imaging.getMinX()
    xmax = xmin + outer_edge_width
    for amp in (8, 16):
        imarr = ccd.segments[amp].image.getArray()
        imarr[:, xmin:xmax] += signal
    #
    # Amps 1 & 9 have roll-off adjacent to the serial overscan:
    #
    xmax = imutils.imaging.getMaxX() + 1
    xmin = xmax - outer_edge_width
    for amp in (1, 9):
        imarr = ccd.segments[amp].image.getArray()
        imarr[:, xmin:xmax] += signal
    #
    # Loop over all amps, set signal in perimeter and around blooming
    # stop.
    #
    ymax = imutils.imaging.getMaxY() + 1
    ymin = ymax - bloom_stop_width
    for i, amp in enumerate(imutils.allAmps):
        image = ccd.segments[amp].image
        imarr = image.getArray()
        #
        # Set signal in row direction along perimeter.
        #
        imarr[0:outer_edge_width, :] += signal
        #
        # Set signal around blooming stop
        #
        imarr[ymin:ymax, :] += signal
    #
    # Write the images of the mask regions to the FITS file.
    #
    ccd.writeto(tmp_mask_image)
    #
    # Use BrightPixels code to detect the mask regions and write the mask file.
    #
    mask = afwImage.MaskU(image.getDimensions())
    mask.addMaskPlane(mask_plane)
    for amp in imutils.allAmps:
        bright_pixels = BrightPixels(tmp_mask_image, amp)
        pixels, columns = bright_pixels.find(gain, ethresh=signal/2.)
        bright_pixels.write_mask(outfile, mask_plane)
        
    if cleanup:
        os.remove(tmp_mask_image)
        
if __name__ == '__main__':
    mask_file = 'CCD250_DEFECTS_mask.fits'
    ccd250_mask(mask_file, cleanup=False)
