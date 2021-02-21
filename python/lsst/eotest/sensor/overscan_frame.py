"""
Function to generate an "overscan" image based on the parallel and
serial overscans to use for overscan correction.
"""
import numpy as np
from astropy.io import fits
import lsst.geom
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD


__all__ = ['make_overscan_frame']


def make_overscan_frame(fits_file, outfile=None, amps=None):
    """
    Use overscan regions to do column-wise, then row-wise overscan
    profile-based overscan frame.
    """
    if amps is None:
        amps = imutils.allAmps(fits_file)
    ccd = MaskedCCD(fits_file)
    par_corner = lsst.geom.Point2I(0, ccd.amp_geom.imaging.getHeight())
    par_extent = lsst.geom.Extent2I(ccd.amp_geom.full_segment.getWidth(),
                                    ccd.amp_geom.parallel_overscan.getHeight())
    par_bbox = lsst.geom.Box2I(par_corner, par_extent)
    ser_corners = ccd.amp_geom.serial_overscan.getCorners()
    ser_extent = lsst.geom.Extent2I(ccd.amp_geom.serial_overscan.getWidth(),
                                    ccd.amp_geom.full_segment.getHeight())
    ser_bbox = lsst.geom.Box2I(ser_corners[0], ser_extent)

    oscan_arrays = dict()
    for amp in ccd:
        image = ccd[amp].getImage()

        parallel = image.Factory(image, par_bbox).array
        p_profile = np.sum(parallel, 0)/parallel.shape[0]

        overscan_array = np.zeros(image.array.shape)

        image.array[:,] -= p_profile
        overscan_array[:,] += p_profile

        serial = image.Factory(image, ser_bbox).array
        s_profile = np.sum(serial, 1)/serial.shape[1]

        for col in range(ccd.amp_geom.full_segment.getWidth()):
            image.array[:serial.shape[0], col] -= s_profile
            overscan_array[:serial.shape[0], col] += s_profile
        oscan_arrays[amp] = overscan_array

    if outfile is not None:
        with fits.open(fits_file) as hdus:
            for amp in ccd:
                hdus[amp].data = oscan_arrays[amp]
            hdus.writeto(outfile, overwrite=True)

    return ccd
