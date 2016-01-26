"""
@brief For pairs of flats obtain for a range of exposures, compute the
photon transfer curve and compute and write out the full well.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import numpy as np
import pyfits
from lsst.eotest.pyfitsTools import pyfitsTableFactory, pyfitsWriteto
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def find_flat2(flat1):
    pattern = flat1.split('flat1')[0] + 'flat2*.fits'
    try:
        flat2 = glob.glob(pattern)[0]
        return flat2
    except IndexError:
        return flat1

exptime = lambda x : imutils.Metadata(x, 1).get('EXPTIME')

def glob_flats(full_path, outfile='ptc_flats.txt'):
    flats = glob.glob(os.path.join(full_path, '*_flat?.fits'))
    output = open(outfile, 'w')
    for item in flats:
        output.write('%s\n' % item)
    output.close()

def find_flats(args):
    files = args.files(args.flats, args.flats_file_list)
    file1s = sorted([item.strip() for item in files
                     if item.find('flat1') != -1])
    return [(f1, find_flat2(f1)) for f1 in file1s]

class FlatPairStats(object):
    def __init__(self, fmean, fvar):
        self.flat_mean = fmean
        self.flat_var = fvar

def flat_pair_stats(ccd1, ccd2, amp, mask_files=(), bias_frame=None):
    if ccd1.md.get('EXPTIME') != ccd2.md.get('EXPTIME'):
        raise RuntimeError("Exposure times for files %s, %s do not match"
                           % (ccd1.imfile, ccd2.imfile))
    #
    # Mean and variance calculations that account for masks (via
    # ccd1.stat_ctrl, which is the same for both MaskedImages).
    #
    mean = lambda im : afwMath.makeStatistics(im, afwMath.MEAN,
                                              ccd1.stat_ctrl).getValue()
    var = lambda im : afwMath.makeStatistics(im, afwMath.VARIANCE,
                                             ccd1.stat_ctrl).getValue()
    #
    # Extract imaging region for segments of both CCDs.
    #
    image1 = ccd1.unbiased_and_trimmed_image(amp)
    image2 = ccd2.unbiased_and_trimmed_image(amp)
    #
    # Use serial overscan for bias region.
    #
    b1 = ccd1[amp].Factory(ccd1[amp], ccd1.amp_geom.serial_overscan)
    b2 = ccd2[amp].Factory(ccd2[amp], ccd2.amp_geom.serial_overscan)
    if ccd1.imfile == ccd2.imfile:
        # Don't have pairs of flats, so estimate noise and gain
        # from a single frame, ignoring FPN.
        fmean = mean(image1)
        fvar = var(image1)
    else:
        #
        # Make a deep copy since otherwise the pixel values in image1
        # would be altered in the ratio calculation.
        #
        fratio_im = afwImage.MaskedImageF(image1, True)
        fratio_im /= image2
        fratio = mean(fratio_im)
        image2 *= fratio
        fmean = (mean(image1) + mean(image2))/2.

        fdiff = afwImage.MaskedImageF(image1, True)
        fdiff -= image2
        fvar = var(fdiff)/2.

    return FlatPairStats(fmean, fvar)

class PtcConfig(pexConfig.Config):
    """Configuration for ptc task"""
    output_dir = pexConfig.Field("Output directory", str, default='.')
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class PtcTask(pipeBase.Task):
    """Task to compute photon transfer curve from flat pair dataset"""
    ConfigClass = PtcConfig
    _DefaultName = "PtcTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, infiles, mask_files, gains, binsize=1, 
            bias_frame=None):
        outfile = os.path.join(self.config.output_dir, 
                               '%s_ptc.fits' % sensor_id)
        ptc_stats = dict([(amp, ([], [])) for amp in imutils.allAmps])
        exposure = []
        file1s = sorted([item for item in infiles if item.find('flat1')  != -1])
        for flat1 in file1s:
            flat2 = find_flat2(flat1)
            if self.config.verbose:
                self.log.info("processing %s" % flat1)
            exposure.append(exptime(flat1))
            ccd1 = MaskedCCD(flat1, mask_files=mask_files,
                             bias_frame=bias_frame)
            ccd2 = MaskedCCD(flat2, mask_files=mask_files,
                             bias_frame=bias_frame)
            for amp in imutils.allAmps:
                results = flat_pair_stats(ccd1, ccd2, amp,
                                          mask_files=mask_files,
                                          bias_frame=bias_frame)
                ptc_stats[amp][0].append(results.flat_mean)
                ptc_stats[amp][1].append(results.flat_var)
        output = pyfits.HDUList()
        output.append(pyfits.PrimaryHDU())
        colnames = ['EXPOSURE']
        units = ['seconds']
        columns = [np.array(exposure, dtype=np.float)]
        for amp in imutils.allAmps:
            colnames.extend(['AMP%02i_MEAN' % amp, 'AMP%02i_VAR' % amp])
            units.extend(['ADU', 'ADU**2'])
            columns.extend([np.array(ptc_stats[amp][0], dtype=np.float),
                            np.array(ptc_stats[amp][1], dtype=np.float)])
        formats = 'E'*len(colnames)
        fits_cols = [pyfits.Column(name=colnames[i], format=formats[i],
                                   unit=units[i], array=columns[i])
                     for i in range(len(columns))]
        output.append(pyfitsTableFactory(fits_cols))
        output[-1].name = 'PTC_STATS'
        pyfitsWriteto(output, outfile, clobber=True)
