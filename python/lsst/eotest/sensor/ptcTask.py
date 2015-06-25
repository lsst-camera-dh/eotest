"""
@brief For pairs of flats obtain for a range of exposures, compute the
photon transfer curve (using pair_stats.py) and compute and write out
the full well.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import numpy as np
import pyfits
from lsst.eotest.pyfitsTools import pyfitsTableFactory, pyfitsWriteto
import lsst.eotest.image_utils as imutils
from pair_stats import pair_stats

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
            for amp in imutils.allAmps:
                results = pair_stats(flat1, flat2, amp, mask_files=mask_files,
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
