"""
@brief For pairs of flats obtain for a range of exposures, compute the
photon transfer curve (using pair_stats.py) and compute and write out
the full well.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import lsst.eotest.image_utils as imutils
from pair_stats import pair_stats

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def find_flat2(flat1):
    pattern = flat1.split('flat1')[0] + 'flat2*.fits'
    flat2 = glob.glob(pattern)[0]
    return flat2

exptime = lambda x : afwImage.readMetadata(x, 1).get('EXPTIME')

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
    def run(self, sensor_id, infiles, mask_files, gains, binsize=1):
        outfile = os.path.join(self.config.output_dir, '%s_ptc.txt' % sensor_id)
        output = open(outfile, 'w')
        file1s = sorted([item for item in infiles if item.find('flat1')  != -1])
        for flat1 in file1s:
            flat2 = flat1.replace('flat1', 'flat2')
            if self.config.verbose:
                self.log.info("processing %s" % flat1)
            exposure = exptime(flat1)
            output.write('%12.4e' % exposure)
            for amp in imutils.allAmps:
                results = pair_stats(flat1, flat2, amp, mask_files=mask_files)
                output.write('  %12.4e  %12.4e' % (results.flat_mean,
                                                   results.flat_var))
            output.write('\n')
            output.flush()
        output.close()
