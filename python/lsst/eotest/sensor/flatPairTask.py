"""
@brief Compute detector response vs incident flux for pairs of flats.
These data are to be used for linearity and full-well measurments.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from DetectorResponse import DetectorResponse

import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def pair_mean(flat1, flat2, amp):
    """
    Compute the mean pixel value for a given segment, averaged between
    two exposures (assumed to have the same exposure time). The inputs
    flat1 and flat2 are MaskedCCD objects.  Any masks are assumed to
    have been enabled in the stat_ctrl attribute.
    """
    #
    # Remove bias using fit to overscan region and trim to imaging region.
    #
    im1 = flat1.unbiased_and_trimmed_image(amp)
    im2 = flat2.unbiased_and_trimmed_image(amp)
    #
    # Compute the mean DN of the images.
    #
    stats1 = afwMath.makeStatistics(im1, afwMath.MEAN, flat1.stat_ctrl)
    stats2 = afwMath.makeStatistics(im2, afwMath.MEAN, flat2.stat_ctrl)

    avg_mean_value = (stats1.getValue(afwMath.MEAN) +
                      stats1.getValue(afwMath.MEAN))/2.
    return avg_mean_value

def find_flat2(flat1):
    pattern = flat1.split('flat1')[0] + 'flat2*.fits'
    flat2 = glob.glob(pattern)[0]
    return flat2

class FlatPairConfig(pexConfig.Config):
    """Configuration for flat pair task"""
    output_dir = pexConfig.Field("Output directory", str, default='.')
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class FlatPairTask(pipeBase.Task):
    """Task to compute detector response vs incident flux for pairs of flats."""
    ConfigClass = FlatPairConfig
    _DefaultName = "FlatPairTask"
    
    @pipeBase.timeMethod
    def run(self, sensor_id, infiles, mask_files, gains, detrespfile=None):
        self.sensor_id = sensor_id
        self.infiles = infiles
        self.mask_files = mask_files
        self.gains = gains
        if detrespfile is None:
            #
            # Compute detector response from flat pair files.
            #
            detrespfile = self.extract_det_response()
        #
        # Perform full well and linearity analyses.
        #
        detresp = DetectorResponse(detrespfile)

        outfile = os.path.join(self.config.output_dir,
                               '%s_flat_pair_results.txt' % self.sensor_id)
        output = open(outfile, 'w')
        if self.config.verbose:
            self.log.info("Segment    full well (e-/pixel)   max. frac. dev.")
        for amp in imutils.allAmps:
            try:
                full_well, fp = detresp.full_well(amp)
            except RuntimeError:
                full_well, fp = detresp.full_well(amp, frac_offset=0.05)
            maxdev, fit_pars = detresp.linearity(amp)
            if self.config.verbose:
                self.log.info('%s            %.1f             %12.4e' 
                              % (imutils.channelIds[amp], full_well, maxdev))
            output.write("%i  %12.4e  %12.4e  %12.4e  %12.4e\n" 
                         % (amp, full_well, maxdev, fit_pars[0], fit_pars[1]))
    def extract_det_response(self):
        outfile = os.path.join(self.config.output_dir, 
                               '%s_det_response.txt' % self.sensor_id)
        file1s = sorted([item for item in self.infiles 
                         if item.find('flat1')  != -1])
        if self.config.verbose:
            self.log.info("writing to %s" % outfile)
        output = open(outfile, 'w')
        for file1 in file1s:
            if self.config.verbose:
                self.log.info("processing %s" % file1)
            file2 = find_flat2(file1)
    
            flat1 = MaskedCCD(file1, mask_files=self.mask_files)
            flat2 = MaskedCCD(file2, mask_files=self.mask_files)

            if flat1.md.get('EXPTIME') != flat2.md.get('EXPTIME'):
                raise RuntimeError("Exposure times do not match for:\n%s\n%s\n"
                                   % (file1, file2))
    
            flux = abs(flat1.md.get('EXPTIME')*flat1.md.get('MONDIODE') +
                       flat2.md.get('EXPTIME')*flat2.md.get('MONDIODE'))/2.

            output.write('%12.4e' % flux)
            #
            # Convert to e- and write out for each segment.
            #
            for amp in imutils.allAmps:
                output.write('  %12.4e' % (pair_mean(flat1, flat2, amp)
                                           *self.gains[amp]))
            output.write('\n')
            output.flush()
        output.close()
        return outfile

if __name__ == '__main__':
    parser = TaskParser('Compute detector response vs incident flux')
    parser.add_argument('-f', '--flats', type=str,
                        help='flat pairs file pattern')
    parser.add_argument('-F', '--flats_file_list', type=str,
                        help='list of flat pairs')
    args = parser.parse_args()
    sensor = args.sensor()
    sensor_id = args.sensor_id
