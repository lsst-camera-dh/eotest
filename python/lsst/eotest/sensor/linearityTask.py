"""
@brief Compute detector response vs incident flux for pairs of flats.
These data are to be used for linearity measurments.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from EOTestResults import EOTestResults
from DetectorResponse import DetectorResponse

import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from flatPairTask import pair_mean, find_flat2

class LinearityConfig(pexConfig.Config):
    """Configuration for flat pair task"""
    output_dir = pexConfig.Field("Output directory", str, default='.')
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class LinearityTask(pipeBase.Task):
    """Task to compute detector response vs incident flux from
       flat pair dataset."""
    ConfigClass = LinearityConfig
    _DefaultName = "LinearityTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, infiles, mask_files, gains, detrespfile=None,
            bias_frame=None, median_stack=None):
        self.sensor_id = sensor_id
        self.infiles = infiles
        self.mask_files = mask_files
        self.gains = gains
        self.bias_frame = bias_frame
        self.median_stack = median_stack
        if detrespfile is None:
            #
            # Compute detector response from flat pair files.
            #
            detrespfile = self.extract_det_response()
        #
        # Perform linearity analyses.
        #
        detresp = DetectorResponse(detrespfile)

        outfile = self.config.eotest_results_file
        if outfile is None:
            outfile = os.path.join(self.config.output_dir,
                                   '%s_eotest_results.fits' % self.sensor_id)
        all_amps = imutils.allAmps(detrespfile)
        output = EOTestResults(outfile, namps=len(all_amps))
        if self.config.verbose:
            self.log.info("Amp        max. frac. dev.")
        for amp in all_amps:
            try:
                maxdev, fit_pars, Ne, flux = detresp.linearity(amp)
            except:
                maxdev = None
            if self.config.verbose:
                self.log.info('%2i            %s' % (amp, maxdev))
            if maxdev is not None:
                output.add_seg_result(amp, 'MAX_FRAC_DEV', float(maxdev))
        output.write()
    def _create_detresp_fits_output(self, nrows):
        self.output = fits.HDUList()
        self.output.append(fits.PrimaryHDU())
        all_amps = imutils.allAmps()
        colnames = ['flux'] + ['AMP%02i_SIGNAL' % i for i in all_amps]
        formats = 'E'*len(colnames)
        units = ['None'] + ['e-']*len(all_amps)
        columns = [np.zeros(nrows, dtype=np.float) for fmt in formats]
        fits_cols = [fits.Column(name=colnames[i], format=formats[i],
                                 unit=units[i], array=columns[i])
                     for i in range(len(units))]
        hdu = fitsTableFactory(fits_cols)
        hdu.name = 'DETECTOR_RESPONSE'
        self.output.append(hdu)
    def extract_det_response(self):
        outfile = os.path.join(self.config.output_dir,
                               '%s_det_response_linearity.fits' % self.sensor_id)
        file1s = sorted([item for item in self.infiles 
                         if item.find('flat1')  != -1 or 
                         item.find('linearity_flat') != -1])
        if self.config.verbose:
            self.log.info("writing to %s" % outfile)
        self._create_detresp_fits_output(len(file1s))
        for row, file1 in enumerate(file1s):
            if self.config.verbose:
                self.log.info("processing %s" % file1)
            try:
                file2 = find_flat2(file1)
            except IndexError:
                # Just use flat1 again since only average is taken and
                # FPN subtraction isn't needed.
                file2 = file1

            flat1 = MaskedCCD(file1, mask_files=self.mask_files,
                              bias_frame=self.bias_frame)
            flat2 = MaskedCCD(file2, mask_files=self.mask_files,
                              bias_frame=self.bias_frame)

            if flat1.md.get('EXPTIME') != flat2.md.get('EXPTIME'):
                raise RuntimeError("Exposure times do not match for:\n%s\n%s\n"
                                   % (file1, file2))

            if (flat1.md.get('MONDIODE') != 0 and
                flat2.md.get('MONDIODE') != 0):
                flux = abs(flat1.md.get('EXPTIME')*flat1.md.get('MONDIODE') +
                           flat2.md.get('EXPTIME')*flat2.md.get('MONDIODE'))/2.
            else:
                flux = flat1.md.get('EXPTIME')

            self.output[-1].data.field('FLUX')[row] = flux
            for amp in flat1:
                # Convert to e- and write out for each segment.
		signal = pair_mean(flat1, flat2, amp, self.median_stack)*self.gains[amp]
                self.output[-1].data.field('AMP%02i_SIGNAL' % amp)[row] = signal
        self.output[0].header['NAMPS'] = len(flat1)
        fitsWriteto(self.output, outfile, clobber=True)
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
