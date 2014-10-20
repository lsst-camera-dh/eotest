#!/usr/bin/env python
"""
@brief Translator application to convert RTS2 sensor FITS image files to 
conforming FITS files for analysis by the eotest package.
"""
import os
import sys
import numpy as np
import pyfits
from lsst.eotest.pyfitsTools import pyfitsWriteto
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest

class RTS2_FITS_translator(object):
    def __init__(self, luts, geom, verbose=True, args=None):
        self.luts = luts
        self.verbose = verbose
        self.args = args
        amp_loc = sensorTest.amp_loc[geom['vendor']]
        self.geom = sensorTest.AmplifierGeometry(prescan=geom['prescan'],
                                                 nx=geom['nx'], ny=geom['ny'],
                                                 amp_loc=amp_loc)
        self._sequence_number = -1
    def __call__(self, infile, outfile, clobber=True):
        if self.verbose:
            print "processing", infile

        # Recompute the amplifier geometry using NAXIS[12] from
        # first image extension of input file.
        self.geom.compute_geometry(fitsfile=infile)
        
        self.input = pyfits.open(infile)
        self.output = pyfits.open(infile, do_not_scale_image_data=True)
        prototypes = sensorTest.fits_headers.fits_headers()

        # Primary HDU
        try:
            self.output[0].header.set('MJD', self.output[0].header['JD']-2400000.5)
        except KeyError:
            from astropy.time import Time
            t=Time(self.output[0].header['DATE-OBS'])
            self.output[0].header.set('MJD', t.jd-2400000.5)
        self.output[0].header.set('FILENAME', os.path.basename(infile))
        self.output[0].header.set('SEQNUM', self._seqnum(infile))
        self.output[0].header.set('HEADVER', self._headver())
        if self.args is not None:
            self.output[0].header.set('LSST_NUM', self.args.sensor_id)
            self.output[0].header.set('CCD_MANU', self.args.vendor)
        self._update_keywords(0)

        # TEST_COND and CCD_COND extensions
        self._update_extension('TEST_COND', prototypes['TEST_COND'])
        self._update_extension('CCD_COND', prototypes['CCD_COND'])

        # Image extensions
        self._update_amp_keywords()

        # Special handling for monitoring diode current in BNL data.
        try:
            self.luts[0]['MONDIODE']
        except KeyError:
            self._set_bnl_mondiode_keyword_value()

        pyfitsWriteto(self.output, outfile, clobber=clobber, checksum=True,
                      output_verify='fix')
    def _seqnum(self, infile):
        """This assumes the sequence number is the penultimate token
        in the base filename when split by the '_' delimiter."""
        try:
            tokens = os.path.basename(infile).split('_')
            return tokens[-2]
        except IndexError:
            self._sequence_number += 1
            return "%04i" % self._sequence_number
    def _headver(self):
        """Increment the header version by 1. If it is missing, assume
        there has only been one version."""
        try:
            headver = int(self.output[0].header['HEADVER']) + 1
        except KeyError:
            headver = 2
        return headver
    def _update_keywords(self, ext):
        unresolved_keywords = []
        for key, source in self.luts[ext].items():
            try:
                value = self.input[0].header[source]
            except KeyError:
                unresolved_keywords.append(source)
                value = '' # write an empty string for the missing value
            self.output[ext].header.set(key, value)
        if unresolved_keywords and self.verbose:
            sys.stdout.write("HDU %s: " % ext)
            print "unresolved keywords in source primary hdu:"
            for item in unresolved_keywords:
                print "  %s" % item
    def _update_extension(self, extname, prototype):
        try:
            self.output[extname]
        except KeyError:
            # No header by that name, so add it along with required keys.
            self.output.append(pyfits.ImageHDU())
            self.output[-1].name = extname
            for keyword in prototype:
                self.output[-1].header.set(keyword, prototype[keyword])

        # Set the values from the primary hdu.
        self._update_keywords(extname)
    def _update_amp_keywords(self):
        self.output[0].header.set('DETSIZE', self.geom.DETSIZE)
        for amp in imutils.allAmps:
            self.output[amp].header.set('EXTNAME', imutils.hdu_dict[amp])
            for key in self.geom[amp].keys():
                self.output[amp].header.set(key, self.geom[amp][key])
    def _set_bnl_mondiode_keyword_value(self):
        try:
            self.input['AMP0.MEAS_TIMES']
        except KeyError:
            try:
                self.input['AMP1.MEAS_TIMES']
            except KeyError:
                # Extension does not exist so assume this is not a flat and
                # set to empty string.
                self.output[0].header.set('MONDIODE', '')
                return
        try:
            mean, stdev = self._bnl_mondiode_current()
        except (RuntimeError, ValueError):
            # Transitions in mondiode current data were not detected,
            # so a reliable current value cannot be computed using the
            # present algorithm.  Set the keyword value to an empty
            # string so that downstream tasks will fail if they try to
            # use this value in a calculation.
            mean = ''
        if self.verbose:
            print "Setting MONDIODE to", mean
        self.output[0].header.set('MONDIODE', mean)
    def _bnl_mondiode_current(self):
        try:
            data = self.input['AMP0.MEAS_TIMES'].data
        except KeyError:
            data = self.input['AMP1.MEAS_TIMES'].data
        y_pA, x_t = data.field('AMP0_A_CURRENT'), data.field('AMP0_MEAS_TIMES')
        # The following code has been lifted directly from JohnK's
        # xlatfits.py script at http://git.kuzew.net/lsst/xlatfits.git/
        i = 0;
        cpnts = [];
        downflg = 0;
        upflg   = 0;

        #normalize data
        norm = y_pA/np.max(np.abs(y_pA));

        while (i < len(norm) - 2):
        #check thresholds 
            if (norm[i] <= -0.8 and downflg == 0):
                #make sure it's trending properly
                if (np.sum(np.diff(y_pA[i-2:i+2])) < 0.0):
                    downflg = 1;
                    #print "Found transition at t=", x_t[i],  y_pA[i], i
                    cpnts.append(i);
            elif (norm[i] >= -0.8 and upflg == 0 and downflg == 1):
                if (np.sum(np.diff(y_pA[i-2:i+2])) > 0.0):
                    upflg = 1;
                    #print "Found transition at t=", x_t[i],  y_pA[i], i 
                    cpnts.append(i);
                    break;
            i += 1;

        if (len(cpnts) > 0):
            x1 = cpnts[0];
            if (upflg == 0):
                x2 = len(y_pA) - 1;
            else:
                x2 = cpnts[len(cpnts)-1];

            # Convert from pA to nA
            return np.mean(y_pA[x1:x2])/1e3, np.std(y_pA[x1:x2])/1e3

        raise RuntimeError("Could not compute monitoring photodiode current")

if __name__ == '__main__':
    import sys
    import glob
    import argparse

    parser = argparse.ArgumentParser(description='RTS2 FITS file translator')
    parser.add_argument('inputs', help="File pattern for input files")
    parser.add_argument('-o', '--output_dir', type=str, default='.',
                        help='output directory')
    parser.add_argument('-s', '--sensor_id', type=str, help='sensor id')
    parser.add_argument('-V', '--vendor', type=str, default='E2V',
                        help='Vendor (E2V or ITL)')
    parser.add_argument('-l', '--lab', type=str, default='BNL',
                        help='lab (BNL or Harvard)')
    parser.add_argument('-p', '--policy', type=str, default=None,
                        help='policy file for mapping lab-specific RTS2 keywords to eotest keywords')
    parser.add_argument('-v', '--verbose', action='store_true',
                        default=False, help='verbosity flag')
    
    args = parser.parse_args()

    if args.policy is None:
        sys.path.append(os.path.join(os.environ['EOTEST_DIR'], 'policy'))
        from RTS2_FITS_LUTs import *
    else:
        sys.path.insert(0, os.path.split(args.policy)[0])
        exec("from %s import *" % os.path.basename(args.policy).strip('.py'))

    rts2_translator = RTS2_FITS_translator(RTS2_FITS_LUTs[args.lab],
                                           sensor_geom[args.vendor],
                                           verbose=args.verbose,
                                           args=args)
    infiles = glob.glob(args.inputs)
    for infile in infiles:
        outfile = os.path.join(args.output_dir, os.path.basename(infile))
        rts2_translator(infile, outfile)
