#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import numpy as np
import astropy.io.fits as fits
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('System gain characterization from Fe55 catalog')
parser.add_argument('fe55_catalog', type=str,
                    help='catalog of Fe55 charge clusters')
parser.add_argument('-c', '--chiprob_min', type=float, default=0.1,
                    help='Mininum chi-square probability for cluster fit')
parser.add_argument('-p', '--plot', action='store_true', default=False,
                    help='Plot distribution and fit')

args = parser.parse_args()

catalog = fits.open(args.fe55_catalog)

gains = {}
for amp in imutils.allAmps:
    chiprob = catalog[amp].data.field('CHIPROB')
    index = np.where(chiprob > args.chiprob_min)
    dn = catalog[amp].data.field('DN')[index]
    if args.verbose:
        sys.stdout.write("amp = %s, len(DN) = %i, " % (amp, len(dn)))

    if len(dn) > 2:
        plot_filename = "Fe55_dist_%s_amp%02i.png" % (args.sensor_id, amp)
        try:
            gains[amp], peak, sigma = \
                sensorTest.fe55_gain_fitter(dn, make_plot=args.plot,
                                            title='Amp %i' % amp,
                                            plot_filename=plot_filename)
        except RuntimeError as e:
            print(e)
            continue

    if args.verbose:
        try:
            print("gain = %.2f" % (gains[amp],))
        except KeyError:
            print()

results_file = args.results_file
if results_file is None:
    results_file = os.path.join(args.output_dir,
                                '%s_eotest_results.fits' % args.sensor_id)

if args.verbose:
    print("Writing to %s" % results_file)

results = sensorTest.EOTestResults(results_file)
for amp in gains:
    results.add_seg_result(amp, 'GAIN', gains[amp])
results.write(overwrite=True)
