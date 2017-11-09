"""
@brief Compute linearity by fitting mean(DN)/lamp current vs exposure time.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
from builtins import zip
import os
import sys
import numpy as np
import glob
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

import image_utils as imUtils
import pipeline.pipeline_utils as pipeUtils


def glob_flats(pattern, outfile='flats.txt'):
    flats = glob.glob(pattern)
    flats.sort()
    output = open(outfile, 'w')
    for item in flats:
        output.write('%s\n' % item)
    output.close()


class MeanSignal(imUtils.SubRegionSampler):
    def __init__(self, dx=100, dy=100, nsamp=1000, imaging=imUtils.imaging):
        imUtils.SubRegionSampler.__init__(self, dx, dy, nsamp, imaging)

    def __call__(self, infile):
        signals = {}
        for amp in imUtils.allAmps:
            samples = []
#            im = imUtils.unbias_and_trim(afwImage.ImageF(infile,
#                                                         imUtils.dm_hdu(amp)))
            im = imUtils.trim(afwImage.ImageF(infile, imUtils.dm_hdu(amp)))
            for x, y in zip(self.xarr, self.yarr):
                subim = self.subim(im, x, y)
                samples.append(imUtils.mean(subim))
            signals[amp] = imUtils.median(samples)
        return signals


def compute_mean_signal(flat_list, outfile='linearity_results.txt',
                        verbose=True):
    signal_estimator = MeanSignal()
    output = open(outfile, 'w')
    for flat in open(flat_list):
        infile = flat.strip()
        if infile[-11:] != '_flat1.fits':
            continue
        if verbose:
            print("processing", infile)
        md = afwImage.readMetadata(infile, 1)
        exptime = md.get('EXPTIME')
        kphot = np.abs(md.get('MONDIODE'))
        output.write('%12.4e  %12.4e' % (exptime, kphot))
        mean_signals = signal_estimator(infile)
        for amp in imUtils.allAmps:
            output.write('  %12.4e' % (mean_signals[amp]*gains[amp]))
        output.write('\n')
        output.flush()
    output.close()


if __name__ == '__main__':
    if len(sys.argv) >= 3:
        flat_pattern = sys.argv[1].replace('\\', '')
        outputdir = sys.argv[2]
        flat_list = 'linearity_flats.txt'
        glob_flats(flat_pattern, outfile=flat_list)
    else:
        try:
            sensor_id = os.environ['SENSOR_ID']
            flat_list = '%s_FLAT.txt' % sensor_id
            print(flat_list)
            outputdir = os.environ['OUTPUTDIR']
        except KeyError:
            print("usage: python linearity_task.py <flats pattern> <output directory> [<gains>=5.5]")
            sys.exit(1)

    gains, sensor = pipeUtils.setup(sys.argv, 3)

    try:
        os.makedirs(outputdir)
    except OSError:
        pass

    linearity_file = '%s_linearity_points.txt' % sensor_id.replace('-', '_')
    linearity_file = os.path.join(outputdir, linearity_file)

    compute_mean_signal(flat_list, outfile=linearity_file)
    #
    # Read in the data and fit
    #
    data = np.recfromtxt(linearity_file).transpose()
    exposure = data[0]
    lamp_current = data[1]
    print("Segment    max. frac. deviation")
    maxdevs = []
    for amp in imUtils.allAmps:
        indx = np.where((data[amp+1] > 100.) & (data[amp+1] < 9e4))
        signal = data[amp+1]/lamp_current
        results = np.polyfit(exposure[indx], signal[indx], 1)
        sensor.add_seg_result(amp, 'linefit_Slope', results[0])
        sensor.add_seg_result(amp, 'linefit_Intercept', results[1])
        f = np.poly1d(results)
        fvals = f(exposure[indx])
        maxDeviation = max(np.abs(fvals - signal[indx])/fvals)
        maxdevs.append(maxDeviation)
        sensor.add_seg_result(amp, 'maxDeviation', maxDeviation)
        print("%s         %.4f" % (imUtils.channelIds[amp], maxDeviation))
    sensor.add_ccd_result('maxDeviation', max(maxdevs))
