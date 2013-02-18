"""
@brief Compute linearity by fitting mean(DN)/lamp current vs exposure time.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import glob
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

from image_utils import dm_hdu, trim, unbias_and_trim, imaging, allAmps, \
    SubRegionSampler, mean, median, channelIds

def glob_flats(pattern, outfile='flats.txt'):
    flats = glob.glob(pattern)
    flats.sort()
    output = open(outfile, 'w')
    for item in flats:
        output.write('%s\n' % item)
    output.close()

class MeanSignal(SubRegionSampler):
    def __init__(self, dx=100, dy=100, nsamp=1000, imaging=imaging):
        SubRegionSampler.__init__(self, dx, dy, nsamp, imaging)
    def __call__(self, infile):
        signals = {}
        for amp in allAmps:
            samples = []
#            im = unbias_and_trim(afwImage.ImageF(infile, dm_hdu(amp)))
            im = trim(afwImage.ImageF(infile, dm_hdu(amp)))
            for x, y in zip(self.xarr, self.yarr):
                bbox = afwGeom.Box2I(afwGeom.Point2I(int(x), int(y)),
                                     afwGeom.Extent2I(self.dx, self.dy))
                subim = im.Factory(im, bbox)
                samples.append(mean(subim))
            signals[amp] = median(samples)
        return signals

def compute_mean_signal(flat_list, outfile='linearity_results.txt',
                        verbose=True):
    signal_estimator = MeanSignal()
    output = open(outfile, 'w')
    for flat in open(flat_list):
        infile = flat.strip()
        if verbose:
            print "processing", infile
        md = afwImage.readMetadata(infile, 1)
        exptime = md.get('EXPTIME')
        kphot = np.abs(md.get('K_PHOT_CURRENT'))
        output.write('%12.4e  %12.4e' % (exptime, kphot))
        mean_signals = signal_estimator(infile)
        for amp in allAmps:
            output.write('  %12.4e' % (mean_signals[amp]*gains[amp]))
        output.write('\n')
        output.flush()
    output.close()

if __name__ == '__main__':
    import os
    import sys
    from database.SensorDb import SensorDb, NullDbObject
    from database.SensorGains import SensorGains

    if len(sys.argv) >= 3:
        flat_pattern = sys.argv[1].replace('\\', '')
        linearity_file = sys.argv[2]
        try:
            gains = SensorGains(float(sys.argv[3]))
        except IndexError:
            print "Setting system gain to 5.5 e-/DN for all segments."
            gains = SensorGains(5.5)
        sensor = NullDbObject()
        flat_list = 'linearity_flats.txt'
        glob_flats(flat_pattern, outfile=flat_list)
    else:
        try:
            flat_list = os.environ['FLAT_LIST']
            linearity_file = os.environ['LINEARITY_OUTFILE']
            sensor_id = os.environ['SENSOR_ID']
            vendor = os.environ['CCD_VENDOR']
            sensorDb = SensorDb(os.environ['DB_CREDENTIALS'])
            sensor = sensorDb.getSensor(vendor, sensor_id)
            gains = SensorGains(vendor=vendor, vendorId=sensor_id)
        except KeyError:
            print "usage: python linearity_task.py <flats pattern> <linearity output file> [<gains>=5.5]"
            sys.exit(1)
    
    compute_mean_signal(flat_list, outfile=linearity_file)
    #
    # Read in the data and fit
    #
    data = np.recfromtxt(linearity_file).transpose()
    exposure = data[0]
    lamp_current = data[1]
    print "Segment    max. frac. deviation"
    for amp in allAmps:
        indx = np.where((data[amp+1] > 100.) & (data[amp+1] < 9e4))
        signal = data[amp+1]/lamp_current
        results = np.polyfit(exposure[indx], signal[indx], 1)
        sensor.add_seg_result(amp, 'ptcSlope', results[0])
        sensor.add_seg_result(amp, 'ptcIntercept', results[1])
        f = np.poly1d(results)
        fvals = f(exposure)
        maxDeviation = max(np.abs(fvals - signal)/fvals)
        sensor.add_seg_result(amp, 'maxDeviation', maxDeviation)
        print "%s         %.2f" % (channelIds[amp], maxDeviation)
