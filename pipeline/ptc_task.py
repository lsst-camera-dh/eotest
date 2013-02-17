"""
@brief For pairs of flats obtain for a range of exposures, compute the
photon transfer curve (using pair_stats.py) and write as an output
file for use by later tasks, such as full_well_task.py,
linearity_task.py, etc..

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import pyfits
from pair_stats import pair_stats
import lsst.afw.math as afwMath

mean = lambda x : afwMath.makeStatistics(x, afwMath.MEAN).getValue()

exptime = lambda x : pyfits.open(x)[0].header['EXPTIME']

def glob_flats(full_path, outfile='ptc_flats.txt'):
    flats = glob.glob(os.path.join(full_path, '*_flat?.fits'))
    output = open(outfile, 'w')
    for item in flats:
        output.write('%s\n' % item)
    output.close()

def find_flats_from_file(infile):
    flat1s, flat2s = [], []
    for line in open(infile):
        filename = line.strip()
        if filename[-11:] == '_flat1.fits':
            flat1s.append(filename)
        elif filename[-11:] == '_flat2.fits':
            flat2s.append(filename)
    flat1s.sort()
    flats = []
    for flat1 in flat1s:
        flat2 = flat1.replace('_flat1.fits', '_flat2.fits')
        if flat2 in flat2s:
            flats.append((flat1, flat2))
    return flats

def find_flats(full_path, outfile='ptc_flats.txt'):
    glob_flats(full_path, outfile=outfile)
    return find_flats_from_file(outfile)

def accumulate_stats(flats, outfile='ptc_results.txt', verbose=True):
    """Run pair_stats.py to find mean and variance (in units of DN) as a
    function of exposure time."""
    output = open(outfile, 'w')
    for file1, file2 in flats:
        if verbose:
            print "processing", file1
        exposure = exptime(file1)
        output.write('%12.4e' % exposure)
        for hdu in range(16):
            results, b1, b2 = pair_stats(file1, file2, hdu+2)
            output.write('  %12.4e  %12.4e'%(results.flat_mean,
                                             results.flat_var))
        output.write('\n')
        output.flush()
    output.close()

if __name__ == '__main__':
    import sys
    from full_well import full_well
    from database.SensorDb import SensorDb, NullDbObject
    
    try:
        vendor = os.environ['CCD_VENDOR']
    except KeyError:
        vendor = 'e2v'

    flat_list = 'ptc_flats.txt'

    if len(sys.argv) == 4:
        full_path = sys.argv[1]
        ptcfile = sys.argv[2]
        sensor_id = sys.argv[3]
        glob_flats(full_path, outfile=flat_list)
        pipeline_task = False
    else:
        #
        # Pipeline input
        #
        flat_list = os.environ['PTC_FLAT_LIST']
        ptcfile = os.environ['PTC_OUTFILE']
        sensor_id = os.environ['SENSOR_ID']
        pipeline_task = True

    if pipeline_task:
        sensorDb = SensorDb(os.environ["DB_CREDENTIALS"])
        sensor = sensorDb.getSensor(vendor, sensor_id)
    else:
        sensor = NullDbObject(vendor, sensor_id)

    flats = find_flats_from_file(flat_list)
    accumulate_stats(flats, outfile=ptcfile)

    #
    # Full well calculations.
    #
    gains = SensorGains(vendor=vendor, vendorId=sensor_id)
    full_well_values = []
    for segment in range(16):
        full_well_est = full_well(ptcfile, segment, gain=)
        full_well_values.append(full_well_est)
        print '%02o  %i' % (segment, full_well_est)
        sensor.add_seg_result(segment, 'fullWell', full_well_est)
    sensor.add_ccd_result('fullWellMean', mean(full_well_values))
