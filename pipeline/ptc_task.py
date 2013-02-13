import os
import glob
import pyfits
from pair_stats import pair_stats

exptime = lambda x : pyfits.open(x)[0].header['EXPTIME']

sensor_id = '000-00'
datadir = '/nfs/farm/g/lsst/u1/testData/SIMData/%s/flats' % sensor_id

#files = glob.glob(os.path.join(datadir, 'flat_*_0.fits'))
#files.sort()

files = glob.glob('/nfs/farm/g/lsst/u1/testData/HarvardData/112-02/ptc/logain/112_02_ptc_logain_*_flat1.fits')
files.sort()

output = open('ptc_results.txt', 'w')
for file1 in files:
    print file1
    file2 = file1.replace('_flat1.fits', '_flat2.fits')
    if not os.path.isfile(file2):
        continue
    exposure = exptime(file1)
    output.write('%12.4e' % exposure)
    for hdu in range(16):
        results, b1, b2 = pair_stats(file1, file2, hdu+2)
        output.write('  %12.4e  %12.4e'%(results.flat_mean, results.flat_var))
    output.write('\n')
    output.flush()
output.close()
