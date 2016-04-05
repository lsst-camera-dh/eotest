"""
Script to translate vendor data from e2v to operationally compliant
FITS files for analysis with eotest.
"""
import os
import glob
import subprocess
import astropy.io.fits as fits
import lsst.eotest.sensor as sensorTest

#
# test_types: fe55 dark flat lambda trap sflat_nnn spot
# image_types: bias dark fe55
# filenames: <sensor_id>_<test_type>_<image_type>_<seqno>_<time_stamp>.fits
#

def translate(infile, test_type, image_type, seqno, time_stamp='000'):
    pass
    foo = fits.open(infile)
    detxsize = 8*foo[1].header['NAXIS1']
    detysize = 2*foo[1].header['NAXIS2']
    ampGeom = sensorTest.AmplifierGeometry(detxsize=detxsize, detysize=detysize)
    sensor_id = foo[0].header['DEV_ID']
    exptime = foo[0].header['EXPOSURE']
    foo[0].header['EXPTIME'] = exptime
    foo[0].header['MONOWL'] = foo[0].header['WAVELEN']
    foo[0].header['MONDIODE'] = foo[0].header['LIGHTPOW']
    foo[0].header['CCDTEMP'] = foo[0].header['TEMP_MEA']
    foo[0].header['DETSIZE'] = ampGeom.DETSIZE
    for hdu in range(1, 17):
        amp = foo[hdu].header['AMPNO']
        foo[hdu].header['DETSIZE'] = ampGeom[amp]['DETSIZE']
        foo[hdu].header['DATASEC'] = ampGeom[amp]['DATASEC']
        foo[hdu].header['DETSEC'] = ampGeom[amp]['DETSEC']
    outfile = "%(sensor_id)s_%(test_type)s_%(image_type)s_%(seqno)s_%(time_stamp)s.fits" % locals()
    outdir = os.path.join(test_type, time_stamp)
    try:
        os.makedirs(outdir)
    except OSError:
        pass
    outfile = os.path.join(outdir, outfile)
    foo.writeto(outfile, clobber=True, checksum=True)

trr_files = lambda x : sorted(glob.glob(os.path.join('Final_TRR_Data_Set', x)))
#
# Fe55
#
infiles = trr_files('Xray Gain and PSF/11093*.fits')
for iframe, infile in enumerate(infiles):
    print "processing", os.path.basename(infile)
    seqno = '%03i' % iframe
    translate(infile, 'fe55', 'fe55', seqno)

#
# bias
#
infiles = trr_files('Noise - Zero frames/11093*.fits')
for iframe, infile in enumerate(infiles):
    print "processing", os.path.basename(infile)
    seqno = '%03i' % iframe
    translate(infile, 'fe55', 'bias', seqno)

#
# dark
#
infiles = trr_files('Dark 5 images/11093*.fits')
for iframe, infile in enumerate(infiles):
    print "processing", os.path.basename(infile)
    seqno = '%03i' % iframe
    translate(infile, 'dark', 'dark', seqno)

#
# flat
#
infiles = trr_files('satlin - multi/11093*.fits')
exptime = lambda x : fits.open(x)[0].header['EXPOSURE']
for iframe, infile in enumerate(infiles):
    print "processing", os.path.basename(infile)
    seqno = '%03i_flat1' % exptime(infile)
    translate(infile, 'flat', 'flat', seqno)

#
# lambda
#
infiles = trr_files('QE and PRNU/11093*qe*.fits')
wl = lambda x: fits.open(x)[0].header['WAVELEN']
for infile in infiles:
    print "processing", os.path.basename(infile)
    seqno = "%04i" % wl(infile)
    translate(infile, 'lambda', 'flat', seqno)
#
# trap
#
infiles = trr_files('Traps/11093*cycl*fits')
for iframe, infile in enumerate(infiles):
    print "processing", os.path.basename(infile)
    seqno = '%03i' % iframe
    translate(infile, 'trap', 'flat', seqno)
#
# sflat_nnn
#
infiles = trr_files('superflat low/11093*.fits')
for iframe, infile in enumerate(infiles):
    print "processing", os.path.basename(infile)
    seqno = '%03i' % iframe
    translate(infile, 'sflat_500', 'flat', seqno)
#
# spot
#
infiles = trr_files('Crosstalk/11093*.fits')
for iframe, infile in enumerate(infiles):
    print "processing", os.path.basename(infile)
    seqno = '%03i' % iframe
    translate(infile, 'spot', 'flat', seqno)
