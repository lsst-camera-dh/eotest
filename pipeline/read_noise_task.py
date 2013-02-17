"""
@brief Compute read noise distributions for a sample of images.  Bias
and system readout noise exposures, the latter for determining the
noise contribution from the electronics, must be provided.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import sys
import glob
    
import numpy as np
import pyfits

from read_noise import NoiseDists, median, stdev
from pipeline.file_handling import get_file_list, export_file_list
from database.SensorDb import SensorDb, NullDbObject
from database.SensorGains import SensorGains

def write_read_noise_dists(outfile, Nread, Nsys, gains, bias, sysnoise):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header.update('BIASFILE', bias)
    output[0].header.update('SYSNFILE', sysnoise)
    for hdu, sigread, sigsys in zip(range(len(Nread)), Nread, Nsys):
        nread_col = pyfits.Column(name="CCD_READ_NOISE", format="E",
                                  unit="e- rms", array=sigread)
        nsys_col = pyfits.Column(name="SYSTEM_NOISE", format="E",
                                 unit="e- rms", array=sigsys)
        output.append(pyfits.new_table((nread_col, nsys_col)))
        output[hdu+1].name = "AMP%02o" % hdu
        output[0].header.update("GAIN%02i" % hdu, gains[hdu])
        output[0].header.update("RNOISE%02i" % hdu, median(sigread))
        output[0].header.update("RNSTDV%02i" % hdu, stdev(sigread))
    output.writeto(outfile, clobber=True)

def get_input_files(sensordir):
    bias_files = glob.glob(os.path.join(sensordir, 'Fe55', 'Fe55_bias_*'))
    bias_files.sort()

    system_noise_files = glob.glob(os.path.join(sensordir, 'system_noise',
                                                 'system_noise_*'))
    system_noise_files.sort()
    return Fe55_files, bias_files, system_noise_files

if __name__ == '__main__':
    if len(sys.argv) >= 4:
        bias_pattern = sys.argv[1]
        sysnoise_pattern = sys.argv[2]
        sensor_id = sys.argv[2]
        outdir = sys.argv[3]
        try:
            gains = SensorGains(float(sys.argv[4]))
        except IndexError:
            print "Setting system gain to 5.5 e-/DN for all segments."
            gains = SensorGains(5.5)
        bias_files = glob.glob(bias_pattern)
        system_noise_files = glob_glob(sysnoise_pattern)
        pipeline_task = False
    else:
        try:
            bias_files = get_file_list('BIAS')
            system_noise_files = get_file_list('SYSNOISE')
            sensor_id = os.environ['SENSOR_ID']
            vendor = os.environ['CCD_VENDOR']
            outdir = os.environ['OUTPUTDIR']
            gains = SensorGains(vendorId=sensor_id, vendor=vendor)
            pipeline_task = True
        except:
            print "usage: python read_noise_task.py <sensordir> <outputdir> <gain>"
            sys.exit(1)

    if not os.path.isdir(outdir):
        try:
            os.mkdirs(outdir)
        except OSError:
            pass

    if pipeline_task:
        sensorDb = SensorDb(os.environ["DB_CREDENTIALS"])
        sensor = sensorDb.getSensor(vendor, sensor_id)
    else:
        vendor = 'e2v'
        sensor = NullDbObject(vendor, sensor_id)

    nhdu = 16

    outfiles = []
    Nread_dists = [[] for x in range(nhdu)]
    for i, bias, sysnoise in zip(range(len(bias_files)),
                                 bias_files, system_noise_files):
        outfile = "ccd_read_noise_%s_%02i.fits" % (sensor_id, i)
        outfile = os.path.join(outdir, outfile)
        outfiles.append(outfile)
        
        print "Processing", bias, sysnoise, "->", outfile 
        noise_dists = NoiseDists(gains)
        
        Ntot = noise_dists(bias)
        Nsys = noise_dists(sysnoise)
        
        Nread = np.sqrt(Ntot*Ntot - Nsys*Nsys)
        write_read_noise_dists(outfile, Nread, Nsys, gains, bias, sysnoise)

        for hdu in range(nhdu):
            Nread_dists[hdu].extend(Nread[hdu])

    for hdu, dist in enumerate(Nread_dists):
        sensor.add_seg_result(hdu, 'readNoise', median(dist))
            
    if pipeline_task:
        export_file_list(outfiles, "READNOISE")
