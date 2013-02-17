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

from image_utils import allAmps, channelIds
from read_noise import NoiseDists, median, stdev
from pipeline.file_handling import get_file_list, export_file_list
from database.SensorDb import SensorDb, NullDbObject
from database.SensorGains import SensorGains

def write_read_noise_dists(outfile, Nread, Nsys, gains, bias, sysnoise):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header.update('BIASFILE', bias)
    output[0].header.update('SYSNFILE', sysnoise)
    for amp, sigread, sigsys in zip(allAmps, Nread, Nsys):
        nread_col = pyfits.Column(name="CCD_READ_NOISE", format="E",
                                  unit="e- rms", array=sigread)
        nsys_col = pyfits.Column(name="SYSTEM_NOISE", format="E",
                                 unit="e- rms", array=sigsys)
        output.append(pyfits.new_table((nread_col, nsys_col)))
        output[amp].name = "AMP%s" % channelIds[amp]
        output[0].header.update("GAIN%s" % channelIds[amp], gains[amp])
        output[0].header.update("RNOISE%s" % channelIds[amp], median(sigread))
        output[0].header.update("RNSTDV%s" % channelIds[amp], stdev(sigread))
    output.writeto(outfile, clobber=True)

if __name__ == '__main__':
    if len(sys.argv) >= 5:
        bias_pattern = sys.argv[1].replace('\\', '')
        sysnoise_pattern = sys.argv[2].replace('\\', '')
        sensor_id = sys.argv[3]
        outdir = sys.argv[4]
        try:
            gains = SensorGains(float(sys.argv[5]))
        except IndexError:
            print "Setting system gain to 5.5 e-/DN for all segments."
            gains = SensorGains(5.5)
        bias_files = glob.glob(bias_pattern)
        system_noise_files = glob.glob(sysnoise_pattern)
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
        except KeyError:
            print "usage: python read_noise_task.py <bias file pattern> <sysnoise file pattern> <sensor id> <output dir> [<gain>=5.5]"
            sys.exit(1)

    if not os.path.isdir(outdir):
        try:
            os.makedirs(outdir)
        except OSError:
            pass

    if pipeline_task:
        sensorDb = SensorDb(os.environ["DB_CREDENTIALS"])
        sensor = sensorDb.getSensor(vendor, sensor_id)
    else:
        vendor = 'e2v'
        sensor = NullDbObject(vendor, sensor_id)

    outfiles = []
    Nread_dists = dict([(amp, []) for amp in allAmps])
    for amp, bias, sysnoise in zip(allAmps, bias_files, system_noise_files):
        outfile = "ccd_read_noise_%s_%s.fits" % (sensor_id, channelIds[amp])
        outfile = os.path.join(outdir, outfile)
        outfiles.append(outfile)
        
        print "Processing", bias, sysnoise, "->", outfile 
        noise_dists = NoiseDists(gains)
        
        Ntot = noise_dists(bias)
        Nsys = noise_dists(sysnoise)
        
        Nread = np.sqrt(Ntot*Ntot - Nsys*Nsys)
        write_read_noise_dists(outfile, Nread, Nsys, gains, bias, sysnoise)

        for amp in allAmps:
            Nread_dists[amp].extend(Nread[amp-1])

    print "Segment    read noise"
    for amp in allAmps:
        med_dist = median(Nread_dists[amp])
        sensor.add_seg_result(amp, 'readNoise', med_dist)
        print "%s         %.4f" % (channelIds[amp], med_dist)
            
    if pipeline_task:
        export_file_list(outfiles, "READNOISE")
