"""
@brief Compute read noise distributions for a sample of images.  The
system gain for each channel is determined from Fe55 data, and the
corresponding bias frame is used to characterize the total system
noise.  Readout noise frames are used for determining the noise
contribution from the electronics.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import sys
import glob
    
import numpy as np
import pyfits

from read_noise import NoiseDists, median, stdev
from xray_gain import hdu_gains
from file_handling import get_file_list, export_file_list
from database.SensorDb import SensorDb, NullDbObject

def write_read_noise_dists(outfile, Nread, Nsys, gains, fe55, bias, sysnoise):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header.update('FE55FILE', fe55)
    output[0].header.update('BIASFILE', bias)
    output[0].header.update('SYSNFILE', sysnoise)
    for hdu, sigread, sigsys in zip(range(len(Nread)), Nread, Nsys):
        nread_col = pyfits.Column(name="CCD_READ_NOISE", format="E",
                                  unit="e- rms", array=sigread)
        nsys_col = pyfits.Column(name="SYSTEM_NOISE", format="E",
                                 unit="e- rms", array=sigsys)
        output.append(pyfits.new_table((nread_col, nsys_col)))
        output[hdu+1].name = "AMP%02i" % hdu
        output[0].header.update("GAIN%02i" % hdu, gains[hdu])
        output[0].header.update("RNOISE%02i" % hdu, median(sigread))
        output[0].header.update("RNSTDV%02i" % hdu, stdev(sigread))
    output.writeto(outfile, clobber=True)

def get_input_files(sensordir):
    Fe55_files = glob.glob(os.path.join(sensordir, 'Fe55', 'Fe55_exp_*'))
    Fe55_files.sort()

    bias_files = glob.glob(os.path.join(sensordir, 'Fe55', 'Fe55_bias_*'))
    bias_files.sort()

    system_noise_files = glob.glob(os.path.join(sensordir, 'system_noise',
                                                 'system_noise_*'))
    system_noise_files.sort()
    return Fe55_files, bias_files, system_noise_files

if __name__ == '__main__':
    if len(sys.argv) == 3:
        sensordir = sys.argv[1]
        outdir = sys.argv[2]
        sensor_id = os.path.basename(sensordir)
        Fe55_files, bias_files, system_noise_files = get_input_files(sensordir)
        pipeline_task = False
    else:
        try:
            Fe55_files = get_file_list('FE55')
            bias_files = get_file_list('BIAS')
            system_noise_files = get_file_list('SYSNOISE')
            sensor_id = os.environ['SENSOR_ID']
            outdir = os.environ['OUTPUTDIR']
            pipeline_task = True
        except:
            print "usage: python read_noise_task.py <sensordir> <outputdir>"
            sys.exit(1)

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    #
    # Get vendor from environment, otherwise assume "e2v".
    #
    try:
        vendor = os.environ['VENDOR']
    except KeyError:
        vendor = 'e2v'

    try:
        sensorDb = SensorDb(os.environ["DB_CREDENTIALS"])
        sensor = sensorDb.getSensor(vendor, sensor_id)
    except:
        print "using NullDbObject"
        sensor = NullDbObject(vendor, sensor_id)

    nhdu = 16

    outfiles = []
    Nread_dists = [[] for x in range(nhdu)]
    gain_dists = [[] for x in range(nhdu)]
    for i, fe55, bias, sysnoise in zip(range(len(Fe55_files)), Fe55_files,
                                       bias_files, system_noise_files):
        outfile = "ccd_read_noise_%s_%02i.fits" % (sensor_id, i)
        outfile = os.path.join(outdir, outfile)
        outfiles.append(outfile)
        
        print "Processing", fe55, bias, sysnoise, "->", outfile 
        gains = hdu_gains(fe55)
        noise_dists = NoiseDists(gains)
        
        Ntot = noise_dists(bias)
        Nsys = noise_dists(sysnoise)
        
        Nread = np.sqrt(Ntot*Ntot - Nsys*Nsys)
        for hdu in range(nhdu):
            Nread_dists[hdu].extend(Nread[hdu])
            gain_dists[hdu].append(gains[hdu])
    
        write_read_noise_dists(outfile, Nread, Nsys, gains,
                               fe55, bias, sysnoise)

    seg_gains = [median(gain_dists[hdu]) for hdu in range(nhdu)]
    sensor.add_ccd_result('gainMedian', median(seg_gains))
    for hdu, dist in enumerate(Nread_dists):
        sensor.add_seg_result(hdu, 'readNoise', median(dist))
        sensor.add_seg_result(hdu, 'gain', seg_gains[hdu])
            
    if pipeline_task:
        export_file_list(outfiles, "READNOISE")
