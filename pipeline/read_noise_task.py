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

from read_noise import NoiseDists, median
from xray_gain import hdu_gains

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
        output[hdu+1].header.update("GAINFE55", gains[hdu])
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

def get_file_list(prefix):
    numfiles = int(os.environ["NUM%sFILES" % prefix])
    my_files = []
    for i in range(numfiles):
        my_files.append(os.environ["%s_%02i" % (prefix, i)])
    return my_files

def export_file_list(files, prefix):
    import pipeline
    pipeline.setVariable("NUM%sFILES" % prefix, "%s" % len(files))
    for i, item in enumerate(files):
        pipeline.setVariable("%s_%02i" % (prefix, i), item)

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
            print "usage: python get_read_noise.py <sensordir> <outputdir>"
            sys.exit(1)

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    outfiles = []
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
    
        write_read_noise_dists(outfile, Nread, Nsys, gains,
                               fe55, bias, sysnoise)
    if pipeline_task:
        export_file_list(outfiles, "READNOISE")
