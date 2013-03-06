"""
@brief Compute read noise distributions for a sample of images.  Bias
and system readout noise exposures, the latter for determining the
noise contribution from the electronics, must be provided.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import sys, traceback
import glob
import argparse
import numpy as np
import pyfits
import image_utils as imUtils
import pipeline.pipeline_utils as pipeUtils
from database.SensorDb import NullDbObject
from database.SensorGains import SensorGains

from read_noise import NoiseDists, median, stdev

def write_read_noise_dists(outfile, Ntot, Nsys, gains, bias, sysnoise):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header.update('BIASFILE', bias)
    output[0].header.update('SYSNFILE', sysnoise)
    for amp, sigtot, sigsys in zip(imUtils.allAmps, Ntot, Nsys):
        nread_col = pyfits.Column(name="TOTAL_NOISE", format="E",
                                  unit="e- rms", array=sigtot)
        nsys_col = pyfits.Column(name="SYSTEM_NOISE", format="E",
                                 unit="e- rms", array=sigsys)
        output.append(pyfits.new_table((nread_col, nsys_col)))
        output[amp].name = "AMP%s" % imUtils.channelIds[amp]
        output[0].header.update("GAIN%s" % imUtils.channelIds[amp], gains[amp])
        output[0].header.update("SIGTOT%s" % imUtils.channelIds[amp],
                                median(sigtot))
        output[0].header.update("SIGSYS%s" % imUtils.channelIds[amp],
                                median(sigsys))
    output.writeto(outfile, clobber=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(\
                      description='Compute Read Noise')
    parser.add_argument('-b', '--bias', help="input bias file pattern", type=str)
    parser.add_argument('-n', '--noise', help="input noise file pattern", type=str)
    parser.add_argument('-g', '--gain', help="input gain file", type=str)
    parser.add_argument('-i', '--id', help="sensor id", type=str)
    parser.add_argument('-o', '--outdir', help="output directory", type=str)
    parser.add_argument('-v', '--verbose', help="turn verbosity on", \
             action='store_true', default=False)
    args = parser.parse_args()

    if len(sys.argv) >= 5:
        bias_pattern = args.bias.replace('\\', '')
        sysnoise_pattern = args.noise.replace('\\', '')
        sensor_id = args.id
        outdir = args.outdir
        bias_files = glob.glob(bias_pattern)
        bias_files.sort()
        system_noise_files = glob.glob(sysnoise_pattern)
        system_noise_files.sort()
        gain_file = args.gain
        pipeline_task = False
    else:
        try:
            sensor_id = os.environ['SENSOR_ID']
            bias_files = pipeUtils.get_file_list('BIAS', sensor_id)
            system_noise_files = pipeUtils.get_file_list('SYSNOISE', sensor_id)
            outdir = os.environ['OUTPUTDIR']
            pipeline_task = True
        except KeyError:
            print "usage: python read_noise_task.py <bias file pattern> <sysnoise file pattern> <sensor id> <output dir> [<gain>=5.5]"
            sys.exit(1)

    if pipeline_task:
        gains, sensor = pipeUtils.setup(sys.argv, 5)
    else:
        g=[]
        sensor = NullDbObject()
        try:
            f = pyfits.open(gain_file)
            hdr = f[0].header
            for amp in imUtils.allAmps:
                g.append(hdr["GAIN%s" % imUtils.channelIds[amp]])
            gains = SensorGains(g)
        except:
            print 'Failed to access gains from: ', gain_file
            traceback.print_exc(file=sys.stdout)
            sys.exit()

    if not os.path.isdir(outdir):
        try:
            os.makedirs(outdir)
        except OSError:
            pass

    outfiles = []
    Nread_dists = dict([(amp, []) for amp in imUtils.allAmps])
    for i, bias, sysnoise in zip(range(len(bias_files)), bias_files,
                                 system_noise_files):
        outfile = "%s_read_noise_%03i.fits" % (sensor_id.replace('-', '_'), i)
        outfile = os.path.join(outdir, outfile)
        outfiles.append(outfile)
        
        print "Processing", bias, sysnoise, "->", outfile 
        noise_dists = NoiseDists(gains)
        
        Ntot = noise_dists(bias)
        Nsys = noise_dists(sysnoise)
        
        write_read_noise_dists(outfile, Ntot, Nsys, gains, bias, sysnoise)

    print "Segment    read noise"
    for amp in imUtils.allAmps:
        var = median(Ntot[amp-1])**2 - median(Nsys[amp-1])**2
        if var >= 0:
            Nread = np.sqrt(var)
        else:
            Nread = -1
        sensor.add_seg_result(amp, 'readNoise', Nread)
        print "%s         %.4f" % (imUtils.channelIds[amp], Nread)
