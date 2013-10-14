#!/usr/bin/env python

"""
@brief Compute read noise distributions for a sample of images.  Bias
and system readout noise exposures, the latter for determining the
noise contribution from the electronics, must be provided.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import pyfits
import lsst.test_scripts.image_utils as imutils
import lsst.test_scripts.sensor as sensorTest

def _write_read_noise_dists(outfile, Ntot, Nsys, gains, bias, sysnoise):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header.update('BIASFILE', bias)
    output[0].header.update('SYSNFILE', sysnoise)
    for amp in imutils.allAmps:
        sigtot, sigsys = Ntot[amp], Nsys[amp]
        nread_col = pyfits.Column(name="TOTAL_NOISE", format="E",
                                  unit="e- rms", array=sigtot)
        nsys_col = pyfits.Column(name="SYSTEM_NOISE", format="E",
                                 unit="e- rms", array=sigsys)
        output.append(pyfits.new_table((nread_col, nsys_col)))
        output[amp].name = "AMP%s" % imutils.channelIds[amp]
        output[0].header.update("GAIN%s" % imutils.channelIds[amp], gains[amp])
        output[0].header.update("SIGTOT%s" % imutils.channelIds[amp],
                                imutils.median(sigtot))
        output[0].header.update("SIGSYS%s" % imutils.channelIds[amp],
                                imutils.median(sigsys))
    output.writeto(outfile, clobber=True)

parser = sensorTest.TaskParser('Compute Read Noise')
parser.add_argument('-b', '--bias', type=str,
                    help='bias file pattern')
parser.add_argument('-B', '--bias_file_list', type=str,
                    help='list of bias files')
parser.add_argument('-n', '--noise', type=str, 
                    help='system noise file pattern')
parser.add_argument('-N', '--noise_file_list', type=str,
                    help='list of system noise files')
parser.add_argument('-x', '--dx', default=100, type=int,
                    help='subregion size in pixels along x-direction')
parser.add_argument('-y', '--dy', default=100, type=int,
                    help='subregion size in pixels along y-direction')
parser.add_argument('-S', '--nsamp', default=1000, type=int,
                    help='number of subregions to sample')
args = parser.parse_args()

print args.bias
print args.noise

bias_files = args.files(args.bias, args.bias_file_list)
print bias_files
system_noise_files = args.files(args.noise, args.noise_file_list)
print system_noise_files

sensor_id = args.sensor_id
sensor = args.sensor()
gains = args.system_gains()
mask_files = args.mask_files()

outfiles = []
Nread_dists = dict([(amp, []) for amp in imutils.allAmps])
for i, bias, sysnoise in zip(range(len(bias_files)), bias_files,
                             system_noise_files):
    outfile = "%s_read_noise_%03i.fits" % (sensor_id.replace('-', '_'), i)
    outfile = os.path.join(args.output_dir, outfile)
    outfiles.append(outfile)
    
    if args.verbose:
        print "Processing", bias, sysnoise, "->", outfile 

    #
    # Determine the nominal imaging region from the bias file.
    #
    ccd = sensorTest.MaskedCCD(bias)
    imaging = ccd.seg_regions[imutils.allAmps[0]].imaging
    sampler = imutils.SubRegionSampler(args.dx, args.dy, args.nsamp,
                                       imaging=imaging)

    Ntot = sensorTest.noise_dists(bias, gains, sampler, mask_files=mask_files)
    Nsys = sensorTest.noise_dists(sysnoise, gains, sampler,
                                  mask_files=mask_files)

    _write_read_noise_dists(outfile, Ntot, Nsys, gains, bias, sysnoise)

print "Segment    read noise"
for amp in imutils.allAmps:
    var = imutils.median(Ntot[amp])**2 - imutils.median(Nsys[amp])**2
    if var >= 0:
        Nread = np.sqrt(var)
    else:
        Nread = -1
    sensor.add_seg_result(amp, 'readNoise', Nread)
    print "%s         %.4f" % (imutils.channelIds[amp], Nread)
