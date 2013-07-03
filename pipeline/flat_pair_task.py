"""
@brief Compute detector response vs incident flux for pairs of flats.
These data are to be used for linearity and full-well measurments.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
from DetectorResponse import DetectorResponse
import lsst.afw.math as afwMath
import image_utils as imutils
from MaskedCCD import MaskedCCD
from pipeline.TaskParser import TaskParser

def pair_mean(flat1, flat2, amp):
    """
    Compute the mean pixel value for a given segment, averaged between
    two exposures (assumed to have the same exposure time). The inputs
    flat1 and flat2 are MaskedCCD objects.  Any masks are assumed to
    have been enabled in the stat_ctrl attribute.
    """
    #
    # Remove bias using fit to overscan region and trim to imaging region.
    #
    im1 = imutils.unbias_and_trim(flat1[amp])
    im2 = imutils.unbias_and_trim(flat2[amp])
    #
    # Compute the mean DN of the images.
    #
    stats1 = afwMath.makeStatistics(im1, afwMath.MEAN, flat1.stat_ctrl)
    stats2 = afwMath.makeStatistics(im2, afwMath.MEAN, flat2.stat_ctrl)

    avg_mean_value = (stats1.getValue(afwMath.MEAN) +
                      stats1.getValue(afwMath.MEAN))/2.
    return avg_mean_value

def extract_det_response(args, outfile):
    files = args.files(args.flats, args.flats_file_list)
    file1s = sorted([item.strip() for item in files
                     if item.find('flat1')  != -1])
    gains = args.system_gains()
    mask_files = args.mask_files()

    print "writing to", outfile
    output = open(outfile, 'w')
    for file1 in file1s:
        print "processing", file1
        file2 = file1.replace('flat1', 'flat2')
    
        flat1 = MaskedCCD(file1, mask_files=mask_files)
        flat2 = MaskedCCD(file2, mask_files=mask_files)

        if flat1.md.get('EXPTIME') != flat2.md.get('EXPTIME'):
            raise RuntimeError("Exposure times do not match for:\n%s\n%s\n"
                               % (file1, file2))
    
        flux = abs(flat1.md.get('EXPTIME')*flat1.md.get('MONDIODE') +
                   flat2.md.get('EXPTIME')*flat2.md.get('MONDIODE'))/2.

        output.write('%12.4e' % flux)
        #
        # Convert to e- and write out for each segment.
        #
        for amp in imutils.allAmps:
            output.write('  %12.4e'%(pair_mean(flat1, flat2, amp)*gains[amp]))
        output.write('\n')
        output.flush()
    output.close()
    return outfile

if __name__ == '__main__':
    parser = TaskParser('Compute detector response vs incident flux')
    parser.add_argument('-f', '--flats', type=str,
                        help='flat pairs file pattern')
    parser.add_argument('-F', '--flats_file_list', type=str,
                        help='list of flat pairs')
    args = parser.parse_args()
    sensor = args.sensor()
    sensor_id = args.sensor_id
    outfile = os.path.join(args.output_dir, '%s_det_response.txt' % sensor_id)
    #
    # Compute detector response from flat pair files.
    #
    extract_det_response(args, outfile)
    #
    # Perform full well and linearity analyses.
    #
    detresp = DetectorResponse(outfile)

    print "Segment    full well (e-/pixel)"
    full_well_values = []
    for amp in imutils.allAmps:
        try:
            full_well = detresp.full_well(amp)
        except RuntimeError:
            full_well = detresp.full_well(amp, frac_offset=0.05)
        print '%s            %.1f' % (imutils.channelIds[amp], full_well)
        full_well_values.append(full_well)
        sensor.add_seg_result(amp, 'fullWell', full_well)
    sensor.add_ccd_result('fullWellMean', imutils.mean(full_well_values))

    print "Segment    max. frac. deviation"
    maxdevs = []
    for amp in imutils.allAmps:
        maxdev, fit_pars = detresp.linearity(amp)
        maxdevs.append(maxdev)
        sensor.add_seg_result(amp, 'maxDeviation', maxdev)
        sensor.add_seg_result(amp, 'linefit_Slope', fit_pars[0])
        sensor.add_seg_result(amp, 'linefit_Intercept', fit_pars[1])
        print "%s             %.4f" % (imutils.channelIds[amp], maxdev)
    sensor.add_ccd_result('maxDeviation', max(maxdevs))
