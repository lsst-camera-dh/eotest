"""
@brief Compute charge transfer (in)efficiency.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import pyfits
import image_utils as imutils
from pipeline.TaskParser import TaskParser
from eperTask import EPERTask

def superflat_sum(files, outfile='summed.fits'):
    summed = pyfits.open(files[0])
    for infile in files[1:]:
        addend = pyfits.open(infile)
        for amp in imutils.allAmps:
            summed[amp].data += addend[amp].data
    summed.writeto(outfile, clobber=True)
    return outfile

if __name__ == '__main__':
    parser = TaskParser('Compute charge transfer efficiency in parallel and serial directions using extended pixel edge response technique')
    parser.add_argument('-f', '--superflat_pattern', type=str,
                        help='superflat dataset file pattern')
    parser.add_argument('-F', '--superflat_file_list', type=str,
                        help='list of superflat files')
    args = parser.parse_args()
    sensor = args.sensor()
    sensor_id = args.sensor_id
    outfile = os.path.join(args.output_dir, '%s_cti_values.txt' % sensor_id)

    files = args.files(args.superflat_pattern, args.superflat_file_list)
    summed_file = superflat_sum(files)

    overscans = 3
    #
    # Compute serial CTE
    #
    s_task = EPERTask()
    s_task.config.direction = 's'
    s_task.config.verbose = False
    s_task.config.cti = True
    scti = s_task.run(summed_file, imutils.allAmps, overscans)
    #
    # Compute parallel CTE
    #
    p_task = EPERTask()
    p_task.config.direction = 'p'
    p_task.config.verbose = False
    p_task.config.cti = True
    pcti = p_task.run(summed_file, imutils.allAmps, overscans)
    #
    # Write results to the db and text file
    #
    sensor.add_ccd_result('ctiParallelMean', imutils.mean(pcti.values()))
    sensor.add_ccd_result('ctiSerialMean', imutils.mean(scti.values()))
    output = open(outfile, 'w')
    output.write('amp  parallel_cti  serial_cti\n')
    for amp in imutils.allAmps:
        sensor.add_seg_result(amp, 'ctiParallel', pcti[amp])
        sensor.add_seg_result(amp, 'ctiSerial', scti[amp])
        output.write('%02i  %12.4e  %12.4e\n' % (amp, pcti[amp], scti[amp]))
    output.close()
