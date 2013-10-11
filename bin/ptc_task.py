"""
@brief For pairs of flats obtain for a range of exposures, compute the
photon transfer curve (using pair_stats.py) and compute and write out
the full well.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import image_utils as imutils
from pipeline.TaskParser import TaskParser
from pair_stats_new import pair_stats

def find_flat2(flat1):
    pattern = flat1.split('flat1')[0] + 'flat2*.fits'
    flat2 = glob.glob(pattern)[0]
    return flat2

exptime = lambda x : afwImage.readMetadata(x, 1).get('EXPTIME')

def glob_flats(full_path, outfile='ptc_flats.txt'):
    flats = glob.glob(os.path.join(full_path, '*_flat?.fits'))
    output = open(outfile, 'w')
    for item in flats:
        output.write('%s\n' % item)
    output.close()

def find_flats(args):
    files = args.files(args.flats, args.flats_file_list)
    file1s = sorted([item.strip() for item in files
                     if item.find('flat1') != -1])
    return [(f1, find_flat2(f1)) for f1 in file1s]

def accumulate_stats(flats, outfile='ptc_results.txt', mask_files=(),
                     verbose=True):
    """
    Run pair_stats.py to find mean and variance (in units of DN) as a
    function of exposure time.
    """
    output = open(outfile, 'w')
    for file1, file2 in flats:
        if verbose:
            print "processing", file1
        exposure = exptime(file1)
        output.write('%12.4e' % exposure)
        for amp in imutils.allAmps:
            results = pair_stats(file1, file2, amp, mask_files=mask_files)
            output.write('  %12.4e  %12.4e' % (results.flat_mean,
                                               results.flat_var))
        output.write('\n')
        output.flush()
    output.close()

if __name__ == '__main__':
    parser = TaskParser('Compute photon transfer curve.')
    parser.add_argument('-f', '--flats', type=str,
                        help='flat pairs file pattern')
    parser.add_argument('-F', '--flats_file_list', type=str,
                        help='list of flat pairs')
    args = parser.parse_args()

    sensor = args.sensor()
    sensor_id = args.sensor_id
    gains = args.system_gains()
    mask_files = args.mask_files()

    outfile = os.path.join(args.output_dir, '%s_ptc.txt' % sensor_id)
    files = find_flats(args)
    accumulate_stats(files, outfile=outfile)
