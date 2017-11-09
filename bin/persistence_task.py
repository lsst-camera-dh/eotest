#!/usr/bin/env python

"""
@brief Persistence test.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Compute image persistence.')
parser.add_argument('--pre_flat_darks', type=str,
                    help='file pattern for pre-flat dark files')
parser.add_argument('--post_flat_darks', type=str,
                    help='file pattern for post-flat dark files')
parser.add_argument('--flat', type=str,
                    help='filename of saturated flat')
args = parser.parse_args()

task = sensorTest.PersistenceTask()
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose

pre_flat_darks = args.files(args.pre_flat_darks, None)
post_flat_darks = args.files(args.post_flat_darks, None)
if args.verbose:
    print 'processing pre-flat files:'
    for item in pre_flat_darks:
        print '  ', item

    print 'processing post-flat files:'
    for item in post_flat_darks:
        print '  ', item

    print 'processing flat file:', args.flat

task.run(args.sensor_id, pre_flat_darks, args.flat,
         post_flat_darks, args.mask_files(pre_flat_darks[0]),
         args.system_gains())
