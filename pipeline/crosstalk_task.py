"""
@brief Task to produce crosstalk matrix from a set of spot image files.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import pylab
from TaskParser import TaskParser
from crosstalk import make_crosstalk_matrix

parser = TaskParser('Compute crosstalk from a set of spot images')
parser.add_argument('-f', '--xtalk_files', type=str,
                    help='file pattern for crosstalk files')
parser.add_argument('-F', '--xtalk_file_list', type=str,
                    help='list of crosstalk files')

args = parser.parse_args()
sensor_id = args.sensor_id
xtalk_files = args.files(args.xtalk_files, args.xtalk_file_list)

xtalk = make_crosstalk_matrix(xtalk_files, mask_files=args.mask_files())
xtalk.plot_matrix('%s' % sensor_id)
pylab.savefig(os.path.join(args.output_dir, '%s_xtalk_matrix.png' % sensor_id))
xtalk.write(os.path.join(args.output_dir, '%s_xtalk_matrix.txt' % sensor_id))
