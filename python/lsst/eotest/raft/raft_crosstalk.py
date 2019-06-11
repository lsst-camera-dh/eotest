from __future__ import absolute_import, print_function

import os
import glob
import lsst.eotest.image_utils as imutils
import itertools
from lsst.eotest.sensor import CrosstalkTask
import argparse
import errno
import multiprocessing as mp
import time

def median_stack_images(main_dir, output_dir, sensor_list, xpos, ypos):

    print(xpos, ypos)
    for sensor_id in sensor_list:

        try:
            os.makedirs(os.path.join(output_dir, sensor_id))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        xtalk_files = glob.glob(os.path.join(main_dir, 'xtalk_{0}_{1}*'.format(xpos, ypos),
                                             '*_{0}.fits'.format(sensor_id)))
        infile = os.path.join(output_dir, sensor_id,
                              '{0}_{1}_{2}_median.fits'.format(sensor_id, xpos, ypos))
        imutils.fits_median_file(xtalk_files, infile, bitpix=-32)

    return True

class CrosstalkManager():

    def __init__(self, main_dir, xtalk_dict, gains_dict, bias_dict, output_dir='./'):

        self.main_dir = main_dir
        self.xtalk_dict = xtalk_dict
        self.gains_dict = gains_dict
        self.bias_dict = bias_dict
        self.output_dir = output_dir

    def run_crosstalk_task(self, sensor_id, sensor_id2):

        print(sensor_id, sensor_id2)

        gains = self.gains_dict[sensor_id]
        gains2 = self.gains_dict[sensor_id2]
        
        bias_frame = self.bias_dict[sensor_id]
        bias_frame2 = self.bias_dict[sensor_id2]

        infiles = []
        infiles2 = []
        for xpos, ypos in self.xtalk_dict[sensor_id]:
            xtalk_files = os.path.join(output_dir, sensor_id, '{0}_{1}_{2}_median.fits'.format(sensor_id, xpos, ypos))
            infiles.append(xtalk_files)

            xtalk_files2 = os.path.join(output_dir, sensor_id2, '{0}_{1}_{2}_median.fits'.format(sensor_id2, xpos, ypos))
            infiles2.append(xtalk_files2)

        print("running {0}, {1} crosstalk".format(sensor_id, sensor_id2))
        
        try:
            xtalktask = CrosstalkTask()
            xtalktask.config.threshold = 60000.
            xtalktask.config.output_dir = self.output_dir
            xtalktask.run(sensor_id, infiles, gains, bias_frame=bias_frame, 
                          infiles2=infiles2, sensor_id2=sensor_id2, gains2=gains2, bias_frame2=bias_frame2)
        except Exception as e:
            pass

def get_central_sensor(xpos, ypos):
    """Calculate the central CCD from given projector x and y positions."""

    intervals = [(x*127.0/3.-317.5, (x+1)*127.0/3.-317.5) for x in range(15)]
    
    for i, interval in enumerate(intervals):
        if interval[0] < float(xpos) <= interval[1]:
            raft_x = i // 3
            sensor_x = i % 3
        if interval[0] < float(ypos) <= interval[1]:
            raft_y = i // 3
            sensor_y = i % 3

    sensor_id = 'R{0}{1}_S{2}{3}'.format(raft_x, raft_y, sensor_x, sensor_y)

    return sensor_id
    
def BOT(main_dir, output_dir='./'):

    ## Get sensor list
    sensor_list = ['R22_S10', 'R22_S11', 'R22_S12', 'R22_S20', 'R22_S21', 'R22_S22',
                   'R10_S10', 'R10_S11', 'R10_S12', 'R10_S20', 'R10_S21', 'R10_S22']

    ## Get bias and gain per sensor
    gains_dict = {}
    bias_dict = {}
    for sensor_id in sensor_list:
        gains = {i : 1.0 for i in range(1, 17)} # replace with function to find relevant eotest
        bias_frame = None # replace with function to find relevant superbias
        gains_dict[sensor_id] = gains
        bias_dict[sensor_id] = None

    ## Construct dictionary of acquisition information
    xtalk_dict = dict()
    directory_list = [x.path for x in os.scandir(main_dir) if os.path.isdir(x.path)]
    pos_list = []
    for acquisition_dir in directory_list:
        basename = os.path.basename(acquisition_dir)
        if "xtalk" not in basename:
            continue
        xpos, ypos = basename.split('_')[-4:-2]
        central_sensor = get_central_sensor(xpos, ypos)
        pos = (xpos, ypos)
        pos_list.append(pos)
        if central_sensor in xtalk_dict:
            xtalk_dict[central_sensor].add(pos)
        else:
            xtalk_dict[central_sensor] = {pos}

    ## Create median stacks for all files
    a = time.time()
    calib_pool = mp.Pool(processes=8)
    calib_results = []
    for xpos, ypos in pos_list:

        calib_results.append(calib_pool.apply_async(median_stack_images, args=(main_dir, output_dir, sensor_list, xpos, ypos)))

    calib_results = [p.get() for p in calib_results]
    b = time.time()
    print(b-a)
        
    ## Run crosstalk for each sensor pair
    c = time.time()
    sensor_pairs = list(itertools.product(xtalk_dict, xtalk_dict))
    xtalk_pool = mp.Pool(processes=8)
    xtalk_results = []

    xtalk_manager = CrosstalkManager(main_dir, xtalk_dict, gains_dict, bias_dict, output_dir=output_dir)
    ## For each sensor pair, calculate crosstalk
    for sensor_id, sensor_id2 in sensor_pairs:

        try:
            os.makedirs(os.path.join(output_dir, sensor_id))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        try:
            os.makedirs(os.path.join(output_dir, sensor_id2))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        xtalk_results.append(xtalk_pool.apply_async(xtalk_manager.run_crosstalk_task, args = (sensor_id, sensor_id2)))

    xtalk_results = [p.get() for p in xtalk_results]
    d = time.time()
    print(d-c)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('eotest_dir', type=str)
    parser.add_argument('--output_dir', '-o', type=str, default='./')
    args = parser.parse_args()

    main_dir = args.eotest_dir
    output_dir = args.output_dir

    BOT(main_dir, output_dir=output_dir)






        
