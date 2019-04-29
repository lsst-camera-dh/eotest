from __future__ import absolute_import, print_function

import os
import glob
import lsst.eotest.image_utils as imutils
import itertools
from lsst.eotest.sensor import CrosstalkTask
import argparse

class CrosstalkData():

    def __init__(self, sensor_id, gains, bias_frame=None):

        self.sensor_id = sensor_id
        self.gains = gains
        self.bias_frame = bias_frame
        self.image_dict = {}

    def add_image(self, central_sensor, infile):

        try:
            self.image_dict[central_sensor].append(infile)
        except KeyError:
            self.image_dict[central_sensor] = [infile]

    def infiles(self, sensor_id):

        return self.image_dict[sensor_id]

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
    ## For BOT testing
    sensor_list = ['R22_S10', 'R22_S11', 'R22_S12', 'R22_S20', 'R22_S21', 'R22_S22',
                   'R10_S10', 'R10_S11', 'R10_S12', 'R10_S20', 'R10_S21', 'R10_S22']

    gains = {i : 1.0 for i in range(1, 17)}
    xtalk_data = {sensor_id : CrosstalkData(sensor_id, gains) for sensor_id in sensor_list}

    ## Get acquisition directories from the main directory
    directory_list = [x.path for x in os.scandir(main_dir) if os.path.isdir(x.path)]
    projector_positions = set()
    for acquisition_dir in directory_list:
        basename = os.path.basename(acquisition_dir)
        if "xtalk" not in basename:
            continue
        xpos, ypos = basename.split('_')[-4:-2]
        projector_positions.add((xpos, ypos))

    ## Find images for each unique position, per sensor
    for xpos, ypos in projector_positions:

        central_sensor = get_central_sensor(xpos, ypos)

        for sensor_id in sensor_list:
            infiles = glob.glob(os.path.join(main_dir, 'xtalk_{0}_{1}*'.format(xpos, ypos),
                                             '*_{0}*.fits'.format(sensor_id)))
            outfile = os.path.join(output_dir, 
                                   '{0}_{1}_{2}_median.fits'.format(sensor_id, xpos, ypos))
            imutils.fits_median_file(infiles, outfile, bitpix=-32)
            xtalk_data[sensor_id].add_image(central_sensor, outfile)

    sensor_pairs = list(itertools.product(xtalk_data, xtalk_data))

    for sensor_id, sensor_id2 in sensor_pairs:

        print(sensor_id, sensor_id2)

        try:
            gains = xtalk_data[sensor_id].gains
            bias_frame = xtalk_data[sensor_id].bias_frame
            infiles = xtalk_data[sensor_id].infiles(sensor_id)

            gains2 = xtalk_data[sensor_id2].gains
            bias_frame2 = xtalk_data[sensor_id2].bias_frame
            infiles2 = xtalk_data[sensor_id2].infiles(sensor_id)

            xtalktask = CrosstalkTask()
            xtalktask.config.threshold = 60000.
            xtalktask.config.output_dir = output_dir
            xtalktask.run(sensor_id, infiles, gains, bias_frame=bias_frame, 
                          infiles2=infiles2, sensor_id2=sensor_id2, gains2=gains2, bias_frame2=bias_frame2)
        except Exception as e:
            print("Error occurred, skipping...")
            print(e)

def TS8(eotest_dir, output_dir='./'):
    ## For BOT testing
    sensor_list = ['S10', 'S11', 'S12', 'S20', 'S21', 'S22',
                   'S10', 'S11', 'S12', 'S20', 'S21', 'S22']

    gains = {i : 1.0 for i in range(1, 17)}
    xtalk_data = {sensor_id : CrosstalkData(sensor_id, gains) for sensor_id in sensor_list}

    ## Get acquisition directories from the main directory
    file_list = [f for f in os.listdir(os.path.join(main_dir, 'S00')) if os.isfile(f)]
    projector_positions = set()
    for filename in file_list:
        xpos, ypos = basename.split('_')[-5:-3]
        projector_positions.add((xpos, ypos))

    ## Find images for each unique position, per sensor
    for xpos, ypos in projector_positions:

        central_sensor = get_central_sensor(xpos, ypos)

        for sensor_id in sensor_list:
            infiles = glob.glob(os.path.join(main_dir, sensor_id, 
                                             '*_{0}_{1}_*.fits}'.format(xpos, ypos)))
            outfile = os.path.join(output_dir, 
                                   '{0}_{1}_{2}_median.fits'.format(sensor_id, xpos, ypos))
            imutils.fits_median_file(infiles, outfile, bitpix=-32)
            xtalk_data[sensor_id].add_image(central_sensor, infile)

    sensor_pairs = list(itertools.product(xtalk_data, xtalk_data))

    for sensor_id, sensor_id2 in sensor_pairs:

        gains = xtalk_data[sensor_id].gains
        bias_frame = xtalk_data[sensor_id].bias_frame
        infiles = sorted(xtalk_data[sensor_id].infiles(sensor_id))

        gains2 = xtalk_data[sensor_id2].gains
        bias_frame2 = xtalk_data[sensor_id2].bias_frame
        infiles2 = sorted(xtalk_data[sensor_id2].infiles(sensor_id))

        xtalktask = CrosstalkTask()
        xtalktaks.config.threshold = 60000.
        xtalktask.config.output_dir = output_dir
        xtalktask.run(sensor_id, infiles, gains, bias_frame=bias_frame, 
                      infiles2=infiles2, sensor_id2=sensor_id2, gains2=gains2, bias_frame2=bias_frame2)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('eotest_dir', type=str)
    parser.add_argument('--output_dir', '-o', type=str, default='./')
    args = parser.parse_args()

    main_dir = args.eotest_dir
    output_dir = args.output_dir

    BOT(main_dir, output_dir=output_dir)






        
