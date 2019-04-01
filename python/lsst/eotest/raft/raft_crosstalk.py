from __future__ import absolute_import, print_function

import os
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor import CrosstalkTask

class CrosstalkData():

    def __init__(self, sensor_id, image_dict, sensor_pos_keys, gains, 
                 bias_frame=None):

        super(CrosstalkData, self).__init__()
        self.sensor_id = sensor_id
        self.sensor_pos_keys = sensor_pos_keys
        self.gains = gains
        self.bias_frame = bias_frame
        self.image_dict = {}
        self.make_dict(image_dict)

    def make_dict(self, image_dict):

        for key in image_dict:

            value = image_dict[key]
            
            if isinstance(value, list):
                outfile = '{0}_{1}_median_stack.fits'.format(self.sensor_id, key)
                imutils.fits_median_file(value)
                self.image_dict[key] = outfile
            else:
                self.image_dict[key] = value

class CrosstalkButler():

    def __init__(self, sensor_list, output_dir='./'):

        self.sensor_dict = {sensor:None for sensor in sensor_list}
        self.output_dir = output_dir
    
    def sensor_ingest(self, sensor_id, sensor_pos_keys, image_dict, 
                      gains, bias_frame=None):

        data = CrosstalkData(sensor_id, image_dict, sensor_pos_keys, gains, 
                             bias_frame)
        self.sensor_dict[sensor_id] = data

    def run_sensor(self, sensor_id, sensor_id2=None):
        """Run crosstalk task for single sensor pair."""

        print("Running crosstalk for {0} x {1}".format(sensor_id, sensor_id2))

        data = self.sensor_dict[sensor_id]
        infiles = [data.image_dict[key] for key in data.sensor_pos_keys]
        crosstalktask = CrosstalkTask()
        crosstalktask.config.output_dir = self.output_dir
        if sensor_id2 is not None:
            data2 = self.sensor_dict[sensor_id2]
            infiles2 = [data2.image_dict[key] for key in data.sensor_pos_keys]
            crosstalktask.run(sensor_id, infiles, data.gains, bias_frame=data.bias_frame,
                              sensor_id2=sensor_id2, infiles2=infiles2, gains2=data2.gains, 
                              bias_frame2=data2.bias_frame)
        else:
            crosstalktask.run(sensor_id, infiles, data.gains, bias_frame=data.bias_frame)

    def run_sensor_all(self, sensor_id):
        """Run crosstalk task for all inter-sensor crosstalk for given aggressor."""

        for sensor_id2 in self.sensor_dict:
            self.run_sensor(sensor_id, sensor_id2)

    def run_all(self):

        for sensor_id in self.sensor_dict:
            self.run_sensor_all(sensor_id)

## Assume 1 exposure per CCD
# Corresponding images for all CCDs (list)
# Name of aggressor CCD
# Loop over all CCDs including aggressor with crosstalk task and output files.
        
        
