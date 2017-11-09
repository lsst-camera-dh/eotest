from __future__ import print_function
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import image_utils as iu
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetect
import lsst.afw.display.ds9 as ds9


import glob
import numpy as np

for agg in range(1, 17):
    #output file name
    outfile = "xtalk_sysdet_output.txt"

    #write column names to output file
    f = open(outfile, "w+")
    f.write('\t'.join(["agg_amp", "vic_01", "vic_02", "vic_03", "vic_04", "vic_05", "vic_06", "vic_07",
                       "vic_08", "vic_09", "vic_10", "vic_11", "vic_12", "vic_13", "vic_14", "vic_15", "vic_16", "\n"]))

    #find appropriate image
    files = glob.glob("/Users/amali/Desktop/900nm*.fits")

    fname = files[0]

    #bias correct aggressor image
    im_a = afwImage.ImageF(fname, agg+1)
    im_a2 = iu.unbias_and_trim(im_a)

    #######The DM stack stuff I couldn't get to work
    #im_a = afwImage.ExposureF(fname, agg+1)
    #im_a_i = im_a.getMaskedImage().getImage()
    #im_a2 = iu.unbias_and_trim(im_a_i)

    #ds9.mtv(im_a2)

    #threshold = afwDetect.Threshold(30000)
    #fs = afwDetect.FootprintSet(im_a.getMaskedImage(), threshold)
    #fs.setMask(im_a.getMaskedImage().getMask(), "DETECTED")

    #ds9.mtv(im_a)

    ##########using numpy instead
    #locate spot pixels = counts > n above median
    n = 10000
    im_a_arr = im_a2.getArray()
    med = np.median(im_a_arr)
    detectedpos = np.where(im_a_arr > med+n)
    detectedcts = im_a_arr[detectedpos]

    #mean of detected counts
    agg_mean = np.mean(detectedcts)
    print(agg_mean)

    #get victim counts and ratio
    ratio = []
    for amp in range(1, 17):
        im_v = afwImage.ImageF(fname, amp+1)
        im_v2 = iu.unbias_and_trim(im_v)

        im_v_arr = im_v2.getArray()
        vic_mean = np.mean(im_v_arr[detectedpos])

        ratio.append(str(vic_mean/agg_mean))

    f.write(str(agg)+'\t' + '\t'.join(ratio) + '\t' + '\n')
f.close()
