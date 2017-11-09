import argparse
import glob
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9
import lsst.afw.math as afwMath
import numpy as np
import matplotlib.pyplot as p


#generate counts vs. exposure time data for a directory of flat fields
def linearstats(directory, infilebase, outfile, amps, x0, y0, boxsize):
    #get list of files
    files = glob.glob(directory+infilebase)

    #write output file header
    f = open(outfile, 'w+')
    f.write('\t'.join(['filename', 'amp', 'exptime', 'photodiode', 'mediancounts', '\n']))

    for filename in files:
        #strip directory from filename
        fname = filename.split('/')[-1]

        #get exposure time from header
        hdr = afwImage.readMetadata(filename, 1)
        exptime = hdr.get('EXPTIME')
        kphot = hdr.get('K_PHOT_CURRENT')

        for amp in amps:
            #define selected region
            box = afwGeom.Box2I(afwGeom.Point2I(x0, y0), afwGeom.Extent2I(boxsize, boxsize))

            #read in selected region of file
            im = afwImage.ExposureF(filename, amp+1, box)

            #get median of region of image
            box_median = afwMath.makeStatistics(im.getMaskedImage(), afwMath.MEDIAN).getValue()

            #write filename, amp, exposure time, photodiode measurement, and region median to file
            f.write('\t'.join([str(fname), str(amp), str(exptime), str(kphot), str(box_median), '\n']))

    f.close()


def linearfit(infile, gain, amps):
    #read in file output from linearstats
    lindata = np.recfromtxt(infile, names=True)

    for amp in amps:
        #position of chosen amp data
        whichamp = np.where(lindata['amp'] == amp)

        #extract exposure time data
        linexp = lindata['exptime'][whichamp]

        #extract pixel median in DN
        linmedian = lindata['mediancounts'][whichamp]

        #convert DN median to electron median
        linmedian *= gain

        #select points from 100 e-/pixel to 90000e-/pixel with exptimes below saturation
        maxpoint = np.where(linmedian == max(linmedian))[0][0]
        print maxpoint

        selected = np.where((linmedian[:maxpoint] > 100) & (linmedian[:maxpoint] < 90000))
        selectedpoints = linmedian[selected]
        exptimes = linexp[selected]

        #fit a line to selected points
        fit = np.polyfit(exptimes, selectedpoints, 1)

        #find points on the fit line
        fitpoints = np.polyval(fit, exptimes)

        #plot figure -- for sanity checking
        p.figure()
        p.plot(linexp, linmedian, '.b')
        p.plot(exptimes, fitpoints, '-r')
        p.show()

        #calculate deviations from fit line as percent of expected value
        deviations = []
        for i in range(len(selectedpoints)):
            deviation = abs((fitpoints[i] - selectedpoints[i])/float(fitpoints[i])*100)
            print deviation
            deviations.append(deviation)

        #find maximum deviation
        print max(deviations)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate linearity data for a set of flat field exposures')
    parser.add_argument('-d', '--direc', default='./', type=str,
                        help="directory of files to work on. Include /")
    parser.add_argument('-f', '--infiles', type=str, default='*.fits',
                        help="file string to search for; default= *.fits")
    parser.add_argument('-o', '--outfile', type=str, default='linearity_results.txt',
                        help="output file name; default=linearity_results.txt")
    parser.add_argument('-x', '--x0', type=int, default=200,
                        help="x0 pixel position for region of interest; default 200")
    parser.add_argument('-y', '--y0', type=int, default=900,
                        help="y0 pixel position for region of interest; default 900")
    parser.add_argument('-s', '--size', type=int, default=100, help="box size in pixels; default 100")
    parser.add_argument('-a', '--amps', help="amps to be analyzed, separated by a space",
                        type=int, nargs='+', default=range(1, 17))
    parser.add_argument('-g', '--gain', help="system gain", type=float, default=5)
    args = parser.parse_args()

    linearstats(args.direc, args.infiles, args.outfile, args.amps, args.x0, args.y0, args.size)
    linearfit(args.outfile, args.gain, args.amps)
