import argparse, glob
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9
import lsst.afw.math as afwMath
import numpy as np

#generate counts vs. exposure time data for a directory of flat fields
def linearity(directory, infilebase, outfile, amps, x0, y0, boxsize):
	#get list of files
	files = glob.glob(directory+infilebase)
	
	#write output file header
	f=open(directory+outfile, 'w+')
	f.write('amp\tx0\ty0\tboxsize\tmedian\texptime\n')
	
	for filename in files:
		#get exposure time from header
		hdr = afwImage.readMetadata(filename, 1)
		exptime = hdr.get('EXPTIME')
		
		for amp in amps:
			#define selected region
			box = afwGeom.Box2I(afwGeom.Point2I(x0, y0), afwGeom.Extent2I(boxsize, boxsize))  
		
			#read in selected region of file
			im = afwImage.ExposureF(filename, amp+1, box)
			
			#get median of region of image
			box_median = afwMath.makeStatistics(im.getMaskedImage(), afwMath.MEDIAN).getValue()
			
			#write amp, box parameters, region median, and exptime to file
			f.write(str(amp) + '\t' + str(x0) + '\t' + str(y0) + '\t' + str(boxsize) + '\t' + str(box_median) + '\t' + str(exptime) + '\n')
	
	f.close()
	


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generate linearity data for a set of flat field exposures')
	parser.add_argument('-d', '--direc', default='./', type=str, help="directory of files to work on. Include /")
	parser.add_argument('-f', '--infiles', type=str, default='*.fits', help="file string to search for; default= *.fits")
	parser.add_argument('-o', '--outfile', type=str, default='linearity_results.txt', help="output file name; default=linearity_results.txt")
	parser.add_argument('-x', '--x0', type=int, default=200, help="x0 pixel position for region of interest; default 200")
	parser.add_argument('-y', '--y0', type=int, default=900, help="y0 pixel position for region of interest; default 900")
	parser.add_argument('-s', '--size', type=int, default=100, help="box size in pixels; default 100")
	parser.add_argument('-a', '--amps', help="amps to be analyzed, separated by a space", type=int, nargs = '+', default = range(1,17))
	args = parser.parse_args()

	linearity(args.direc, args.infiles, args.outfile, args.amps, args.x0, args.y0, args.size)