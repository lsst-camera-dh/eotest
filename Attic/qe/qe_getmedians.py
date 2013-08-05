import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import image_utils as iu
import lsst.afw.math as afwMath
import glob

amps = range(1,17)
filepath = "/Volumes/Datadrive/CCDTesting/112-01/qe/qe/*.fits"
outfilepath = "qe_medians.txt"

###############Take median of each of a list of images

#List relevant images
files = glob.glob(filepath)

#setup output file

o = open(outfilepath, "w+")

o.write('\t'.join(['filename', 'exptime', 'amp','wavelength','median','photodiode', '\n']))

for i,fname in enumerate(files):
	#get wavelength, pd reading, and exptime
	md = afwImage.readMetadata(fname, 1)
	wl = md.get("MONO_WAVELENG")
	photodiode = md.get("K_PHOT_CURRENT")
	exptime = md.get("EXPTIME")

	for amp in amps:
		#open image, trim and debias, and get median
		im = afwImage.ImageF(fname, amp+1)
		im2 = iu.unbias_and_trim(im)
		
		ampmedian = afwMath.makeStatistics(im2, afwMath.MEDIAN).getValue()
		
		o.write('\t'.join([fname, str(exptime), str(amp), str(wl), str(ampmedian), str(photodiode), '\n']))
	
	print  "%i of %i done " %(i, len(files))
