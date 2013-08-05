import glob
import lsst.afw.image as afwImage
import pair_stats as ps
import numpy as np

outfile = 'ptctest.txt'
infiles1 = '/Volumes/Datadrive/CCDTesting/112-02/ptc/logain/*flat1.fits'
infiles2 = '/Volumes/Datadrive/CCDTesting/112-02/ptc/logain/*flat2.fits'
amps = [2,3]


#write names line to file
f=open(outfile, 'w+')
f.write('\t'.join(['amp', 'exptime', 'flat_mean', 'flat_var', 'gain', 'noise', '\n']))


#get lists of relevant file pairs
files1 = glob.glob(infiles1)
files2 = glob.glob(infiles2)


#calculate stats for each pair
for i in range(len(files1)):

	#grab exposure time
	hdr = afwImage.readMetadata(files1[i], 1)
	exptime = hdr.get('EXPTIME')

	for amp in amps:
		stats = ps.pair_stats(files1[i], files2[i], hdu=amp+1)
		flat_mean = stats[0].flat_mean
		flat_var = stats[0].flat_var
		gain = stats[0].gain
		noise = stats[0].noise
		
		f.write('\t'.join([str(amp), str(exptime), str(flat_mean), str(flat_var), str(gain), str(noise), '\n']))
	print str(i) + ' of ' + str(len(files1)) + ' file pairs done\n'
f.close()
