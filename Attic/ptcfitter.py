import glob
import lsst.afw.image as afwImage
import pair_stats as ps
import numpy as np

datafile = 'ptctest.txt'
outfile = 'ptcfit.txt'
amps = [2]
n=50


data = np.recfromtxt(datafile, names=True)

#open output file
f = open(outfile, 'w+')
f.write('\t'.join(['amp', 'exptime', 'flat_mean', 'flat_var', 'ptc_fit', '\n']))


fmean = []
fvar = []
#read in data for specified amp
for amp in amps:
	fmean = data['flat_mean'][np.where(data['amp']==amp)]
	fvar = data['flat_var'][np.where(data['amp']==amp)]


#generate linear fit to first n points
#NOTE: I'm not sure what n should be -- maybe it will have to be changed to flux level or something similar
fit = np.polyfit(fmean[:n], fvar[:n], 1)
print fit

#figure out gain
fitgain = 1/fit[0]
print 'gain = ' + str(fitgain)