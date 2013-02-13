import os
import glob
import pyfits
from pair_stats import pair_stats

exptime = lambda x : pyfits.open(x)[0].header['EXPTIME']

def find_flats(full_path):
    flat1s = glob.glob(os.path.join(full_path, '*_flat1.fits'))
    flat1s.sort()
    flats = []
    for file1 in flat1s:
        file2 = file1.replace('_flat1.fits', '_flat2.fits')
        if os.path.isfile(file1) and os.path.isfile(file2):
            flats.append(file1, file2)
    return flats

def accumulate_stats(flats, outfile='ptc_results.txt'):
    output = open(outfile, 'w')
    for file1, file2 in flats:
        exposure = exptime(file1)
        output.write('%12.4e' % exposure)
        for hdu in range(16):
            results, b1, b2 = pair_stats(file1, file2, hdu+2)
            output.write('  %12.4e  %12.4e'%(results.flat_mean,
                                             results.flat_var))
        output.write('\n')
        output.flush()
    output.close()

if __name__ == '__main__':
    import numpy as np
    import pylab_plotter as plot
    
    full_path = '/nfs/farm/g/lsst/u1/testData/HarvardData/112-02/ptc/logain'
    outfile = 'ptc_results.txt'

    gain = 5
    
#    flats = find_flats(full_path)
#    accumulate_stats(flats, outfile=outfile)
    
    data = np.recfromtxt(outfile)
    data = data.transpose()
    exptime = data[0]
    meanDN = data[1]
    varDN = data[2]
    #
    # Using gain estimate, find exposure where e-=5e4.
    #
    meanNe = gain*meanDN
    npts = 0
    while (meanNe[npts] < 5e4):
        npts += 1
    #
    # Perform linear fit to first npts points and obtain better gain estimate.
    #
    results = np.polyfit(meanDN[:npts], varDN[:npts], 1)
    gain = 1./results[0]
    print "initial fit gain =", gain
    meanNe = gain*meanDN

    #
    # Select region for final linear fit as specified in E/O doc.
    #
#    indx = np.where((meanNe > 3000) & (meanNe < 5e4))
    indx = np.where((meanNe > 100) & (meanNe < 9e4))

    results = np.polyfit(meanDN[indx], varDN[indx], 1)
    
    #
    # update gain
    #
    gain = 1./results[0]
    print "updated gain =", gain

    meanNe = gain*meanDN
    varNe = gain*varDN
    plot.xyplot(meanNe, varNe, xname='mean(e-)', yname='var(e-)')
    plot.xyplot(meanNe[indx], varNe[indx], oplot=1, color='r')

    jndx = np.where(meanNe < 9e4)
    
    f = lambda x : results[0]*x + results[1]
    deviations = (varNe[jndx] - f(meanNe[jndx]))/varNe[jndx]

    x = np.linspace(min(meanNe[jndx]), max(meanNe[jndx]), 100)
    plot.curve(x, f(x), oplot=1, lineStyle=':')

    print "gain = ", gain
    print "max. fractional deviation:", max(np.abs(deviations))

    plot.xyplot(meanNe[jndx], deviations, xname='mean(e-)',
                yname='(var(e-) - f(mean(e-)))/var(e-)')
    plot.hline(0)
