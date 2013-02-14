import numpy as np
import numpy.linalg as linalg
import pylab_plotter as plot

def linear_fit_fixedpt(x, y, x0, y0):
    dx = x - x0
    dy = y - y0
    slope = sum(dx*dy)/sum(dx*dx)
    return slope

def quadratic_fit_fixedpt(x, y, x0, y0):
    A = np.matrix([[sum((x**2 - x0**2)**2), sum((x-x0)*(x**2 - x0**2))],
                   [sum((x-x0)*(x**2 - x0**2)), sum((x-x0)**2)]])
    B = [sum((y - y0)*(x**2 - x0**2)), sum((y - y0)*(x - x0))]
    return linalg.solve(A, B)

ptcfile = 'ptc_results.txt'

data = np.recfromtxt(ptcfile)
data = data.transpose()

hdu = 1

exptime = data[0]
meanDN = data[hdu*2 + 1]
varDN = data[hdu*2 + 2]

#
# Estimate gain
#
nmin = 10
nmax = 100
results = np.polyfit(meanDN[nmin:nmax], varDN[nmin:nmax], 1)
gain = 1./results[0]

meanNe = gain*meanDN
varNe = gain*varDN

fwdev = 0.10

#
# Find the reference e- signal level
#
xref = 1000.
indx = 0
while meanNe[indx] < xref:
    indx += 1

#
# Linear fit local to this point.
#
imin = indx - 3
imax = indx + 2
f = np.poly1d(np.polyfit(meanNe[imin:imax], varNe[imin:imax], 1))
fref = f(xref)

#
# Loop over signal level points greater than this, compute linear and
# quadratic fit, and compute maximum deviation.
#
dvarmax = 0
imax = indx + 10
while dvarmax < fwdev:
#    f1 = np.poly1d(np.polyfit(meanNe[indx:imax], varNe[indx:imax], 1))
#    f2 = np.poly1d(np.polyfit(meanNe[indx:imax], varNe[indx:imax], 2))
    slope = linear_fit_fixedpt(meanNe[indx:imax], varNe[indx:imax], xref, fref)
    f1 = lambda x : slope*(x - xref) + fref

    quadfit = quadratic_fit_fixedpt(meanNe[indx:imax], varNe[indx:imax],
                                    xref, fref)
    f2 = lambda x : quadfit[0]*(x**2 - xref**2) + quadfit[1]*(x - xref) + fref
    
    fracdev = lambda x : np.abs(f1(x) - f2(x))/f1(x)
    x = meanNe[imax]
#    dvarmax = max(np.abs(f1(xref) - f2(xref))/f1(xref),
#                  np.abs(f1(x) - f2(x))/f1(x))
    dvarmax = fracdev(x)
    print imax, x, dvarmax, 1./slope
    imax += 1

plot.xyplot(meanNe, varNe, xname='mean(e-)', yname='var(e-)')
x = np.linspace(xref, meanNe[imax], 100)
plot.curve(x, f1(x), oplot=1, lineStyle=':')
plot.curve(x, f2(x), oplot=1, lineStyle='--')

plot.curve(x, np.abs(f1(x) - f2(x))/f1(x))
plot.hline(fwdev)
