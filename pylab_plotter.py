"""
@brief Interface to pylab routines.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import sys
import numpy as np
import pylab

_windows = {}
_win_id = -1

class Window(object):
    def __init__(self, id=None, subplot=111):
        global _win_id, _windows
        if id is None:
            _win_id += 1
            id = _win_id
        self.fig = pylab.figure(id)
        self.axes = [self.fig.add_subplot(subplot)]
        self.id = id
        self.handles = []
        _windows[id] = self
    def set_title(self, title, iax=None):
        if iax is None:
            iax = -1
        self.axes[iax].set_title(title)
    def select_subplot(self, *args, **kwds):
        self.axes.append(self.fig.add_subplot(*args, **kwds))

def newDataArray(cols, colnames):
    nt = {}
    for label, col in zip(colnames, cols):
        nt[label] = np.array(col)
    return nt

def setAxis(xrange=None, yrange=None):
    axisrange = list(pylab.axis())
    if xrange is not None:
        axisrange[:2] = xrange
    if yrange is not None:
        axisrange[2:] = yrange
    pylab.axis(axisrange)

def contour(x, y, z, levels=None, xname='x', yname='y',
            xrange=None, yrange=None, oplot=0, new_win=True):
    global _win_id, _windows
    if oplot == 0 and new_win:
        win = Window()
    else:
        win = _windows[_win_id]
    handle = pylab.contour(x, y, z, levels=levels)
    pylab.clabel(handle, inline=1, fontsize=10)
    if not oplot:
        pylab.xlabel(xname)
        pylab.ylabel(yname)
    setAxis(xrange, yrange)
    win.handles.append(handle)
    return win

def histogram(x, bin_edges=None, bins=50,
              xname='x', yname='entries / bin', oplot=0, color='k',
              linestyle='solid', weights=None,
              xrange=None, yrange=None, linewidth=1, ylog=False,
              new_win=True):
    global _win_id, _windows
    if oplot == 0 and new_win:
        win = Window()
    else:
        win = _windows[_win_id]
    if weights is None:
        weights = np.ones(len(x))
    if bin_edges is not None:
        bins = bin_edges
    handle = pylab.hist(x, bins=bins, range=xrange, histtype='step',
                        align='mid', color=color, linestyle=linestyle,
                        weights=weights, linewidth=linewidth, log=ylog)
    if not oplot:
        pylab.xlabel(xname)
        pylab.ylabel(yname)
    setAxis(xrange, yrange)
    win.handles.append(handle)
    return win

def xrange_errors(xmin, xmax, y, color='k'):
    for x1, x2, yy in zip(xmin, xmax, y):
        pylab.plot([x1, x2], [yy, yy], '%s-' % color)

def xerrors(x, y, xerr, color='k'):
    for xx, yy, err in zip(x, y, xerr):
        if err == xx == 1:
            err = 0.999 # workaround for plotting bug in pylab for log plots
        pylab.plot([xx-err, xx+err], [yy, yy], '%s-' % color)

def yerrors(x, y, yerr, color='k'):
    for xx, yy, err in zip(x, y, yerr):
        if err == yy == 1:
            err = 0.999 # workaround for plotting bug in pylab for ylog plots
        pylab.plot([xx, xx], [yy-err, yy+err], '%s-' % color)

def XYPlot(nt, xname, yname, xerr=None, yerr=None, xlog=0, ylog=0,
           oplot=0, color='k', xrange=None, yrange=None,
           new_win=True):
    if xerr is not None:
        xerr = nt[xerr]
    if yerr is not None:
        yerr = nt[yerr]
    win = xyplot(nt[xname], nt[yname], xerr=xerr, yerr=yerr,
                 xlog=xlog, ylog=ylog, oplot=oplot, color=color)
    return win

def xyplot(x, y, xerr=None, yerr=None, xlog=0, ylog=0, 
           xname='x', yname='y', oplot=0, color='k',
           xrange=None, yrange=None, new_win=True,
           markersize=3):
    global _win_id, _windows
    if oplot == 0 and new_win:
        win = Window()
    else:
        win = _windows[_win_id]
    marker = '%so' % color
    if xlog and ylog:
        handle = pylab.loglog(x, y, marker, markersize=markersize)
    elif xlog:
        handle = pylab.semilogx(x, y, marker, markersize=markersize)
    elif ylog:
        handle = pylab.semilogy(x, y, marker, markersize=markersize)
    else:
        handle = pylab.plot(x, y, marker, markersize=markersize)
    if xerr is not None:
        xerrors(x, y, xerr, color=color)
    if yerr is not None:
        yerrors(x, y, yerr, color=color)
    if not oplot:
        pylab.xlabel(xname)
        pylab.ylabel(yname)
        setAxis(xrange, yrange)
    win.handles.append(handle)
    return win

def curve(x, y, xlog=0, ylog=0, xname='x', yname='y', 
          oplot=0, color='k', lineStyle='-', linewidth=1,
          xrange=None, yrange=None, new_win=True):
    global _win_id, _windows
    if oplot == 0 and new_win:
        win = Window()
    else:
        win = _windows[_win_id]
    marker = '%s%s' % (color, lineStyle)
    if xlog and ylog:
        handle = pylab.loglog(x, y, marker, markersize=3, linewidth=linewidth)
    elif xlog:
        handle = pylab.semilogx(x, y, marker, markersize=3,
                                linewidth=linewidth)
    elif ylog:
        handle = pylab.semilogy(x, y, marker, markersize=3,
                                linewidth=linewidth)
    else:
        handle = pylab.plot(x, y, marker, markersize=3, linewidth=linewidth)
    if not oplot:
        pylab.xlabel(xname)
        pylab.ylabel(yname)
        setAxis(xrange, yrange)
    win.handles.append(handle)
    return win

def vline(x, color='k', lineStyle=':', linewidth=1):
    xmin, xmax, ymin, ymax = pylab.axis()
    curve([x, x], [ymin, ymax], color=color, lineStyle=lineStyle,
          oplot=1, linewidth=linewidth)

def hline(y, color='k', lineStyle=':', yrange=None, linewidth=1):
    xmin, xmax, ymin, ymax = pylab.axis()
    curve([xmin, xmax], [y, y], color=color, lineStyle=lineStyle, 
          oplot=1, yrange=yrange, linewidth=linewidth)

def clear(win_id=None):
    global _win_id, _windows
    if win_id is not None:
        ids = (win_id,)
    else:
        ids = _windows.keys()
    for id in ids:
        Window(id)
        pylab.clf()
    _win_id = -1

def save(outfile):
    pylab.savefig(outfile)

def show():
    pylab.show()

def prompt(prompt=None):
    if (prompt):
        sys.stderr.write(prompt)
    else:
        sys.stderr.write("Hit return to continue: ")
    x = sys.stdin.readline()
    return x

if __name__ == '__main__':
    x = np.linspace(0, 10)
    y = x**2
    curve(x, y, xname='x values', yname='y values')
    xyplot(x, y, oplot=1, color='r')

