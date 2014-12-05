import numpy as np
import lsst.afw.math as afwMath
import pylab_plotter as plot

def getProcessedImageArray(ccd, amp, nx=10, ny=10):
    """
    Subtract afwMath local background model from an de-biased and
    trimmed image segment and return the underlying ndarray.
    """
    # Define extent of local background region
    bg_ctrl = afwMath.BackgroundControl(nx, ny, ccd.stat_ctrl)
    image = ccd.unbiased_and_trimmed_image(amp)
    bg = afwMath.makeBackground(image, bg_ctrl)
    bg_image = bg.getImageF()
    image -= bg.getImageF()
    return image.getImage().getArray()

class TrapFinder(object):
    """
    Class to find traps using pixel dipole "correlator" description in
    Kotov et al. 2014, NIMA 
    (http://www.sciencedirect.com/science/article/pii/S0168900214012066).
    """
    def __init__(self, ccd, amp, C2_thresh=10, C3_thresh=15,
                 nx=10, ny=10, edge_rolloff=10):
        """
        ccd = MaskedCCD object.
        amp = Amplifier segment to analyze.
        masked_image = afw.MaskedImage of a single amplifier segment.
        C2_thresh = Theshold for C2 correlator value.  Accepted pixels have
                    C2 < -C2_thresh.
        C3_thresh = Threshold for C3 value, where C3 = 2*C2 + A0**2 + A1**2
                    (A0, A1 = bg-subtracted pixel values of dipole. 
                    NB: C2 = A0*A1 < 0 for a candidate dipole.).  Accepted
                    pixels have C3 > C3_thresh
        nx, ny = size of local background region in pixels
        edge_rolloff = 10.  Ignore this number of rows near the sensor edge
                       to avoid edge rolloff regions.
        """
        self.imarr = getProcessedImageArray(ccd, amp, nx, ny).transpose()
        self.C2_thresh = C2_thresh
        self.C3_thresh = C3_thresh
        self.edge_rolloff = edge_rolloff
        self.prescan = ccd.amp_geom.prescan_width
        self.row_max = ccd.amp_geom.imaging.getHeight()+1
    def find(self, regfile='ds9.reg'):
        nx, ny = self.imarr.shape
        my_arrays = 'ix iy c2 c3 a0 a1'.split()
        for item in my_arrays:
            exec('%(item)s = []' % locals())
        for icol in range(nx):
            results = self.process_column(icol)
            for i, item in enumerate(my_arrays):
                exec('%(item)s.extend(results[%(i)i])' % locals())
        for item in my_arrays:
            exec('%(item)s = np.array(%(item)s)' % locals())
        self._write_reg_file(regfile, ix, iy)
        return ix, iy, c2, c3, a0, a1
    def _write_reg_file(self, regfile, ix, iy):
        reg_output = open(regfile, 'w')
        reg_output.write("""# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image
""")
        for coord in zip(ix, iy):
            reg_output.write("point(%i,%i) # point=circle\n" % coord)
        reg_output.close()
    def process_column(self, icol, plot_stats=False, oplot=0):
        """
        Process a single column and return a tuple of
        ix = x-pixel in image coordinates
        iy = y-pixel in image coordinates
        C2 = Correlator value = A0*A1
        C3 = 2*C2 + A0**2 + A1**2
        A0 = Background-subtracted pixel value of first pixel of dipole
        A1 = Background-subtracted pixel value of second pixel of dipole
        """
        col0 = self.imarr[icol][self.edge_rolloff:self.row_max]
        sigma = np.std(col0)
        col = col0/sigma
        C2 = col[1:]*col[:-1]
        C3 = np.abs(2*C2 + col[1:]**2 + col[:-1]**2)
        index = np.where((C2 < -self.C2_thresh) & (C3 < self.C3_thresh))
        iy, a1 = [], []
        for irow in index[0]:
            iy.append(irow + self.edge_rolloff + 1)
            a1.append(col0[irow + 1])
            if col0[irow] < a1[-1]:
                # We have a forward trap so increment y-coordinate by 1
                iy[-1] += 1
        ix = np.ones(len(index[0]))*(icol + self.prescan + 1)
        iy = np.array(iy)
        c2 = C2[index]
        c3 = C3[index]
        a0 = col0[index]
        a1 = np.array(a1)
        if plot_stats:
            self._plot_stats(icol, C2, C3, oplot)
        return ix, iy, c2, c3, a0, a1
    def _plot_stats(self, icol, C2, C3, oplot):
        plot.pylab.ion()
        if oplot == 0:
            rows = range(self.edge_rolloff+1, self.edge_rolloff + len(C2) + 1)
            win0 = plot.curve(rows, C2, xname='row', yname='C2')
            win0.set_title('Column %i' % icol)
            plot.hline(-self.C2_thresh)
        win1 = plot.xyplot(C2, C3, xname='C2', yname='C3', oplot=oplot)
        win1.set_title('Column %i' % icol)
        plot.hline(self.C3_thresh)
        plot.vline(-self.C2_thresh)

if __name__ == '__main__':
    from MaskedCCD import MaskedCCD
    trap_file = '../20141118-184914/114-03_trap_000.fits'
    ccd = MaskedCCD(trap_file)
    amp = 4

    finder = TrapFinder(ccd, amp)
    results = finder.find()
    for row in zip(*results):
        print row

    index = np.where(results[2] == min(results[2]))
    finder.process_column(results[0][index[0][0]]-4, True)
    for icol in range(509):
        finder.process_column(icol, True, oplot=1)
    plot.xyplot(results[2], results[3], oplot=1, color='r')
