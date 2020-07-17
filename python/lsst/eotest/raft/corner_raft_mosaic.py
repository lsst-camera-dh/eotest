import copy
import numpy as np
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor.EOTestPlots import cmap_range
import lsst.eotest.sensor as sensorTest
import lsst.eotest.raft as raftTest

# These are nominal locations of one of the corners of each imaging
# segment for the corresponding amp in the guide and wavefront sensors
# in a coordinate system that is rotated so that the wavefront sensors
# are in the lower right corner of the raft.
slot_data = dict()
slot_data['SG0'] = {amp: (4126, 4827 + (amp - 1)*509) for amp in range(1, 9)}
slot_data['SG0'].update({amp: (125, 8390 - (amp - 9)*509)
                         for amp in range(9, 17)})
slot_data['SG1'] = {amp: (8390 - (amp - 1)*509, 4126) for amp in range(1, 9)}
slot_data['SG1'].update({amp: (4827 + (amp - 9)*509, 125)
                         for amp in range(9, 17)})
slot_data['SW1'] = {amp: (602 + (amp-1)*509, 4126) for amp in range(1, 9)}
slot_data['SW0'] = {amp: (4165 - (amp-1)*509, 125) for amp in range(1, 9)}

class CornerRaftMosaic:
    def __init__(self, fits_files, nx=9000, ny=9000, bias_frames=None):
        self.fits_files = fits_files
        with fits.open(list(fits_files.values())[0]) as hdus:
            self.raft_name = hdus[0].header['RAFTNAME']
            try:
                self.wl = hdus[0].header['MONOWL']
            except KeyError:
                self.wl = None
        self.image_array = np.zeros((nx, ny), dtype=np.float32)
        for slot, filename in fits_files.items():
            if slot not in slot_data:
                continue
            bias_frame = bias_frames[slot] if bias_frames is not None else None
            ccd = sensorTest.MaskedCCD(filename, bias_frame=bias_frame)
            for amp in ccd:
                if slot.startswith('SW'):
                    self._set_sw_segment(slot, ccd, amp)
                elif slot == 'SG0':
                    self._set_sg0_segment(slot, ccd, amp)
                elif slot == 'SG1':
                    self._set_sg1_segment(slot, ccd, amp)

    def _set_sw_segment(self, slot, ccd, amp):
        seg_array = np.array(copy.deepcopy(ccd.unbiased_and_trimmed_image(amp)
                                           .getImage().getArray()),
                             dtype=np.float32)
        xx, yy = slot_data[slot][amp]
        xmin = xx - ccd.amp_geom.nx
        xmax = xx
        if slot == 'SW0':
            # Flip amps in serial direction and parallel directions.
            # This is equivalent to a 180 degree rotation about the
            # segment center.
            seg_array = seg_array[::-1, ::-1]
            ymin = yy
            ymax = yy + ccd.amp_geom.ny
        else:
            ymin = yy - ccd.amp_geom.ny
            ymax = yy
        self.image_array[ymin:ymax, xmin:xmax] = seg_array

    def _set_sg0_segment(self, slot, ccd, amp):
        seg_array = np.array(copy.deepcopy(ccd.unbiased_and_trimmed_image(amp)
                                           .getImage().getArray()),
                             dtype=np.float32)
        seg_array = seg_array.transpose()
        xx, yy = slot_data[slot][amp]
        ymin = yy - ccd.amp_geom.nx
        ymax = yy
        if amp < 9:
            seg_array = seg_array[:, ::-1]
            xmin = xx - ccd.amp_geom.ny
            xmax = xx
        else:
            xmin = xx
            xmax = xx + ccd.amp_geom.ny
        self.image_array[ymin:ymax, xmin:xmax] = seg_array

    def _set_sg1_segment(self, slot, ccd, amp):
        seg_array = np.array(copy.deepcopy(ccd.unbiased_and_trimmed_image(amp)
                                           .getImage().getArray()),
                             dtype=np.float32)
        xx, yy = slot_data[slot][amp]
        seg_array = seg_array[:, ::-1]
        xmin = xx - ccd.amp_geom.nx
        xmax = xx
        if amp < 9:
            seg_array = seg_array[::-1, :]
            ymin = yy - ccd.amp_geom.ny
            ymax = yy
        else:
            ymin = yy
            ymax = yy + ccd.amp_geom.ny
        self.image_array[ymin:ymax, xmin:xmax] = seg_array

    def plot(self, title=None, cmap=plt.cm.hot, nsig=5, figsize=(10, 10),
             binsize=10, flipx=True, textcolor='c', annotation='',
             rotate180=False, vrange=None):
        """
        Render the raft mosaic.

        Parameters
        ----------
        title : str, optional
            The plot title. If None (default), then build the title
            from the RAFTNAME and MONOWL primary header keyword values.
        cmap : matplotlib.colors.Colormap, optional
            The color map to use. Default: matplotlib.pyplot.cm.hot.
        nsig : float, optional
            The n-sigma value for the sigma clipping used to determine
            the pixel value range over which the color map is mapped.
        figsize : (float, float), optional
            The width x height size of the figure in inches. Default: (10, 10).
        binsize : int, optional
            Rebin the plotted image data by binsize*binsize,
            averging over the coarser bin.  Default: 10
        flipx : bool, optional
            Flip full raft mosaic in x so that parity of image matches
            LCA-13381. Default: True
        textcolor : str, optional
            Color of the text for the segment and sensor labeling.
            Default: 'c' (cyan)
        annotation : str, optional
            Description of the plot, e.g., pixel units (ADU or e-),
            gain-corrected, bias-subtracted.  Default: ''
        rotate180 : bool [False]
            Flag to rotate the mosaic by 180 degrees to match the
            orientation of the focalplane mosiacs created for the
            BOT-level plots.
        vrange : (float, float) [None]
            Range of pixel values to plot.  If None, then the cmap_range
            function will be used.
        """
        plt.rcParams['figure.figsize'] = figsize
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        output_array = imutils.rebin_array(self.image_array, binsize,
                                           use_mean=True)
        if flipx:
            output_array = output_array[:, ::-1]
        if rotate180:
            ny, nx = output_array.shape
            rotated_array = np.zeros((nx, ny), dtype=output_array.dtype)
            for j in range(ny):
                rotated_array[:, ny-1-j] = output_array[::-1, j]
            output_array = rotated_array
        image = ax.imshow(output_array, interpolation='nearest', cmap=cmap)
        # Set range and normalization of color map based on sigma-clip
        # of pixel values.
        if vrange is None:
            vmin, vmax = cmap_range(output_array, nsig=nsig)
        else:
            vmin, vmax = vrange
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        image.set_norm(norm)
        if title is None:
            if self.wl is None:
                title = self.raft_name
            else:
                title = "%s, %i nm" % (self.raft_name, self.wl)
        ax.set_title(title)
        fig.colorbar(image)
        # Turn off ticks and tick labels for x- and y-axes.
        plt.tick_params(axis='both', which='both',
                        top='off', bottom='off', left='off', right='off',
                        labelbottom='off', labelleft='off')
#        # Label segments by sensor bay and segment number.
#        for slot in self.fits_files:
#            seg_coords = list(self._amp_coords[slot].values())[-8]
#            xmin, xmax, ymin, ymax = seg_coords
#            xx = float(xmax + xmin)/2./float(self.nx)
#            if flipx:
#                xx = 1 - xx
#            yy = 1. - (float(ymax - ymin)*0.05 + ymin)/float(self.ny)
#            if rotate180:
#                xx = 1 - xx - 7*np.abs(xmax - xmin)/float(self.nx)
#                yy = 1 - yy + 1.9*np.abs(ymax - ymin)/float(self.ny)
#            plt.annotate('%s' % slot,
#                         (xx, yy), xycoords='axes fraction',
#                         size='x-small', horizontalalignment='center',
#                         verticalalignment='center', color=textcolor)
#            for amp, seg_coords in list(self._amp_coords[slot].items()):
#                xmin, xmax, ymin, ymax = seg_coords
#                xx = float(xmax + xmin)/2./float(self.nx)
#                if flipx:
#                    xx = 1. - xx
#                if amp <= 8:
#                    yy = 1. - (float(ymax - ymin)*0.85 + ymin)/float(self.ny)
#                else:
#                    yy = 1. - (float(ymax - ymin)*0.15 + ymin)/float(self.ny)
#                if rotate180:
#                    xx = 1 - xx
#                    yy = 1 - yy
#                plt.annotate('%s' % imutils.channelIds[amp],
#                             (xx, yy), xycoords='axes fraction',
#                             size='x-small', horizontalalignment='center',
#                             verticalalignment='center', color=textcolor)
#        plt.annotate(annotation, (1, -0.1), xycoords='axes fraction',
#                     horizontalalignment='right', verticalalignment='bottom')
        return fig
