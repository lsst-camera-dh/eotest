"""
Module to use PCA modeling of bias frames.   This code is based on a
jupyter notebook from Andrew Bradshaw.
"""
import pickle
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from sklearn.decomposition import PCA
import lsst.geom
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .AmplifierGeometry import makeAmplifierGeometry


__all__ = ['CCD_bias_PCA', 'make_overscan_frame']


def get_amp_stack(fits_files, amp):
    """Get a list of numpy arrays of pixel data for the specified amp.

    Parameters
    ----------
    fits_files: list
        List of FITS filenames.
    amp: int
        Desired amp.

    Returns
    -------
    list of numpy arrays.
    """
    amp_stack = []
    for item in fits_files:
        with fits.open(item) as hdus:
            amp_stack.append(np.array(hdus[amp].data, dtype=float))
    return amp_stack


class CCD_bias_PCA(dict):
    """
    Class to compute mean bias frames and PCA-based models of the overscan
    subtraction derived from an ensemble of bias frames.
    """
    def __init__(self, std_max=20, xstart=2, ystart=0, ncomp_x=6, ncomp_y=8):
        """
        Parameters
        ----------
        std_max: float [15]
            Cutoff for stdev of amp ADU values for inclusion in the PCA
            training set.
        xstart: int [2]
            Starting pixel for the PCA modeling in the serial direction.
        ystart: int [0]
            Starting pixel for the PCA modeling in the parallel direction.
        ncomp_x: int [6]
            Number of PCA components to fit in the serial direction.
        ncomp_y: int [8]
            Number of PCA components to fit in the parallel direction.
        """
        super().__init__()
        self.std_max = std_max
        self.xstart = xstart
        self.ystart = ystart
        self.ncomp_x = ncomp_x
        self.ncomp_y = ncomp_y
        self.x_oscan_corner = None
        self.y_oscan_corner = None

    def compute_pcas(self, fits_files, amps=None, verbose=False,
                     fit_full_segment=True):
        """
        Compute mean bias and PCA models of serial and parallel
        overscans using a list of bias frame FITS files for a
        particular CCD.

        Parameters
        ----------
        fits_files: list
            List of bias frame FITS files for a single CCD.
        amps: list-like [None]
            A list of amps to model. If None, then do all amps in the CCD.
        verbose: bool [False]
            Flag to print the progress of computing PCAs for each amp.
        fit_full_segment: bool [True]
            Use the full amplifier segment in deriving the PCAs.  If False,
            then use the parallel and serial overscan regions.
        """
        amp_geom = makeAmplifierGeometry(fits_files[0])
        self.x_oscan_corner = amp_geom.imaging.getEndX()
        self.y_oscan_corner = amp_geom.imaging.getEndY()
        if amps is None:
            amps = imutils.allAmps(fits_files[0])
        for amp in amps:
            if verbose:
                print(f"amp {amp}")
            amp_stack = get_amp_stack(fits_files, amp)
            self[amp] \
                = self._compute_amp_pcas(amp_stack,
                                         fit_full_segment=fit_full_segment,
                                         verbose=verbose)

    def _compute_amp_pcas(self, amp_stack, fit_full_segment=True,
                          verbose=False, sigma=3):
        # Compute the mean bias image from the stack of amp data.
        mean_amp = np.mean(np.array(amp_stack), axis=0)
        if verbose:
            print("np.std(mean_amp):", np.std(mean_amp))

        # Assemble the training set of mean-subtracted images from the
        # stack of raw amplifier data.  Also subtract the mean of the
        # per-amp overscan corner from each image, and apply a noise
        # cut of self.std_max for inclusion in the training set.
        training_set = list()
        for i, _ in enumerate(amp_stack):
            imarr = _.copy()
            imarr -= mean_amp
            imarr -= self.mean_oscan_corner(imarr)
            sigma = np.std(imarr)
            if sigma > self.std_max:
                print('rejected frame:', i, sigma, self.std_max)
                continue
            # Apply sigma clipping to mask pixel defects.
            training_set.append(sigma_clip(imarr, sigma=sigma))
        if verbose:
            print("training set size:", len(training_set))

        # Create the ensemble of profiles in the serial direction.
        if fit_full_segment:
            # This uses the full data segment, rather than just the parallel
            # overscan.
            x_profs = [np.mean(_, axis=0)[self.xstart:] for _ in training_set]
        else:
            # Use the parallel overscan region instead of full
            # segment.
            x_profs = [np.mean(_[self.y_oscan_corner:, :], axis=0)[self.xstart:]
                       for _ in training_set]

        # Run the PCA fit for the serial direction
        pcax = PCA(self.ncomp_x)
        pcax.fit(x_profs)

        # Use the previous serial direction decomposition to do the
        # fitting in the parallel direction.
        y_profs = []
        for imarr in training_set:
            # Build the serial model, using the pcax basis set, and fit
            # to the full segment data in the serial direction.
            _ = pcax.transform(imarr[self.ystart:, self.xstart:])
            proj = pcax.inverse_transform(_)
            serial_model = np.mean(proj, axis=0)

            # Subtract the serial model from the mean-subtracted training image.
            new_imarr = imarr.copy()[self.ystart:, self.xstart:] - serial_model

            # Add the resulting profile to the y-ensemble
            if fit_full_segment:
                y_profs.append(np.mean(new_imarr, axis=1))
            else:
                # Just use the serial overscan data.
                y_profs.append(np.mean(new_imarr[:, self.x_oscan_corner:],
                                       axis=1))


        # Run the PCA fit for the parallel direction.
        pcay = PCA(self.ncomp_y)
        pcay.fit(y_profs)

        return pcax, pcay, mean_amp

    def mean_oscan_corner(self, imarr, buff=0):
        """
        Compute the mean pixel value of the region common to
        the parallel and serial overscan regions.

        Parameters
        ----------
        imarr: numpy.array
            Array of pixel values from the full segment of the amplifier.
        buff: int [2]
            Buffer pixels to offset from the imaging region for computing the
            mean value, e.g., to avoid trailed charge.

        Returns
        -------
        float
        """
        return np.mean(imarr[self.y_oscan_corner + buff:,
                             self.x_oscan_corner + buff:])

    def to_pickle(self, outfile):
        """
        Write the CCD_bias_PCA object as a pickle object.

        Parameters
        ----------
        outfile: str
            Filename of output pickle file.
        """
        with open(outfile, 'wb') as fd:
            pickle.dump(self, fd)

    @staticmethod
    def read_pickle(infile):
        """
        Read a CCD_bias_PCA object from a pickle file.

        Parameters
        ----------
        infile: str
            Filename of input pickle file.

        Returns
        -------
        CCD_bias_PCA object.
        """
        with open(infile, 'rb') as fd:
            my_instance = pickle.load(fd)
        return my_instance

    def make_bias_frame(self, raw_file, outfile, residuals_file=None):
        """
        Construct the PCA model bias frame for one of the bias files
        and optionally write the bias-subtracted file.

        Parameters
        ----------
        raw_file: str
            Filename of raw single CCD FITS file.
        outfile: str
            Filename of the output bias frame.
        residuals_file: str [None]
            Filename of the output bias-subtracted frame. If None, then
            the file is not written.
        """
        with fits.open(raw_file) as hdus:
            for amp in range(1, 17):
                pcax, pcay, mean_amp = self[amp]
                imarr = hdus[amp].data - mean_amp
                corner_mean = self.mean_oscan_corner(imarr)
                imarr -= corner_mean
                bias_model = mean_amp + corner_mean

                # Build the serial PCA-based model using the parallel overscan
                _ = pcax.transform(imarr[self.y_oscan_corner:, self.xstart:])
                projx = pcax.inverse_transform(_)
                serial_model = np.mean(projx, axis=0)

                # Build the parallel PCA-based model using the serial overscan
                # after subtracting the serial_model
                imarr[self.ystart:, self.xstart:] -= serial_model
                bias_model[self.ystart:, self.xstart:] += serial_model
                _ = pcay.transform(imarr[self.ystart:, self.x_oscan_corner:].T)
                projy = pcay.inverse_transform(_)
                parallel_model = np.mean(projy, axis=0)
                imarr[self.ystart:, self.xstart:] \
                    = (imarr[self.ystart:, self.xstart:].T - parallel_model).T
                bias_model = (bias_model.T + parallel_model).T
                bias_model \
                    += self.mean_oscan_corner(hdus[amp].data - bias_model)
                hdus[amp].data = bias_model
            hdus.writeto(outfile, overwrite=True)

        if residuals_file is not None:
            with fits.open(raw_file) as resids, fits.open(outfile) as bias:
                for amp in range(1, 17):
                    resids[amp].data = (np.array(resids[amp].data, dtype=float)
                                        - bias[amp].data)
                resids.writeto(residuals_file, overwrite=True)

def make_overscan_frame(fits_file, outfile=None):
    """
    Use overscan regions to do column-wise, then row-wise overscan
    profile-based overscan frame.
    """
    ccd = MaskedCCD(fits_file)
    par_corner = lsst.geom.Point2I(0, ccd.amp_geom.imaging.getHeight())
    par_extent = lsst.geom.Extent2I(ccd.amp_geom.full_segment.getWidth(),
                                    ccd.amp_geom.parallel_overscan.getHeight())
    par_bbox = lsst.geom.Box2I(par_corner, par_extent)
    ser_corners = ccd.amp_geom.serial_overscan.getCorners()
    ser_extent = lsst.geom.Extent2I(ccd.amp_geom.serial_overscan.getWidth(),
                                    ccd.amp_geom.full_segment.getHeight())
    ser_bbox = lsst.geom.Box2I(ser_corners[0], ser_extent)

    oscan_arrays = dict()
    for amp in ccd:
        image = ccd[amp].getImage()

        parallel = image.Factory(image, par_bbox).array
        p_profile = np.sum(parallel, 0)/parallel.shape[0]

        overscan_array = np.zeros(image.array.shape)

        image.array[:,] -= p_profile
        overscan_array[:,] += p_profile

        serial = image.Factory(image, ser_bbox).array
        s_profile = np.sum(serial, 1)/serial.shape[1]

        for col in range(ccd.amp_geom.full_segment.getWidth()):
            image.array[:serial.shape[0], col] -= s_profile
            overscan_array[:serial.shape[0], col] += s_profile
        oscan_arrays[amp] = overscan_array

    if outfile is not None:
        with fits.open(fits_file) as hdus:
            for amp in ccd:
                hdus[amp].data = oscan_arrays[amp]
            hdus.writeto(outfile, overwrite=True)

    return ccd
