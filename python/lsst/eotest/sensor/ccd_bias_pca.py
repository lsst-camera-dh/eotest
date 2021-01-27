"""
Module to use PCA modeling of bias frames.   This code is based on a
jupyter notebook from Andrew Bradshaw.
"""
import pickle
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from sklearn.decomposition import PCA
import lsst.eotest.image_utils as imutils
from lsst.eotest.fitsTools import fitsWriteto
from .AmplifierGeometry import makeAmplifierGeometry


__all__ = ['CCD_bias_PCA']


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
        self.pca_bias_file = None

    def compute_pcas(self, fits_files, outfile_prefix, amps=None,
                     verbose=False, fit_full_segment=True, sigma=3,
                     use_median=True):
        """
        Compute mean bias and PCA models of serial and parallel
        overscans using a list of bias frame FITS files for a
        particular CCD.

        Parameters
        ----------
        fits_files: list
            List of bias frame FITS files for a single CCD.
        outfile_prefix: str
            Prefix of output files containing the mean/median bias frame,
            `f'{outfile_prefix}_pca_bias.fits'`, and the pickle file
            containing the PCA models, `f'{outfile_prefix}_pca_bias.pickle'`.
        amps: list-like [None]
            A list of amps to model. If None, then do all amps in the CCD.
        verbose: bool [False]
            Flag to print the progress of computing PCAs for each amp.
        fit_full_segment: bool [True]
            Use the full amplifier segment in deriving the PCAs.  If False,
            then use the parallel and serial overscan regions.
        sigma: int [3]
            Value to use for sigma-clipping the amp-level images that
            are included in the training set.
        use_median: bool [True]
            Compute the median of the stacked images for the mean_amp
            image.  If False, then compute the mean.
        """
        amp_geom = makeAmplifierGeometry(fits_files[0])
        self.x_oscan_corner = amp_geom.imaging.getEndX()
        self.y_oscan_corner = amp_geom.imaging.getEndY()
        if amps is None:
            amps = imutils.allAmps(fits_files[0])
        with fits.open(fits_files[0]) as mean_bias_frame:
            for amp in amps:
                if verbose:
                    print(f"amp {amp}")
                amp_stack = get_amp_stack(fits_files, amp)
                pcax, pcay, mean_amp \
                    = self._compute_amp_pcas(amp_stack,
                                             fit_full_segment=fit_full_segment,
                                             verbose=verbose,
                                             sigma=sigma,
                                             use_median=use_median)
                self[amp] = pcax, pcay
                mean_bias_frame[amp].data = mean_amp
            self.pca_bias_file = f'{outfile_prefix}_pca_bias.fits'
            fitsWriteto(mean_bias_frame, self.pca_bias_file, overwrite=True)
        pickle_file = f'{outfile_prefix}_pca_bias.pickle'
        self.to_pickle(pickle_file)
        return pickle_file, self.pca_bias_file

    def _compute_amp_pcas(self, amp_stack, fit_full_segment=True,
                          verbose=False, sigma=3, use_median=True):
        # Compute the mean bias image from the stack of amp data.
        if use_median:
            mean_amp = np.median(np.array(amp_stack), axis=0)
        else:
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

    @staticmethod
    def read_model(pca_model_file, pca_bias_file):
        """
        Read in the PCA model and associated PCA bias frame for computing
        the bias corrections.

        Parameters
        ----------
        pca_model_file: str
            Pickle file containing the PCA model of the bias correction.
        pca_bias_file: str
            FITS file containing the mean images of each amplifier that
            were used to fit the PCA model.

        Returns
        -------
        CCD_bias_PCA object with the pca_bias_file FITS file explicitly
        set.
        """
        my_instance = CCD_bias_PCA.read_pickle(pca_model_file)
        my_instance.pca_bias_file = pca_bias_file
        return my_instance

    def pca_bias_correction(self, amp, image_array):
        """
        Compute the bias model based on the PCA fit.  This should be
        subtracted from the raw data for the specified amp in order to
        apply an overscan+bias correction.

        Parameters
        ----------
        amp: int
            Amplifier for which to compute the correction.
        image_array: numpy.array
            Array containing the pixel values for the full segment of
            the specified amp.

        Returns
        -------
        numpy.array: Array with the pixel values of the computed correction
            for the full segment.

        """
        pcax, pcay = self[amp]
        mean_amp = fits.getdata(self.pca_bias_file, amp).astype('float')

        imarr = image_array - mean_amp
        corner_mean = self.mean_oscan_corner(imarr)
        imarr -= corner_mean

        # Build the serial PCA-based model using the parallel overscan
        _ = pcax.transform(imarr[self.y_oscan_corner:, self.xstart:])
        projx = pcax.inverse_transform(_)
        serial_model = np.mean(projx, axis=0)

        # Build the parallel PCA-based model using the serial overscan
        # after subtracting the serial_model
        imarr[self.ystart:, self.xstart:] -= serial_model
        _ = pcay.transform(imarr[self.ystart:, self.x_oscan_corner:].T)
        projy = pcay.inverse_transform(_)
        parallel_model = np.mean(projy, axis=0)

        bias_model = mean_amp + corner_mean
        bias_model[self.ystart:, self.xstart:] += serial_model
        bias_model = (bias_model.T + parallel_model).T
        bias_model += self.mean_oscan_corner(image_array - bias_model)

        return bias_model

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
                hdus[amp].data = self.pca_bias_correction(amp, hdus[amp].data)
            fitsWriteto(hdus, outfile, overwrite=True)

        if residuals_file is not None:
            with fits.open(raw_file) as resids, fits.open(outfile) as bias:
                for amp in range(1, 17):
                    resids[amp].data = (np.array(resids[amp].data, dtype=float)
                                        - bias[amp].data)
                resids.writeto(residuals_file, overwrite=True)
