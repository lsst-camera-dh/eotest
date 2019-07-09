import numpy as np
import scipy.special


def cte_matrix(npix, cti):
    """
    Compute the CTE matrix so that the apparent charge q_i in the i-th
    pixel is given by

    q_i = Sum_j cte_matrix_ij q0_j

    where q0_j is the initial charge in j-th pixel.  The corresponding
    python code would be

    >>> cte = cte_matrix(npix, cti)
    >>> qout = numpy.dot(cte, qin)

    Parameters
    ----------
    npix : int
        Total number of pixels in either the serial or parallel
        directions.
    cti : float
        The charge transfer inefficiency.

    Returns
    -------
    numpy.array
        The npix x npix numpy array containing the CTE matrix.

    Notes
    -----
    This implementation is based on
    Janesick, J. R., 2001, "Scientific Charge-Coupled Devices", Chapter 5,
    eqs. 5.2a,b.

    """
    old_settings = np.seterr(invalid='ignore', under='ignore')
    # Fill CTE matrix using the exact expression
    my_matrix = np.zeros((npix, npix), dtype=np.float)
    for i in range(npix):
        jvals = np.arange(i+1)
        my_matrix[i, :i+1] = (scipy.special.binom(i, i - jvals)
                              *cti**(i - jvals)*(1 - cti)**(jvals + 1))

    np.seterr(**old_settings)
    np.seterr(under='ignore')
    # For a large number of transfers, the binomial coefficient
    # diverges while the cti**(1-jvals) factor underflows, and the
    # resulting element can be a nan.  Replace those entries with the
    # Poisson approximation.
    ivals = np.array([np.ones(npix)*i for i in range(npix)])
    jvals = np.array([np.arange(npix) for _ in range(npix)])
    index = np.where(my_matrix != my_matrix)
    if len(index[0]) > 0:
        lamb = jvals[index]*cti
        kval = ivals[index] - jvals[index]
        my_matrix[index] = (1 - cti)*(np.exp(-lamb)*lamb**(kval)
                                      /scipy.special.factorial(kval))
    np.seterr(**old_settings)
    return my_matrix
