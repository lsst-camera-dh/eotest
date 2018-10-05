import numpy as np
import scipy.special


def cte_matrix(npix, cti, ntransfers=20, nexact=30):
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
    ntransfers : int, optional
        Maximum number of transfers to consider as contributing to
        a target pixel.
    nexact : int, optional
        Number of transfers to use exact the binomial distribution
        expression, otherwise use Poisson's approximation.

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
    ntransfers = min(npix, ntransfers)
    nexact = min(nexact, ntransfers)
    my_matrix = np.zeros((npix, npix), dtype=np.float)
    for i in range(1, npix):
        jvals = np.concatenate((np.arange(1, i+1), np.zeros(npix-i)))
        index = np.where(i - nexact < jvals)
        j = jvals[index]
        my_matrix[i-1, :][index] \
            = scipy.special.binom(i, j)*(1 - cti)**i*cti**(i - j)
        if nexact < ntransfers:
            index = np.where((i - nexact >= jvals) & (i - ntransfers < jvals))
            j = jvals[index]
            my_matrix[i-1, :][index] \
                = (j*cti)**(i - j)*np.exp(-j*cti)/scipy.special.factorial(i - j)
    return my_matrix
