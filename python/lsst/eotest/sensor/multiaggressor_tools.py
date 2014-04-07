"""
Data and tools to generate simulation input for multiaggressor crosstalk
pattern.
"""
import lsst.eotest.image_utils as imutils

_xpos = [3856, 3364, 2872, 2380, 1888, 1396,  904,  412,
         260,  792, 1324, 1856, 2388, 2920, 3452, 3984]
_ypos = [1980, 1215,  450, 1880, 1115,  350, 1780, 1015,
         3754, 2324, 3089, 3854, 2424, 3189, 3954, 2524]

class AmpCoords(object):
    """
    Convert pixel coordinates in detector (CCD) geometry to amplifier
    coordinates.
    """
    def __init__(self, nx, ny):
        self.nx, self.ny = nx, ny
    def __call__(self, xx, yy):
        if yy <= self.ny:
            amp = 8 - int((xx - 1)/self.nx)
        else:
            amp = int((xx - 1)/self.nx) + 9
        if amp <= 8:
            x = ((xx - 1) % self.nx) + 1
            y = yy
        else:
            x = self.nx + 1 - (((xx - 1) % self.nx) + 1)
            y = 2*self.ny - yy + 1
        return amp, x, y

def multiaggressor_amplifier_coords(nx, ny, xpos=_xpos, ypos=_ypos):
    """
    Return dictionaries (keyed by amplifier) of x, y pixel locations
    of aggressor images translated from detector coordinates to
    amplifier coordinates.
    """
    amp_coords = AmpCoords(nx, ny)
    x, y = {}, {}
    for amp, xx, yy in zip(imutils.allAmps, xpos, ypos):
        my_amp, x[amp], y[amp] = amp_coords(xx, yy)
        if amp != my_amp:
            raise RuntimeError("multiaggressor amp mismatch")
    return x, y
