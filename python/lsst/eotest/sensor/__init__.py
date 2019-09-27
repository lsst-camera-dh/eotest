"""
Modules to perform sensor characterization from laboratory data.
"""
from __future__ import absolute_import
import warnings
#
# Low-level classes and functions
#
from .MaskedCCD import MaskedCCD, add_mask_files
from .TaskParser import TaskParser
from .EOTestResults import EOTestResults
from .EOTestPlots import EOTestPlots, plot_flat, fe55_zoom
from .EOTestReport import EOTestReport
from .AmplifierGeometry import *
from .fe55_psf import PsfGaussFit
from .fe55_gain_fitter import fe55_gain_fitter
from .Fe55GainFitter import Fe55GainFitter
from .fe55_yield import Fe55Yield
from .Fe55PixelStats import Fe55PixelStats
from .BrightPixels import BrightPixels
from .DarkPixels import DarkPixels
from .read_noise import noise_dists
from .eperTask import EPERTask
from .DetectorResponse import DetectorResponse
from .prnu import prnu
from .Traps import Traps
from .TrapFinder import TrapFinder
from .rolloff_mask import rolloff_mask, pixel_counts
from .generate_mask import generate_mask
from .cte_matrix import cte_matrix
from .total_noise_histograms import *
from .tearing_statistics import *
from .divisidero_tearing import *
#
# Pipe tasks
#
from .fe55Task import Fe55Task
from .brightPixelsTask import BrightPixelsTask
from .darkPixelsTask import DarkPixelsTask
from .readNoiseTask import ReadNoiseTask
from .darkCurrentTask import DarkCurrentTask
from .crosstalkTask import CrosstalkTask
from .cteTask import CteTask, superflat
from .prnuTask import PrnuTask
from .trapTask import TrapTask
from .qeTask import QeTask
from .flatPairTask import FlatPairTask, find_flat2
from .linearityTask import LinearityTask
from .ptcTask import PtcTask
from .persistenceTask import PersistenceTask
from .fe55CteTask import Fe55CteTask
from .BFTask import BFTask
from .overscanTask import OverscanTask
try:
    from .spotTask import SpotTask
except Exception as eobj:
    message = '\nSpotTask import raised a ModuleNotFoundError:\n' + str(eobj)
    warnings.warn(message)
#
# Turn off debug messages emitted by LSST Stack
#
try:
    import lsst.log
except ImportError:
    pass
else:
    try:
        lsst.log.setLevel(lsst.log.getDefaultLoggerName(), lsst.log.INFO)
    except AttributeError:
        lsst.log.setLevel(lsst.log.getDefaultLogger().getName(), lsst.log.INFO)
