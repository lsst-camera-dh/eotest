"""
Characterize CTE using Fe55 cluster asymmetries.
"""
from __future__ import absolute_import, print_function
import matplotlib.pyplot as plt
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .Fe55PixelStats import Fe55PixelStats

class Fe55CteConfig(pexConfig.Config):
    "Configuration for Fe55CteTask"
    output_dir = pexConfig.Field("Output directory", str, default=".")
    direction = pexConfig.Field("Readout direction", str, default="serial")
    verbose = pexConfig.Field("Verbosity flag", bool, default=True)

class Fe55CteTask(pipeBase.Task):
    "Task for characterizing CTE using Fe55 cluster asymmetries."
    ConfigClass = Fe55CteConfig
    _DefaultName = "Fe55CteTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, fe55_files, mask_files):
        "Run the Fe55 pixel asymmetry analysis"
        if self.config.direction == 'serial':
            pixel_coord, pix0, pix1 = 'x', 'p3', 'p5'
        elif self.config.direction == 'parallel':
            pixel_coord, pix0, pix1 = 'y', 'p1', 'p7'
        else:
            raise RuntimeError("Readout direction must be either " +
                               "'serial' or 'parallel'")

        pixel_stats = Fe55PixelStats(fe55_files, mask_files=mask_files,
                                     logger=self.log)

        # Plot the histograms of p3 and p5 (or p1 and p7).
        pixel_stats.pixel_hists(pix0, pix1)
        plt.savefig('%(sensor_id)s_Fe55_%(pix0)s_%(pix1)s_hists.png' % locals())

        # Plot the profiles of p3, p5, and p5-p3 (or p1, p7, and p7-p1)
        # and fit a line to the profile differences.
        fig, results = pixel_stats.pixel_diff_profile(pixel_coord, pix0, pix1)
        plt.savefig('%(sensor_id)s_Fe55_%(pix0)s_%(pix1)s_profiles.png'
                    % locals())

        outfile = '%(sensor_id)s_Fe55_slopes_%(pix0)s_%(pix1)s.txt' % locals()
        with open(outfile, 'w') as output:
            output.write(str(results) + '\n')
        return results
