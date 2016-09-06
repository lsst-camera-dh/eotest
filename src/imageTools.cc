#include <iostream>
#include "ndarray.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/eotest/ImageTools.h"

namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace pexLogging = lsst::pex::logging;

namespace lsst {
namespace eotest {

   afwImage::Image<float>
   ImageTools::rebin(const afwImage::Image<float> & input,
                     unsigned int binsize) {
      if (binsize == 1) {
         return input;
      }

      // Get Array reference to pixel data in image.
      const ndarray::Array<const float, 2, 1> & inarray(input.getArray());

      int nx0(input.getWidth());
      int ny0(input.getHeight());

      // Create new Array object to contain rebinned data.
      int nx(nx0/binsize);
      int ny(ny0/binsize);
      ndarray::Array<float, 2, 1> outarray = 
         ndarray::allocate(ndarray::makeVector(ny, nx));

      // Initialize with zeros
      for (int jout(0); jout < ny; jout++) {
         for (int iout(0); iout < nx; iout++) {
            outarray[jout][iout] = 0;
         }
      }

      // Loop over input pixels, adding content to appropriate bin in
      // output array.
      for (int j(0); j < ny0; j++) {
         int jout(j/binsize);
         if (jout < ny) {
            for (int i(0); i < nx0; i++) {
               int iout(i/binsize);
               if (iout < nx) {
                  outarray[jout][iout] += inarray[j][i];
               }
            }
         }
      }

      // Cast as Image<float> and return.
      return afwImage::Image<float>(outarray);
   }

} // namespace eotest
} // namespace lsst
