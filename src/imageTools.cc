#include <iostream>
#include "ndarray.h"
#include "lsst/afw/image/Image.h"
#include "lsst/eotest/ImageTools.h"

namespace lsst {
namespace eotest {

   lsst::afw::image::Image<float> 
   ImageTools::rebin(const lsst::afw::image::Image<float> & input,
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
      return lsst::afw::image::Image<float>(outarray);
   }
   
   lsst::afw::image::Image<float>
   ImageTools::applyCTI(const lsst::afw::image::Image<float> & input,
                        double cti) {
      if (cti == 0) {
         return input;
      }
      return input;
   }

} // namespace eotest
} // namespace lsst
         
