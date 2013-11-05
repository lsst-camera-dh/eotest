#include "ndarray.h"
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
   
   afwImage::Image<float>
   ImageTools::applyCTI(const afwImage::Image<float> & input,
                        const afwGeom::Box2I & serial_overscan,
                        double pcti, double scti, bool verbose) {
      if (pcti == 0 && scti == 0) {
         return input;
      }

      double pcte(1. - pcti);
      double scte(1. - scti);

      // Prepare a working copy and remove temporarily remove bias.
      bool deep_copy;
      afwImage::Image<float> 
         work_image(afwImage::Image<float>(input, deep_copy=true));
                                           
      const afwImage::Image<float>
         biassec(afwImage::Image<float>(input, serial_overscan,
                                        afwImage::LOCAL, deep_copy=false));
      float bias_med(afwMath::makeStatistics(biassec,
                                             afwMath::MEDIAN).getValue());
      work_image -= bias_med;
      
      int nx(work_image.getWidth());
      int ny(work_image.getHeight());

      // Access a view of the ndarray in the working image.
      ndarray::Array<float, 2, 1> work_array(work_image.getArray());

      // Create a ndarray to store output.
      ndarray::Array<float, 2, 1> outarray = 
         ndarray::allocate(ndarray::makeVector(ny, nx));

      if (pcti != 0) {
         // Perform parallel charge transfer.
         // Loop over rows.
         for (int j(0); j < ny-2; j++) {
            if (j % 100 == 0 and verbose) {
               pexLogging::TTrace<0>("lsst.eotest.ImageTools.applyCTI", 
                                     "  row %d", j);
            }
            // Copy bottom row to output.
            outarray[j].deep() = pcte*work_array[0];
            // Calculate new shifted frame.
            for (int jj(0); jj < ny-2-j; jj++) {
               work_array[jj].deep() = (pcti*work_array[jj] +
                                        pcte*work_array[jj+1]);
            }
            // Last row just has the deferred charge.
            work_array[ny-1].deep() = pcti*work_array[ny-1];
         }
      } else {
         // No PCTI, so just deep copy the bias-subtracted, but
         // otherwise pristine, work_array to the outarray.
         outarray.deep() = work_array;
      }

      if (scti != 0) {
         // Set the working_array to the new outarray, if necessary.
         if (pcti != 0) {
            work_array.deep() = outarray;
         }
         // Perform serial charge transfer.
         // Loop over columns.
         for (int i(0); i < nx-2; i++) {
            if (i % 100 == 0 and verbose) {
               pexLogging::TTrace<0>("lsst.eotest.ImageTools.applyCTI", 
                                     "  column %d", i);
            }
            // Copy last column to output.
            for (int j(0); j < ny; j++) {
               outarray[j][i] = scte*work_array[j][0];
            }
            // Calculate new shifted frame.
            for (int ii(0); ii < nx-2-i; ii++) {
               for (int jj(0); jj < ny; jj++) {
                  work_array[jj][ii] = (scti*work_array[jj][ii]
                                        + scte*work_array[jj][ii+1]);
               }
            }
            // Last column just has the deferred charge.
            for (int j(0); j < ny; j++) {
               work_array[j][nx-1] = scti*work_array[j][nx-1];
            }
         }
      }

      // Output image.
      afwImage::Image<float> output(outarray);

      // Add the readout bias back.
      output += bias_med;

      return output;
   }

} // namespace eotest
} // namespace lsst
         
