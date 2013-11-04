#ifndef lsst_eotest_ImageTools_h
#define lsst_eotest_ImageTools_h

namespace lsst {

   namespace afw {
      namespace image {
         template<typename T> class Image;
      }
   }

   namespace eotest {
   
      class ImageTools {

      public:

#ifndef SWIG         
         ImageTools() {}
#endif
         static lsst::afw::image::Image<float>
         rebin(const lsst::afw::image::Image<float> & input,
               unsigned int binsize);

         static lsst::afw::image::Image<float>
         applyCTI(const lsst::afw::image::Image<float> & input, double cti);
         
      private:

      };

   } // namespace eotest
} // namespace lsst

#endif // lsst_eotest_ImageTools_h
