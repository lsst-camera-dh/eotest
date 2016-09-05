#ifndef lsst_eotest_ImageTools_h
#define lsst_eotest_ImageTools_h

namespace lsst {

   namespace afw {
      namespace image {
         template<typename T> class Image;
      }
      namespace geom {
         class Box2I;
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

      private:

      };

   } // namespace eotest
} // namespace lsst

#endif // lsst_eotest_ImageTools_h
