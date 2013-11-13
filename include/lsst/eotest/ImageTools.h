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

         static lsst::afw::image::Image<float>
         applyCTI(const lsst::afw::image::Image<float> & input,
                  const lsst::afw::geom::Box2I & serial_overscan,
                  double pcti=0, double scti=0, bool verbose=false);

      private:

      };

   } // namespace eotest
} // namespace lsst

#endif // lsst_eotest_ImageTools_h
