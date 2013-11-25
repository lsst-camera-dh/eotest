IDL LSST Data Simulation:

LSST image data are stored in FITS format files. Each FITS format file
has 16 extensions that include image data. Each FITS file extension contains 
the pixel data from a single CCD output gernated from a single amplifier
and containing data from 1 of the 16 LSST CCD outputs.

Files that begin with 'seg_' work on a single 'segment' representing
the output of a single amplifier out of the 16 amplifiers present on
an LSST CCD. These work on a single FITS file extension.

Files that begin with 'lsst_' work on a complete LSST image file, including
all 16 image extensions.

Function/Procedure descriptions:

function seg_valid, infile, exten, WRITE=write, MESSAGE=message

check that the file exists, is readable or writable as required and 
that the extension number specified makes sense.

  
SEG_ROWVAL
Set the value of a single row in a FITS file extension. Optionally 
include only a certain range of columns in the segment.


SEG_POISSON
Add poisson noise to the imaging array pixels in a FITS file extension.

SEG_PIXVAL
Set the value of pixels in a FITS file extension.


SEG_OVERVAL
Set the value of all prescan and overscan area pixels in a FITS 
file extension.

SEG_NOISE
Add noise to the pixels in a FITS file extension. Optionally works on 
image area, parallel overscan, serial prescan, and serial overscan

SEG_IMVAL
Set the value of pixels in a FITS file extension. Optionally works on 
image area, parallel overscan, serial prescan, and serial overscan

SEG_FE55_ADD
Add simulated fe55 x-ray charge packets to an image segment

SEG_EXPOSE
Add simulated flux to an image segment

SEG_CONST_ADD
Add a constant value to every pixel in an image segment

SEG_COLVAL
Set the value of a column in a FITS file extension.

SEG_BIAS_SET
Set the bias offset level in a segment to a given mean value

SEG_APPLY_DEFECTS
Add defects to an image segment using a defect mask image as guide

SEG_ADD_SIGNAL
Add signal to an image segment.

SEG_ADD_DEFECT
Put a defect code in place of a single pixel in a FITS file extension.

LSST_PAIR_CREATE
Create a pair of 16 extension FITS files 

LSST_MASK_IMAGE
Create a 16 extension FITS file defect mask image

LSST_IMAGE_CREATE
Create a 16 extension FITS file 

LSST_FLAT_STACK
Create COUNT pairs of 16 extension FITS file 


LSST_FLAT_PAIR
Create a pair of 16 extension FITS files 

LSST_FLAT_IMAGE
Create a 16 extension FITS file simulated flat field image

LSST_FLAT_IMAGE
Create a 16 extension FITS file simulated flat field image


LSST_DEFECT_IMAGE
Create a 16 extension FITS file defect mask image with defect pixels 
identified via bit mask. If no bits are set the pixel is neither masked 
nor a defect. Bits are defined in the following manner:
0 : masked pixel, masked for any reason, not used in calculations
1 : dead pixel: pixel does not respond to light
2 : dark pixel: pixel has less response to light vs other pixels
3 : dark column : pixel in a column defined as dark due to dark defects
4 : bright pixel : pixel has excess dark current
5 : bright column : pixel in a column defined as bright due to bright defects
6 : trap pixel : a pixel that traps charge in excess of some limit
7 : trap column : pixel in a column that contains a trap
   etc. 

Pixels may fall into multiple categories. Some defect examples:  
Masked pixel : 0x1
Dead Pixel : 0x2
Dark pixel not in a dark column : 0x4
Functional pixel in a dark column : 0x8
Dark pixel in a dark column : 0x4 || 0x8 = 0xc
etc.

LSST_DATASET
Create a set of simulated image data files 

LSST_DARK_IMAGE
Create a 16 extension FITS file simulated flat field image

LSST_BIAS_CREATE
Create a 16 extension FITS file simulated bias image

LSST_APPLY_DEFECTS
Add defects to an image segment using a defect mask image as guide 

