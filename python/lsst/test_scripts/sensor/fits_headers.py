"""
@brief pyfits-based tools for handling FITS header keywords.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import sys
import numpy as np
import pyfits

#import simulation
#
#_module_path = os.path.dirname(simulation.__file__)
_module_path = os.environ['TEST_SCRIPTS_DIR']
template_file = os.path.join(_module_path, 'policy', 'fits_header_template.txt')
template_used_file = os.path.join(_module_path, 'policy',
                                  'fits_header_template_used.txt')

def _cast(value):
    """
    Cast input strings to their 'natural' types.  Strip single quotes
    from strings.
    """
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    if value.strip() == 'T':
        return True
    if value.strip() == 'F':
        return False
    return value.strip("'")

def fits_headers(template=template_file):
    """
    Create a set of CCD-level FITS headers according to the FITS
    template file, which is supposed to implement the FITS standard
    for sensors (LCA-10140).
    """
    headers = []
    hdr = pyfits.header.Header()
    for line in open(template):
        # Skip comments and whitespace lines.
        if line[0] == '#' or len(line.strip()) == 0:
            continue
        if line[:3] == 'END':
            headers.append(hdr)
            hdr = pyfits.header.Header()
            continue
        data = line.split('=')
        key, value = data[0].strip(), '='.join(data[1:]).strip()
        data = value.split('/')
        value, comment = data[0].strip(), '/'.join(data[1:]).strip()
        hdr.update(key, _cast(value), comment=comment)
    return headers

def check_keywords(infile, template=template_file, verbose=True):
    """
    Check that the keywords in a the specified FITS header template
    file are present.  The default file is based on the FITS standard
    document for sensors, LCA-10140.

    @return Dictionary of missing keywords by header extension number.
    """
    prototype_headers = fits_headers(template=template)
    input = pyfits.open(infile)
    report = []
    missing_keys = {}
    for i, prototype, input_hdu in zip(range(len(prototype_headers)),
                                       prototype_headers, input):
        missing_keys[i] = [keyword for keyword in prototype.keys()
                           if not input_hdu.header.has_key(keyword)]
        if missing_keys[i]:
            report.append("Checking HDU #%i, '%s'. Missing required keywords:"
                          % (i, input_hdu.name))
            for key in missing_keys[i]:
                report.append("  %s" % key)
    if verbose:
        if report:
            for line in report:
                print line
        else:
            print "No missing keywords"
    return missing_keys

if __name__ == '__main__':
    import lsst.test_scripts.image_utils as imutils
    
    infile = 'fits_header_template.txt'
    phdr, ihdr = fits_headers(infile)

    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header = phdr.copy()
    for amp in imutils.allAmps:
        output.append(pyfits.ImageHDU())
        output[-1].header = ihdr.copy()
        output[-1].header.update('DETSEC', imutils.detsec(amp))
        output[-1].header.update('CHANNEL', amp)
        output[-1].header.update('EXTNAME',
                                 'SEGMENT%s' % imutils.channelIds[amp])
        output[-1].data = np.ones((2022, 542), dtype=np.float32)*amp
        del output[-1].header['BSCALE']
        del output[-1].header['BZERO']
        
    output.writeto('test_output.fits', clobber=True, output_verify='fix',
                   checksum=True)
