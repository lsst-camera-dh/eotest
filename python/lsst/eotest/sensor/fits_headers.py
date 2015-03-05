"""
@brief pyfits-based tools for handling FITS header keywords.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import sys
import numpy as np
import pyfits
from collections import OrderedDict
import lsst.eotest.image_utils as imutils
from AmplifierGeometry import parse_geom_kwd

_module_path = os.environ['EOTEST_DIR']
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
    headers = OrderedDict()
    hdr = pyfits.header.Header()
    for line in open(template):
        # Skip comments and whitespace lines.
        if line[0] == '#' or len(line.strip()) == 0:
            continue
        if line[:3] == 'END':
            if len(headers) == 0:
                # First hdu must be the Primary HDU.
                headers['PRIMARY'] = hdr
            else:
                # Subsequent ones must be extensions with an EXTNAME
                headers[hdr['EXTNAME']] = hdr
            hdr = pyfits.header.Header()
            continue
        data = line.split('=')
        key, value = data[0].strip(), '='.join(data[1:]).strip()
        data = value.split('/')
        value, comment = data[0].strip(), '/'.join(data[1:]).strip()
        hdr[key] = (_cast(value), comment)
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
    missing_headers = []
    #
    for i, extname in enumerate(prototype_headers):
        prototype = prototype_headers[extname]
        if i < 17:
            # Check the first 17 input headers (PHDU + 16 image
            # extensions) by index i, since EXTNAME is often not set in
            # the image extensions.
            try:
                input_hdu = input[i]
            except IndexError:
                missing_headers.append(extname)
                continue
        else:
            # Check for remaining prototype headers by extension name.
            try:
                input_hdu = input[extname]
            except KeyError:
                missing_headers.append(extname)
                continue
        # Check for required keywords.
        missing_keys[extname] = [keyword for keyword in prototype.keys()
                                 if keyword not in input_hdu.header.keys()]
        if missing_keys[extname]:
            report.append("Checking HDU #%i, '%s'. Missing keywords:"
                          % (i, input_hdu.name))
            for key in missing_keys[extname]:
                report.append("  %s" % key)
    if missing_headers:
        report.append("Missing headers:")
        for item in missing_headers:
            report.append("  %s" % item)
    if verbose:
        if report:
            for line in report:
                print line
        else:
            print "No missing keywords or extensions"
    return missing_keys

def check_noao_keywords(infile, verbose=True):
    defects = []
    input = pyfits.open(infile)
    geom_kwd = lambda x : parse_geom_kwd(input[0].header[x])
    # Sanity check for DETSIZE in PHDU
    try:
        pdetsize = geom_kwd('DETSIZE')
        if (pdetsize['xmin'] != 1 or pdetsize['ymin'] != 1 or
            pdetsize['xmax'] <= pdetsize['xmin'] or
            pdetsize['ymax'] <= pdetsize['ymin']):
            defects.append("Primary HDU DETSIZE fails sanity check: %s" 
                           % input[0].header['DETSIZE'])
    except KeyError:
        pdetsize = None
        defects.append("Primary HDU: missing DETSIZE keyword")
    for extnum in imutils.allAmps:
        geom_kwd = lambda x : parse_geom_kwd(input[extnum].header[x])
        try:
            detsize = geom_kwd('DETSIZE')
        except KeyError:
            detsize = None
            defects.append("HDU %i: missing DETSIZE keyword" % extnum)
        try:
            detsec = geom_kwd('DETSEC')
        except KeyError:
            detsec = None
            defects.append("HDU %i: missing DETSEC keyword" % extnum)
        try:
            datasec = geom_kwd('DATASEC')
        except KeyError:
            datasec = None
            defects.append("HDU %i: missing DATASEC keyword" % extnum)

        if (detsize is not None and 
            (detsize['xmin'] != 1 or detsize['ymin'] != 1 or 
             detsize['xmax'] <= detsize['xmin'] or
             detsize['ymax'] <= detsize['ymin'])):
            defects.append("HDU %i DETSIZE fails sanity check: %s" 
                           % (extnum, input[extnum].header['DETSIZE']))
        if (pdetsize is not None and detsize is not None and
            (detsize['xmax'] >= pdetsize['xmax'] or
             detsize['ymax'] >= pdetsize['ymax'])):
            defects.append("HDU %i DETSIZE too big relative to PHDU values: %s"
                           % (extnum, input[extnum].header['DETSIZE']))
        if datasec is not None and datasec['xmin'] <= 1:
            defects.append("HDU %i, No prescan pixels implied by DATASEC: %s"
                           % (extnum, input[extnum].header['DATASEC']))
        if (datasec is not None and detsec is not None and 
            datasec['ymax'] != detsize['ymax']/2.):
            defects.append("HDU %i, Inconsistent parallel DATASEC: %s"
                           % (extnum, input[extnum].header['DATASEC']))
        if (datasec is not None and detsec is not None and 
            (abs(datasec['xmin'] - datasec['xmax'])
             != abs(detsec['xmin'] - detsec['xmax']) or
             abs(datasec['ymin'] - datasec['ymax'])
             != abs(detsec['ymin'] - detsec['ymax']))):
            defects.append("HDU %i, inconsistent DETSEC and DATASEC sizes: %s %s" 
                           % (extnum, input[extnum].header['DETSEC'],
                              input[extnum].header['DATASEC']))
    if verbose and defects:
        for item in defects:
            print item
    return defects
