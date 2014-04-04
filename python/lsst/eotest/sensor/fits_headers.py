"""
@brief pyfits-based tools for handling FITS header keywords.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import sys
import numpy as np
import pyfits
from collections import OrderedDict

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
    missing_headers = []
    for i, extname in enumerate(prototype_headers):
        prototype = prototype_headers[extname]
        try:
            input_hdu = input[i]
        except IndexError:
            missing_headers.append(prototype['EXTNAME'])
            continue
        missing_keys[i] = [keyword for keyword in prototype.keys()
                           if not input_hdu.header.has_key(keyword)]
        if missing_keys[i]:
            report.append("Checking HDU #%i, '%s'. Missing required keywords:"
                          % (i, input_hdu.name))
            for key in missing_keys[i]:
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
            print "No missing keywords"
    return missing_keys
