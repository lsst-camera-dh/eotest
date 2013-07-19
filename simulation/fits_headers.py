import os
import numpy as np
import pyfits

import simulation

_module_path = os.path.dirname(simulation.__file__)
_tplfile = os.path.join(_module_path, 'fits_header_template.txt')

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

def fits_headers(template=_tplfile):
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

if __name__ == '__main__':
    import image_utils as imutils
    
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
