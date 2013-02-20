#script to check all FITS headers in a directory 

import sys, traceback
import glob, warnings
import pyfits as pyf
import re, os, fnmatch
import argparse

headerVersion = 0

def fixHeader(dir, r):
    if r==True:
        #recursively find FITS files in specified directory
        files = []
        for root, dirnames, filenames in os.walk(dir):
          for filename in fnmatch.filter(filenames, '*.fits'):
              files.append(os.path.join(root, filename))
    
    else:
        #find FITS files in this directory only
        files = glob.glob(os.path.join(dir,'*.fits'))
    
    for filename in files:
        try:
            #open file in update mode
            f = pyf.open(filename, memmap=True, ignore_missing_end=True, mode='update')
        
            #reorder DETSEC values
            a = re.compile('\[([0-9]+):([0-9]+),([0-9]+):([0-9]+)]')
            for x in f[1:17]:
                (x1,x2,y1,y2) = a.match(x.header['DETSEC']).groups()
                x.header['DETSEC'] = '[{0}:{1},{2}:{3}]'.format(x2,x1,y1,y2)
        
            #grab primary header
            hdr = f[0].header

            # Check to see if HeaderVersion keyword already exists
            # If so, and version is sufficiently updated, exit
            # Otherwise, update the header version number
            if 'HDRVER' in hdr.keys():
                ver = hdr['HDRVER']
                if ver >= headerVersion:
                    print 'Skipping file: ', filename, ' Header at Version: ', ver
                    continue                   
                else:
                    hdr['HDRVER'] = headerVersion
            # If the header version keyword does not yet exist, add it
            else:
                hdr.update('HDRVER', headerVersion)
        
            #truncate keithley idn because HIERARCH isn't compatible with CONTINUE
            if 'K_PHOT.IDN' in hdr.keys():
                idn1 = hdr['K_PHOT.IDN']
                hdr['K_PHOT.IDN'] = idn1[0:50]
                idn1 = hdr['K_BIAS.IDN']
                hdr['K_BIAS.IDN'] = idn1[0:50]
        
        
            #make an empty header
            newhdr = hdr.copy()
            newhdr.clear()
            #print hdr
            newkeys = hdr.keys()
            #Generate new keywords
            for i in range(len(hdr.keys())):
                #if keyword contains ., replace it with _
                if '.' in hdr.keys()[i]:
                    newkwddot = newkeys[i].replace('.', '_')
                    #print newkwddot
                else:
                    newkwddot = newkeys[i]
                #if in a hierarch, add an extra space before the keyword
                if len(str(newkwddot)) > 8:
                    newkwd = " " + str(newkwddot)
                else:
                    newkwd = newkwddot

                newkeys[i] = newkwd
                #print newkeys[i]
                # write new keywords to new header
                if len(str(hdr.values()[i])) < 80:
                    warnings.resetwarnings()
                    warnings.filterwarnings('ignore', category=UserWarning, append=True)
                    newhdr.append((newkeys[i], hdr.values()[i], hdr.comments[i]))
                    warnings.resetwarnings()
                    warnings.filterwarnings('always', category=UserWarning, append=True)
                else:
                    warnings.resetwarnings()
                    warnings.filterwarnings('ignore', category=UserWarning, append=True)
                    newhdr.append((newkeys[i], str(hdr.values()[i])[0:79], hdr.comments[i]))
                    warnings.resetwarnings()
                    warnings.filterwarnings('always', category=UserWarning, append=True)
        except: 
            print "Failed to update file: ", filename
            f.close()
            traceback.print_exc(file=sys.stdout)
            continue
            
        #update file's actual header
        f[0].header = newhdr
        f.close()
        print filename + " done."
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser( \
       description='Fix FITS headers to comply with LSST DM Stack.')
    parser.add_argument('-d', '--dir', help="directory of files to fix")
    parser.add_argument('-r', '--recursive', help= \
       "whether to search recursively through the directory. Default = False", \
       action='store_true', default = False)
    args = parser.parse_args()

    fixHeader(args.dir, args.recursive)
