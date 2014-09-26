#script to check all FITS headers in a directory 

import re
import sys, traceback
import glob, warnings
import pyfits as pyf
import re, os, fnmatch
import argparse

headerVersion = 1
testTypeList = ['DARK', 'FE55', 'FLAT', 'LAMBDA', 'NOISE', 'PPUMP', 'SFLAT', 'SPOT', 'XRAY', 'XTALK']

def get_input_files(dir, r):
    """ Create list of input files given directory and recursive flag. """

    if r==True:
        #recursively find FITS files in specified directory
        files = []
        for root, dirnames, filenames in os.walk(dir):
          for filename in fnmatch.filter(filenames, '*.fits'):
              files.append(os.path.join(root, filename))
    
    else:
        #find FITS files in this directory only
        files = glob.glob(os.path.join(dir,'*.fits'))

    return files

def update_version(hdr):

    # Check to see if HeaderVersion keyword already exists
    # If so, and version is sufficiently updated, exit
    # Otherwise, update the header version number
    if 'HDRVER' in hdr.keys():
        ver = hdr['HDRVER']
        if ver >= headerVersion:
            print 'Skipping file, Header at Version: ', ver
            return                   
        else:
            hdr['HDRVER'] = headerVersion
            # If the header version keyword does not yet exist, add it
    else:
        hdr['HDRVER'] = (headerVersion, 'Header Version Number')


def extract_path(filename):
    """ Extract the filename and test type from a path """
    
    ret_type = 'UNKNOWN'
    sensor_id = 'UNKNOWN'

    # Normalize the path
    normfilename = os.path.abspath(filename)
    orgpath, fname = os.path.split(normfilename)

    subdirs = orgpath.split(os.sep)
    upperSubDirs = [subStr.upper() for subStr in subdirs]
    rplSubDirs = [subStr.replace('SUPERFLAT','SFLAT') for subStr in upperSubDirs]
    
    print upperSubDirs
    print rplSubDirs

    for dir in rplSubDirs:
        for type in testTypeList:
            if type in dir:
                ret_type = type         

        idMatch = re.match('[0-9][0-9][0-9][-|_][0-9][0-9]', dir)
        if idMatch:
            sensor_id = idMatch.group() 

    print fname, ret_type, sensor_id
    return fname, ret_type, sensor_id

def get_id_and_type(filename):
    """ Extract sensor id and type of file from the filename. """  

    imgTypeList = ['BIAS', 'DARK', 'FLAT', 'PPUMP', 'QE', 'SFLAT']

    imgType = 'UNKNOWN'
    testType = 'UNKKNOWN'
    seqNum = 0
    sensor_id = 'UNKNOWN'

    try:
        fname, type, sensor_id = extract_path(filename)

        # Reset FE55 to XRAY if necessary
        if type.upper() in testTypeList:
            testType = type.upper()
            if testType == 'FE55':
                testType = 'XRAY'

        subList = fname.split('_')
        upperCaseList = [subStr.upper() for subStr in subList]

        print upperCaseList

        imgMatch = []
        endOfFile = '.FITS'
        # Find the image type and seq number from the filename 
        for str in upperCaseList:
            # We have the ending file string, where seq number might be
            if endOfFile in str:
                seqMatch = re.match('[0-9][0-9][0-9]', str)
                if seqMatch:
                    # Extract the sequence number
                    seqStr = seqMatch.group()
                    seqNum = int(seqStr)

            for img in imgTypeList:
                if img in str:
                    imgMatch.append(img)

        # if we find more than 1 image type, prefer BIAS and QE over others
        if len(imgMatch) > 1:
            if 'BIAS' in imgMatch:
                imgType = 'BIAS'
            if 'QE' in imgMatch:
                imgType = 'QE'
        elif len(imgMatch) > 0:
            imgType = imgMatch[0]

    except:
        print "failed to extract id and type from filename for file: ", filename
        traceback.print_exc(file=sys.stdout)

    return sensor_id, imgType, testType, seqNum
     

def extendHeader(files):
    """ Add new draft keywords. """

    for filename in files:
        try:
            print filename
            #open file in update mode
            f = pyf.open(filename, memmap=True, ignore_missing_end=True, mode='update')
            #grab primary header
            hdr = f[0].header

            update_version(hdr)

        except:
            print 'Failed to open file: ', filename
            traceback.print_exc(file=sys.stdout)
            continue 
            
        # Assuming this is an env var
        if 'CCD_MANU' not in hdr.keys():
            try:
                vendor = os.environ['CCD_MANU']
                hdr['CCD_MANU'] = vendor
            except: 
                print 'Failed to update CCD_MANU for file: ', filename

        try:
            sensor_id, imgType, testType, seqNum = get_id_and_type(filename)

            if 'LSST_NUM' not in hdr.keys():
                hdr['LSST_NUM'] = sensor_id
                print "setting LSST_NUM ", sensor_id

            if 'IMGTYPE' not in hdr.keys():
                hdr['IMGTYPE'] = imgType
                print "setting IMGTYPE: ", imgType

            if 'TESTTYPE' not in hdr.keys():
                hdr['TESTTYPE'] = testType
                print "setting TESTTYPE: ", testType

            if 'SEQNUM' not in hdr.keys():
                hdr['SEQNUM'] = seqNum
                print "setting SEQNUM: ", seqNum
        except:
            traceback.print_exc(file=sys.stdout)
            continue

        f.close()
        print filename + " done."

def fixHeader(files):
    
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

            update_version(hdr)
        
            #truncate keithley idn because HIERARCH isn't compatible with CONTINUE
            if 'K_PHOT.IDN' in hdr.keys():
                idn1 = hdr['K_PHOT.IDN']
                hdr['K_PHOT.IDN'] = idn1[0:50]
            if 'K_BIAS.IDN' in hdr.keys():
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
    parser.add_argument('-e', '--extend', help= \
       "Add additional useful keywords", action='store_true', default=False)
    args = parser.parse_args()

    files = get_input_files(args.dir, args.recursive)

    if args.extend:
        extendHeader(files)
    else:
        fixHeader(files)
  
