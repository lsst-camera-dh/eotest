#QE Calibration

import numpy as np
import pylab as p
import matplotlib as mpl
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import glob

####################################Setup
########User Editable Parameters

#which amps to measure
amps = range(1, 16)

#list of gain measured for each amp, from Fe55 measurement (just a placeholder here)
gain = np.ones(16)*5.2

#########Output File

#file with qe for each monochromator setting
outfile_all = "qe_output_all.txt"

#file with qe per LSST filter band
outfile_bands = "qe_output_bands.txt"
#########Input Files

#name of file output from qe_getmedians.py
imdatafile = "qe_medians.txt"

#NIST photodiode calibration data
f1 = "OD142.csv"
f2 = "OD143.csv"

#Harvard setup calibration file
f3 = "WLscan.txt"

########Constants
h = 6.626e-34
c = 299792458

#window transmission - will go away once we recalibrate
windowtr = 0.932

#photodiode surface area
pdarea = 3.14159*((1.13/2)*10**(-2))**2

#pixel area
pxarea = 10e-6**2

#######

#e2v measurements
e2vwl = [400, 500, 600, 800, 900, 1000]
e2vqe = [47.5, 82.9, 93.1, 89.8, 75.5, 23.5]


##########
#LSST Specs
lsstwl = [375, 475, 625, 750, 875]
lsstqe = [40, 80, 88, 85, 80]


############Get auxiliary data

#Get NIST calibration data
od142all = np.recfromtxt(f1, delimiter=",", skip_header=5, names="wl, qe, pcterror, none, none2")
od143all = np.recfromtxt(f2, delimiter=",", skip_header=5, names="wl, qe, pcterror, none, none2")

#Get relevant info from NIST calibration files
wl = od142all["wl"]
od142qe = od142all["qe"]
od142err = od142all["pcterror"]

od143qe = od143all["qe"]
od143err = od143all["pcterror"]

#Get relevant info from Harvard setup calib file
alldata = np.recfromtxt(f3, names=True)
wls = alldata["wl"]
intsphere = alldata["intsphere"]
ccdpos = alldata["ccdpos"]


#Calibrate IS to Detector data by dividing out QE of photodiodes
for i in range(len(intsphere)):
    intsphere[i] = intsphere[i]/od143qe[np.where(wl == wls[i])]

for i in range(len(ccdpos)):
    ccdpos[i] = ccdpos[i]/od142qe[np.where(wl == wls[i])]


#Get calibrated fraction of light received at detector position
fractionatccd = ccdpos/intsphere


#Read in image data and photodiode readings
imdata = np.recfromtxt(imdatafile, names=True)
amp = imdata['amp']
imwl = imdata['wavelength']
median = imdata['median']
pd = imdata['photodiode']*-1.0e-12
exptime = imdata['exptime']

#Calibrate PD reading
pdiscal = []
for j in range(len(pd)):
    pdiscal.append(pd[j]/od143qe[np.where(abs(wls - imwl[j]) <= 0.1)])

pdiscal = np.array(pdiscal)


#########Finally calculate QE

f4 = open(outfile_all, "w+")
f4.write("amp\twl\tqe\n")

#Note: only calculates 400-1000nm because that's what we have calib data for at the moment
for a in amps:
    for w in imwl[np.where((amp == a) & (imwl < 1001))]:
        qe = 100 * (median[np.where((amp == a) & (imwl == w))]*h*c*pdarea) / ((1.0/gain[a])*exptime[np.where(imwl == w)]
                                                                              [0]*pdiscal[np.where(imwl == w)][0]*fractionatccd[np.where(abs(wls-w) <= 0.1)]*pxarea*(w*1e-9))
        f4.write('\t'.join([str(a), str(w), str(qe[0]), '\n']))

f4.close()

##########Calculate QE per band

f5 = open(outfile_bands, "w+")
f5.write('\t'.join(['amp', 'u', 'g', 'r', 'i', 'z', 'y', '\n']))

qedata = np.recfromtxt(outfile_all, names=True)
qe_all = qedata["qe"]
qe_amp = qedata["amp"]
qe_wl = qedata["wl"]


for a in amps:
    #u
    if len(qe_all[np.where((qe_amp == a) & (qe_wl >= 321) & (qe_wl <= 391))]) > 0:
        qeu = np.mean(qe_all[np.where((qe_amp == a) & (qe_wl >= 321) & (qe_wl <= 391))])
    else:
        qeu = 'no data'
    #g
    if len(qe_all[np.where((qe_amp == a) & (qe_wl >= 402) & (qe_wl <= 552))]) > 0:
        qeg = np.mean(qe_all[np.where((qe_amp == a) & (qe_wl >= 402) & (qe_wl <= 552))])
    else:
        qeg = 'no data'
    #r
    if len(qe_all[np.where((qe_amp == a) & (qe_wl >= 552) & (qe_wl <= 691))]) > 0:
        qer = np.mean(qe_all[np.where((qe_amp == a) & (qe_wl >= 552) & (qe_wl <= 691))])
    else:
        qer = 'no data'
    #i
    if len(qe_all[np.where((qe_amp == a) & (qe_wl >= 691) & (qe_wl <= 818))]) > 0:
        qei = np.mean(qe_all[np.where((qe_amp == a) & (qe_wl >= 691) & (qe_wl <= 818))])
    else:
        qei = 'no data'
    #z
    if len(qe_all[np.where((qe_amp == a) & (qe_wl >= 818) & (qe_wl <= 922))]) > 0:
        qez = np.mean(qe_all[np.where((qe_amp == a) & (qe_wl >= 818) & (qe_wl <= 922))])
    else:
        qez = 'no data'
    #y
    if len(qe_all[np.where((qe_amp == a) & (qe_wl >= 930) & (qe_wl <= 1070))]) > 0:
        qey = np.mean(qe_all[np.where((qe_amp == a) & (qe_wl >= 930) & (qe_wl <= 1070))])
    else:
        qey = 'no data'
    f5.write('\t'.join([str(a), str(qeu), str(qeg), str(qer), str(qei), str(qez), str(qey), '\n']))
