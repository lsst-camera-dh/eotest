import numpy as np
from lsst.eotest.sensor.PhotodiodeResponse import PhotodiodeResponse

def convert_Harvard_files(wlscan_file, ccd_cal_file, sph_cal_file,
                          outfile):
    ccd_pd = PhotodiodeResponse(ccd_cal_file)
    sph_pd = PhotodiodeResponse(sph_cal_file)
    wlscan = np.recfromtxt(wlscan_file, skip_header=1,
                           names='wl, sph_pos, ccd_pos')
    output = open(outfile, 'w')
    output.write('wl(nm)  cal_pd_sens(A/W)  cal_QE  pd_ratio  mon_pd_charge(pC) cal_pd_charge(pC)\n')
    for wl, sph_current, ccd_current in zip(wlscan['wl'], wlscan['sph_pos'],
                                            wlscan['ccd_pos']):
        pd_sens = sph_pd(wl)
        pd_ratio = sph_current/ccd_current
        cal_qe = 0
        mon_pd_charge = pd_ratio
        cal_pd_charge = 1
        output.write('%(wl).1f      %(pd_sens).5f        %(cal_qe).3f   %(pd_ratio).6f   %(mon_pd_charge)12.4e     %(cal_pd_charge)12.4e\n' % locals())
    output.close()

if __name__ == '__main__':
    wlscan_file = 'WLscan.txt'
    ccd_cal_file = 'OD143.csv'
    sph_cal_file = 'OD142.csv'
    outfile = 'pd_Cal_Harvard.txt'
    convert_harvard_files(wlscan_file, ccd_cal_file, sph_cal_file, outfile)
