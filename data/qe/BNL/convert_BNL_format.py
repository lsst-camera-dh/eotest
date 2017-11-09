import numpy as np


def convert_BNL_format(old_pd_cal_file, outfile):
    old_data = np.recfromtxt(old_pd_cal_file, skip_header=1,
                             names='truewl, sens, eff, monowl, ccdfrac, foo')

    wl_vals = old_data['monowl']
    pd_sens_vals = old_data['sens']/1e3  # convert to A/W
    pd_ratio_vals = old_data['ccdfrac']

    output = open(outfile, 'w')
    output.write('wl(nm)  cal_pd_sens(A/W)  cal_QE  pd_ratio  mon_pd_charge(pC) cal_pd_charge(pC)\n')
    for wl, pd_sens, pd_ratio in zip(wl_vals, pd_sens_vals, pd_ratio_vals):
        cal_qe = 0
        mon_pd_charge = pd_ratio
        cal_pd_charge = 1
        output.write(
            '%(wl).1f      %(pd_sens).5f        %(cal_qe).3f   %(pd_ratio).6f   %(mon_pd_charge)12.4e     %(cal_pd_charge)12.4e\n' % locals())
    output.close()


if __name__ == '__main__':
    old_pd_cal_file = 'pd_Cal_mar2013.txt'
    outfile = 'pd_Cal_mar2013_v1.txt'
    convert_BNL_format(old_pd_cal_file, outfile)
