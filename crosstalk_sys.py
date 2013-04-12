import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import image_utils as iu
import lsst.afw.math as afwMath
import glob
import numpy as np

#output file name
outfile = "xtalk_sys_output.txt"

#write output file column names
f = open(outfile, "w+")
f.write('\t'.join(["agg_amp", "vic_01", "vic_02", "vic_03", "vic_04", "vic_05", "vic_06", "vic_07", "vic_08", "vic_09", "vic_10", "vic_11", "vic_12", "vic_13", "vic_14", "vic_15", "vic_16", "\n" ]))

#for each aggressor amp
for agg in range(1,17):
	if agg<10:
		aggstr = "0" + str(agg)
	else:
		aggstr = str(agg)
	
	#get list of relevant files
	files = glob.glob("/Volumes/Datadrive/CCDTesting/112-03/xtalk/system/*ccd60k*seg" + aggstr+ ".fits")
	
	#for each file
	for fname in files:
		#get aggressor column
		im_a = afwImage.ImageF(fname, agg+1)
		im_a2 =iu.unbias_and_trim(im_a)
		im_a_arr = im_a2.getArray()
		
		aggcol = np.where(im_a_arr == np.max(np.max(im_a_arr, 0)))[1][0]
		aggsignal = np.mean(im_a_arr[:,aggcol])
		
		#print "aggsignal=%f" %aggsignal
		
		#get victim signal and take the ratio of the two
		ratio = []
		for amp in range(1,17):
			
			im_v = afwImage.ImageF(fname, amp+1)
			im_v2 = iu.unbias_and_trim(im_v)
			
			im_v_arr = im_v2.getArray()
			
			vicsignal = np.mean(im_v_arr[:,aggcol])
			ratio.append(str(vicsignal/aggsignal))
			
			#print "vicsignal=%f" %vicsignal
			#print "ratio = %s" %ratio[amp-1]
		
		#write to output file
		f.write(str(agg)+'\t' + '\t'.join(ratio) + '\t' + '\n')
#close output file
f.close()