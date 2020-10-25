'''Email: zhuangrui0452@sina.com'''
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian2DKernel, Tophat2DKernel,Gaussian1DKernel
from astropy.modeling.models import Gaussian2D
from astropy.io import fits

from mpl_toolkits import mplot3d
import math

from photutils import HanningWindow, TukeyWindow, CosineBellWindow,SplitCosineBellWindow, TopHatWindow
from photutils import create_matching_kernel
from photutils.psf.matching import resize_psf


def folder_file_name(file_dir):
	n=[]
	for root,dirs,files in os.walk(file_dir):
		n=files
	return n


def new_path(img_name,prefix):
	path_name=img_name.split('/')
	img_name_last=path_name[-1]
	path_name.pop()
	path_prefix=''
	for m in path_name:
		path_prefix+=m+'/'
	full_path=path_prefix+prefix+img_name_last	
	#print(full_path)
	return full_path

def Calibrate_x_y(folder_name,BAO_Fit_Name,MAG_INDEX=None,index_l='MAG'):
	m=folder_file_name(folder_name)
	os.system('rm -rf IPHAS_IOS_Cal')
	os.system('mkdir IPHAS_ISO_Cal')
	for z in m:
		img_name=z
		location=folder_name+'/'+img_name
		if MAG_INDEX:
			fluxscaleindex=MAG_INDEX
			ISO_Calibrate(location,fluxscaleindex)
		else:
			fluxscaleindex=Flux_Scale_Index(IPHA=location,BAO=BAO_Fit_Name,index=index_l)
			print('#########index###############')
			print(fluxscaleindex)
			if index_l=='MAG':
				ISO_Calibrate(location,fluxscaleindex)
			elif index_l=='FLUX':
				ISO_FLUX_Calibrate(location,fluxscaleindex)
		full_name=new_path(img_name,'ISO_Cal_')
		ISO_Cal_Location=folder_name+'/'+full_name
		os.system('mv %s IPHAS_ISO_Cal' %ISO_Cal_Location)
		
		
def Flux_Scale_Index(IPHA='r414157-3.fits',BAO='L20160919_04629_001901+6501_180S_S2_178471.fits',index='MAG'):
	os.system('cp %s IPHA_Sextractor_param' %BAO)
	BAO=BAO.split('/')[1]
	os.system('cd IPHA_Sextractor_param && sex %s ' %BAO)
	os.system('cd IPHA_Sextractor_param && mv %s %s' %('test.cat','test2.cat'))
	os.system('cd IPHA_Sextractor_param && rm %s' %BAO)
	
	os.system('cp %s IPHA_Sextractor_param' %IPHA)
	IPHA=IPHA.split('/')[1]
	os.system('cd IPHA_Sextractor_param && sex %s ' %IPHA)
	os.system('cd IPHA_Sextractor_param && rm %s' %IPHA)
			
	import format2csv
	difference_cat=format2csv.Get_Same_Sex_x_y(IPHAS_Cat='IPHA_Sextractor_param/test.cat',BAO_Cat='IPHA_Sextractor_param/test2.cat',result_file='IPHA_Sextractor_param/jg.csv')
	if index=='MAG':
		FluxScaleIndex=format2csv.Flux_Scale_Index_From_Cat_MAG(difference_cat)
	elif index=='FWHM':
		FluxScaleIndex=format2csv.Convolve_Index_From_Cat(difference_cat)
	elif index=='FLUX':
		FluxScaleIndex=format2csv.Flux_Scale_Index_From_Cat_FLUX(difference_cat)
	elif index=='Alpha':
		FluxScaleIndex=format2csv.Alpha_diff(difference_cat)
	elif index=='Beta':
		FluxScaleIndex=format2csv.Beta_diff(difference_cat)
	os.system('cd IPHA_Sextractor_param && rm test.cat')
	os.system('cd IPHA_Sextractor_param && rm test2.cat')
	os.system('cd IPHA_Sextractor_param && rm test.csv')
	os.system('cd IPHA_Sextractor_param && rm test2.csv')
	os.system('cd IPHA_Sextractor_param && rm jg.csv')
	return float(FluxScaleIndex)


Calibrate_x_y('after','before/before_scamp_1000_1000.fits',index_l='Alpha')
#Calibrate_x_y('after','before/before_scamp_1000_1000.fits',index_l='Beta')
