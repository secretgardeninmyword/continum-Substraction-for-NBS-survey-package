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
	
	
def conv_image(img_name,FWHM_orign,FWHM_desinate):
	# FWHM_orign=3.01
	# FWHM_desinate=20.36
	kernel_name = Gaussian2DKernel(math.sqrt((FWHM_desinate/FWHM_orign)**2-1)*FWHM_orign/(2*math.sqrt(2*math.log(2))))
	hdul=fits.open(img_name)
	try:
		hdr=hdul[0].header
		img_data=hdul[0].data
		convolved_kernel = convolve(img_data, kernel_name)
		hdul[0].data=convolved_kernel	
	except:
		hdr=hdul[1].header
		img_data=hdul[1].data
		convolved_kernel = convolve(img_data, kernel_name)
		hdul[1].data=convolved_kernel
	hdul.writeto('Conv_'+img_name,overwrite=True)

def Conv_Image_With_Kernel(img_name,kernel_name):
	# FWHM_orign=3.01
	# FWHM_desinate=20.36
	#kernel_name = Gaussian2DKernel(math.sqrt((FWHM_desinate/FWHM_orign)**2-1)*FWHM_orign/(2*math.sqrt(2*math.log(2))))
	hdul=fits.open(img_name)
	try:
		hdr=hdul[0].header
		img_data=hdul[0].data
		convolved_kernel = convolve(img_data, kernel_name)
		hdul[0].data=convolved_kernel	
	except:
		hdr=hdul[1].header
		img_data=hdul[1].data
		convolved_kernel = convolve(img_data, kernel_name)
		hdul[1].data=convolved_kernel
	full_name=new_path(img_name,'FWHM_Changed_')
	hdul.writeto(full_name,overwrite=True)	
	return full_name

def Conv_Image_With_Kernel_List(folder_name,kernel_name):
	m=folder_file_name(folder_name)
	os.system('mkdir IPHAS_Convolve')
	for z in m:
		img_name=z
		location=folder_name+'/'+img_name
		Conv_Image_With_Kernel(location,kernel_name)
		full_name=new_path(img_name,'FWHM_Changed_')
		os.system('mv %s/%s IPHAS_Convolve' %(folder_name,full_name))	
			
def just_conv(img_name,sigma=None,FWHM_orign=None,FWHM_desinate=None,one_flame='yes'):
	if FWHM_orign and FWHM_desinate:
		print('calculating kernel')
		sigma=math.sqrt((FWHM_desinate/FWHM_orign)**2-1)*FWHM_orign/(2*math.sqrt(2*math.log(2)))
		print(sigma)
		kernel_name_= Gaussian2DKernel(sigma)
	elif sigma:
		print('just using sigma')
		kernel_name_=Gaussian2DKernel(sigma)
	hdul=fits.open(img_name)
	try:
		hdr=hdul[0].header
		img_data=hdul[0].data
		convolved_kernel = convolve(img_data, kernel_name_)
		hdul[0].data=convolved_kernel	
	except:
		hdr=hdul[1].header
		img_data=hdul[1].data
		convolved_kernel = convolve(img_data, kernel_name_)
		hdul[1].data=convolved_kernel
	
	full_name=new_path(img_name,'Conv_')
	hdul.writeto(full_name,overwrite=True)
    # if one_flame=='yes':
	    # os.system('rm -rf IPHAS_PSFex_Convolve')
	    # os.system('mkdir %s' %'IPHAS_PSFex_Convolve')
	    # os.system('mv IPHA/Conv_r414564-2.fits IPHAS_PSFex_Convolve')
	return full_name

def PSFEx_Conv(IPHA_fits,BAO_fits):
	IPHA_model=PSFEx_Model(IPHA_fits)
	BAO_model=PSFEx_Model(BAO_fits)
	kernel_name=FWHM_Calibrate(IPHA_model,BAO_model,Save_Kernel='no')
	kernel_name=kernel_name/kernel_name.sum()
	
	hdul=fits.open(IPHA_fits)
	try:
		hdr=hdul[0].header
		img_data=hdul[0].data
		convolved_kernel = convolve(img_data, kernel_name)
		hdul[0].data=convolved_kernel	
	except:
		hdr=hdul[1].header
		img_data=hdul[1].data
		convolved_kernel = convolve(img_data, kernel_name)
		hdul[1].data=convolved_kernel
	full_name=new_path(IPHA_fits,'Conv_')
	hdul.writeto(full_name,overwrite=True)
	os.system('rm %s' % IPHA_model)	
	os.system('rm %s' % BAO_model)		
	return full_name
	
def PSFEx_Model(img_name):
	os.system('cp %s PSFEx_param' %img_name)
	img_name=img_name.split('/')[1]
	os.system('cd PSFEx_param && sex %s'  %img_name)
	os.system('cd PSFEx_param && psfex %s'  %'test.cat')
	os.system('cd PSFEx_param && rm  psfex.xml')
	os.system('cd PSFEx_param && rm  test.psf')
	os.system('cd PSFEx_param && rm  chi_test.fits')
	os.system('cd PSFEx_param && rm  snap_test.fits')
	os.system('cd PSFEx_param && rm  samp_test.fits')
	os.system('cd PSFEx_param && rm  resi_test.fits')
	os.system('cd PSFEx_param && rm  %s' %img_name)	
	os.system('cd PSFEx_param && rm  test.cat')	
	model=cut_image('PSFEx_param/proto_test.fits',0,25,0,25)
	proto_fits_name='proto_'+img_name
	os.system('cd PSFEx_param && mv Cutted_proto_test.fits %s' %proto_fits_name)
	os.system('cd PSFEx_param && rm  proto_test.fits')
	return 'PSFEx_param/'+proto_fits_name



def Conv_Image_With_Kernel_List_PSF_Model(folder_name,BAO_fits):
	m=folder_file_name(folder_name)
	Convolved_folder_name='IPHAS_PSFex_Convolve'
	os.system('mkdir %s' %Convolved_folder_name)
	for z in m:
		img_name=z
		location=folder_name+'/'+img_name
		full_path=PSFEx_Conv(location,BAO_fits)
		#full_name=new_path(img_name,'FWHM_Changed_')
		os.system('mv %s %s' %(full_path,Convolved_folder_name))	



def cut_image(img_name,x_min=None,x_max=None,y_min=None,y_max=None):
	hdul=fits.open(img_name)
	if x_min==None or x_max==None or y_min==None or y_max==None:
		x_min=0
		x_max=39
		y_min=0
		y_max=39
		# x_min=37-18
		# x_max=37+19
		# y_min=32-18
		# y_max=32+19
	else:
		pass
	try:
		img_data=hdul[0].data
		#L=img_data[239:273,580:622]
		L=img_data[y_min:y_max,x_min:x_max]
		hdul[0].data=L
	except:
		img_data=hdul[1].data
		#L=img_data[1149:1259,1183:1292]
		L=img_data[y_min:y_max,x_min:x_max]
		hdul[1].data=L
	full_name=new_path(img_name,'Cutted_')
	hdul.writeto(full_name,overwrite=True)

def reverse_image(img_name):
	hdul=fits.open(img_name)
	try:
		hdr=hdul[0].header
		img_data=hdul[0].data
		L=img_data*-1
		hdul[0].data=L	
	except:
		hdr=hdul[1].header
		img_data=hdul[1].data
		L=img_data*-1
		hdul[1].data=L
	full_name=new_path(img_name,'Reversed_')
	hdul.writeto(full_name,overwrite=True)
	return full_name

def Reverse_Image_List(folder_name):
	m=folder_file_name(folder_name)
	os.system('rm -rf IPHAS_Reverse')
	os.system('mkdir IPHAS_Reverse')
	for z in m:
		img_name=z
		location=folder_name+'/'+img_name
		reverse_image(location)
		full_name=new_path(img_name,'Reversed_')
		Reverse_Cal_Location=folder_name+'/'+full_name		
		os.system('mv %s IPHAS_Reverse' %Reverse_Cal_Location)
	
def ISO_Calibrate(img_name,Flux_Scale_Index):
	coefficient=1/((100**0.2)**(Flux_Scale_Index))
	#coefficient=Flux_Scale_Index
	hdul=fits.open(img_name)
	try:
		hdr=hdul[0].header
		img_data=hdul[0].data
		L=img_data*coefficient
		hdul[0].data=L	
	except:
		hdr=hdul[1].header
		img_data=hdul[1].data
		L=img_data*coefficient
		hdul[1].data=L
	full_name=new_path(img_name,'ISO_Cal_')
	hdul.writeto(full_name,overwrite=True)
	return full_name

def ISO_FLUX_Calibrate(img_name,Flux_Scale_Index):
	coefficient=Flux_Scale_Index
	hdul=fits.open(img_name)
	try:
		hdr=hdul[0].header
		img_data=hdul[0].data
		L=img_data*coefficient
		hdul[0].data=L	
	except:
		hdr=hdul[1].header
		img_data=hdul[1].data
		L=img_data*coefficient
		hdul[1].data=L
	full_name=new_path(img_name,'ISO_Cal_')
	hdul.writeto(full_name,overwrite=True)
	return full_name

def ISO_Calibrate_List(folder_name,BAO_Fit_Name,MAG_INDEX=None,index_l='MAG'):
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
	difference_cat=format2csv.Get_Same_Sex(IPHAS_Cat='IPHA_Sextractor_param/test.cat',BAO_Cat='IPHA_Sextractor_param/test2.cat',result_file='IPHA_Sextractor_param/jg.csv')
	if index=='MAG':
		FluxScaleIndex=format2csv.Flux_Scale_Index_From_Cat_MAG(difference_cat)
	elif index=='FWHM':
		FluxScaleIndex=format2csv.Convolve_Index_From_Cat(difference_cat)
	elif index=='FLUX':
		FluxScaleIndex=format2csv.Flux_Scale_Index_From_Cat_FLUX(difference_cat)
	os.system('cd IPHA_Sextractor_param && rm test.cat')
	os.system('cd IPHA_Sextractor_param && rm test2.cat')
	os.system('cd IPHA_Sextractor_param && rm test.csv')
	os.system('cd IPHA_Sextractor_param && rm test2.csv')
	os.system('cd IPHA_Sextractor_param && rm jg.csv')
	return float(FluxScaleIndex)
	

def FWHM_Calibrate_List_Sex(folder_name,BAO_Fit_Name):
	m=folder_file_name(folder_name)
	os.system('rm -rf IPHAS_FWHM_Cal')
	os.system('mkdir IPHAS_FWHM_Cal')
	for z in m:
		img_name=z
		location=folder_name+'/'+img_name

		FWHMscaleindex=Flux_Scale_Index(IPHA=location,BAO=BAO_Fit_Name,index='FWHM')
		print('#########index###############')
		print(FWHMscaleindex)
		full_name=just_conv(location,sigma=FWHMscaleindex,one_flame='no')
		print(full_name)
		os.system('mv %s IPHAS_FWHM_Cal' %full_name)


def Normalize_IPHAS_Fits_Header(file_name):
	hdul = fits.open(file_name)
	newdata = hdul[1].data
	new_hdr=hdul[1].header                   
	#hdu = fits.PrimaryHDU(data=newdata,header=new_hdr)
	hdu = fits.PrimaryHDU(data=newdata)
	hdul = fits.HDUList([hdu])
	hdul.writeto('normal_'+file_name,overwrite=True)  
	hdul.close()

def Pixel_Scale_Calibrate(img_name,input_pixel_scale, output_pixel_scale):
	hdul=fits.open(img_name)
	hdr=hdul[0].header
	img_data=hdul[0].data
	img_data=resize_psf(img_data,input_pixel_scale, output_pixel_scale)
	hdul[0].data=img_data
	full_name=new_path(img_name,'Pixel_Scale_')
	hdul.writeto(full_name,overwrite=True) 
	
	
def FWHM_Calibrate(img_1,img_2,Save_Kernel='yes'):                    # img_1 is the fits which FWHM is smaler than FWHM of img_2
	#window = TopHatWindow(0.35)
	#window=HanningWindow()
	#window=TukeyWindow(alpha=0.4)
	window=CosineBellWindow(alpha=0.235)
	#window=SplitCosineBellWindow(alpha=0.4, beta=0.3)
	try:
		hdul1=fits.open(img_1)
		img_data1=hdul1[0].data
	except:
		hdul1=fits.open(img_1)
		img_data1=hdul1[1].data
	hdul2=fits.open(img_2)
	img_data2=hdul2[0].data
	# print(img_data1.shape)
	# print(img_data2.shape)
	#input()
	kernel = create_matching_kernel(img_data1,img_data2, window=window)
	print(kernel.shape)
	#print(kernel.shape)
	final_img_data1=convolve(img_data1,kernel)
	hdul1[0].data=final_img_data1
	full_name=new_path(img_1,'Kernel_')
	if Save_Kernel=='yes':
		hdul1.writeto(full_name,overwrite=True)
	else:
		pass
	return kernel

def Swarp_Batch_Process(BAO_Fit_Name,Folder_IPHAS):
	Folder_IPHAS_List=folder_file_name(Folder_IPHAS)
	Folder_Basis=BAO_Fit_Name.split('/')[0]
	img_name_basis=BAO_Fit_Name.split('/')[-1]
	os.system('mkdir IPHAS_Swarp')

	for r in Folder_IPHAS_List:
		img_name_IPHAS=r
		IPHAS_Location=Folder_IPHAS+'/'+img_name_IPHAS
		#print(img_name_basis[:-5])
		#print(img_name_IPHAS[-14:])
		os.system('cd '+Folder_Basis+' && swarp '+img_name_basis+' ../'+IPHAS_Location)
		os.system('mv '+Folder_Basis+'/coadd.fits '+'IPHAS_Swarp/Swarped_'+img_name_basis[:-5]+img_name_IPHAS[-14:])
		os.system('cd '+Folder_Basis+' && rm coadd.weight.fits')
		os.system('cd '+Folder_Basis+' && rm swarp.xml')

def Batch_Process(BAO_Source_Folder_Name,IPHA_Source_Folder_Name):
	BAO_Folder_List=folder_file_name(BAO_Source_Folder_Name)
	for Bao_Fit in BAO_Folder_List:
		Bao_Fit_Full_Path=BAO_Source_Folder_Name+'/'+Bao_Fit
		ISO_Calibrate_List(IPHA_Source_Folder_Name,Bao_Fit_Full_Path)
		Reverse_Image_List('IPHAS_ISO_Cal')
		Conv_Image_With_Kernel_List('IPHAS_Reverse',kernel_name)
		Swarp_Batch_Process(Bao_Fit_Full_Path,'IPHAS_Convolve')
		os.system('rm -rf IPHAS_Convolve')
		os.system('rm -rf IPHAS_ISO_Cal')
		os.system('rm -rf IPHAS_Reverse')		


'''psfex mode Conv'''



ISO_Calibrate_List(folder_name='IPHA',BAO_Fit_Name='BAO/l2016.fits',index_l='FLUX')
Conv_Image_With_Kernel_List_PSF_Model(folder_name='IPHAS_ISO_Cal',BAO_fits='BAO/l2016.fits')
Reverse_Image_List(folder_name='IPHAS_PSFex_Convolve')
Swarp_Batch_Process(BAO_Fit_Name='BAO/l2016.fits',Folder_IPHAS='IPHAS_Reverse')
