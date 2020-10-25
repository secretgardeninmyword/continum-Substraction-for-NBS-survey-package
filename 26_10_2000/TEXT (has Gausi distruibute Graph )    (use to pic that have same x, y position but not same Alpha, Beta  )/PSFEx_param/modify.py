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
	


def cut_image(img_name):
	hdul=fits.open(img_name)
	x_min=0
	x_max=25
	y_min=0
	y_max=25
	# x_min=37-18
	# x_max=37+19
	# y_min=32-18
	# y_max=32+19
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

def Read_psf_hdul(file_name):
	hdul = fits.open(file_name)
	#newdata = hdul[0].data
	new_hdr=hdul[1].header
	print(new_hdr)  	

#cut_image('proto_test.fits')
Read_psf_hdul('test.psf')
