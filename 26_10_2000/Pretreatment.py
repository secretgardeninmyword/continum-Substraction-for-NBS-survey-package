import os
import datetime

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import numpy as np


def folder_file_name(file_dir):
	n=[]
	for root,dirs,files in os.walk(file_dir):
		n=files
	return n

	
def get_list_index(list_name,n):
	total_data=[]
	for m in list_name:
		unit_data=m[n:n+1,]
		total_data.append(unit_data)
	return total_data	

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
	
def get_above_overscan(i,img_data):
	index_main=i*1340
	black_one=img_data[0:5280,(i+1)*1340-20:(i+1)*1340]
	black_two=img_data[5280:5300,i*1340:(i+1)*1340]
	tran=np.transpose(black_two)
	fi=np.vstack((tran,black_one))
	middle_black=np.median(fi)
	return middle_black


def get_belove_overscan(i,img_data):
	index_main=i*1340
	black_one=img_data[5320:10600,(i+1)*1340-20:(i+1)*1340]
	black_two=img_data[5300:5320,i*1340:(i+1)*1340]
	tran=np.transpose(black_two)
	fi=np.vstack((tran,black_one))
	middle_black=np.median(fi)
	return middle_black


def get_above_pic(i,index_main,img_data):
	middle_black=get_above_overscan(i,img_data)
	main_above=img_data[0:5280,index_main:index_main+1320]
	final=main_above-middle_black
	return final


def get_belove_pic(i,index_main,img_data):
	middle_black=get_belove_overscan(i,img_data)
	main_belove=img_data[5320:10600,index_main:index_main+1320]
	final_belove=main_belove-middle_black
	return final_belove

	
def combine_above_pic(n,list_one):
	n_above=list_one[0]	
	for m in range(n):
		n_above=np.hstack((n_above,list_one[m+1]))
	return n_above


def combine_belove_pic(n,list_two):
	n_belove=list_two[0]
	for m in range(n):
		n_belove=np.hstack((n_belove,list_two[m+1]))
	return n_belove

	
def main_get_rid_of_overscan(file_name,write_to_name):     #get rid of overscan and make the size of picture to 10560*10560 '''
	hdul=fits.open(file_name)
	img_data=hdul[0].data
	list_one=[]
	list_two=[]
	'''above pic'''
	for i in range(8):
		index_main=i*1340
		final=get_above_pic(i,index_main,img_data)
		list_one.append(final)
	'''belove pic'''
	for s in range(8):
		index_main=s*1340
		final_belove=get_belove_pic(s,index_main,img_data)
		list_two.append(final_belove)
	n_final=np.vstack((combine_above_pic(7,list_one),combine_belove_pic(7,list_two)))
	#plt.imshow(n_final)
	plt.show()
	hdul[0].data=n_final
	hdul.writeto(write_to_name)


def main_get_rid_of_list_pic_of_overscan(get_folder_name,write_to_folder_name):
	m=folder_file_name(get_folder_name)
	list_1=[]
	for z in m:
		name=z
		location=get_folder_name+'/'+z
		write_to_location=write_to_folder_name+'/'+z
		main_get_rid_of_overscan(location,write_to_location)


	
def get_pic_data(file_name):                           					#get the picture info data and return the matrix'''
	hdul=fits.open(file_name)
	img_data=hdul[0].data
	return img_data


def get_average_back_noise(file_one,file_two,file_three):
	average_data=(get_pic_data(file_one)+get_pic_data(file_two)+get_pic_data(file_three))/3
	return average_data
	
	
def get_rid_of_ground(pic_name,bias_one,bias_two,bias_three):								#combine the bias(background noise) and picture'''
	hdul=fits.open(pic_name)
	img_data=hdul[0].data
	back_ground=get_average_back_noise(bias_one,bias_two,bias_three)
	final_pic=img_data-back_ground
	hdul[0].data=final_pic
	hdul.writeto('Get_Rid_Of_bias_'+pic_name)

def get_rid_of_ground_list(folder_name,write_to_folder,bias_one,bias_two,bias_three):
	m=folder_file_name(folder_name)
	for z in m:
		name=z
		location=folder_name+'/'+z
		write_to_location=write_to_folder+'/'+z
		get_rid_of_ground_new(location,write_to_location,bias_one,bias_two,bias_three)
	
	
def get_rid_of_ground_new(source_folder,write_to_folder,bias_one,bias_two,bias_three):								#combine the bias(background noise) and picture'''
	hdul=fits.open(source_folder)
	img_data=hdul[0].data
	back_ground=get_average_back_noise(bias_one,bias_two,bias_three)
	final_pic=img_data-back_ground
	hdul[0].data=final_pic
	hdul.writeto(write_to_folder)


def combine_five_pic_to_one(source_pic_one,source_pic_two,source_pic_three,source_pic_four,source_pic_five,named_pic):    # combine five pic to one (get the mid-value)
	hdul=fits.open(source_pic_one)
	img_data=hdul[0].data
	pic_1=get_pic_data(source_pic_one)
	pic_2=get_pic_data(source_pic_two)
	pic_3=get_pic_data(source_pic_three)
	pic_4=get_pic_data(source_pic_four)
	pic_5=get_pic_data(source_pic_five)
	m=np.arange(10560)
	o=np.arange(10560)
	p=np.arange(10560)
	q=np.arange(10560)
	for n in range(2640):
		one=pic_1[n:n+1,]
		two=pic_2[n:n+1,]
		three=pic_3[n:n+1,]
		four=pic_4[n:n+1,]
		five=pic_5[n:n+1,]
		final=np.vstack((one,two,three,four,five))
		average_line=np.median(final,axis=0)
		m=np.vstack((m,average_line))
	m=m[1:2641,]
	for n in range(2640,5280):
		one=pic_1[n:n+1,]
		two=pic_2[n:n+1,]
		three=pic_3[n:n+1,]
		four=pic_4[n:n+1,]
		five=pic_5[n:n+1,]
		final=np.vstack((one,two,three,four,five))
		average_line=np.median(final,axis=0)
		o=np.vstack((o,average_line))
	o=o[1:2641,]
	for n in range(5280,7920):
		one=pic_1[n:n+1,]
		two=pic_2[n:n+1,]
		three=pic_3[n:n+1,]
		four=pic_4[n:n+1,]
		five=pic_5[n:n+1,]
		final=np.vstack((one,two,three,four,five))
		average_line=np.median(final,axis=0)
		p=np.vstack((p,average_line))
	p=p[1:2641,]
	for n in range(7920,10560):
		one=pic_1[n:n+1,]
		two=pic_2[n:n+1,]
		three=pic_3[n:n+1,]
		four=pic_4[n:n+1,]
		five=pic_5[n:n+1,]
		final=np.vstack((one,two,three,four,five))
		average_line=np.median(final,axis=0)
		q=np.vstack((q,average_line))
	q=q[1:2641,]
	img_data_final=np.vstack((m,o,p,q))
	hdul[0].data=img_data_final
	hdul.writeto(named_pic)


def combine_all_list_pic_to_one(folder_name,opened_file_name):			#like the above fuction combine_five_pic_to_one
	m=folder_file_name(folder_name)
	list_1=[]
	for z in m:
		name="'"+z+"'"
		location=folder_name+'/'+z
		data=get_pic_data(location)
		list_1.append(data)
	
	hdul=fits.open(opened_file_name)
	img_data=hdul[0].data
	m=np.arange(10560)
	o=np.arange(10560)
	p=np.arange(10560)
	q=np.arange(10560)
	for n in range(2640):
		list_2=get_list_index(list_1,n)
		final=np.vstack(list_2)
		average_line=np.median(final,axis=0)
		m=np.vstack((m,average_line))
		print(n)
	m=m[1:2641,]
	for n in range(2640,5280):
		list_2=get_list_index(list_1,n)
		final=np.vstack(list_2)
		average_line=np.median(final,axis=0)
		o=np.vstack((o,average_line))
		print(n)
	o=o[1:2641,]
	for n in range(5280,7920):
		list_2=get_list_index(list_1,n)
		final=np.vstack(list_2)
		average_line=np.median(final,axis=0)
		p=np.vstack((p,average_line))
		print(n)
	p=p[1:2641,]
	for n in range(7920,10560):
		list_2=get_list_index(list_1,n)
		final=np.vstack(list_2)
		average_line=np.median(final,axis=0)
		q=np.vstack((q,average_line))
		print(n)
	q=q[1:2641,]
	img_data_final=np.vstack((m,o,p,q))
	hdul[0].data=img_data_final
	hdul.writeto(folder_name+'.FITS')


def normalization_combined_pic(source_pic,normalization_type):                               # normalize the pic (in the belove case, we normalize the mid-value pic) 
	hdul=fits.open(source_pic)
	img_data=hdul[0].data
	if normalization_type=='max':
		max_value=img_data.max()
		max_value_matrix=np.ones((10560,10560))*max_value
		normaliztion_matrix=img_data/max_value_matrix
		hdul[0].data=normaliztion_matrix
		hdul.writeto('normalization_max_value_'+source_pic)

	if normalization_type=='mean':
		mean_value=img_data.mean()
		mean_value_matrix=np.ones((10560,10560))*mean_value
		normaliztion_matrix=img_data/mean_value_matrix
		hdul[0].data=normaliztion_matrix
		hdul.writeto('normalization_mean_value_'+source_pic)

	if normalization_type=='min':
		min_value=img_data.min()
		min_value_matrix=np.ones((10560,10560))*min_value
		normaliztion_matrix=img_data/min_value_matrix
		hdul[0].data=normaliztion_matrix
		hdul.writeto('normalization_min_value_'+source_pic)

	if normalization_type=='mid':
		mid_value=np.median((img_data))
		mid_value_matrix=np.ones((10560,10560))*mid_value
		normaliztion_matrix=img_data/mid_value_matrix
		hdul[0].data=normaliztion_matrix
		hdul.writeto('normalization_mid_value_'+source_pic)


def get_rid_of_negative_plate(source_pic,negative_plate): 	
	negative_plate=get_pic_data(negative_plate)	
	hdul=fits.open(source_pic)
	img_data=hdul[0].data
	get_rid_of_negative_plate=img_data/negative_plate
	hdul[0].data=get_rid_of_negative_plate
	hdul.writeto('Get_Rid_Of_Negative_Plate'+source_pic)

def pre_get_rid_of_negative_plate(source_pic,negative_plate,write_to_name,write_to_folder): 	
	negative_plate=get_pic_data(negative_plate)	
	hdul=fits.open(source_pic)
	img_data=hdul[0].data
	get_rid_of_negative_plate=img_data/negative_plate
	hdul[0].data=get_rid_of_negative_plate
	hdul.writeto(write_to_folder+'/'+'Get_Rid_Of_Negative_Plate'+write_to_name)

def get_rid_of_negative_plate_list(folder_name,negative_plate,write_to_folder):
	m=folder_file_name(folder_name)
	for z in m:
		name=z
		location=folder_name+'/'+z
		pre_get_rid_of_negative_plate(location,negative_plate,name,write_to_folder)

def added_astrometry_wcs_header(wcs_file,img_name):
	hdul_wcs=fits.open(wcs_file)
	hdul_img=fits.open(img_name)
	hdul_img[0].header=hdul_wcs[0].header
	hdul_img.writeto(img_name,overwrite=True)
	#print(hdul[0].header)

# def Change_HDUL_fmt(img_name,model):
	# hdul_img=fits.open(img_name)
	# model_hdul=fits.open(model)
	# #model[0].header=hdul_img[1].header
	# model[0].data=hdul_img[1].data
	
	# #hdul_img[0].header=hdul_img[1].header
	# new_hdul.writeto('Changed'+img_name)
		
def cut_image(img_name):
	hdul=fits.open(img_name)
	try:
		img_data=hdul[0].data
		#L=img_data[1760:5280,1760:5280]
		#L=img_data[2860:4180,2860:4180]
		#L=img_data[2360:4680,2360:4680]
		#L=img_data[2700:3520,2700:3520]
		L=img_data[897:972,572:604]
		hdul[0].data=L
	except:
		img_data=hdul[1].data
		L=img_data[1149:1259,1183:1292]
		hdul[1].data=L
	hdul.writeto('Cutted_'+img_name,overwrite=True)

def cut_image_with_wcs(img_name,position,size):
	from astropy.nddata import Cutout2D
	from astropy.wcs import WCS
	
	# hdul=fits.open(file_name)
	# hdr=hdul[0].header
	hdu=fits.open(img_name)[0]
	wcs=WCS(hdu.header)
	
	cutout=Cutout2D(hdu.data,position=position,size=size,wcs=wcs)
	hdu.data=cutout.data
	hdu.header.update(cutout.wcs.to_header())
	full_name=new_path(img_name,'Cutted_')
	hdu.writeto(full_name,overwrite=True)

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
	hdul.writeto('Reversed_'+img_name)
	
def manipulate_image(img_name,number):
	hdul=fits.open(img_name)
	try:
		hdr=hdul[0].header
		img_data=hdul[0].data
		L=img_data*number
		hdul[0].data=L	
	except:
		hdr=hdul[1].header
		img_data=hdul[1].data
		L=img_data*number
		hdul[1].data=L
	hdul.writeto(img_name,overwrite=True)
	
def calibrate_image(img_name):
	hdul=fits.open(img_name)
	img_data=hdul[0].data
	for x in range(len(img_data)):
		for y in range(len(img_data)):
			if img_data[x,y]-0<0.1 or img_data[x,y]-0>-0.1 :
				img_data[x,y]=0

	hdul[0].data=img_data
	hdul.writeto('Calibrated'+img_name)

def pre_cut_image_list(source_pic,img_name):
	hdul=fits.open(source_pic)
	img_data=hdul[0].data
	L=img_data[1760:8800,1760:8800]
	hdul[0].data=L
	hdul.writeto(img_name)

def cut_image_list(folder_name,write_to_folder):
	m=folder_file_name(folder_name)
	for z in m:
		name=z
		location=folder_name+'/'+z
		write_to_location=write_to_folder+'/'+z
		pre_cut_image_list(location,write_to_location)

def put_the_pic_name_to_sh_script(script_name,folder_name):  #folder name is the folder which contain the pic
	with open(script_name,'w') as f:
		m=folder_file_name(folder_name)
		for file_name in m:
			prefix='./solve-field ../pic_un_treated/'
			volume=prefix+file_name+'\n'
			f.writelines(volume)
	
def HDUL_Edit(file_name):
	curr_time=datetime.datetime.now()
	hdul=fits.open(file_name)
	hdr=hdul[0].header
	hdr['COMMENT']='############LAST CHANGE TIME############'
	hdr['COMMENT']=str(curr_time.date())
	hdr['OBSERVER']=('jg','wasd')
	hdul[0].header=hdr
	hdul.writeto('HDUL_CHANGED_'+file_name)

def Read_HDUList_head(file_name):
	hdul=fits.open('file_name')
	hdr=hdul[0].header
	print(repr(hdr))


def HDUL_Add_Header(file_name):
	'''require pixel frame is 5280*5280'''
	#this fuction is useful for scamp, since scamp requires fits that contain below info in their hdr, this info such as CTYPE1, CTYPE2, CUNIT1
	hdul=fits.open(file_name)
	hdr=hdul[0].header
	hdr['CTYPE1']='RA---TAN'
	hdr['CTYPE2']='DEC--TAN'
	hdr['CUNIT1']='deg'
	hdr['CUNIT2']='deg'
	hdr['CRPIX1']=hdr['NAXIS1']/2
	hdr['CRPIX2']=hdr['NAXIS2']/2
	hdr['CRVAL1']=float(hdr.comments['RA'])
	hdr['CRVAL2']=float(hdr.comments['DEC'])
	hdr['CD1_1']=2.854528829e-04
	hdr['CD1_2']=1.123165721e-06
	hdr['CD2_1']=-1.163116820e-06
	hdr['CD2_2']=2.854487758e-04
	hdul[0].header=hdr
	hdul.writeto('HDUL_Changed_'+file_name)
	
def HDUL_Add_Header_1000(file_name):
	'''require pixel frame is 1000*1000'''
	#this fuction is useful for scamp, since scamp requires fits that contain below info in their hdr, this info such as CTYPE1, CTYPE2, CUNIT1
	hdul=fits.open(file_name)
	hdr=hdul[0].header
	hdr['CTYPE1']='RA---TAN'
	hdr['CTYPE2']='DEC--TAN'
	hdr['CUNIT1']='deg'
	hdr['CUNIT2']='deg'
	hdr['CRPIX1']=hdr['NAXIS1']/2
	hdr['CRPIX2']=hdr['NAXIS2']/2
	hdr['CRVAL1']=float(hdr.comments['RA'])
	hdr['CRVAL2']=float(hdr.comments['DEC'])
	hdr['CD1_1']=0.000285432368127
	hdr['CD1_2']=-8.33631732182E-07
	hdr['CD2_1']=1.14874111805E-06
	hdr['CD2_2']=0.00028598603957
	hdul[0].header=hdr
	hdul.writeto('HDUL_Changed_'+file_name)
	
def write_ahead_file(file_name):
	hdul=fits.open(file_name)
	hdr=hdul[0].header
	print(hdr)
	#get the file name, which doesn't contain suffix
	i=''
	for m in file_name:
		if m=='.':
			break
		i=i+m
	#write the hdr content to the file which the suffix of it is .ahead
	print('write to '+i+'.ahead')
	w = open(i+'.ahead', 'w')
	for m in hdr:
		sad=str(m)+'='+str(hdr[m])+'\n'
		w.writelines(sad)	
		
def Add_WCS_Header(file_name):
	#put the hdr info which scamp produces(named test.head) into fits file
	hdul=fits.open(file_name)
	hdr=hdul[0].header
	del hdr['CTYPE1']
	del hdr['CTYPE2']
	del hdr['CUNIT1']
	del hdr['CUNIT2']
	del hdr['CRPIX1']
	del hdr['CRPIX2']
	del hdr['CRVAL1']
	del hdr['CRVAL2']
	del hdr['CD1_1']
	del hdr['CD1_2']
	del hdr['CD2_1']
	del hdr['CD2_2']
	m=fits.Header.fromtextfile('test.head')
	hdr.extend(m)
	hdul[0].header=hdr
	del hdr['FLXSCALE']
	hdul.writeto('WCS_Added_'+file_name,overwrite=True)	
''' 5 pic condition'''
# main_get_rid_of_overscan('/home/rui/Astro/Astro_Pic/B20160919_00000_183702-1715_0S_CL_178354.FITS','B20160919_00000_183702-1715_0S_CL_178354.FITS')	
# get_rid_of_ground('L20160919_04625_185235+0458_180S_S2_178365.FITS','B20160919_00000_183629-1715_0S_CL_178352.FITS','B20160919_00000_183647-1715_0S_CL_178353.FITS','B20160919_00000_183702-1715_0S_CL_178354.FITS')	
# combine_five_pic_to_one('Get_Rid_Of_bias_L20160919_04625_185212+0502_180S_S2_178361.FITS','Get_Rid_Of_bias_L20160919_04625_185217+0454_180S_S2_178362.FITS',
	# 'Get_Rid_Of_bias_L20160919_04625_185218+0459_180S_S2_178364.FITS','Get_Rid_Of_bias_L20160919_04625_185230+0502_180S_S2_178363.FITS',
	# 'Get_Rid_Of_bias_L20160919_04625_185235+0458_180S_S2_178365.FITS','mid_value_matrix.fits')
# normalization_combined_pic('mid_value_matrix.fits','mid')
# get_rid_of_negative_plate('Get_Rid_Of_bias_L20160919_04625_185230+0502_180S_S2_178363.FITS','normalization_mid_value_mid_value_matrix.fits')
#cut_image('Get_Rid_Of_Negative_PlateL20160919_04626_191832+0858_180S_S2_178382.FITS')
'''list pic condition'''
#main_get_rid_of_list_pic_of_overscan('overscan_data','unoverscan_data')
#main_get_rid_of_list_pic_of_overscan('123','234')
# get_rid_of_ground_list('unoverscan_data','get_rid_of_overground','B20160919_00000_183629-1715_0S_CL_178352.FITS',
	# 'B20160919_00000_183647-1715_0S_CL_178353.FITS','B20160919_00000_183702-1715_0S_CL_178354.FITS')

#combine_all_list_pic_to_one('12','1.FITS')
#normalization_combined_pic('get_rid_of_overground.FITS','mid')
#get_rid_of_negative_plate_list('get_rid_of_overground','normalization_mid_value_get_rid_of_overground.FITS','get_rid_of_negative_ground')
#cut_image_list('get_rid_of_negative_ground','cuted frame')
#put_the_pic_name_to_sh_script('program.sh','cuted frame')
#HDUL_Edit('Get_Rid_Of_Negative_PlateL20160919_04625_185230+0502_180S_S2_178363.FITS')





'''      scamp     Change Header          '''
#write_ahead_file('Get_Rid_Of_Negative_PlateL20160919_04626_191658+1101_180S_S2_178390.FITS')    #write fits headr info to .ahead file 



# HDUL_Add_Header('jg.fits')                   #some hdr info must be writed in to fits file, such as CDi_j
#Add_WCS_Header('scamp_178471.fits')
'''astrometry add wcs'''
#added_astrometry_wcs_header('wcs.fits','Cutted_WCS_edited_CuttedAfter_calibrated_Get_Rid_Of_Negative_PlateL20160919_04629_001901+6501_180S_S2_178471.fits')

'''inverse, cut or calibrate fits'''
#cut_image('r414564-2.fits')
#reverse_image('r414564-2.fits')
#manipulate_image('Reversed_r414564-2.fits',1/(100**0.2)**(16.3108-13.6765))
#manipulate_image('Reversed_r414564-2.fits',1.5)

#calibrate_image('coadd.fits')




#cut_image_with_wcs(img_name='WCS_Added_scamp_178471.fits',position=(3380,3380),size=(1000,1000))
