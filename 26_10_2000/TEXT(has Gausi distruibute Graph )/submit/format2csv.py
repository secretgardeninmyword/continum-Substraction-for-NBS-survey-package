'''
You need to put 'default.conv', 'default.param', 'default.sex' in folder 
whose name is IPHA_Sextractor_param
####The request for default.param of sextractor
1 Number
3 FLUX_APER
4 FLUXERR_APER
5 MAG_APER
10 ALPHA_J2000
11 DELTA_J2000
19 FWHM_IMAGE
####The request for default.sex of sextractor
CATALOG_TYPE     ASCII 
'''

import os
import csv
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from scipy.stats import norm
from scipy.optimize import curve_fit 

def alter(file_,old_str,new_str):
    file_data = ""
    with open(file_, "r") as f:
        for line in f:
            if old_str in line:
                print(line)
                line = line.replace(old_str,new_str)
				#print(line)
            file_data += line
    with open(file_,"w") as f:
        f.write(file_data)
        f.close()



def Get_Same_Sex(IPHAS_Cat='test.cat',BAO_Cat='test2.cat',result_file='jg.csv'):
	IPHAS_Csv=IPHAS_Cat.split('.')[0]+'.csv'
	BAO_Csv=BAO_Cat.split('.')[0]+'.csv'
	os.system("sed 's/^[ ]*//g' %s | sed 's/ \+/,/g' > %s" %(IPHAS_Cat,IPHAS_Csv))
	os.system("sed 's/^[ ]*//g' %s | sed 's/ \+/,/g' > %s" %(BAO_Cat,BAO_Csv))
	f_result=open(result_file,'w')
	f1= open(IPHAS_Csv,'r') 
	reader1=csv.reader(f1)
	
	f_result.writelines('NUMBER ,'+'FLUX_APER ,'+'FLUXERR_APER ,'+'MAG_APER ,'
		+'ALPHA_J2000 ,'+'DELTA_J2000 ,'+'FWHM_IMAGE ,'+'MAG_ISO ,'+'NAME   \n')
	for i in reader1:
		ALPHA_J2000_1=float(i[9])
		DELTA_J2000_1=float(i[10])
		f2=open(BAO_Csv,'r')
		reader2=csv.reader(f2)
		for j in reader2:
			ALPHA_J2000_2=float(j[9])
			DELTA_J2000_2=float(j[10])	
			if abs(ALPHA_J2000_1-ALPHA_J2000_2)<=15/3600*0.20 and abs(DELTA_J2000_1-DELTA_J2000_2)<=1/3600*0.20:
				f_result.writelines(i[0]+' ,'+i[2]+' ,'+i[3]+' ,'+i[4]+' ,'+i[9]+' ,'+i[10]+' ,'+i[18]+' ,'+i[20]+' ,IPHAS'+'\n')	
				f_result.writelines(j[0]+' ,'+j[2]+' ,'+j[3]+' ,'+j[4]+' ,'+j[9]+' ,'+j[10]+' ,'+j[18]+' ,'+j[20]+' ,BAO'+'\n')
				#f_result.writelines('\n')
	return result_file

def Flux_Scale_Index_From_Cat_MAG(difference_cat):	
	list_BAO=[]
	list_IPHA=[]
	difference_mag=[]
	result_file=difference_cat
	f_Same_Sex=open(result_file,'r')
	reader_flux_scale=csv.reader(f_Same_Sex)
	index_n=1
	for i in reader_flux_scale:	
		try:
			if index_n%2==0:
				list_IPHA.append(float(i[3]))
			elif index_n%2==1:
				list_BAO.append(float(i[3]))
		except:
			pass
		index_n+=1
	for j in range(len(list_BAO)):
		difference_mag.append(list_BAO[j]-list_IPHA[j])
	difference_mag=np.array(difference_mag)
	print(np.median(difference_mag))
	return np.median(difference_mag)

def Flux_Scale_Index_From_Cat_FLUX(difference_cat):	
	list_BAO=[]
	list_IPHA=[]
	difference_mag=[]
	result_file=difference_cat
	f_Same_Sex=open(result_file,'r')
	reader_flux_scale=csv.reader(f_Same_Sex)
	index_n=1
	for i in reader_flux_scale:	
		try:
			if index_n%2==0:
				list_IPHA.append(float(i[1]))
			elif index_n%2==1:
				list_BAO.append(float(i[1]))
		except:
			pass
		index_n+=1
	for j in range(len(list_BAO)):
		difference_mag.append(list_BAO[j]/list_IPHA[j])
	difference_mag=np.array(difference_mag)
	median_value=plot_Histogram(difference_mag)
	return median_value
	'''
	return np.median(difference_mag)
    '''
def plot_Histogram_pre(array_data):
	DIF_value=15
	Histo_num=[]
	Histo_value=[]
	median_value=np.median(array_data)
	min_DIF=(np.median(array_data)-array_data.min())/DIF_value
	max_DIF=(array_data.max()-np.median(array_data))/DIF_value
	#DIF=(array_data.max()-array_data.min())/5
	if min_DIF>max_DIF:
		DIF=max_DIF
	else:
		DIF=min_DIF
	#print(array_data)
	#print('#############')
	print('mean value :',median_value)
	print('DIF: ',DIF)
	print('std: ',np.std(array_data))
	for n in range(DIF_value*-1,DIF_value):
		#print(n+9,n+1+9)
		w=array_data[(n*DIF+median_value<array_data) & (array_data<(n+1)*DIF+median_value)]
		Histo_num.append(len(w))
		Histo_value.append((n+0.5)*DIF+median_value)
		print('#',n*DIF+median_value,'#')
		# print('##')
		# print((n+0.5)*DIF+median_value)
		# print(n*DIF+median_value,(n+1)*DIF+median_value)
		# print('###')
		#print(len(w))
	#plt.bar(Histo_value,Histo_num)
	plt.scatter(Histo_value,Histo_num)
	plt.plot(Histo_value,Histo_num)


	# plt.scatter(Histo_value, g(Histo_value), label='Gaussian')

	
	plt.show()

def func2(x, a, b, c,d):                    
    return a*np.exp(-(x-b)**2/(2*c**2))+d

def plot_Histogram(array_data):	
	upper_limit=array_data.max()
	lower_limit=array_data.min()
	std_deviation=3
	mean=np.mean(array_data)
	sigma=np.std(array_data)
	median=np.median(array_data)
	print('mean(unfilter): ',mean)
	print('median(unfilter): ',median)
	available_data=array_data[(mean-std_deviation*sigma<array_data) & (array_data<mean+std_deviation*sigma)]
	mean=np.mean(available_data)
	sigma=np.std(available_data)
	median=np.median(available_data)
	print('mean: ',mean)
	print('median: ',median)


	n,bins,patchs=plt.hist(available_data,50)
	#plt.cla()
	width=bins[1]-bins[0]
	mid_bins=bins[1:]-width/2
	plt.scatter(mid_bins,n,color='r')
	
	''' just use paramater mean sigma, and not to fit data
	# y=norm.pdf(bins,mean,sigma)
	# print(y.mean())
	# plt.plot(bins,y,'--')
	'''
	x = np.linspace(lower_limit, upper_limit, 1000)
	P0=[36,mean,sigma]
	popt, _ = curve_fit(func2, mid_bins, n)
	a1=popt[0]
	b1=popt[1]
	c1=popt[2]
	d1=popt[3]
	print('a1: ',a1)
	print('b1: ',b1)
	print('c1: ',c1)
	print('d1: ',d1)
	plt.plot(x,func2(x,a1,b1,c1,d1))
	print('mean value:',b1)
	plt.pause(3)
	return b1






	
def Convolve_Index_From_Cat(difference_cat):	
	list_BAO=[]
	list_IPHA=[]
	difference_FWHM=[]
	result_file=difference_cat
	f_Same_Sex=open(result_file,'r')
	reader_flux_scale=csv.reader(f_Same_Sex)
	index_n=1
	for i in reader_flux_scale:	
		try:
			if index_n%2==0:
				list_IPHA.append(float(i[6]))
			elif index_n%2==1:
				list_BAO.append(float(i[6]))
		except:
			pass
		index_n+=1
	for j in range(len(list_BAO)):
		FWHM_desinate=list_BAO[j]
		FWHM_orign=list_IPHA[j]
		#print('#############')
		#print(FWHM_desinate)
		#print(FWHM_orign)
		try:
			sigma=math.sqrt((FWHM_desinate/FWHM_orign)**2-1)*FWHM_orign/(2*math.sqrt(2*math.log(2)))
			difference_FWHM.append(sigma)
			#print(difference_FWHM)
		except:
			pass
	difference_FWHM=np.array(difference_FWHM)
	mean_value=plot_Histogram(difference_FWHM)
	'''
	return np.median(difference_FWHM)
	'''
	return mean_value
