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
		+'ALPHA_J2000 ,'+'DELTA_J2000 ,'+'FWHM_IMAGE  \n')
	for i in reader1:
		ALPHA_J2000_1=float(i[9])
		DELTA_J2000_1=float(i[10])
		f2=open(BAO_Csv,'r')
		reader2=csv.reader(f2)
		for j in reader2:
			ALPHA_J2000_2=float(j[9])
			DELTA_J2000_2=float(j[10])	
			if abs(ALPHA_J2000_1-ALPHA_J2000_2)<=15/3600*0.20 and abs(DELTA_J2000_1-DELTA_J2000_2)<=1/3600*0.20:
				f_result.writelines(i[0]+' ,'+i[2]+' ,'+i[3]+' ,'+i[4]+' ,'+i[9]+' ,'+i[10]+' ,'+i[18]+'\n')	
				f_result.writelines(j[0]+' ,'+j[2]+' ,'+j[3]+' ,'+j[4]+' ,'+j[9]+' ,'+j[10]+' ,'+i[18]+'\n')
				#f_result.writelines('\n')
	return result_file

def Flux_Scale_Index_From_Cat(difference_cat):	
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
	# print('#############stop#################')
	# input()
	return np.median(difference_mag)

# difference_cat=Get_Same_Sex(IPHAS_Cat='IPHAS.cat',BAO_Cat='BAO.cat',result_file='jg.csv')
#Flux_Scale_Index_From_Cat('jg.csv')
