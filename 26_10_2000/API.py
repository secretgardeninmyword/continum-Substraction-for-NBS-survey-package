from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy import wcs
import numpy as np
import requests
import json
import random
from urllib.request import urlopen, Request
from urllib.request import urlopen, Request
from urllib.parse import urlencode, quote
import time
import os
import wget

def folder_file_name(file_dir):
	n=[]
	for root,dirs,files in os.walk(file_dir):
		n=files
	return n

def json2python(data):
    try:
        return json.loads(data)
    except:
        pass
    return None
python2json = json.dumps

def Get_file_args(fn):
	if fn is not None:
		try:
			f = open(fn, 'rb')
			file_args = (fn, f.read())
			return file_args
		except IOError:
			print('File %s does not exist' % fn)
			raise

def up_load_file(fn,session,json):
	file_args=Get_file_args(fn)
	url='http://nova.astrometry.net/api/upload'
	# asd={"allow_commercial_use": "d", "allow_modifications": "d", 
		# "publicly_visible": "y", "session": session,"apikey": "wsjkkufuuwxvjovl"}
	# json1=json.dumps(asd)
	boundary_key = ''.join([random.choice('0123456789') for i in range(19)])
	boundary = '===============%s==' % boundary_key
	headers = {'Content-Type':
			   'multipart/form-data; boundary="%s"' % boundary}
	data_pre = (
		'--' + boundary + '\n' +
		'Content-Type: text/plain'+'\r\n' +
	
		'MIME-Version: 1.0'+'\r\n' +
		'Content-disposition: form-data; name="request-json"'+'\r\n' +
		'\r\n' +
		json + '\n' +
		'--' + boundary + '\n' +
		'Content-Type: application/octet-stream'+'\r\n' +
		'MIME-Version: 1.0'+'\r\n' +
		'Content-disposition: form-data; name="file"; filename="%s"' % file_args[0] +
	    '\r\n' + '\r\n') 
	data_post = (
		'\n' + '--' + boundary + '--\n')
	data = data_pre.encode()+file_args[1]  + data_post.encode()
	print('#######start to request#######')
	request = Request(url=url, headers=headers, data=data)
	f = urlopen(request)
	txt = f.read()
	str1=str(txt,encoding='utf-8')
	result=eval(str1)
	#print('txt type',type(txt))
	#print('Got json:', txt)
	print('Got result',result)
	stat = result['status']
	print('Got status:', stat)
	return result
	
def upload_info(result,json):	
	print('#########start to uploading#########')
	sub_id=result['subid']
	url_return= 'http://nova.astrometry.net/api/'+'submissions/%s' % sub_id	
	data2 = {'request-json': json}
	#print('Sending form data:', data2)
	data2 = urlencode(data2)
	data2 = data2.encode('utf-8')
	print('Sending data:', data2)
	headers = {}
	request = Request(url=url_return, headers=headers, data=data2)
	f = urlopen(request)
	txt2 = f.read()
	#print('Got json:', txt2)
	str2=str(txt2,encoding='utf-8')
	result2=eval(str2)
	print('Got result:', result2)
	return result2,sub_id,url_return

'''check if the job has done and the return id'''
def sub_status(url_return):    
	print('######getting status#######')
	data2 = {'request-json': json}
	print('Sending form data:', data2)
	data2 = urlencode(data2)
	data2 = data2.encode('utf-8')
	#print('Sending data:', data2)
	headers = {}	
	request = Request(url=url_return, headers=headers, data=data2)
	json_duplicate={'jobs':[]}
	f = urlopen(request)
	txt2 = f.read()
	#print('txt',txt2)
	str2=str(txt2,encoding='utf-8')
	try:
		result2=eval(str2)
		print('Got result:', result2)
		return result2
	except:
		pass
	return json_duplicate

def get_job_status(solved_id):
	print('#######getting job status######')
	url_return='http://nova.astrometry.net/api/jobs/'+str(solved_id)+'/'
	data2 = {'request-json': json}
	print('Sending form data:', data2)
	data2 = urlencode(data2)
	data2 = data2.encode('utf-8')
	#print('Sending data:', data2)
	headers = {}
	request = Request(url=url_return, headers=headers, data=data2)
	f = urlopen(request)
	txt = f.read()
	#print('txt',txt)
	str2=str(txt,encoding='utf-8')
	result2=eval(str2)
	print('Got result:', result2)
	return result2
	
def get_solved_id(url_return):
	while True:
		stat3 = sub_status(url_return)
		jobs = stat3['jobs']
		if len(jobs):
			for j in jobs:
				if j is not None:
					break
			if j is not None:
				#print('Selecting job id', j)
				solved_id = j
				print('##########solved id##########')
				print(solved_id)
				print('####################')
				break
		time.sleep(5)	
	return solved_id

def check_job(solved_id):
	while True:
		stat = get_job_status(solved_id)
		print('Got job status:', stat)
		if stat.get('status','') in ['success']:
			success = (stat['status'] == 'success')
			break
		time.sleep(5)
	
def Down_loading_fits(solved_id,write_to_name1):
	print('#########Downloading fits#########')
	url_final='http://nova.astrometry.net/new_fits_file/'+str(solved_id)+'/'
	f3 = urlopen(url_final)
	txt3 = f3.read()
	w = open('After_calibrated_'+write_to_name, 'wb')
	w.write(txt3)
	w.close()
	print('Wrote to', write_to_name1)

def Down_loading_fits_wget(solved_id,write_to_name1):
	print('#########Downloading fits#########')
	url_final='http://nova.astrometry.net/new_fits_file/'+str(solved_id)+'/'
	write_name='After_calibrated_'+write_to_name1
	wget.download(url_final,out=write_name)
	print('Wrote to', write_name)
####################################################

'''necessary varinate '''
r = requests.post('http://nova.astrometry.net/api/login', data={'request-json': json.dumps({"apikey": "wsjkkufuuwxvjovl"})})
m=r.text
dic=json.loads(m)
session=dic['session']

request_json={"allow_commercial_use": "d", "allow_modifications": "d", 
	"publicly_visible": "y", "session": session,"apikey": "wsjkkufuuwxvjovl"}
json=json.dumps(request_json)


'''main fanction'''



file_name_dic=folder_file_name('Pic_ready_to_upload')
for fn in file_name_dic:
	location='Pic_ready_to_upload'+'/'+fn
	result=up_load_file(location,session,json)
	write_to_name=fn
	result2,sub_id,url_return=upload_info(result,json)
	solved_id=get_solved_id(url_return)	
	check_job(solved_id)
	Down_loading_fits_wget(solved_id,write_to_name)
	









