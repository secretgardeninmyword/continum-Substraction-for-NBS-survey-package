import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

import astropy.units as u
import math
import astropy.coordinates as coord


class data_list():
	def __init__(self,ra,dec):
		self.ra=ra
		self.dec=dec

def write_position_to_txt(ra,dec,file_name):
	sk5_coord=coord.SkyCoord(ra*u.degree,dec*u.degree)
	with open(file_name,'wt') as f:
		for m in sk5_coord:
			coordinate=m.to_string('hmsdms')
			f.writelines(coordinate+'\n')
			
def convert_nparray_to_list(nparray):
	m=[]
	for unit1 in nparray:
		m.append(unit1)
	return m

'''fmt like L20171202_06367_023246+6128_180S_S2_257608.FITS'''
def extract_source_data(file_name):
	data=[]
	RA_t=[]
	DEC_t=[]
	
	with open(file_name,'r') as f:
		for line in f:
			#print(line)
			data.append(line)
	
	for n in data:
		RA_t.append(n[16:22])
		DEC_t.append(n[22:27])
	return data_list(RA_t,DEC_t)
'''fmt like "XXhXXmXX.XXXXs +XXdXXmXX.XXXXs"'''
def extract_data_format_hms_dms(file_name):
	data=[]
	RA_t=[]
	DEC_t=[]
	with open(file_name,'r') as f:
		for line in f:
			#print(line)
			data.append(line)		
	for n in data:
		RA_t.append(n[0:14])
		DEC_t.append(n[15:29])
	return data_list(RA_t,DEC_t)
'''fmt like "NAME  XX.XXXXXX  XX.XXXXXX"'''	
def extract_data_format_d_d(file_name):
	data=[]
	RA_t=[]
	DEC_t=[]
	with open(file_name,'r') as f:
		for line in f:
			#print(line)
			data.append(line)	
	m=0	
	for n in data:	
		RA_t.append(float(n[5:15]))
		DEC_t.append(float(n[17:26]))
		m+=1
		#if m==50:
			#break
	return data_list(RA_t,DEC_t)
def Boundary_position(ra,dec):
	try:
		fk5_coord=coord.SkyCoord(ra*u.degree,dec*u.degree)
	except:
		fk5_coord=coord.SkyCoord(ra,dec)
	gltc_coord=fk5_coord.transform_to('galactic')
	restricted_ra=[]
	restricted_dec=[]
	for m in gltc_coord:

		print('calculating whether in boundery.......')
		out_of_place='false'
			
		#if m.b.degree>=6.5 or m.b.degree<=-6.5:          # 3*3degree region
		if m.b.degree>=6 or m.b.degree<=-6:					#2*2degree region
			out_of_place='yes'

		if out_of_place=='false':
			restricted_dec.append(m.b.degree)
			restricted_ra.append(m.l.degree)
		gltc_coord=np.delete(gltc_coord,0)
	result_data_glt=coord.SkyCoord(restricted_ra*u.degree, restricted_dec*u.degree,frame='galactic')
	result_data_fk5=result_data_glt.transform_to('fk5')
	RA=convert_nparray_to_list(result_data_fk5.ra.degree)
	DEC=convert_nparray_to_list(result_data_fk5.dec.degree)
	return data_list(RA,DEC)

def calculate_valid_position(ra_data,dec_data,whether_write_to_file='no',file_name='valid_position.txt'):
	fram1=coord.SkyCoord(ra_data*u.degree, dec_data*u.degree,frame='fk5')
	restricted_ra=[]
	restricted_dec=[]
	for m in fram1:
		ttl_nbr=len(fram1)-1
		print('calculating.......')
		#print(fram1)
		same_place='false'
		if ttl_nbr==0:
			pass
		else:
			for n in range(ttl_nbr):	
				sep=m.separation(fram1[n+1])
				if sep.degree<=1:
					same_place='yes'
					break
				else :
					pass
		if same_place=='false':
			restricted_dec.append(m.dec.degree)
			restricted_ra.append(m.ra.degree)
		fram1=np.delete(fram1,0)
	result_data=coord.SkyCoord(restricted_ra*u.degree, restricted_dec*u.degree,frame='fk5')
	RA=convert_nparray_to_list(result_data.ra.degree)
	DEC=convert_nparray_to_list(result_data.dec.degree)
	if whether_write_to_file=='yes':
		write_position_to_txt(RA,DEC,file_name)
	return data_list(RA,DEC)

def converg_to_RA_DEC_coord(RA_t,DEC_t):
	RA=[]
	DEC=[]
	
	for m in RA_t:  
		ra = int(m[0:2])*15+int(m[2:4])*15/60+int(m[4:6])*15/3600
		RA.append(ra)
	for m in DEC_t: 
		if m[0]=='+':
			dec = int(m[1:3])+int(m[3:5])/60*1
		elif m[0]=='-':
			dec = (int(m[1:3])+int(m[3:5])/60)*(-1)
		DEC.append(dec)
	return data_list(RA,DEC)

def scale_the_field_to_3degree_multiple_3degree(RA,DEC,full_filled='no'):    # 3 degree or 2 degree
	RA_expend=[]
	DEC_expend=[]
	
	''' 3*3degree region'''	
	# for offset in range(len(RA)):
		# for n in range(31):
			# DEC_expend.append(DEC[offset]-1.5+0.1*n)
			# RA_expend.append(RA[offset]-1.5)
			# DEC_expend.append(DEC[offset]-1.5+0.1*n)
			# RA_expend.append(RA[offset]+1.5)
			
			# DEC_expend.append(DEC[offset]-1.5)
			# RA_expend.append(RA[offset]-1.5+0.1*n)
			# DEC_expend.append(DEC[offset]+1.5)
			# RA_expend.append(RA[offset]-1.5+0.1*n)
	''' 2*2degree region'''	
	# if full_filled=='yes':	
		# for offset in range(len(RA)):               it is used to full-fill the ractangle, but with the use of fill_between, we don't need it anymore
			# for n in range(21):
				# DEC_expend.append(DEC[offset]-1+0.1*n)
				# RA_expend.append(RA[offset]-1)
				# DEC_expend.append(DEC[offset]-1+0.1*n)
				# RA_expend.append(RA[offset]-0.5)
	
				# DEC_expend.append(DEC[offset]-1+0.1*n)
				# RA_expend.append(RA[offset]-0)
				# DEC_expend.append(DEC[offset]-1+0.1*n)
				# RA_expend.append(RA[offset]+0.5)
				# DEC_expend.append(DEC[offset]-1+0.1*n)
				# RA_expend.append(RA[offset]+1)
	
	
							
				# DEC_expend.append(DEC[offset]-1)
				# RA_expend.append(RA[offset]-1+0.1*n)
				# DEC_expend.append(DEC[offset]-0.5)
				# RA_expend.append(RA[offset]-1+0.1*n)	
				# DEC_expend.append(DEC[offset]-0)
				# RA_expend.append(RA[offset]-1+0.1*n)
				# DEC_expend.append(DEC[offset]+0.5)
				# RA_expend.append(RA[offset]-1+0.1*n)
				# DEC_expend.append(DEC[offset]+1)
				# RA_expend.append(RA[offset]-1+0.1*n)
	
	for offset in range(len(RA)):
		for n in range(21):
			DEC_expend.append(DEC[offset]-1+0.1*n)
			RA_expend.append(RA[offset]-1)
			DEC_expend.append(DEC[offset]-1+0.1*n)
			RA_expend.append(RA[offset]+1)
			
			DEC_expend.append(DEC[offset]-1)
			RA_expend.append(RA[offset]-1+0.1*n)
			DEC_expend.append(DEC[offset]+1)
			RA_expend.append(RA[offset]-1+0.1*n)
							
	return data_list(RA_expend,DEC_expend)

def print_the_celestial_coord(RA_expend,DEC_expend,coord_frame='fk5',plot_milki_way_boundry='yes',input_coords_type='fk5',combine_order='none',axes_extra='none',color='blue',valid_point='null'):
	nra=np.array(RA_expend)
	ndec=np.array(DEC_expend)


	if 	combine_order=='none' or combine_order=='first':
		fig, axes = plt.subplots(subplot_kw={'projection': 'aitoff'})
		
		axes.set_title(coord_frame)
		
	else:
		pass
	if axes_extra=='none':
		pass
	else:
		axes=axes_extra
	if input_coords_type=='fk5'	:
		if coord_frame=='fk5':
			icrs=coord.SkyCoord(nra*u.degree, ndec*u.degree,frame=coord_frame)
			try:
				axes.plot(icrs.ra.wrap_at(180*u.deg).radian,icrs.dec.radian,linestyle='none',
					marker='.',markersize=0.1,color=color,alpha=1)
			except:
				axes.plot(icrs.ra.wrap_at(180*u.deg).radian,icrs.dec.radian,linestyle='none',
					marker='.',markersize=1)				
			# if color=='yellow':
				# axes.plot(icrs.ra.wrap_at(180*u.deg).radian,icrs.dec.radian,linestyle='none',
					# marker='.',markersize=1,color='y')
			# if color=='green':
				# axes.plot(icrs.ra.wrap_at(180*u.deg).radian,icrs.dec.radian,linestyle='none',
					# marker='.',markersize=1,color='green')
			# else:
				# axes.plot(icrs.ra.wrap_at(180*u.deg).radian,icrs.dec.radian,linestyle='none',
					# marker='.',markersize=1)			
		if coord_frame=='galactic':
			icrs_fk5=coord.SkyCoord(nra*u.degree, ndec*u.degree,frame='fk5')
			icrs=icrs_fk5.transform_to('galactic')

			axes.plot(icrs.l.wrap_at(180*u.deg).radian,icrs.b.radian,linestyle='none',
				marker='.',markersize=1)		
	if input_coords_type=='galactic':
		if coord_frame=='galactic':
			icrs=coord.SkyCoord(nra*u.degree, ndec*u.degree,frame='galactic')
			

			axes.plot(icrs.l.wrap_at(180*u.deg).radian,icrs.b.radian,linestyle='none',
				marker='.',markersize=1)		
		if coord_frame=='fk5':
			icrs_fk5=coord.SkyCoord(nra*u.degree, ndec*u.degree,frame='galactic')
			icrs=icrs_fk5.transform_to('fk5')

			axes.plot(icrs.ra.wrap_at(180*u.deg).radian,icrs.dec.radian,linestyle='none',
				marker='.',markersize=1)	
				
	if 	combine_order=='none' or combine_order=='last':
		pass
	else:
		plot_milki_way_boundry='no'
				
	if plot_milki_way_boundry=='no':
		if combine_order=='none':
			plt.show()
		if combine_order=='last':
			plt.show()
		else:
			pass
	fill_2degree2_region(axes,valid_point,coord_frame,color)		
	if plot_milki_way_boundry=='yes':
		print_the_celestial_coord_boundry(axes,coord_frame=coord_frame)
	return axes

def fill_2degree2_region(axes,valid_point,coord_frame,color):
	if valid_point=='null':
		pass
	else:
		np_point_ra=np.array(valid_point.ra)
		np_point_dec=np.array(valid_point.dec)
		icrs=coord.SkyCoord(np_point_ra*u.degree, np_point_dec*u.degree,frame=coord_frame)
		icrs_unit=coord.SkyCoord(1*u.degree, 1*u.degree,frame=coord_frame)
		for n in icrs:
			n_ra=n.ra.wrap_at(180*u.deg).radian
			n_dec=n.dec.wrap_at(180*u.deg).radian
			unit=icrs_unit.ra.wrap_at(180*u.deg).radian
			try:
				axes.fill_between([n_ra-unit,n_ra+unit],n_dec-unit,n_dec+unit,facecolor=color,alpha=1)	
			except:
				axes.fill_between([n_ra-unit,n_ra+unit],n_dec-unit,n_dec+unit,alpha=1)	
			
def print_the_celestial_coord_boundry(axes,coord_frame):
	milki_way_above_l=np.arange(29,215,0.1)
	milki_way_above_b=np.ones(1860)*5
	milki_way_belove_l=np.arange(29,215,0.1)
	milki_way_belove_b=np.ones(1860)*(-5)
	milki_way_l=np.arange(29,215,1)
	milki_way_b=np.ones(186)*0
	gltc_above=coord.SkyCoord(milki_way_above_l*u.degree,
		milki_way_above_b*u.degree,
		frame='galactic')
	gltc_belove=coord.SkyCoord(milki_way_belove_l*u.degree,
		milki_way_belove_b*u.degree,
		frame='galactic')
	gltc=coord.SkyCoord(milki_way_l*u.degree,
		milki_way_b*u.degree,
		frame='galactic')
	if coord_frame=='fk5':
		print('fk5')
		fk5_above=gltc_above.transform_to('fk5')
		fk5_belove=gltc_belove.transform_to('fk5')
		fk5=gltc.transform_to('fk5')
		#fig, axes = plt.subplots(subplot_kw={'projection': 'aitoff'})
		axes.plot(fk5_above.ra.wrap_at(180*u.deg).radian,fk5_above.dec.radian,linestyle='none',
			marker='.',markersize=1,color='r')
		axes.plot(fk5_belove.ra.wrap_at(180*u.deg).radian,fk5_belove.dec.radian,linestyle='none',
			marker='.',markersize=1,color='r')
		axes.plot(fk5.ra.wrap_at(180*u.deg).radian,fk5.dec.radian,linestyle='dashed',
			marker='.',markersize=1,color='g')
		
		
	if coord_frame=='galactic':
		#fig, axes = plt.subplots(subplot_kw={'projection': 'aitoff'})
		print('galactic')
		axes.plot(gltc_above.l.wrap_at(180*u.deg).radian,gltc_above.b.radian,
			marker='.',markersize=0.1,color='r')
		axes.plot(gltc_belove.l.wrap_at(180*u.deg).radian,gltc_belove.b.radian,
			marker='.',markersize=0.1,color='r')
		axes.plot(gltc.l.wrap_at(180*u.deg).radian,gltc.b.radian,linestyle='dashed',
			marker='.',markersize=0.1,color='g')
	plt.savefig('fig_2.eps',format='eps')
	#plt.savefig('fig_1.png',format='png',dpi=1000,figsize=(10,4))	
	plt.show()
	
def print_centrial_galactic_region(scale_to_ractangle='yes'):
	list_longitutde=[]
	list_latitude=[]
	for n in range(29,215,2):
		for m in range(-5,5,2):
			list_longitutde.append(n)
			list_latitude.append(m)
	if scale_to_ractangle=='yes':
		Expended_Position=scale_the_field_to_3degree_multiple_3degree(list_longitutde,list_latitude)
		#print(Expended_Position.ra)
		m=converge_cellestial_type(Expended_Position.ra,Expended_Position.dec,input_type='galactic',output_type='fk5')
		return data_list(m.ra,m.dec)
	if scale_to_ractangle=='no':
		m=converge_cellestial_type(list_longitutde,list_latitude,input_type='galactic',output_type='fk5')
		return data_list(m.ra,m.dec)
def refine_print_centrial_galactic_region(RA,DEC,already_exist_ra,already_exist_dec,whether_write_to_file='yes',file_name='galactic_free_position.txt'):
	fram1=coord.SkyCoord(RA*u.degree, DEC*u.degree,frame='fk5')
	fram_already=coord.SkyCoord(already_exist_ra*u.degree, already_exist_dec*u.degree,frame='fk5')
	restricted_ra=[]
	restricted_dec=[]
	for m in fram1:
		ttl_nbr=len(fram1)-1
		print('calculating.......')
		#print(fram1)
		same_place='false'
		if ttl_nbr==0:
			pass
		else:
			for n in fram_already:	
				sep=m.separation(n)
				#if sep.degree<=3:          #choose this when format is 3*3
				if sep.degree<=2:			#or choose this when format is 2*2
					same_place='yes'
					break
				else :
					pass
		if same_place=='false':
			restricted_dec.append(m.dec.degree)
			restricted_ra.append(m.ra.degree)
		fram1=np.delete(fram1,0)
	result_data=coord.SkyCoord(restricted_ra*u.degree, restricted_dec*u.degree,frame='fk5')
	RA=convert_nparray_to_list(result_data.ra.degree)
	DEC=convert_nparray_to_list(result_data.dec.degree)
	if whether_write_to_file=='yes':
		write_position_to_txt(RA,DEC,file_name)
	return data_list(RA,DEC)
	

def converge_cellestial_type(RA,DEC,input_type,output_type):
	try:
		icrs=coord.SkyCoord(RA*u.degree, DEC*u.degree,frame=input_type)
	except:
		icrs=coord.SkyCoord(RA, DEC,frame=input_type)
	icrs=icrs.transform_to(output_type)
	try:
		return data_list(icrs.ra,icrs.dec)
	except:
		return data_list(icrs.l,icrs.b)
	


	
	
''' different expressing, internal class can transfer form hour/min/sec
to degree or inverse
#fk5=coord.FK5('00h19m01.920s','+65d01m48.000s')
w=coord.SkyCoord('00:19:01.920 +65:01:48',unit=(u.hourangle,u.deg),frame='fk5')
#fk5=coord.FK5(ra=4.758*u.degree,dec=65.03*u.degree)
c = coord.SkyCoord(ra=4.758*u.degree,dec=65.03*u.degree, frame='fk5')
print(c)
#print(fk5)
print(w)
'''


'''draw the valid position of target, all the target have been detected 
'''
# Source_data=extract_source_data('S2_NBS_filename.txt')
# Formatted_Source_data=converg_to_RA_DEC_coord(Source_data.ra,Source_data.dec)
# Position_in_Boundry=Boundary_position(Formatted_Source_data.ra,Formatted_Source_data.dec)
# Valid_position_in_degree=calculate_valid_position(np.array(Position_in_Boundry.ra),np.array(Position_in_Boundry.dec),whether_write_to_file='yes',file_name='S2_NBS_valid_number.txt')        # you must use numpy typed list

# Expended_Position_With_Ractangle=scale_the_field_to_3degree_multiple_3degree(Valid_position_in_degree.ra,Valid_position_in_degree.dec)
# print_the_celestial_coord(Expended_Position_With_Ractangle.ra,Expended_Position_With_Ractangle.dec,coord_frame='fk5',plot_milki_way_boundry='yes')


'''draw the undetected position whitin galaxy'''
''' by the help of outgrid.txt, we no longer need to product the position 
point by ourself
# Formatted_Source_data=extract_data_format_hms_dms('S2_NBS_valid_number.txt')

# #######used to check detacted region whether in region required, but now because extract already-fixed data, it is unneccesary
# ####Position_in_Boundry=Boundary_position(Formatted_Source_data.ra,Formatted_Source_data.dec)     
# ####Valid_position_in_degree=calculate_valid_position(np.array(Position_in_Boundry.ra),np.array(Position_in_Boundry.dec),whether_write_to_file='no')   # you must use numpy typed list
# ###################
'''


# alreay_exit_position=converge_cellestial_type(Formatted_Source_data.ra,Formatted_Source_data.dec,input_type='fk5',output_type='fk5')
# Formatted_Source_data=extract_data_format_d_d('outgrid.txt')
# centrial_region=Boundary_position(Formatted_Source_data.ra,Formatted_Source_data.dec)
# refined_regin=refine_print_centrial_galactic_region(np.array(centrial_region.ra),np.array(centrial_region.dec),np.array(alreay_exit_position.ra),np.array(alreay_exit_position.dec),
#	file_name='galactic_free_position.txt')                                                                             ########you must use numpy typed list
# rectangled_refine_position=scale_the_field_to_3degree_multiple_3degree(refined_regin.ra,refined_regin.dec,full_filled='no')
# print_the_celestial_coord(rectangled_refine_position.ra,rectangled_refine_position.dec,coord_frame='fk5',plot_milki_way_boundry='yes')






'''refine and plot detectable region'''
'''refine and plot detectable region'''
'''refine and plot detectable region'''
'''refine and plot detectable region'''
# Formatted_Source_data=extract_data_format_hms_dms('avalible_position.txt')
# Valid_position_in_degree=Boundary_position(Formatted_Source_data.ra,Formatted_Source_data.dec)
# Expended_Position_With_Ractangle=scale_the_field_to_3degree_multiple_3degree(Valid_position_in_degree.ra,Valid_position_in_degree.dec,full_filled='no')
# print_the_celestial_coord(Expended_Position_With_Ractangle.ra,Expended_Position_With_Ractangle.dec,coord_frame='fk5',plot_milki_way_boundry='yes')



'''combine already-detected frame, all-galactic frame, detectable fram to one pic'''
'''combine already-detected frame, all-galactic frame, detectable fram to one pic'''
'''combine already-detected frame, all-galactic frame, detectable fram to one pic'''
'''combine already-detected frame, all-galactic frame, detectable fram to one pic'''


################# centrial_region=print_centrial_galactic_region('yes')
################ axes=print_the_celestial_coord(centrial_region.ra,centrial_region.dec,coord_frame='fk5',combine_order='first')

Formatted_Source_data=extract_data_format_d_d('outgrid.txt')
Valid_position_in_degree=Boundary_position(Formatted_Source_data.ra,Formatted_Source_data.dec)
Expended_Position_With_Ractangle=scale_the_field_to_3degree_multiple_3degree(Valid_position_in_degree.ra,Valid_position_in_degree.dec,full_filled='no')
axes=print_the_celestial_coord(Expended_Position_With_Ractangle.ra,Expended_Position_With_Ractangle.dec,coord_frame='fk5',combine_order='first',color='royalblue')


Formatted_Source_data=extract_data_format_hms_dms('avalible_position.txt')
Valid_position_in_degree=Boundary_position(Formatted_Source_data.ra,Formatted_Source_data.dec)
Expended_Position_With_Ractangle=scale_the_field_to_3degree_multiple_3degree(Valid_position_in_degree.ra,Valid_position_in_degree.dec)
axes=print_the_celestial_coord(Expended_Position_With_Ractangle.ra,Expended_Position_With_Ractangle.dec,coord_frame='fk5',combine_order='middle',axes_extra=axes,
	color='yellow',valid_point=Valid_position_in_degree)





Formatted_Source_data=extract_data_format_hms_dms('S2_NBS_valid_number.txt')
Valid_position_in_degree=Boundary_position(Formatted_Source_data.ra,Formatted_Source_data.dec)
Expended_Position_With_Ractangle=scale_the_field_to_3degree_multiple_3degree(Valid_position_in_degree.ra,Valid_position_in_degree.dec,full_filled='yes')
print_the_celestial_coord(Expended_Position_With_Ractangle.ra,Expended_Position_With_Ractangle.dec,coord_frame='fk5',combine_order='last',plot_milki_way_boundry='yes',
	axes_extra=axes,valid_point=Valid_position_in_degree,color='red')



