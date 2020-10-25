import ephem 
import astropy.units as u
import math
import astropy.coordinates as coord



class data_list():
	def __init__(self,ra,dec,name='none',avalible_time='0'):
		self.ra=ra
		self.dec=dec
		self.name=name
		self.time=avalible_time

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

def sideral_time(date):
	gwd=ephem.Observer()
	gwd.long='117.5'
	gwd.lat='40.39'
	gwd.elevation=0
	gwd.date='2020/'+date
	LST=gwd.sidereal_time()
	fmt_LST=''
	for n in str(LST):
		if n==':':
			fmt_LST+=' '
		else:
			fmt_LST+=n
	LST_coord=coord.SkyCoord(fmt_LST+' +00 00 00',unit=(u.hourangle, u.deg),frame='fk5')
	return LST_coord




def get_avalible_region(ra,dec):
	avalible_region_ra=[]
	avalible_region_dec=[]
	avalible_region_date=[]
	transit_time_list=[]
	lunar_date={'05/05':['12h22m','19h42m'],'05/06':['12h23m','19h41m'],'05/07':['12h24m','19h40m'],'05/08':['12h25m','19h39m'],'05/09':['12h26m','19h38m']}
	UT_8=coord.SkyCoord('08h00m00s','00d00m',frame='fk5')
	UT_1=coord.SkyCoord('01h00m00s','00d00m',frame='fk5')
	UT_24=coord.SkyCoord('24h00m00s','00d00m',frame='fk5')
	UT_extra=coord.SkyCoord('01h00m00s','00d00m',frame='fk5')
	fram1=coord.SkyCoord(ra,dec,frame='fk5')
	
	for date in lunar_date.keys():
		twiling_down_time_UTC=coord.SkyCoord(lunar_date[date][0],'00d00m',frame='fk5')
		twiling_up_time_UTC=coord.SkyCoord(lunar_date[date][1],'00d00m',frame='fk5')
		#twiling_down_time_local=twiling_down_time_UTC.ra+UT_8.ra-UT_extra.ra
		#twiling_up_time_local=twiling_up_time_UTC.ra+UT_8.ra-UT_24.ra+UT_extra.ra
		twiling_down_time_local=twiling_down_time_UTC.ra+UT_8.ra
		twiling_up_time_local=twiling_up_time_UTC.ra+UT_8.ra-UT_24.ra
		#print(twiling_down_time_local.hms)
		#print(twiling_up_time_local.hms,'\n')		
		for region in fram1:
			LST_coord=sideral_time(date)
			transit_time=region.ra-LST_coord.ra+UT_8.ra

			if transit_time<0:
				#print('today')
				transit_time=transit_time+coord.SkyCoord('24h00m00s','00d00m',frame='fk5').ra
				dt=date
				#print(transit_time.hms)
			else:
				#print('tommorow')
				day=date[3:5]
				dt=date[0:3]+str(int(day)+1)


			#print(transit_time.hms)
			
			#if twiling_down_time_local-UT_1.ra*0.5<transit_time or transit_time<twiling_up_time_local+UT_1.ra*0.5:
			if twiling_down_time_local-UT_1.ra<transit_time or transit_time<twiling_up_time_local+UT_1.ra: 
			#if twiling_down_time_local<transit_time or transit_time<twiling_up_time_local:                                    #make sure the transit time is between twil time, or at least near the
				avalible_region_ra.append(region.ra)																		   #tran sit time(meas it's altitude is high between twil)
				avalible_region_dec.append(region.dec)
				avalible_region_date.append(dt)
				transit_time_list.append(transit_time.hms)
				#print(transit_time.hms)
				#print(twiling_down_time_local.hms)
				#print(twiling_up_time_local.hms,'\n')
	write_position_to_txt(avalible_region_ra,avalible_region_dec,avalible_region_date,transit_time_list)
	return data_list(avalible_region_ra,avalible_region_dec,avalible_time=avalible_region_date)
				
def write_position_to_txt(ra,dec,date,transit_time,file_name='avalible_position.txt'):
	sk5_coord=coord.SkyCoord(ra*u.degree,dec*u.degree)
	with open(file_name,'wt') as f:
		n=0
		for m in sk5_coord:
			coordinate=m.to_string('hmsdms')
			transittime=str(transit_time[n])
			hour=transittime[12:15]+'h'
			minute=transittime[19:23]+'m'
			second=transittime[27:32]+'s'
			transittime=hour+minute+second

			f.writelines(coordinate+' '+date[n]+' '+transittime+'\n')
			#f.writelines(coordinate+'\n')
			n+=1


m=extract_data_format_hms_dms('galactic_free_position.txt')
n=get_avalible_region(m.ra,m.dec)


# city = ephem.city('Beijing')
# #print('%s %s' % (city.lat, city.lon))

# sitka = ephem.Observer()
# sitka.date = '2020/05/10'
# sitka.lat = '40.39'
# sitka.lon = '117.5'
# m = ephem.Mars()
# print(sitka.next_transit(m))
