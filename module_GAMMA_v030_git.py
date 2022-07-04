# coding=utf8

### Module, with general functions, in purpose for GAMMA processing script.

## Author: Zlatomir Dimitrov (zlatomir.dimitrov@space.bas.bg), 14-NOV-2018
## 	 License: MIT License, see LICENSES #
##   Friedrich-Schiller-Universität - JENA, Lehrstuhl für Fernerkundung (Uni-JENA), Space Research and Technology Institute (SRTI-BAS)

### Ziel: PROCESS TANDEM-X DATA ###

### всички команди в GAMMA са нужни!

# ------------------------- DISPLAY ------------------------------------- #

def Display(datasets,width,code):	

	'''	  DISPLAY DATA WITH GAMMA  '''

	### DICTIONARY:
	  # code = 1: DEM display - "hgt" 						: command - dis2hgt
	  # code = 2: Two TDX imageries - "slc", "mli" 			: command - dis2pwr
	  # code = 3: Raster imagery - interferogram, LUT, etc.	: command - dismph
	import os
	from datetime import datetime

	path_proc='/homes3/geoinf/ye68hah/GAMMA/Gamma2018_Training/L5_Phaseheight_TDX/proc/'
	lf=os.path.join(path_proc,"Display-cmd.log.sh")
	#
	#exist=check_exists("disp_cmd.log",path_proc)
	#if exist==True:
	#	os.system("rm -f "+lf)	
	#
	t = datetime.now()
	f=open(lf,'a')
	s="\n*** "+str(t)+" ***\n"
	f.write(s)

	# init, dict:  #'code':'command',
	Disp={0:'rasmph',\
			1:'dismph',     \
			2:'dis2pwr',    \
			3:'dis2hgt',    \
			4:'dis_linear', \
			5:'dishgt',     \
			6:'disrmg',		\
			7:'dis2mph',    
	}

	if code==0:
		# rasmph
		cmd=Disp[0]+' '+datasets+' '+str(width)+' - - - - - - - '+datasets+'.bmp'+' 0'
		#print(cmd)

	elif code==1:
		# dismph 
		cmd=Disp[1]+' '+datasets+' '+str(width)

	elif code==2:
		# dis2pwr
		cmd=Disp[2]+' '+datasets[0]+' '+datasets[1]+' '+str(width[0])+' '+str(width[1])

	elif code==7:
		# dis2mph
		cmd=Disp[7]+' '+datasets[0]+' '+datasets[1]+' '+str(width)+' '+str(width)
		
	elif code==3:
		# dis2hgt
		start='-'; nlines='-'; roff='-'; azoff='-'; meters_per_cycle = '60'
		cmd=Disp[3]+' '+datasets[0]+' '+datasets[1]+' '+str(width[0])+' '+str(width[1])+' '+start+' '+nlines+' '+roff+' '+azoff+' '+meters_per_cycle
	
	elif code==4:
		# dis_linear
		cmd=Disp[4]+' '+datasets+' '+str(width)
	
	elif code==5:
		# dishgt
		cmd=Disp[5]+' '+datasets[0]+' '+datasets[1]+' '+str(width)+' - - - 40'
	
	elif code==6:
		# disrmg
		cmd=Disp[6]+' '+datasets[0]+' '+datasets[1]+' '+str(width)+' 1 1 - 1'        


	else:
		exit('Display> None of correct parameters provided! \n   ...Exit.')


	# EXECUTE:	
	output = exec_shell_cmd(cmd)
	#
	f.write(cmd+"\n")

# ------------------------- Parse Dictionary, for exec GAMMA commands ------------------------------------- #
def parse_dict(cmd,params):

	'''  Parse Dictionary, for execution of GAMMA commands, finally provide a string   '''
	
	cmd+=' '
	###	
	for val in params.itervalues():	
		cmd+=str(val)+' '
	####
	return cmd

# ------------------------- Find files needed, using pyroSAR ------------------------------------- #

def find_files(datadir,file,ext):

	'''	FIND FILES : with pyroSAR   '''
	
	from pyroSAR.ancillary import finder

	FILES=[]

	if len(ext)!=0: # search by File-Extention: #here, for DEMs ; WARNING, тук използвам - file за DEM!
		# pyroSAR:
		founded = finder(datadir, ["*."+ext], recursive=True) 
		if(not(len(founded))):
			err("Input Data NOT Found!",file,1)
		else:		
			for val in founded:
				#append founded file names in to return array
				FILES.append(val.replace(datadir,'').strip('/'))

				# DEM:
				if(file=='DEMs'):
					if (val.find("dem")!=-1 or val.find("srtm")!=-1 or val.find("height")!=-1 or val.find("dem")!=-1):			
						print '# Found, DEM file: "'+val+'"'
					else:
						err("DEM file not found! \n> For TanDEM-X needed DEM files are: dem.tif, dem.tif, knn_height.tif.",ft,1)
				elif(file=='TDXs'):
					pass
					#
					#
					#

	else:	# search by file name:	; WARNING, тук използвам - file за ИМЕТО НА ФАЙЛА! 		
		# ext - НЕ СЕ ИЗПОЛЗВА ТУК!
		founded = finder(datadir, [file], recursive=True) 

		for val in founded:
			#append founded file names in to return array
			FILES.append(val.replace(datadir,'').strip('/'))   # връща само един елемент, на масив.

	return FILES

# ------------------------- File handling for Int.parameters, miro ------------------------------------- #

def handle_file(fn,tdx,**kwargs):
	#
	# open file
	# write file
	#	
	with open(fn,'a') as f:
		f.write("*** Interferometric parameters of the TanDEM-X pair ***\n* Master: %s, Slave: %s *\n\n"%(tdx[0],tdx[1]))
		f.write("# Ground Range resolution, will be: ")

		for key,value in kwargs.items():	# iterate on Dict! Miro

			u = "" if (key=='kz') else "[m]"	## Inline IF assign to var! Bravo,Miro!
			s = " %s = %f %s \n"%(key,value,u)
			f.write(s) 

	# end!

# ------------------------- Parse PAR file Metadata, miro ------------------------------------- #

def parsePAR(parf,partoparse):

	''' PARSE PAR METADATA FILE : miro version '''

	parameter=[]

	with open(parf, 'r') as f: 
		for w in f.readlines():
			#print w
			#if not(w.find(partoparse)): width=float((((' '.join(w.splitlines())).strip(partoparse)).strip(':')).strip(' ')) 	–– хубаво, ама, работи само, ако e една стойност; ако е низ от стойности - НЕ! Затова, в array!
			
			if not(w.find(partoparse)): 
				sp = w.strip(' ').split() 
				for i in range(len(sp)):
					s=sp[i]
					#print s
					if(i>0): parameter.append(s)
				 

	if not len(parameter): parameter.append('Parameter <'+partoparse+'> Not found!')

	return parameter

# ------------------------- Parse PAR file Metadata, pyroSAR ------------------------------------- #

def parsePAR_pyroSAR(path,parfile,parameters):

	'''	PARSE PAR METADATA  : with pyroSAR - "ISPPar" class   '''

	import os
	from pyroSAR.gamma import ISPPar # import class
		
	readed_par=[]

	# ИЗГРАЖДАНЕ НА DICTIONARY: 
	resdict=ISPPar(os.path.join(path,parfile))
	# асоцииране на СВОЙСТВАТА на Dictionary - resdict:
	parlist = vars(resdict)

	for p in parameters:
		par = parlist[p]
		readed_par.append(par)

	#print readed_par

	return readed_par

# ------------------------- Parse Log-file to extract Bperp  ------------------------------------- #

def parse_log_extract_Bperp(log):

	''' Parse LOG-file from base_perp to extract Bper - component. Could also the Bpara...	'''

	with open(log,'r') as f:
		#
		data=f.readlines()
		i=0
		for d in data:
			i+=1
			#print ("%i: %s"%(i,d))
			if i==45: Bperp=float(d.strip('baseline perpendicular component (m): '))
		##

	#print ("Bperp = %f"%(Bperp))

	return Bperp

# ------------------------- Check existing files  ------------------------------------- #

def check_exists(file,path):

	'''	Check, rather the files we while to create are already creataed? 	'''
	
	exist=False

	# check:
	Files = find_files(path,file,"")

	#print (Files)

	if len(Files)>0: 
		
		try:
			f=open(Files[0],'r')
			if f: exist=True
		except Exception as e:
			exist=False
		
	return exist

# ------------------------- Height difference calculation inbetween of DEM / dem ------------------------------------- #

def exec_shell_cmd(cmd):

	''' EXECUTE BASH SHELL COMMANDS '''

	import subprocess	
	
	e=0
	err=''

	#system(cmd)
	#p = subprocess.Popen([cmd], cwd=datadir)
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	(output, err) = p.communicate()	# взима изхода от stderr i stdout
	p_status = p.wait()	# obrabotva res
	# печат на двете променливи:
	#print "Command output : ", output
	#print "Command exit status/return code : ", p_status

	## Check for Errors:
	e=output.find('ERROR')
	if (e==0 or err is not None):
		print("#")*28
		print output.find('ERROR')
		print ' §§§ '
		print type(err)
		print("<STD-ERR> There was an ERROR with import! \n> Check out the output Errors, if critical - <Ctr>+C to terminate.")
		raw_input("> Continue executing script?")
		print("#")*28		

	# screiben log:
	f=open('log.txt','a')
	s='*\n'+('*')*28+' '+cmd+' '+('*')*28+'\n*\n'
	f.write(s)
	f.write(output)
	f.close()

	return output

# ------------------------- Height difference calculation inbetween of DEM / dem ------------------------------------- #

def calc_height_diff(mode,width):

	''' CALCULATE DEM-difference in height, between DEM & dem   '''

	roff='-'
	loff='-'
	zflg=1
	c0=1
	
	if(mode==0):
		print '#\n# Summ Geoid with Geoid-offset, to produce dem (dem) :'
		# GAMMA command:
		cmd = "float_math "+"dem_geoid "+"geoid_offset "+"dem "+str(width)+' '+str(mode)+' '+roff+' '+loff+' '+str(zflg)+' '+str(c0)

	elif(mode==1):			
		print '#\n# Substract dem from DEM, to calculate DEM - Height Differce (hgtdiff) :'	
		# GAMMA command:
		cmd = "float_math "+"srtm "+"dem "+"hgtdiff "+str(width)+' '+str(mode)+' '+roff+' '+loff+' '+str(zflg)+' '+str(c0)

	# EXECUTE:	
	output = exec_shell_cmd(cmd)

	return output

# ------------------------- import data in GAMMA ------------------------------------- #

def import_in_GAMMA(datadir,dems):

	'''	IMPORT SAR & DEM DATA IN GAMMA   '''

	# промяна, защото ще се изпълни в текущата директория на скрипта.
	#outpath='output/'
	#	
	
	for i in range(len(dems)):
		output='' 	#clear output
		demf = dems[i-1].replace('.tif','')
		parf = dems[i-1]+'.par'
		space=' '
		geoid='-'		#(input) geoid offset relative to the WGS84 datum: gflg:2 GeoTIFF file containing geoid offset in the WGS84 datum 

		if(demf=='srtm'):

			#exists=check_exists(demf,datadir)
			#if exists==True: return None

			print '#\n# import:',dems[i-1]
			gflg=2		#geoid offset correction flag: 2: add interpolated geoid offset relative to the WGS84 datum
			geoid_dem=' geoid_offset'	#(output) resampled geoid offset map in the DEM grid (gflg: 2)
			# GAMMA command:
			cmd = "srtm2dem " + dems[i-1] + space + demf + space + parf + space + str(gflg) + space + str(geoid) + geoid_dem						
			

		elif(demf=='dem'):
			
			#exists=check_exists(demf,datadir)
			#if exists==True: return None

			print '#\n# import:',dems[i-1]
			gflg=0
			demf+='_geoid'
			geoid_dem=''
			# GAMMA command:
			cmd = "srtm2dem " + dems[i-1] + space + demf + space + parf + space + str(gflg) + space + str(geoid) + geoid_dem

			### Read width of DEM, PAR-file:
			width=float(''.join(parsePAR(parf,'width')))	#ВНИМАНИЕ: ТУК връща array!
			#print '\n*** width = %8.4f'% (width)


		elif(demf=='knn_height'):	
			
			#exists=check_exists(demf,datadir)
			#if exists==True: return None

			print '#\n# import:',dems[i-1]	
			# GAMMA command:
			cmd = "par_data_geo " + dems[i-1] + space + parf + space + demf			
		
		# EXECUTE:	
		output = exec_shell_cmd(cmd)
		#
		#raw_input("#\n# press Enter to continue...")
		
	print ("\n### Import Finished!")
	#
	# END, cycle

	## Calculate DEM-difference height:
	#
	#print '\n*** width = %8.4f'% (width)
	#raw_input("#\n# press Enter to continue...")
	#
	# first: dtem_geoid + geoid_offset = dem !
	mode=0
	calc_height_diff(mode,width)
	print '## OK!'

	# second: srtm - dem = hgtdiff !
	mode=1
	calc_height_diff(mode,width)
	print '## OK!'

	## return imported DEMs:
	#return...

	## return imported TDX-SAR-imagery:
	#return...

# ------------------------- Calculate Multi-look factors ------------------------------------- #

def calc_MultiLooks(rr,ra,thi,ml_az):

	'''   CALCULATE MULTI-LOOK FACTORS   '''
	##
	#def multilook(infile, outfile, targetres):
	#https://pyrosar.readthedocs.io/en/latest/_modules/pyroSAR/gamma/util.html#multilook
	# inputed function arguments: 'range_samples','range_pixel_spacing','azimuth_pixel_spacing','incidence_angle

	import math

	thi_r = thi*math.pi/180.0

	rg = rr/math.sin(thi_r)

	raml = ra*ml_az

	mlrg = round(raml/rg)

	#return mlrg,rg
	return rg

# ------------------------- Multi-Look : miro ------------------------------------------------------- #

def Multi_Look_m(datasets,ml_rg,ml_az):

	'''	MULTI-LOOK OF TDX-DATA   '''
	
	# reqiures: rg = calc_MultiLooks(...)

	slc='_slc'; mli='.mli'; par='.par'
	for tdx in datasets:
		# GAMMA command:
		cmd = "multi_look "+tdx+slc+' '+tdx+slc+par+' '+tdx+mli+' '+tdx+mli+par+' '+ml_rg+' '+ml_az
		# EXECUTE:	
		output = exec_shell_cmd(cmd)

	return datasets[0]+mli, datasets[1]+mli

# ------------------------- Multi-Look : pyroSAR ------------------------------------------------------- #

def Multi_Look(datasets,workdir,lgfn,targetres):

	'''	MULTI-LOOK OF TDX-DATA, by using GAMMA-API from pyroSAR   '''
	
	# reqiures: rg = calc_MultiLooks(...)

	import os
	from pyroSAR.gamma.api import isp
	from pyroSAR.gamma import multilook
	from collections import OrderedDict

	slc='_slc'; mli='.mli'; par='.par' 
	orig=workdir[0]; proc=workdir[1]; odata=[]

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' Multi_Look '+('#')*28+'\n')

	#GAMMA:
	for tdx in datasets:
		tdx_mli = os.path.join(proc, tdx+mli)
		tdx_mli_par = os.path.join(proc, tdx+mli+par)
 		
 		exists=check_exists(tdx+mli,proc)
 		if exists==True: 			
 			print("#> Skip Multi-Look. File: %s exists . " %(tdx_mli))
 			#print("#> Skip Multi-Look. Data exists. ")
 			#return None
 		else:
 			tdx_slc = os.path.join(orig, tdx+slc )
 			#pyroSAR:
			ml_rg,ml_az = multilook(tdx_slc, tdx_mli, targetres)

			# build parameters dict:
			#parameters = {'SLC': os.path.join(orig, tdx+slc ),			\
			#			  'SLC_par': os.path.join(orig, tdx+slc+par ),	\
			#			  'MLI': tdx_mli,								\
			#			  'MLI_par': tdx_mli_par,						\
			#			  'rlks': ml_rg,								\
			#			  'azlks': ml_az								\
			#			  #'logpath': proc								\
			#			  }
			#
			# за да си напиша Log.txt.sh:
			parameters = OrderedDict([('SLC', os.path.join(orig, tdx+slc )),			\
						  ('SLC_par', os.path.join(orig, tdx+slc+par)),	\
						  ('MLI', tdx_mli),								\
						  ('MLI_par', tdx_mli_par),						\
						  ('rlks', int(ml_rg)),							\
						  ('azlks', int(ml_az))							\
						  #'logpath': proc								\
						  ])
			#
			#log:
			c=parse_dict('multi_look',parameters)
			lgf.write(c+"\n")
			#
			# EXECUTE, pyroSAR GAMMA-API:
			#isp.multi_look(**parameters)
			
		odata.append(tdx_mli)
		### end!

	lgf.close()

	return odata[0], odata[1], ml_rg,ml_az	# return, datafiles -MLI, with absolute path!

	## multi_look TSX_110604_VV_slc TSX_110604_VV_slc.par TSX_110604_VV.mli TSX_110604_VV.mli.par $ml_rg $ml_az #
	## multi_look TDX_110604_VV_slc TDX_110604_VV_slc.par TDX_110604_VV.mli TDX_110604_VV.mli.par $ml_rg $ml_az #


# ------------------------- Radiometric Calibration : miro ------------------------------------------------------- #

def Rad_Cal_m(datasets):

	'''	RADIOMETRIC CALIBRATION OF MULTI-LOOKED TDX-DATA   '''

	mli='.mli'; par='.par'; cmli='.cmli'

	OFF_par = '-' # (input) ISP offset/interferogram parameter file (enter - for images in MLI geometry)

	for tdx in datasets:
		# GAMMA command:
		cmd = "radcal_MLI "+tdx+mli+' '+tdx+mli+par+' '+tdx+mli+' '+OFF_par+' '+tdx+cmli
		# EXECUTE:	
		output = exec_shell_cmd(cmd)	

	return datasets[0]+cmli, datasets[1]+cmli	# return datafile names -CMLI

# ------------------------- Radiometric Calibration : pyroSAR ------------------------------------------------------- #

def Rad_Cal(datasets,workdir,lgfn):

	'''	RADIOMETRIC CALIBRATION OF MULTI-LOOKED TDX-DATA, by using GAMMA-API from pyroSAR   '''

	import os
	from pyroSAR.gamma.api import isp
	from collections import OrderedDict


	mli='.mli'; par='.par'; cmli='.cmli'; odata=[]
	
	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' Rad_Cal '+('#')*28+'\n')

	for tdx in datasets:
		
		tdx_cmli=os.path.join(workdir, tdx+cmli )

 		exists=check_exists(tdx+cmli,workdir)

 		if exists==True: 			
 			#print("#> Skip - Radiometric Calibration. File: %s exists . " % (tdx_cmli))
 			print("#> Skip - Radiometric Calibration! Data exists. ")
 		else:
			# build parameters dict:
			#parameters = {'MLI': os.path.join(workdir, tdx+mli ),			\
			#			  'MLI_PAR': os.path.join(workdir, tdx+mli+par ),	\
			#			  'OFF_par':'-',									\
			#			  'CMLI': tdx_cmli									\
			#			  #'logpath': proc
			#			  }
			parameters = OrderedDict([('MLI', os.path.join(workdir, tdx+mli )),			\
									  ('MLI_PAR', os.path.join(workdir, tdx+mli+par)),	\
									  ('OFF_par','-'),									\
									  ('CMLI', tdx_cmli)								\
									  #'logpath': proc
						  ])

			#log:
			c=parse_dict('radcal_MLI',parameters)
			lgf.write(c+"\n")
			

			# EXECUTE, pyroSAR GAMMA-API:
			isp.radcal_MLI(**parameters)

		odata.append(os.path.join(workdir,tdx+cmli))
	### end!

	lgf.close()

	return odata[0], odata[1]	# return datafile names -CMLI, with absolute path!	

	## radcal_MLI TSX_110604_VV.mli TSX_110604_VV.mli.par - TSX_110604_VV.cmli #
	## radcal_MLI TDX_110604_VV.mli TDX_110604_VV.mli.par - TDX_110604_VV.cmli #


# ------------------------- LOOK-UP-TABLE ------------------------------------------------------- #

def LUT(tdxm,dem,workdir,lgfn,targetres):

	'''	CALCULATE LOOK-UP-TABLE FOR GEOCODING   '''
	### gc_map - MISSING in GAMMA-API!
	
	import os
		### MISSING IN API!
	from collections import OrderedDict
	from pyroSAR.gamma import ovs

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' LUT '+('#')*28+'\n')

	par='.par'
	# calculate OVS (over sampling factors, for current DEM)
	lat_ovr, lon_ovr = ovs(dem + '.par', targetres)	#Miro, 13-12-2018!
	#
	print('# Over-Sampling-Factors, are: lat_ovr=%f,lon_ovr=%f'%(lat_ovr, lon_ovr))

	# GAMMA:
	OFF_par='-'; sim_sar='tdx_sim_sar'; u='-'; v='-'; inc='inc'; psi='-'; pix='-'; ls_map='ls_map'; frame='-'; ls_mode=2; r_ovr='-'
	lt=os.path.join(workdir,'LT'); DEM_seg=os.path.join(workdir,'gc_dem')
 		
 	exists1=check_exists('LT',workdir)
 	exists2=check_exists('gc_dem',workdir)
 	#
 	if exists1==True and exists2==True: 			
 		#print("#> Skip - LUT-calculation. Files: %s and %s exists . " % (lt,DEM_seg))
 		print("#> Skip - LUT-calculation. Data exists. ")		

 	else:
		# build parameters dict:
		p = OrderedDict([('MLI_par', os.path.join(workdir, tdxm+par )), 	\
						('OFF_par', OFF_par ),								\
						('DEM_par', os.path.join(workdir, dem+par )),		\
						('DEM', os.path.join(workdir, dem )),				\
						('DEM_seg_par',os.path.join(workdir,DEM_seg+par)),	\
						('DEM_seg',  DEM_seg),		                        \
						('lookup_table',  lt),		                        \
						('lat_ovr', lat_ovr),								\
						('lon_ovr', lon_ovr),								\
						('sim_sar', os.path.join(workdir, sim_sar)),		\
						('u', u),											\
						('v', v),											\
						('inc', os.path.join(workdir, inc)),				\
						('psi', psi),										\
						('pix', pix),										\
						('ls_map', os.path.join(workdir, ls_map)),			\
						('frame', frame),									\
						('ls_mode', ls_mode),								\
						('r_ovr', r_ovr)									\
		])
		#
		#small workaround:
		cmd = parse_dict('gc_map',p)
		#
		# EXECUTE:	
		output = exec_shell_cmd(cmd)

		#log:
		lgf.write(cmd+'\n')		
	## end!

	lgf.close()

	return lt, DEM_seg	# return full path&name of LT and gc_dem!

	## 	gc_map TDX_110604_VV.mli.par - dem.par dem gc_dem.par gc_dem lt 5 5 - slp - inc - - ls_map - 2 - #


# ------------------------- GEOCODE : pyroSAR ------------------------------------------------------- #

def geocode(lt,gc_dem,map_width,mli_width,mli_lines,workdir,lgfn):

	'''	GEOCODE DATA, TRANSFORMATION DATA INTO RADAR GEOMETRY BY LUT, by using GAMMA-API from pyroSAR  '''
	
	from pyroSAR.gamma.api import diff
	import os
	from collections import OrderedDict	

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' geocode '+('#')*28+'\n')

	dem='dem'; rdc='.rdc'

	#GAMMA:
	data_out=os.path.join(workdir, dem+rdc)

 	exists=check_exists(dem+rdc,workdir)
 	#
 	if exists==True: 			
 		#print("#> Skip - Transformation to Radar geometry (SLR). Files: %s exists . " % (data_out))
 		print("#> Skip - Transformation to Radar geometry (SLR). Data exists. ")

 	else:
		#parameters = {	'gc_map': lt,			\
		#				'data_in': gc_dem,		\
		#				'width_in': map_width,	\
		#				'data_out': data_out,	\
		#				'width_out': map_width,	\
		#				'nlines_out': mli_lines	\
		#}
		#
		parameters = OrderedDict([('gc_map', lt),				\
									('data_in', gc_dem),		\
									('width_in', map_width),	\
									('data_out', data_out),		\
									('width_out', mli_width),	\
									('nlines_out', mli_lines)	\
					])
		#log:
		c=parse_dict('geocode',parameters)
		lgf.write(c+"\n")

		# EXECUTE, pyroSAR GAMMA-API:
		diff.geocode(**parameters)
	### end!

	lgf.close()
	
	return os.path.join(workdir,dem+rdc)	# return dem.rdc, with absolute path!	

	## geocode lt gc_dem $map_width dem.rdc $mli_width $mli_lines  #


# ------------------------- GEOCODE : pyroSAR ------------------------------------------------------- #

def create_offset(datasets,workdir,ml_rg,ml_az,lgfn):

	'''	Create and update ISP offset and interferogram parameter files 	  '''
							
	import os
	from pyroSAR.gamma.api import isp
	from collections import OrderedDict	

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' create_offset '+('#')*28+'\n')

	slc='_slc'; par='.par'
	tdx1=datasets[0];tdx2=datasets[1]
	orig=workdir[0]; proc=workdir[1]

	# GAMMA:
	OFF_par=os.path.join(proc,'OFF_par'); algorithm=1; iflg=0

 	exists=check_exists('OFF_par',proc)
 	#
 	if exists==True: 			
 		#print("#> Skip - Create Offset for interferogram parameter files. Files - %s exists . " % (OFF_par))
 		print("#> Skip - Create Offset for interferogram parameter files. Data exists. ")
 	else: 		
		# build, parameters dictionary:
		#parameters = {  'SLC1_par': os.path.join(orig, tdx1+slc+par),\
		#				'SLC2_par': os.path.join(orig, tdx2+slc+par),\
		#				'OFF_par': OFF_par,                             \
		#				'algorithm': algorithm,                         \
		#				'rlks': ml_rg,                                  \
		#				'azlks': ml_az,                                 \
		#				'iflg': iflg                                    \
		#}
		parameters = OrderedDict([('SLC1_par', os.path.join(orig, tdx1+slc+par)),\
									('SLC2_par', os.path.join(orig, tdx2+slc+par)),\
									('OFF_par', OFF_par),                           \
									('algorithm', algorithm),                       \
									('rlks', int(ml_rg)),                           \
									('azlks', int(ml_az)),                          \
									('iflg', iflg)                                  \
					])

		#log:
		c=parse_dict('create_offset',parameters)
		lgf.write(c+"\n")

		# EXECUTE, pyroSAR GAMMA-API:
		isp.create_offset(**parameters)
		#
		### end!

	lgf.close()
	
	return OFF_par    # return OFF_par, with absolute path!   

	## create_offset TDX_110604_VV_slc.par TSX_110604_VV_slc.par off_par 1 $ml_rg $ml_az 0 #


# ------------------------- Calculate Height of Ambiguity (HOA), and Vertical wave number (kz) ------------------------------------------------------- #

def calc_HOA_kz(Bperp,inc_angle,tdx1,path_orig):

	'''	Calculate Height of Ambiguity (HOA), and Vertical wave number (kz) '''

	### BE AWARE:
	#	tdx1_slc.par: incidence_angle = 41.4657 [deg], 
	#	tdx2_slc.par: incidence_angle = 41.4963 [deg], where delta=0.03059999999999974 [deg].
	#		
	#	tdx1_slc.par: center_range_slc = 671083.9687 [m],
	#	tdx2_slc.par: center_range_slc = 671032.1409 [m], where delta=51.82779999997001 [m] //> ПОЧТИ, колкото HOA!
	#
	### ИЛИ, разликата е много малка и е в рамките на грешката! Следователно - може да се използва коя и да е от стойностите, за сметките!
	#
	import math
	## from scipy:
	# speed of light in vacuum - speed_of_light = 299792458.0 [m.s^-1]
	speed_of_light = 299792458.0	# m/s

	slc="_slc"; par=".par"

	print ("#\n# parse MLI-par metadata, for 'center_range_slc', 'radar_frequency':")

	parameters=['center_range_slc','radar_frequency']
	val=[]
	val=parsePAR_pyroSAR(path_orig,tdx1+slc+par,parameters)
	#
	#for i in range(len(parameters)):
	#	print ('# # %s = %f (%.2f[m]) ) ',parameters[i],'=',val[i],	)
	##	
	
	R = val[0]; f = val[1] ; lamb = speed_of_light / f # 0.03106657595854922 [m]

	print ('# # %s = %f (%.3f[km]), %s = %f (%.3f[GHz]), %s = %f (%.3f[cm])' % (parameters[0],R,R/1000, parameters[1],f,f/1000000000, 'wavelength',lamb,lamb*100))

	inc_angle_r = inc_angle*math.pi/180.0
	
	HOA = (lamb*R*math.sin(inc_angle_r))/Bperp
	if HOA<0: HOA*=-1 #!

	kz = 2*math.pi/HOA
	#HOA = 2*math.pi/kz

	print("# # calculated HOA = %f, kz = %f . #"%(HOA,kz))

	return HOA,kz

# ------------------------- Calculate Baseline (perpendicular) : pyroSAR ------------------------------------------------------- #

def calc_Baseline(datasets,workdir,off_par,lgfn):

	''' Calculate Interferometric Baseline, perpendicular component to the Radar-LOS from Orbit State Vectors   '''


	# Сега, тъй, като дава разчлични резултати при изпълнението на:  
	# - [base_orbit] - тук, изхода през pyroSAR е един 
	#

	import os
	from pyroSAR.gamma.api import isp   
	from collections import OrderedDict	

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' calc_Baseline '+('#')*28+'\n')

	slc='_slc'; par='.par'
	orig=workdir[0]; proc=workdir[1]
	tdx1=datasets[0]; tdx2=datasets[1]
	##
	log=str(os.path.join(proc,"base_orbit.miro.log"))

	#GAMMA:
	baseline=os.path.join(proc,"tdx_baseline")
	##Bperp=os.path.join(proc,"tdx_Perp_Baseline")  # командата не връща файл!

 	exists=check_exists("tdx_baseline",proc)
 	#
 	if exists==True: 	 			
 		#print("#> Skip - calculation of Interferometric Baseline. Files: %s exists ." % (baseline))
 		print("#> Skip - calculation of Interferometric Baseline. Data exists.")
 	else: 
		# build, parameters dictionary:
		#parameters = {  'SLC1_par': os.path.join(orig, tdx1+slc+par),\
		#				'SLC2_par': os.path.join(orig, tdx2+slc+par),\
		#				'baseline': baseline
		#} 
		#
		# EXECUTE, pyroSAR GAMMA-API:
		#isp.base_orbit(**parameters)

		### NOW, Calculate exalctly the Perpendicular (& parallel) component:
		#parameters.clear()
		#parameters = {  'baseline': baseline,						\
		#				'SLC1_par': os.path.join(orig, tdx1+slc+par),\
		#				'OFF_par': off_par,						 	  \
		#				'time_rev': 1						
		#} 	
		#
		# EXECUTE, pyroSAR GAMMA-API:
		#isp.base_perp(proc,**parameters)	

		### NONE using pyroSAR, because I need the shell-output from [base_orbit].
			# workaround, to create the output file:
		
		cmd = "base_orbit "+os.path.join(orig, tdx1+slc+par)+" "+os.path.join(orig, tdx2+slc+par)+" "+baseline+" > "+log
		# EXECUTE:	
		os.system(cmd)	
		#

		#log:
		lgf.write(cmd+'\n')

	### end!
	Bperp=parse_log_extract_Bperp(log)	# връща - float!

	lgf.close()
	
	return baseline,Bperp    # return Baseline, with absolute path!   

	## base_orbit TDX_110604_VV_slc.par TSX_110604_VV_slc.par base_orb #
	### miro: 
	# base_perp /homes3/geoinf/ye68hah/GAMMA/Gamma2018_Training/L5_Phaseheight_TDX/proc/tdx_baseline /homes3/geoinf/ye68hah/GAMMA/Gamma2018_Training/L5_Phaseheight_TDX/orig/TSX_110604_VV_slc.par /homes3/geoinf/ye68hah/GAMMA/Gamma2018_Training/L5_Phaseheight_TDX/proc/OFF_par - /homes3/geoinf/ye68hah/GAMMA/Gamma2018_Training/L5_Phaseheight_TDX/proc 
	### base_perp.miro.log


# ------------------------- Calculate Raw INTERFEROGRAM : pyroSAR ------------------------------------------------------- #

def raw_Interferogram(datasets,workdir,ml_rg,ml_az,OFF_par,lgfn):

	''' Calculate Raw-Interferogram (unflattened, unwrapped, with movement along RLOS, with whole phase-noise sources)    '''
							
	import os
		### MISSING IN API!
	from collections import OrderedDict

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' raw_Interferogram '+('#')*28+'\n')

	slc='_slc'; par='.par'
	orig=workdir[0]; proc=workdir[1]
	tdx1=datasets[0]; tdx2=datasets[1]

	#GAMMA: #defaults! #except: rp1_flg AND rp2_flg /0: nearest approach (zero-Doppler) phase/
	loff=0; nlines='-'; spsf=1; azf=1; rp1f=0; rp2f=0; slc_1s='-'; slc_2r='-'; slc_1s_par='-'; slc_2r_par='-'; az_b=2.120
	interf=os.path.join(proc, "tdx_raw_Int")

 	exists=check_exists("tdx_raw_Int",proc)
 	#
 	if exists==True: 			
 		#print("#> Skip - calculation of the Raw Interferogram. Files: %s exists ." % (interf))
 		print("#> Skip - calculation of the Raw Interferogram. Data exists." )
 	else:  		
		# build parameters dict:
		p = OrderedDict([('SLC-1', os.path.join(orig, tdx1+slc )),       \
						('SLC-2R', os.path.join(orig, tdx2+slc )),       \
						('SLC1_par', os.path.join(orig, tdx1+slc+par )), \
						('SLC2R_par',os.path.join(orig, tdx2+slc+par )), \
						('OFF_par', OFF_par),                            \
						('interf', interf),                              \
						('rlks', int(ml_rg)),                            \
						('azlks', int(ml_az))                             \
						#('loff', loff),                                  \
						#('nlines', nlines),                              \
						#('sps_flg', spsf),                               \
						#('azf_flg', azf),                                \
						#('rp1_flg', rp1f),                               \
						#('rp2_flg', rp2f),                               \
						#('SLC-1s]', slc_1s),                             \
						#('SLC-2Rs', slc_2r),                             \
						#('SLC-1s_par', slc_1s_par),                      \
						#('SLC-2Rs_par', slc_2r_par),                     \
						#('az_beta', az_b)                                   
		])
		#
		#small workaround:
		cmd = parse_dict('SLC_intf',p)
		#
		#log:
		lgf.write(cmd+'\n')	

		# EXECUTE:  
		output = exec_shell_cmd(cmd)

	lgf.close()
	
	return interf # return raw Interferogram, full path&name !

	## SLC_intf TDX_110604_VV_slc TSX_110604_VV_slc TDX_110604_VV_slc.par TSX_110604_VV_slc.par off_par TSX_TDX_110604_int $ml_rg $ml_az #


# ------------------------- Simulate unwrapped interferometric phase using DEM height : pyroSAR --------------------------------------- #

def phase_simulation(datasets,workdir,OFF_par,dem_r,lgfn):

	''' Simulate unwrapped interferometric phase using DEM height and deformation rate using orbit state Vectors    '''

	import os
	from pyroSAR.gamma.api import diff
	from collections import OrderedDict	

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' phase_simulation '+('#')*28+'\n')

	slc='_slc'; par='.par'
	orig=workdir[0]; proc=workdir[1]
	tdx1=datasets[0]; tdx2=datasets[1]

	#GAMMA: 
	sim_unw=os.path.join(proc,'tdx_sim_unw_Phase'); int_mode=0  # single-pass mode (Tandem-X)
 	
 	exists=check_exists('tdx_sim_unw_Phase',proc)
 	#
 	if exists==True:			
 		#print("#> Skip - Phase Simulation using DEM. Files: %s exists ." % (sim_unw))
 		print("#> Skip - Phase Simulation using DEM. Data exists.")
 	else:  
		# build, parameters dictionary:
		#parameters = {  'SLC1_par': os.path.join(orig, tdx1+slc+par),\
		#				'SLC2R_par': os.path.join(orig, tdx2+slc+par),\
		#				'OFF_par': OFF_par,                            \
		#				'hgt' : dem_r,                                  \
		#				'sim_unw' : sim_unw,                             \
		#				'int_mode' : int_mode,                            \
		#				'ph_mode' : 0     # absolute phase (default)
		#}
		parameters = OrderedDict([ ('SLC1_par', os.path.join(orig, tdx1+slc+par)),\
									('SLC2R_par', os.path.join(orig, tdx2+slc+par)),\
									('OFF_par', OFF_par),                            \
									('hgt' , dem_r),                                  \
									('sim_unw' , sim_unw),                             \
									('SLC_ref_par' , '-'),								\
									('drm', '-'),										 \
									('delta_t','-'),									  \
									('int_mode' , int_mode),                               \
									('ph_mode' , 0)     			# absolute phase (default)
								]) 		 
		#
		# EXECUTE, pyroSAR GAMMA-API:
		diff.phase_sim_orb(**parameters)	

		#log:
		c=parse_dict('phase_sim_orb',parameters)
		lgf.write(c+"\n")

	### end!

	lgf.close()
	
	return sim_unw    # return Baseline, with absolute path!   

	## phase_sim_orb TDX_110604_VV_slc.par TSX_110604_VV_slc.par off_par dem.rdc TSX_TDX_110604_sim_unw_phase - - - 0 0 #


# ------------------------- Calculate Differential INTERFEROGRAM : pyroSAR --------------------------------------- #

def Diff_Interferogram(datasets,workdir,ml_rg,ml_az,OFF_par,sim_unw_phase,lgfn):

	''' Differential interferogram generation from co-registered SLCs and a simulated interferogram    '''

	import os
	from pyroSAR.gamma.api import diff
	from collections import OrderedDict	

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' Diff_Interferogram '+('#')*28+'\n')

	slc='_slc'; par='.par'
	orig=workdir[0]; proc=workdir[1]
	tdx1=datasets[0]; tdx2=datasets[1]

	#GAMMA: 
	diff_int=os.path.join(proc,"tdx_diff_Int")

 	exists=check_exists("tdx_diff_Int",proc)
 	#
 	if exists==True: 			
 		#print("#> Skip - Differential interferogram generation. Files: %s exists . " % (diff_int))
 		print("#> Skip - Differential interferogram generation. Data exists. ")
 	else:  
		# build, parameters dictionary:
		#parameters = {  'SLC_1': os.path.join(orig, tdx1+slc),         \
		#				'SLC_2R': os.path.join(orig, tdx2+slc),        \
		#				'SLC1_par': os.path.join(orig, tdx1+slc+par),  \
		#				'SLC2R_par' : os.path.join(orig, tdx2+slc+par),\
		#				'OFF_par' : OFF_par,                           \
		#				'sim_unw' : sim_unw_phase,                     \
		#				'diff_int' : diff_int,                         \
		#				'rlks' : ml_rg,                                \
		#				'azlks' : ml_az                                         
		#} 
		parameters = OrderedDict([	('SLC_1', os.path.join(orig, tdx1+slc)),         \
									('SLC_2R', os.path.join(orig, tdx2+slc)),        \
									('SLC1_par', os.path.join(orig, tdx1+slc+par)),  \
									('SLC2R_par' , os.path.join(orig, tdx2+slc+par)),\
									('OFF_par' , OFF_par),                           \
									('sim_unw' , sim_unw_phase),                     \
									('diff_int' , diff_int),                         \
									('rlks' , int(ml_rg)),                           \
									('azlks' , int(ml_az))                                         
								]) 

		#
		# EXECUTE, pyroSAR GAMMA-API:
		diff.SLC_diff_intf(**parameters)

		#log:
		c=parse_dict('SLC_diff_intf',parameters)
		lgf.write(c+"\n")

	### end!

	lgf.close()
	
	return diff_int   # return Baseline, with absolute path!  
   
	## SLC_diff_intf TDX_110604_VV_slc TSX_110604_VV_slc TDX_110604_VV_slc.par TSX_110604_VV_slc.par off_par TSX_TDX_110604_sim_unw_phase TSX_TDX_110604_diff $ml_rg $ml_az #


# ------------------------- Calculate Coherence : pyroSAR --------------------------------------- #

def Coherence_adp(tdx_Diff_Int,tdx1_cm,tdx2_cm,workdir,mli_width,lgfn):

	''' (Adaptive) Interferometric Coherence estimation, with consideration of phase slope and texture    '''

	import os
	from pyroSAR.gamma.api import lat
	from collections import OrderedDict	

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' Coherence_adp '+('#')*28+'\n')

	#GAMMA:
	coh=os.path.join(workdir,"tdx_Coh")

 	exists=check_exists("tdx_Coh",workdir)
 	#
 	if exists==True: 				
 		#print("#> Skip - Interferometric Coherence estimation. Files: %s exists ." % (coh))
 		print("#> Skip - Interferometric Coherence estimation. Data exists.")
 	else:  
		# build, parameters dictionary:
		#parameters = {  'interf': tdx_Diff_Int, \
		#				'pwr1': tdx1_cm,        \
		#				'pwr2': tdx2_cm,        \
		#				'slope' : '-',			\
		#				'texture' : '-',		\
		#				'cc_ad' : coh,          \
		#				'width' : mli_width,	\
		#				'box_min' : 3,		    \
		#				'box_max' : 9
		#}
		parameters = OrderedDict([	('interf', tdx_Diff_Int), 	\
									('pwr1', tdx1_cm),        	\
									('pwr2', tdx2_cm),        	\
									('slope' , '-'),			\
									('texture' , '-'),			\
									('cc_ad' , coh),          	\
									('width' , mli_width),		\
									('box_min' , 3),		    \
									('box_max' , 9)
								]) 		                    
		#
		# EXECUTE, pyroSAR GAMMA-API:
		lat.cc_ad(**parameters)	

		#log:
		c=parse_dict('cc_ad',parameters)
		lgf.write(c+"\n")

	### end!

	lgf.close()
	
	return coh   # return Baseline, with absolute path!  
   
	## cc_ad TSX_TDX_110604_diff TSX_110604_VV.cmli TDX_110604_VV.cmli - - TSX_TDX_110604_cc_ad $mli_width 3 9 #


# ------------------------- Adaptive Filtering : pyroSAR --------------------------------------- #

def adaptive_filtering(tdx_Diff_Int,workdir,mli_width,lgfn):

	''' (Adaptive) Filtering , of Interferogram, Coherence  '''
   
	import os
	from pyroSAR.gamma.api import isp
	from collections import OrderedDict	

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' adaptive_filtering '+('#')*28+'\n')

	#GAMMA:
	nfft=8 
	coh_filt=os.path.join(workdir,"tdx_Coh_filt"); int_filt=os.path.join(workdir,"tdx_Diff_Int_filt")

 	exists1=check_exists("tdx_Coh_filt",workdir)
 	exists2=check_exists("tdx_Diff_Int_filt",workdir)
 	#
 	if exists1==True and exists2==True: 			
 		#print("#> Skip - Adaptive filtering. Files: %s and: %s exists ." % (coh_filt,int_filt))
 		print("#> Skip - Adaptive filtering. Data exists ")
 	else: 
		# build, parameters dictionary:
		#parameters = {  'interf': tdx_Diff_Int, \
		#				'sm': int_filt,         \
		#				'cc': coh_filt,         \
		#				'width' : mli_width,    \
		#				'nfft' : nfft   #filtering FFT window size
		#} 
		parameters = OrderedDict([  ('interf', tdx_Diff_Int), 	\
									('sm', int_filt),         	\
									('cc', coh_filt),         	\
									('width' , mli_width),    	\
									('alpha','-'),				\
									('nfft', nfft),				\
									('cc_win','-')   			#filtering FFT window size	alpha='-', nfft='-', cc_win='-'
								])  		 
		#
		# EXECUTE, pyroSAR GAMMA-API:
		isp.adf(**parameters)	

		#log:
		c=parse_dict('adf',parameters)
		lgf.write(c+"\n")

	### end!

	lgf.close()
	
	return int_filt, coh_filt    # return Filtered Interferogram & Filtered Coherence, with absolute path!  
   
	## adf TSX_TDX_110604_diff TSX_TDX_110604_diff_sm TSX_TDX_110604_diff_cc $mli_width - 8 - #


# ------------------------- Phase-to-Height : pyroSAR --------------------------------------- #
def Diff_Height(datasets,workdir,off_par,dem_r,tdx_diffint_flt_ph,lgfn):

	''' Calculate delta-height, from differential interferometric phase using OSV (also, for baseline calculation)  '''

	import os
	from pyroSAR.gamma.api import diff
	from collections import OrderedDict	

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' geocode '+('#')*28+'\n')
	
	slc="_slc";par=".par"
	orig=workdir[0]; proc=workdir[1]
	tdx1_par=os.path.join(orig, datasets[0]+slc+par) 
	tdx2_par=os.path.join(orig, datasets[1]+slc+par)

	#GAMMA:
	dpdh="-"; int_mode=0 # interferometric acquisition mode: 0 = single-pass mode (Tandem-X) !!!
	dh=os.path.join(proc,"tdx_diffint_flt_ph.hgt")

 	exists=check_exists("tdx_diffint_flt_ph.hgt",proc)
 	#
 	if exists==True: 			
 		#print("#> Skip - Calculation of Differential Height. Files: %s exists ." % (dh))
 		print("#> Skip - Calculation of Differential Height. Data exists.")
 	else: 
		# build, parameters dictionary:
		#parameters = {  'SLC1_par': tdx1_par,       \
		#				'SLC2R_par': tdx2_par,      \
		#				'OFF_par': off_par,         \
		#				'hgt' : dem_r,              \
		#				'dp' : tdx_diffint_flt_ph,  \
		#				'dpdh' : dpdh,              \
		#				'dh' : dh,                  \
		#				'SLC_ref_par' : tdx1_par,   \
		#				'int_mode' : int_mode
		#}
		parameters = OrderedDict([  ('SLC1_par', tdx1_par),       \
									('SLC2R_par', tdx2_par),      \
									('OFF_par', off_par),         \
									('hgt' , dem_r),              \
									('dp' , tdx_diffint_flt_ph),  \
									('dpdh' , dpdh),              \
									('dh' , dh),                  \
									('SLC_ref_par' , tdx1_par),   \
									('int_mode' , int_mode)
								]) 

		#
		# EXECUTE, pyroSAR GAMMA-API:
		diff.dh_map_orb(**parameters)	

		#log:
		c=parse_dict('dh_map_orb',parameters)
		lgf.write(c+"\n")

	### end!

	lgf.close()
	
	return dh    # return Filtered Interferogram & Filtered Coherence, with absolute path!  
   
	## dh_map_orb TDX_110604_VV_slc.par TSX_110604_VV_slc.par off_par dem.rdc TSX_TDX_110604_phase - TSX_TDX_110604_phase.hgt TDX_110604_VV_slc.par 0 #


# ------------------------- GEOCODE, by LUT : pyroSAR --------------------------------------- #

def Geocode_back(dataset,workdir,lt,mli_width,map_width,lgfn):

	''' Geocoding of image data using lookup table values   '''

	import os
	from pyroSAR.gamma.api import diff
	from collections import OrderedDict	

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' geocode '+('#')*28+'\n')

	#GAMMA:
	o_data=os.path.join(workdir, dataset+"_geo")

 	exists=check_exists(dataset+"_geo",workdir)
 	#
 	if exists==True:  			
 		#print("#> Skip Geocoding. Files: %s exists ." % (o_data))
 		print("#> Skip Geocoding. Data exists.")
 	else: 
		# build, parameters dictionary:
		parameters = {  'data_in': dataset,         \
						'width_in': mli_width,      \
						'gc_map': lt,               \
						'data_out' : o_data,        \
						'width_out' : map_width,    \
						'dtype' : 0
		}  
		#dtype:
		#    input/output data type (enter - for default)
		#        * 0: FLOAT (default)
		#        * 1: FCOMPLEX
		#        * 2: SUN/BMP/TIFF 8 or 24-bit raster image
		#        * 3: UNSIGNED CHAR
		#        * 4: SHORT
		#        * 5: DOUBLE
		###
		# EXECUTE, pyroSAR GAMMA-API:
		diff.geocode_back(**parameters)		
	### end!

	lgf.close()
	
	return o_data    # return output Geocoded data, with absolute path!  

	## geocode_back TSX_TDX_110604_cc_ad $mli_width lt TSX_TDX_110604_geo_cc_ad $map_width - - 0
	## geocode_back TSX_TDX_110604_phase $mli_width lt TSX_TDX_110604_geo_phase $map_width - - 0
	## geocode_back TSX_TDX_110604_phase.hgt $mli_width lt TSX_TDX_110604_geo_hgt $map_width - - 0
	## geocode_back TSX_110604_VV.mli $mli_width lt TSX_110604_geo_mli $map_width - - 0


# ------------------------- Convert Complex-to-Real : pyroSAR --------------------------------------- #

def conv__complex_to_real(int_cplxph,workdir,mli_width,lgfn):

	''' Calculate real part, imaginary part, intensity, magnitude, or phase of FCOMPLEX data  '''

	import os
	from pyroSAR.gamma.api import disp
	from collections import OrderedDict	

	#log:
	lgf=open(lgfn,'a')
	lgf.write('\n'+('#')*17+' geocode '+('#')*28+'\n')

	#GAMMA:
	int_realph=os.path.join(workdir, "tdx_diffint_flt_ph")

 	exists=check_exists("tdx_diffint_flt_ph",workdir)
 	#
 	if exists==True: 			
 		#print("#> Skip calculation of the Real part from the Complex image. Files: %s exists ." % (int_realph))
 		print("#> Skip calculation of the Real part from the Complex image. Data exists.")
 	else: 
		# build, parameters dictionary:
		#parameters = {  'cpx': int_cplxph,  \
		#				'real': int_realph, \
		#				'width': mli_width, \
		#				'type' : 4
		#}  
		parameters = OrderedDict([ ('cpx', int_cplxph),  	\
									('real', int_realph), 	\
									('width', mli_width), 	\
									('type' , 4)
								])  		
		#output data type:
		#        * 0: real part
		#        * 1: imaginary part
		#        * 2: intensity (re\*re + im\*im)
		#        * 3: magnitude (sqrt(re\*re + im\*im))
		#        * 4: phase (atan2(im, re))
		###

		#log:
		c=parse_dict('dh_map_orb',parameters)
		lgf.write(c+"\n")

		# EXECUTE, pyroSAR GAMMA-API:
		disp.cpx_to_real(**parameters)
	### end!

	lgf.close()
	
	return int_realph    # return Real part of the filtered Interferometric phase, with absolute path!  
   
	## cpx_to_real TSX_TDX_110604_diff_sm TSX_TDX_110604_phase $mli_width 4 #


# ------------------------- Convert to GeoTIFF : pyroSAR --------------------------------------- #

def conv__to_GeoTIFF(dpar,dataset,workdir,lgfn):

	''' Convert geocoded data with DEM parameter file to GeoTIFF format '''

	import os
	from pyroSAR.gamma.api import disp

	#GAMMA:
	geotiff = dataset+".tif"

	# Quatsch! е да пропуснеш финалната стъпка!
	#
	# build, parameters dictionary:
	parameters = {  'DEM_par': dpar,    \
					'data': dataset,    \
					'type': 2,          \
					'GeoTIFF' : geotiff
	}  
	#data type:
		#     0: RASTER 8 or 24 bit uncompressed raster image, SUN (*.ras), BMP:(*.bmp), or TIFF: (*.tif)
		#     1: SHORT integer (2 bytes/value)
		#     2: FLOAT (4 bytes/value)
		#     3: SCOMPLEX (short complex, 4 bytes/value)
		#     4: FCOMPLEX (float complex, 8 bytes/value)
		#     5: BYTE
	###
	# EXECUTE, pyroSAR GAMMA-API:
	disp.data2geotiff(**parameters)	
	### end!

	return geotiff    # return GeoTIFF, with absolute path!  
   

## data2geotiff gc_dem.par TSX_TDX_110604_geo_cc_ad 2 TSX_TDX_110604_geo_cc_ad.tif
## data2geotiff gc_dem.par TSX_TDX_110604_geo_phase 2 TSX_TDX_110604_geo_phase.tif
## data2geotiff gc_dem.par TSX_TDX_110604_geo_hgt 2 TSX_TDX_110604_geo_hgt.tif 
## data2geotiff gc_dem.par TSX_110604_geo_mli 2 TSX_110604_geo_mli.tif



### ------------------------- END! Finished, on: 17-NOV-2018,T=20:28 --------------------------------------- #
