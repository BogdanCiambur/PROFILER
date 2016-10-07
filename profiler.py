import scipy.special as sp
from Tkinter import *
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
from lmfit import Parameters, minimize, fit_report, Parameter
from matplotlib import pyplot as plt
from matplotlib import rc
import construct_model
import bulgemag
from numpy import random, linspace, sqrt, genfromtxt, log10, vstack, array, linspace
from time import time as ti
import PSF_Convolution
import PSF_Vector
import mbh
import total_magnitude
import add_component
import datetime
import pylab
from pylab import rcParams
from indices import datamin, datamax, fitmin, fitmax
from scipy.interpolate import interp1d
import integrator
from format_logax import log_axis

rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def readEllipseOutput(data):
	data = open(data)
	datatab = genfromtxt(data, dtype={'names': ('sma', 'intens', 'intens_err', 'pix_var', 'rms', 'ellip', 'ellip_err', 'PA', 'PA_err', 'X0', 'X0_ERR', 'Y0', 'Y0_ERR', 'GRAD', 'GRAD_ERR', 'GRAD_R_ERR', 'RSMA', 'MAG', 'MAG_LERR', 'MAG_UERR', 'TFLUX_E', 'TFLUX_C', 'TMAG_E', 'TMAG_C', 'NPIX_E', 'NPIX_C', 'A3', 'A3_ERR', 'B3', 'B3_ERR', 'A4', 'A4_ERR', 'B4', 'B4_ERR'),
	'formats': ('f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, skiprows=0, missing_values=('INDEF'), filling_values=0.0)
	data.close()
	return datatab

def readEllipseOutput3col(data):
	data = open(data)
	datatab = genfromtxt(data, dtype={'names': ('sma', 'intens', 'ellip'),'formats': ('i4', 'f4', 'f4')}, missing_values=('INDEF'), filling_values=0.0)
	data.close()
	return datatab
	
def read2cols(data):
	data = open(data)
	datatab = genfromtxt(data, dtype={'names': ('sma', 'intens'),'formats': ('f4', 'f4')}, missing_values=('INDEF'), filling_values=0.0)
	data.close()
	return datatab

def readPSFFile(data):
	data = open(data)
	datatab = genfromtxt(data, dtype={'names': ('sma', 'intens'),'formats': ('f4', 'f4')}, missing_values=('INDEF'), filling_values=0.0)
	data.close()
	return datatab
	
def readtab_ellipsi(data):
	data = open(data)
	datatab2 = genfromtxt(data, dtype={'names': ('sma', 'intens', 'intens_err', 'pix_var', 'rms', 'ellip',
	 											'ellip_err', 'PA', 'PA_err', 'X0', 'X0_ERR', 'Y0', 'Y0_ERR',
	 											'GRAD', 'GRAD_ERR', 'GRAD_R_ERR', 'RSMA', 'MAG', 'MAG_LERR', 
	 											'MAG_UERR', 'TFLUX_E', 'TFLUX_C', 'TMAG_E', 'TMAG_C', 'NPIX_E', 
	 											'NPIX_C', 'A3', 'A3_ERR', 'B3', 'B3_ERR', 'A4', 'A4_ERR', 'B4', 
	 											'B4_ERR', 'NDATA', 'NFLAG', 'NITER', 'STOP', 'A_BIG', 'SAREA', 
	 											'AI2', 'AI2_ERR', 'BI2', 'BI2_ERR', 'AI3', 'AI3_ERR',
	 											'BI3', 'BI3_ERR', 'AI4', 'AI4_ERR', 'BI4', 'BI4_ERR'),
									  'formats': ('f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
									              'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
									              'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
									              'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
									              'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
									              'f4', 'f4')}, 
									  skiprows=0, missing_values=('INDEF'), filling_values=0.0)
	data.close()
	return datatab2
	
def readtab_ellipsi_ext(data):
	data = open(data)
	datatab2 = genfromtxt(data, dtype={'names': ('sma', 'intens', 'intens_err', 'pix_var', 'rms', 'ellip',
	 											'ellip_err', 'PA', 'PA_err', 'X0', 'X0_ERR', 'Y0', 'Y0_ERR',
	 											'GRAD', 'GRAD_ERR', 'GRAD_R_ERR', 'RSMA', 'MAG', 'MAG_LERR', 
	 											'MAG_UERR', 'TFLUX_E', 'TFLUX_C', 'TMAG_E', 'TMAG_C', 'NPIX_E', 
	 											'NPIX_C', 'A3', 'A3_ERR', 'B3', 'B3_ERR', 'A4', 'A4_ERR', 'B4', 
	 											'B4_ERR', 'NDATA', 'NFLAG', 'NITER', 'STOP', 'A_BIG', 'SAREA', 
	 											'AI2', 'AI2_ERR', 'BI2', 'BI2_ERR', 'AI3', 'AI3_ERR',
	 											'BI3', 'BI3_ERR', 'AI4', 'AI4_ERR', 'BI4', 'BI4_ERR', 'AI8',
	 											'AI8_ERR', 'BI8', 'BI8_ERR', 'AI10', 'AI10_ERR', 'BI10',
	 											'BI10_ERR', 'AI12', 'AI12_ERR', 'BI12', 'BI12_ERR'),
									  'formats': ('f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
									              'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
									              'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
									              'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
									              'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
									              'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
									              'f4', 'f4', 'f4', 'f4')}, 
									  skiprows=0, missing_values=('INDEF'), filling_values=0.0)
	data.close()
	return datatab2

global rowv
global ns,ne,ng,np,nf,nc,nsi,ntd,ns_a,ne_a,ng_a,np_a,nf_a,nc_a,nid_a,ntd_a #the '_a' stands for 'active' 
global sersic_entry,exp_entry,gauss_entry,psf_entry,ferrer_entry,corsic_entry,sinh_entry
sersic_entry=[0.0]*8 #structure == var_sersic,mue,re,n,colour,varymue,varyre,varyn
exp_entry=[0.0]*6 #structure == var_exp,h,mu0
id_entry=[0.0]*7 #structure == var_id,z,mu0,zero_coord
gauss_entry=[0.0]*8 #structure == var_gauss,r0r,mu0r,fwhmr
psf_entry=[0.0]*4 #structure == var_psf, mu0p
ferrer_entry=[0.0]*10 #structure == var_ferrer, r_out, mu_0f, alpha_F, beta_F
corsic_entry=[0.0]*14
td_entry=[0.0]*10
ns=0
ne=0
ng=0
np=0
nf=0
nc=0
nid=0
ntd=0

#----------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------widget definitions--------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
def fitcomm():
	global ns_a,ne_a,ng_a,np_a,nf_a,nc_a,nid_a,ntd_a
	ns_a=0
	ne_a=0
	ng_a=0
	np_a=0
	nf_a=0
	nc_a=0
	nid_a=0
	ntd_a=0

	fit_parameters = Parameters()

	if (ns == 1):

		if (sersic_entry[0].get() == 1): 
			add_component.add_s(ns-1,sersic_entry[1].get(),sersic_entry[2].get(),sersic_entry[3].get(),sersic_entry[5].get(),sersic_entry[6].get(),sersic_entry[7].get(),fit_parameters)
			ns_a=ns_a+1
	else: 
		for i in range(0,ns):
			if (sersic_entry[i,0].get() == 1): 
				add_component.add_s(i,sersic_entry[i,1].get(),sersic_entry[i,2].get(),sersic_entry[i,3].get(),sersic_entry[i,5].get(),sersic_entry[i,6].get(),sersic_entry[i,7].get(),fit_parameters)
				ns_a = ns_a + 1
				
	if (nc == 1):

		if (corsic_entry[0].get() == 1): 
			add_component.add_c(nc-1,corsic_entry[1].get(),corsic_entry[2].get(),corsic_entry[3].get(),corsic_entry[4].get(),corsic_entry[5].get(),corsic_entry[6].get(),corsic_entry[8].get(),corsic_entry[9].get(),corsic_entry[10].get(),corsic_entry[11].get(),corsic_entry[12].get(),corsic_entry[13].get(),fit_parameters)
			nc_a=nc_a+1
	else: 
		for i in range(0,nc):
			if (corsic_entry[i,0].get() == 1): 
				add_component.add_c(i,corsic_entry[i,1].get(),corsic_entry[i,2].get(),corsic_entry[i,3].get(),corsic_entry[i,4].get(),corsic_entry[i,5].get(),corsic_entry[i,6].get(),corsic_entry[i,8].get(),corsic_entry[i,9].get(),corsic_entry[i,10].get(),corsic_entry[i,11].get(),corsic_entry[i,12].get(),corsic_entry[i,13].get(),fit_parameters)
				nc_a = nc_a + 1

	if (ne == 1):

		if (exp_entry[0].get() == 1): 
			add_component.add_e(ne-1,exp_entry[1].get(),exp_entry[2].get(),exp_entry[4].get(), exp_entry[5].get(),fit_parameters)
			ne_a=ne_a+1
	else: 
		for i in range(0,ne):
			if (exp_entry[i,0].get() == 1): 
				add_component.add_e(i,exp_entry[i,1].get(),exp_entry[i,2].get(),exp_entry[i,4].get(), exp_entry[i,5].get(),fit_parameters)
				ne_a=ne_a+1
				
	if (nid == 1):

		if (id_entry[0].get() == 1): 
			add_component.add_id(nid-1,id_entry[1].get(),id_entry[2].get(),id_entry[4].get(),id_entry[5].get(),fit_parameters)
			nid_a=nid_a+1
	else: 
		for i in range(0,nid):
			if (id_entry[i,0].get() == 1): 
				add_component.add_id(i,id_entry[i,1].get(),id_entry[i,2].get(),id_entry[i,4].get(),id_entry[i,5].get(),fit_parameters)
				nid_a=nid_a+1
	
	if (ng == 1):

		if (gauss_entry[0].get() == 1): 
			add_component.add_g(ng-1,gauss_entry[1].get(),gauss_entry[2].get(),gauss_entry[3].get(),gauss_entry[5].get(),gauss_entry[6].get(),gauss_entry[7].get(),fit_parameters)
			ng_a=ng_a+1
	else: 
		for i in range(0,ng):
			if (gauss_entry[i,0].get() == 1): 
				add_component.add_g(i,gauss_entry[i,1].get(),gauss_entry[i,2].get(),gauss_entry[i,3].get(),gauss_entry[i,5].get(),gauss_entry[i,6].get(),gauss_entry[i,7].get(),fit_parameters)
				ng_a=ng_a+1
			
	if (np == 1):

		if (psf_entry[0].get() == 1): 
			add_component.add_p(np-1,psf_entry[1].get(),psf_entry[3].get(),fit_parameters)
			np_a=np_a+1
	else: 
		for i in range(0,np):
			if (psf_entry[i,0].get() == 1): 
				add_component.add_p(i,psf_entry[i,1].get(),psf_entry[i,3].get(),fit_parameters)
				np_a=np_a+1
			
	if (nf == 1):

		if (ferrer_entry[0].get() == 1): 
			add_component.add_f(nf-1,ferrer_entry[1].get(),ferrer_entry[2].get(),ferrer_entry[3].get(),ferrer_entry[4].get(),ferrer_entry[6].get(),ferrer_entry[7].get(),ferrer_entry[8].get(),ferrer_entry[9].get(),fit_parameters)
			nf_a=nf_a+1
	else: 
		for i in range(0,nf):
			if (ferrer_entry[i,0].get() == 1): 
				add_component.add_f(i,ferrer_entry[i,1].get(),ferrer_entry[i,2].get(),ferrer_entry[i,3].get(),ferrer_entry[i,4].get(),ferrer_entry[i,6].get(),ferrer_entry[i,7].get(),ferrer_entry[i,8].get(),ferrer_entry[i,9].get(),fit_parameters)
				nf_a=nf_a+1
				
	if (ntd == 1):

		if (td_entry[0].get() == 1): 
			add_component.add_td(ntd-1,td_entry[1].get(),td_entry[2].get(),td_entry[3].get(),td_entry[4].get(),td_entry[6].get(),td_entry[7].get(),td_entry[8].get(),td_entry[9].get(),fit_parameters)
			ntd_a=ntd_a+1
	else: 
		for i in range(0,ntd):
			if (td_entry[i,0].get() == 1): 
				add_component.add_td(i,td_entry[i,1].get(),td_entry[i,2].get(),td_entry[i,3].get(),td_entry[i,4].get(),td_entry[i,6].get(),td_entry[i,7].get(),td_entry[i,8].get(),td_entry[i,9].get(),fit_parameters)
				ntd_a=ntd_a+1

	ellipticity=float(ellipentry.get())
	if (axistype.get() == 'equ'): ellipticity = 0.
	pxscale=float(PSC.get())
	sky_m = float(Sky.get())
	mzp=float(mag0p.get())
	xmin= float(rlo.get())
	xmax= float(rhi.get())
	fit_mini = (float(rflo.get()))
	fit_maxi = (float(rfhi.get()))
	
	itt = input_table_type.get()
	
	if itt=='ell': table=readEllipseOutput(IP.get())
	if itt=='iso': table=readtab_ellipsi(IP.get())
	if itt=='ir': table=read2cols(IP.get())
	
	sizearr=len(table['sma'])
	x=[0.]*sizearr
	xt=float(PSC.get())*table['sma']
	if (axistype.get()=='equ'):
		if (itt != 'ir'):
			print 'Data R column should be the semi-major axis (R_maj). Converting to R_eq...' 
			ell=table['ellip']
			xt[1:len(x)] = xt[1:len(x)]*sqrt(1.-ell[1:len(x)])
		else: print 'Cannot generate equivalent axis without ellipticity column. Assuming the data R column is R_eq.'
	x=x+xt
	data=-2.5*log10(table['intens']/(pxscale**2.)) + mzp
	x_original = array(xt)
	data_original = array(data)
	sampling = 'linear'
	
	
	xmino = datamin(xmin,x_original)
	xmaxo = datamax(xmax,x_original)
	fit_mino = fitmin(fit_mini,x_original)
	fit_maxo = fitmax(fit_maxi,x_original)
	if (fit_mino < xmino): fit_mino=xmino
	if (fit_maxo > xmaxo): fit_maxo=xmaxo
	if (fit_mino < 0): fit_mino = 0
	x_original = x_original[xmino:xmaxo]
	data_original = data_original[xmino:xmaxo]
	
	#print (x[4]-x[3]), (x[3]-x[2]), (x[2]-x[1]), (x[1]-x[0])
	if ((x[2]-x[1]) != (x[1]-x[0]) and (x[2]-x[1]) != (x[3]-x[2])):
		print 'The radial sampling step is not linear (uniform). Interpolating...'
		f = interp1d(x,data,kind='slinear')
		xn = linspace(0.,max(x),(max(x)-min(x))/((x[2]-x[1])))		
		yn = f(xn)
		x = array(xn)
		data = array(yn)
		sampling = 'nonlinear'
		
	print 'The sampling of the data is :', sampling
 
	xmin = datamin(xmin,x)
	if (xmin != 0.): print 'WARNING: Data does not start at R=0. Convolution might not be correct!'	
	xmax = datamax(xmax,x)
	fit_min = fitmin(fit_mini,x)
	fit_max = fitmax(fit_maxi,x)
	
	if (fit_min < xmin): fit_min=xmin
	if (fit_max > xmax): fit_max=xmax
	if (fit_min < 0): fit_min = 0
	
	x=x[xmin:xmax]
	if (x[len(x)-1] != x_original[len(x_original)-1]): x[len(x)-1] = x_original[len(x_original)-1]
	if (x[0] != x_original[0]): x[0] = x_original[0]
	nxt = 2
	extsize = nxt*len(x)
	xx = array(linspace(x[0],1.*nxt*x[len(x)-1],extsize))
	data=data[xmin:xmax]
	PSF_fft_ext = 0.
	PSF_w = 0.
	elle = [0.]*extsize

	if (itt=='ell' or itt=='iso'):
		ell=table['ellip']
		pan=table['PA']
		if itt=='iso': b4=table['BI4']
		if itt=='ell': b4=table['B4']
		grad=table['GRAD']
		sma=table['sma']
		b4=-1.*b4/grad/sma
		if (axistype.get() == 'equ'):
			b4=b4*sqrt(1. - ell)	
		b4=b4[xmino:xmaxo]
		ell=ell[xmino:xmaxo]
		pan=pan[xmino:xmaxo]
		elle[0:len(ell)] = ell[0:len(ell)]
	else:
		ell=[0.]*len(data)
		b4 =[0.]*len(data)
		pan=[0.]*len(data)
	
	n_all=[ns,ne,ng,np,nf,nc,nid,ntd]
	
	psf_type=ps.get()
	PSF_ext=[1e-99]*extsize
	PSF_o = 0.
	if psf_type == 'moffat':
		psf_fwhm=float(PMFV.get())*pxscale
		psf_beta=float(PMBV.get())
		psf_alpha=psf_fwhm/(2.*sqrt((2.**(1./psf_beta))-1.))
		PSF = ['moffat',psf_fwhm,psf_alpha,psf_beta]
		PSF_o = PSF
		PSF_ext = PSF
	elif psf_type == 'gaussian':
		psf_fwhm=float(PGFV.get())*pxscale
		PSF = ['gaussian',psf_fwhm]
		PSF_o = PSF
		PSF_ext = PSF
	#:::::::::::::::::handle numerical PSF here:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	elif psf_type == 'numerical':
		table=readPSFFile(PNV.get())
		sizepsf=len(table['sma'])
		xp=[0.]*sizearr
		xt=float(pxscale)*table['sma']
		yp=log10(table['intens']/pxscale**2)
		yp = 10.**(yp)
		
		interp_psf = interp1d(xt,yp,kind='slinear')
		PSF=interp_psf(x)					#used for convolution
		PSF_o = interp_psf(x_original)		#used to build nuclear component and compute residual on original data vector
		#PSF_ext=[1e-99]*extsize				#used to build nuclear component in the extended, nice-looking profile
		#PSF_ext = interp_psf(x)
		PSF_ext[0:len(x)] = PSF[0:len(x)]
		PSF = PSF/max(PSF)
		PSF_ext = PSF_ext/max(PSF_ext)
		PSF_o = PSF_o/max(PSF_o)
	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	out = minimize(construct_model.residual, fit_parameters, args=(x,n_all,sersic_entry,exp_entry,gauss_entry,psf_entry,ferrer_entry,corsic_entry,id_entry,td_entry,PSF, fit_min,fit_max,ellipticity), kws={'data':data})
	freepar = fit_max-fit_min-out.nfree
	modelf = open('modelf_'+IP.get(),'w')
	modelf.write ('#SMA'+ ' '+ 'INTENS'+' '+'INTEND_data'+'\n')
	sfact = 2.5*log10(pxscale**2)
	if sampling=='linear':
		fit = construct_model.residual(fit_parameters, x, n_all, sersic_entry,exp_entry,gauss_entry,psf_entry,ferrer_entry,corsic_entry,id_entry,td_entry,PSF,fit_min,fit_max,ellipticity)
		for i in range(0,len(x)): modelf.write (str(round(x[i]/pxscale,2))+ '  '+str(10.**((mzp - fit[i] + sfact)/2.5))+ '  '+str(10.**((mzp - data_original[i] + sfact)/2.5))+'\n')
		modelf.close()
		residu=data-fit
		resexx=residu[fit_min:fit_max]
		deltarms=sqrt(sum((resexx)**2)/(len(resexx)-freepar+1))
	elif sampling=='nonlinear':
		fit = construct_model.residual(fit_parameters, x, n_all, sersic_entry,exp_entry,gauss_entry,psf_entry,ferrer_entry,corsic_entry,id_entry,td_entry,PSF,fit_min,fit_max,ellipticity)
		if max(x) < max(x_original): x[len(x)-1] = max(x_original)
		modelint = interp1d(x,fit,kind='slinear')
		fit_original = modelint(x_original)
		for i in range(0,len(x_original)): 
			modelf.write (str(round(x_original[i]/pxscale,2))+ '  '+str(10.**((mzp - fit_original[i] + sfact)/2.5)) + '  '+str(10.**((mzp - data_original[i] + sfact)/2.5))+ '\n')
		residu = data_original - fit_original
        resexx=residu[fit_mino:fit_maxo]
        deltarms=sqrt(sum((resexx)**2)/(len(resexx)-freepar+1))
        modelf.close()
	fit_nice = construct_model.residual(fit_parameters, xx, n_all, sersic_entry,exp_entry,gauss_entry,psf_entry,ferrer_entry,corsic_entry,id_entry,td_entry,PSF_ext,fit_min,fit_max,ellipticity)
	print fit_report(fit_parameters, show_correl=False)
	#----------------------------------------------------------
	logf = open('logfile_'+IP.get(),'a+')
	logf.write ('-------------------------------- Fit results: ------------------------------------\n')
	logf.write ('----------------------------------------------------------------------------------\n')
	tstamp=str(datetime.datetime.now())
	logf.write ('Time stamp of fit: ' + tstamp + '\n')
	logf.write ('Input file: '+IP.get()+'\n')
	if (psf_type != 'numerical'):
		logf.write ('PSF: ' + str(PSF) + '\n')
	else: logf.write ('PSF: ' + PNV.get()+'\n')
	logf.write ('Central ellipticity: '+str(ellipticity)+'\n')
	logf.write ('Zero-point magnitude: '+str(mzp)+'\n')
	logf.write ('Pixel scale ["/px]: '+str(pxscale)+'\n')
	logf.write ('Radial axis: '+str(axistype.get())+'\n')
	logf.write ('\n')
	logf.write ('Components [generated, active]: \n')
	logf.write ('Sersic: '+str(ns)+' created, '+str(ns_a)+' active \n')
	logf.write ('Exponential: '+str(ne)+' created, '+str(ne_a)+' active \n')
	logf.write ('Inclined disc: '+str(nid)+' created, '+str(nid_a)+' active \n')
	logf.write ('Truncated disc: '+str(ntd)+' created, '+str(ntd_a)+' active \n')
	logf.write ('Gaussian: '+str(ng)+' created, '+str(ng_a)+' active \n')
	logf.write ('PSF: '+str(np)+' created, '+str(np_a)+' active \n')
	logf.write ('Ferrer: '+str(nf)+' created, '+str(nf_a)+' active \n')
	logf.write ('\n')
	logf.write ('Free parameters = '+str(freepar)+'\n')
	logf.write ('Delta rms = '+str(deltarms))
	logf.write ('\n')
	logf.write ('-----------------------------Detailed fit report:-------------------------------- \n')
	logf.write ('----------------------------------------------------------------------------------\n')
	logf.write ('\n')
	
	if (itt!='ir'): n_panels=5+2*(pan_el.get()+pan_pa.get()+pan_b4.get())
	else: n_panels=5
	rcParams['figure.figsize'] = 7.5, 0.78*n_panels+1.08
	fig=plt.figure()
	ax1 = plt.subplot2grid((n_panels,3), (0, 0), colspan=3, rowspan=4)
	#----------------------------------------------------------	
	if (fit_min > xmin): ax1.plot(x_original[xmino:fit_mino],data_original[xmino:fit_mino], marker='o', color='white', markersize=7., linestyle='None')
	if (fit_max < xmax): ax1.plot(x_original[fit_maxo:xmaxo],data_original[fit_maxo:xmaxo], marker='o', color='white', markersize=7., linestyle='None')
	ax1.plot(x_original[fit_mino:fit_maxo],data_original[fit_mino:fit_maxo], marker='o', color='red', markersize=7., linestyle='None', label='Data', markeredgecolor='darkred')
	plt.text(0.65,0.95, TitlLabEn.get(), fontsize=19, horizontalalignment='center', verticalalignment='top', transform = ax1.transAxes)
	#----------------------------------------------------------	
	
	if (ns == 1):
		type = 'sersic'
		i=0
		comp=construct_model.component(fit_parameters, xx, ns, i, type, sersic_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
		plt.plot(xx,comp, color=sersic_entry[4].get(), linewidth=1.6)
	
		if (sersic_entry[0].get()==1):
			#vals = fit_parameters.valuesdict()
			#mu_et = "{:10.2f}".format(vals['mu_e0'])
			#r_et =  "{:10.2f}".format(vals['r_e0'])
			#n_Sert = "{:10.2f}".format(vals['n_Ser0'])
			#plt.text(0.75,0.74, r'$\;n_{ }\:\;\,$ = '+str(n_Sert), fontsize=20, horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)
			#plt.text(0.75,0.9, r'$R_{\rm e}\,$ = '+str(r_et), fontsize=20, horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)
			#plt.text(0.75,0.82, r'$\mu_{\rm e}\:$ = '+str(mu_et), fontsize=20, horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)
			magtot=integrator.magnitude(xx,comp,mzp)
			logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
			logf.write ('------------------------------------------------\n')
			logf.write ('\n')
	elif (ns > 1):
		type = 'sersic'
		for i in range(0,ns):
			comp=construct_model.component(fit_parameters, xx, ns, i, type, sersic_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
			plt.plot(xx,comp, color=sersic_entry[i,4].get(), linewidth=1.6)
		
			if (sersic_entry[i,0].get()==1):
				#if (sersic_entry[0,0].get()==1):
					#vals = fit_parameters.valuesdict()
					#mu_et = "{:10.2f}".format(vals['mu_e0'])
					#r_et =  "{:10.2f}".format(vals['r_e0'])
					#n_Sert = "{:10.2f}".format(vals['n_Ser0'])
					#plt.text(0.75,0.74, r'$\;n_{ }\:\;\,$ = '+str(n_Sert), fontsize=20, horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)
					#plt.text(0.75,0.9, r'$R_{\rm e}\,$ = '+str(r_et), fontsize=20, horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)
					#plt.text(0.75,0.82, r'$\mu_{\rm e}\:$ = '+str(mu_et), fontsize=20, horizontalalignment='left', verticalalignment='top', transform = ax1.transAxes)
				magtot=integrator.magnitude(xx,comp,mzp)
				logf.write ('+ Integrated total magnitude mag = '+str(magtot)+'\n')
				logf.write ('------------------------------------------------\n')
				logf.write ('\n')
	if (nc == 1):
		type = 'corsic'
		i=0
		comp=construct_model.component(fit_parameters, xx, nc, i, type, corsic_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
		plt.plot(xx,comp, color=corsic_entry[7].get(), linewidth=1.6)
		if (corsic_entry[0].get()==1):
			magtot=total_magnitude.compute(xx,comp,elle)
			logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
			logf.write ('------------------------------------------------\n')
			logf.write ('\n')
	elif (nc > 1):
		type = 'corsic'
		for i in range(0,ns):
			comp=construct_model.component(fit_parameters, xx, nc, i, type, corsic_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
			plt.plot(xx,comp, color=corsic_entry[i,7].get(), linewidth=1.6)
			if (corsic_entry[i,0].get()==1):
				magtot=total_magnitude.compute(xx,comp,elle)
				logf.write ('+ Integrated total magnitude mag = '+str(magtot)+'\n')
				logf.write ('------------------------------------------------\n')
				logf.write ('\n')
	if (ne == 1):
		type = 'exponential'
		i=0
		comp=construct_model.component(fit_parameters, xx, ne, i, type, exp_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
		plt.plot(xx,comp, color=exp_entry[3].get(), linewidth=1.6)
		vals = fit_parameters.valuesdict()
		
		if (exp_entry[0].get()==1):
			magtot=integrator.magnitude(xx,comp,mzp)
			logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
			logf.write ('------------------------------------------------\n')
			logf.write ('\n')
	elif (ne > 1):
		type = 'exponential'
		for i in range(0,ne):
			comp=construct_model.component(fit_parameters, xx, ne, i, type, exp_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
			plt.plot(xx,comp, color=exp_entry[i,3].get(), linewidth=1.6)
			if (exp_entry[i,0].get()==1):
				magtot=integrator.magnitude(xx,comp,mzp)
				logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
				logf.write ('------------------------------------------------\n')
				logf.write ('\n')
	if (nid == 1):
		type = 'id'
		i=0
		comp=construct_model.component(fit_parameters, xx, nid, i, type, id_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
		plt.plot(xx,comp, color=id_entry[3].get(), linewidth=1.6)
		vals = fit_parameters.valuesdict()
		

	elif (nid > 1):
		type = 'id'
		for i in range(0,nid):
			comp=construct_model.component(fit_parameters, xx, nid, i, type, id_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
			plt.plot(xx,comp, color=id_entry[i,3].get(), linewidth=1.6)
			
	if (ng == 1):
		type = 'gaussian'
		i=0
		comp=construct_model.component(fit_parameters, xx, ng, i, type, gauss_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
		plt.plot(xx,comp, color=gauss_entry[4].get(), linewidth=1.6)
		if (gauss_entry[0].get()==1):
		       	magtot=total_magnitude.compute(xx,comp,elle)
		       	logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
			logf.write ('------------------------------------------------\n')
			logf.write ('\n')
	elif (ng > 1):
		type = 'gaussian'
		for i in range(0,ng):
			comp=construct_model.component(fit_parameters, xx, ng, i, type, gauss_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
			plt.plot(xx,comp, color=gauss_entry[i,4].get(), linewidth=1.6)
			if (gauss_entry[i,0].get()==1):
				magtot=total_magnitude.compute(xx,comp,elle)
				logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
				logf.write ('------------------------------------------------\n')
				logf.write ('\n')
	if (np == 1):
		type = 'psf'
		i=0
		comp=construct_model.component(fit_parameters, xx, np, i, type, psf_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
		plt.plot(xx,comp, color=psf_entry[2].get(), linewidth=1.6)
		if (psf_entry[0].get()==1):
		       	magtot=total_magnitude.compute(xx,comp,elle)
		       	logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
		       	logf.write ('------------------------------------------------\n')
		       	logf.write ('\n')
	elif (np > 1):
		type = 'psf'
		for i in range(0,np):
			comp=construct_model.component(fit_parameters, xx, np, i, type, psf_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
			plt.plot(xx,comp, color=psf_entry[i,2].get(), linewidth=1.6)
			if (psf_entry[i,0].get()==1):
				magtot=total_magnitude.compute(xx,comp,elle)
				logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
				logf.write ('------------------------------------------------\n')
				logf.write ('\n')
	if (nf == 1):
		type = 'ferrer'
		i=0
		comp=construct_model.component(fit_parameters, xx, nf, i, type, ferrer_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
		plt.plot(xx,comp, color=ferrer_entry[5].get(), linewidth=1.6)
		if (ferrer_entry[0].get()==1):
		       	magtot=total_magnitude.compute(xx,comp,elle)
		       	logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
		       	logf.write ('------------------------------------------------\n')
		       	logf.write ('\n')
	elif (nf > 1):
		type = 'ferrer'
		for i in range(0,nf):
			comp=construct_model.component(fit_parameters, xx, nf, i, type, ferrer_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
			plt.plot(xx,comp, color=ferrer_entry[i,5].get(), linewidth=1.6)
			if (ferrer_entry[i,0].get()==1):
				magtot=total_magnitude.compute(xx,comp,elle)
				logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
				logf.write ('------------------------------------------------\n')
				logf.write ('\n')
	if (ntd == 1):
		type = 'td'
		i=0
		comp=construct_model.component(fit_parameters, xx, ntd, i, type, td_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
		plt.plot(xx,comp, color=td_entry[5].get(), linewidth=1.6)
		if (td_entry[0].get()==1):
		       	magtot=total_magnitude.compute(xx,comp,elle)
		       	logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
		       	logf.write ('------------------------------------------------\n')
		       	logf.write ('\n')
	elif (ntd > 1):
		type = 'td'
		for i in range(0,ntd):
			comp=construct_model.component(fit_parameters, xx, ntd, i, type, td_entry, PSF_ext, PSF_fft_ext, PSF_w, logf, elle,ellipticity)
			plt.plot(xx,comp, color=td_entry[i,5].get(), linewidth=1.6)
			if (td_entry[i,0].get()==1):
				magtot=total_magnitude.compute(xx,comp,elle)
				logf.write ('+ Integrated total magnitude = '+str(magtot)+'\n')
				logf.write ('------------------------------------------------\n')
				logf.write ('\n')
	
	logf.write ('Raw fit report (all parameters[numbered from 0-n] and their correlations): \n')
	logf.write ('------------------------------------------------------------------------------------\n')
	logf.write (fit_report(fit_parameters))
	logf.write ('\n----------------------------------------------------------------------------------\n')
	logf.write ('-----------------------------------------END----------------------------------------\n')
	logf.write ('\n')
	logf.close()
	
	ax1.plot(xx,fit_nice, color='black', label='Model', linewidth=1.6)
	plt.gca().invert_yaxis()
	plt.ylabel('$\mu$ [mag arcsec$^{-2}$]', fontsize=20, labelpad=28)
	plt.tick_params(axis='x', labelbottom='off')
	plt.tick_params(axis='y', labelsize=19)
	plt.minorticks_on()
	plt.tick_params('both', length=5, width=1.2, which='major')
	plt.tick_params('both', length=3, width=1, which='minor')
	if (plottype.get() == 'log'): plt.xscale('log')
	
	y_max=max(data)+0.9
	y_min=min(data)-0.9
	
	if (y_max - int(y_max) < 0.2): y_max = y_max+0.3
	if (y_min - int(y_min) < 0.2): y_min = y_min-0.3
	
	if (1 + int(y_max) - y_max < 0.2): y_max = y_max-0.3
	if (1 + int(y_min) - y_min < 0.2): y_min = y_min+0.3

	plt.ylim([y_max,y_min])
	
	if ((plottype.get() == 'log')):
		if min(x)==0.: xminplt=0.75*x[1]
		else: xminplt = 0.75*min(x) 
		plt.xlim(xminplt,1.25*max(x))
		xsky=linspace(min(x[1:len(x)])-0.5*min(x[1:len(x)]),max(x)+0.3,len(x))
		skyline=[sky_m]*len(xsky)
		plt.plot(xsky,skyline, linewidth=1., color='black', ls = ':')
	else: 
		plt.xlim(min(x)-0.25*max(x),1.09*max(x))
		xsky=linspace(min(x)-10.05,1.5*max(x),len(x))
		skyline=[sky_m]*len(xsky)
		plt.plot(xsky,skyline, linewidth=1., color='black', ls = ':')
	
	ax2 = plt.subplot2grid((n_panels,3), (4, 0), colspan=3, sharex=ax1)
	resex=residu
	ax2.plot(x_original,resex, marker='o', color='black', markersize=5.5, linestyle='None', label='Delta')
	ax2.plot(x_original[xmino:fit_mino],resex[xmino:fit_mino], marker='o', color='white', markersize=5.5, linestyle='None')
	ax2.plot(x_original[fit_maxo:xmaxo],resex[fit_maxo:xmaxo], marker='o', color='white', markersize=5.5, linestyle='None')
	ax2.plot(x_original[fit_mino:fit_maxo],resex[fit_mino:fit_maxo], marker='o', color='black', markersize=5.5, linestyle='None')
	plt.ylabel('$\Delta\mu$', fontsize=20, labelpad=7)
	if (n_panels-4-1 > 0): plt.tick_params(axis='x', labelbottom='off')
	plt.tick_params(axis='y', labelsize=19)
	plt.yticks([-0.1,0,0.1])
	if (plottype.get() == 'log'): plt.xscale('log')
	plt.ylim(0.18,-0.18)
	if ((plottype.get() == 'log')): 
		xz=linspace(1e-11,max(x)+0.5*max(x),len(x))
		zer=[0.0]*len(xz)
		plt.plot(xz,zer, linewidth=1.5, color='red')
	else: 
		plt.xlim(min(x)-0.05*max(x),1.09*max(x))
		xz=linspace(min(x)-0.5*max(x),max(x)+0.5*max(x),len(x))
		zer=[0.0]*len(xz)
		plt.plot(xz,zer, linewidth=1.5, color='red')
	deltarms="{:10.4f}".format(deltarms)	
	plt.text(0.05,0.0, r'$\Delta_{rms} = $' + deltarms, fontsize=19, horizontalalignment='left', verticalalignment='bottom', transform = ax2.transAxes)
	plt.minorticks_on()
	plt.tick_params('both', length=5, width=1.2, which='major')
	plt.tick_params('both', length=3, width=1, which='minor')
	
	npan_extra = 5
	
	if pan_el.get()>0 and itt!='ir':
		ax3 = plt.subplot2grid((n_panels,3), (npan_extra, 0), colspan=3, rowspan=2, sharex = ax1)
		ax3.plot(x_original,ell, marker='o', color='white', markersize=5.5, linestyle='None')
		ax3.plot(x_original[fit_mino:fit_maxo],ell[fit_mino:fit_maxo], marker='o', color='black', markersize=5.5, linestyle='None')
		plt.ylabel('$\epsilon$', fontsize=21, labelpad=20)
		if (n_panels-npan_extra-2 > 0): plt.tick_params(axis='x', labelbottom='off')
		plt.tick_params(axis='y', labelsize=19)
		plt.ylim([-0.04,0.92])
		plt.yticks([0.,0.2, 0.4, 0.6, 0.8])
		plt.minorticks_on()
		plt.tick_params('both', length=5, width=1.2, which='major')
		plt.tick_params('both', length=3, width=1, which='minor')
		npan_extra = npan_extra+2
		
	if pan_pa.get()>0 and itt!='ir':
		ax4 = plt.subplot2grid((n_panels,3), (npan_extra, 0), colspan=3, rowspan=2, sharex=ax1)
		ax4.plot(x_original,pan, marker='o', color='white', markersize=6, linestyle='None', label='Delta')
		ax4.plot(x_original[fit_mino:fit_maxo],pan[fit_mino:fit_maxo], marker='o', color='black', markersize=6, linestyle='None', label='Delta')
		plt.ylabel('$PA$ [deg]', fontsize=20)
		if (n_panels-npan_extra-2 > 0): plt.tick_params(axis='x', labelbottom='off')
		plt.tick_params(axis='y', labelsize=19)
		ax4.set_ylim([-90.,90.])
		ax4.set_yticks([-50,0,50])
		plt.minorticks_on()
		plt.tick_params('both', length=5, width=1.2, which='major')
		plt.tick_params('both', length=3, width=1, which='minor')
		npan_extra=npan_extra+2
	
	if pan_b4.get()>0 and itt!='ir':
		ax5 = plt.subplot2grid((n_panels,3), (npan_extra, 0), colspan=3, rowspan=2, sharex=ax1)
		ax5.plot(x_original,b4, marker='o', color='white', markersize=5.5, linestyle='None', label='Delta')
		ax5.plot(x_original[fit_mino:fit_maxo],b4[fit_mino:fit_maxo], marker='o', color='black', markersize=5.5, linestyle='None', label='Delta')
		zer=[0.0]*len(xz)
		ax5.plot(xz,zer, linewidth=1.5, color='red')
		ax5.set_ylabel('$B_{4}$', fontsize=20, labelpad=-10)
		ax5.set_ylim([-0.14,0.14])
		ax5.set_yticks([-0.1,-0.05,0.,0.05,0.1])
		plt.tick_params(axis='y', labelsize=19)
		npan_extra = npan_extra+2
		
	if (axistype.get() == 'maj'): plt.xlabel(r'$R_{\rm maj}$ [arcsec]', fontsize=20)
	if (axistype.get() == 'equ'): plt.xlabel(r'$R_{\rm eq}$ [arcsec]', fontsize=20)
	plt.tick_params(axis='x', labelsize=19)	
	if (plottype.get() == 'log'): plt.xscale('log')
	if ((plottype.get() == 'log')):
		if min(x_original)==0.: xminplt=0.45*x_original[1]
		else: xminplt = 0.45*min(x_original) 
		plt.xlim(xminplt,1.25*max(x))
		ax1.xaxis.set_major_formatter(FuncFormatter(log_axis))
	else: 
			plt.xlim(min(x)-0.05*max(x),1.09*max(x))
	plt.minorticks_on()
	plt.tick_params('both', length=5, width=1.2, which='major')
	plt.tick_params('both', length=3, width=1, which='minor')
	
	fig.subplots_adjust(hspace=0,left=0.18, bottom=0.09*(9./n_panels), right=0.97, top=(1.-0.03*(9./n_panels)))
	#fig.tight_layout()
	plt.show()
	#pylab.savefig('Plot.pdf')
	
def plotcomm():

	print 'Reading input information...'
	pxscale=float(PSC.get())
	mzp=float(mag0p.get())
	xmin= float(rlo.get())
	xmax= float(rhi.get())
	fit_mini = (float(rflo.get()))
	fit_maxi = (float(rfhi.get()))

	itt = input_table_type.get()
	
	if itt=='ell': table=readEllipseOutput(IP.get())
	if itt=='iso': table=readtab_ellipsi(IP.get())
	if itt=='ir': table=read2cols(IP.get())
	
	sizearr=len(table['sma'])
	x=[0.]*sizearr
	xt=float(PSC.get())*table['sma']
	if (axistype.get()=='equ'):
		if (itt != 'ir'):
			print 'Data R column should be the semi-major axis (R_maj). Converting to R_eq...' 
			ell=table['ellip']
			xt[1:len(x)] = xt[1:len(x)]*sqrt(1.-ell[1:len(x)])
		else: print 'Cannot generate equivalent axis without ellipticity column. Assuming the data R column is R_eq.'
	x=x+xt
	xmin = datamin(xmin,x)	
	xmax = datamax(xmax,x)
	
	fit_min = fitmin(fit_mini,x)
	fit_max = fitmax(fit_maxi,x)
	
	if (fit_min < xmin): fit_min=xmin
	if (fit_max > xmax): fit_max=xmax
	if (fit_min < 0): fit_min = 0
	
	x=x[xmin:xmax]
	
	if (itt=='ell' or itt=='iso'):
		ell=table['ellip']
		pan=table['PA']
		if itt=='iso': b4=table['BI4']
		if itt=='ell': b4=table['B4']
		grad=table['GRAD']
		sma=table['sma']
		b4=-1.*b4/grad/sma
		if (axistype.get() == 'equ'):
			b4=b4*sqrt(1. - ell)	
		b4=b4[xmin:xmax]
		ell=ell[xmin:xmax]
		pan=pan[xmin:xmax]
	
	data=-2.5*log10(table['intens']/(pxscale**2.)) + mzp
	data=data[xmin:xmax]
	
	print 'Done'
	print 'Plotting data...'

	if (itt!='ir'): n_panels=4+2*(pan_el.get()+pan_pa.get()+pan_b4.get())
	else: n_panels=4
	
	rcParams['figure.figsize'] = 7.5, 0.78*n_panels+1.08
	fig=plt.figure()
	ax1 = plt.subplot2grid((n_panels,3), (0, 0), colspan=3, rowspan=4)
	#----------------------------------------------------------	
	if (fit_min > xmin): ax1.plot(x[xmin:fit_min],data[xmin:fit_min], marker='o', color='white', markersize=7., linestyle='None')
	if (fit_max < xmax): ax1.plot(x[fit_max:xmax],data[fit_max:xmax], marker='o', color='white', markersize=7., linestyle='None')
	ax1.plot(x[fit_min:fit_max],data[fit_min:fit_max], marker='o', color='red', markersize=7., linestyle='None', label='Data')
	plt.text(0.65,0.95, TitlLabEn.get(), fontsize=19, horizontalalignment='center', verticalalignment='top', transform = ax1.transAxes)
	#----------------------------------------------------------	
	plt.gca().invert_yaxis()
	plt.ylabel('$\mu$ [mag arcsec$^{-2}$]', fontsize=20, labelpad=28)
	plt.tick_params(axis='x', labelsize=19)
	plt.tick_params(axis='y', labelsize=19)
	plt.minorticks_on()
	plt.tick_params('both', length=5, width=1.2, which='major')
	plt.tick_params('both', length=3, width=1, which='minor')
	if (plottype.get() == 'log'): plt.xscale('log')
	plt.ylim([max(data)+0.9,min(data)-0.9])
	
	if ((plottype.get() == 'log')): 
		xz=linspace(1e-11,max(x)+0.5*max(x),len(x))
		zer=[0.0]*len(xz)
		plt.plot(xz,zer, linewidth=1.5, color='red')
	else: 
		plt.xlim(min(x)-0.05*max(x),1.09*max(x))
		xz=linspace(min(x)-0.5*max(x),max(x)+0.5*max(x),len(x))
		zer=[0.0]*len(xz)
		plt.plot(xz,zer, linewidth=1.5, color='red')
	
	npan_extra = 4
	
	if pan_el.get()>0 and itt!='ir':
		ax3 = plt.subplot2grid((n_panels,3), (npan_extra, 0), colspan=3, rowspan=2, sharex = ax1)
		ax3.plot(x,ell, marker='o', color='black', markersize=5.5, linestyle='None')
		plt.ylabel('$\epsilon$', fontsize=21, labelpad=20)
		if (n_panels-npan_extra-2 > 0): plt.tick_params(axis='x', labelbottom='off')
		plt.tick_params(axis='y', labelsize=19)
		plt.ylim([-0.04,0.92])
		plt.yticks([0.,0.2, 0.4, 0.6, 0.8])
		plt.minorticks_on()
		plt.tick_params('both', length=5, width=1.2, which='major')
		plt.tick_params('both', length=3, width=1, which='minor')
		npan_extra = npan_extra+2
		
	if pan_pa.get()>0 and itt!='ir':
		ax4 = plt.subplot2grid((n_panels,3), (npan_extra, 0), colspan=3, rowspan=2, sharex=ax1)
		ax4.plot(x,pan, marker='o', color='black', markersize=6, linestyle='None', label='Delta')
		plt.ylabel('$PA$ [deg]', fontsize=20)
		if (n_panels-npan_extra-2 > 0): plt.tick_params(axis='x', labelbottom='off')
		plt.tick_params(axis='y', labelsize=19)
		ax4.set_ylim([-90.,90.])
		ax4.set_yticks([-50,0,50])
		plt.minorticks_on()
		plt.tick_params('both', length=5, width=1.2, which='major')
		plt.tick_params('both', length=3, width=1, which='minor')
		npan_extra=npan_extra+2
	
	if pan_b4.get()>0 and itt!='ir':
		ax5 = plt.subplot2grid((n_panels,3), (npan_extra, 0), colspan=3, rowspan=2, sharex=ax1)
		ax5.plot(x,b4, marker='o', color='black', markersize=5.5, linestyle='None', label='Delta')
		zer=[0.0]*len(xz)
		ax5.plot(xz,zer, linewidth=1.5, color='red')
		ax5.set_ylabel('$B_{4}$', fontsize=20, labelpad=-10)
		ax5.set_ylim([-0.14,0.14])
		ax5.set_yticks([-0.1,-0.05,0.,0.05,0.1])
		plt.tick_params(axis='y', labelsize=19)
		npan_extra = npan_extra+2
		
	if (axistype.get() == 'maj'): plt.xlabel(r'$R_{\rm maj}$ [arcsec]', fontsize=20)
	if (axistype.get() == 'equ'): plt.xlabel(r'$R_{\rm eq}$ [arcsec]', fontsize=20)
	plt.tick_params(axis='x', labelsize=19)	
	if (plottype.get() == 'log'): plt.xscale('log')
	if ((plottype.get() == 'log')):
		if min(x)==0.: xminplt=0.45*x[1]
		else: xminplt = 0.45*min(x) 
		plt.xlim(xminplt,1.25*max(x))
		ax1.xaxis.set_major_formatter(FuncFormatter(log_axis))
	else: 
			plt.xlim(min(x)-0.05*max(x),1.09*max(x))
	plt.minorticks_on()
	plt.tick_params('both', length=5, width=1.2, which='major')
	plt.tick_params('both', length=3, width=1, which='minor')
	
	fig.subplots_adjust(hspace=0,left=0.18, bottom=0.09*(9./n_panels), right=0.97, top=(1.-0.03*(9./n_panels)))
	plt.show() 
	

def quit():
	mainW.destroy()

def add_sersic():
	global rowv, ns
	global sersic_entry
	se=[0.0]*8
	r=rowv
	se[0]=IntVar()
	se[0].set(1)
	CS=Checkbutton(mainW, text="Sersic", variable=se[0]).grid(row=r, sticky=W)
	r=r+1
	LS12 = Label(mainW, text="Sersic n = ")
	LS12.grid(row=r, column=0,sticky=E)
	se[1] = Entry(mainW, bd = 3, width=11)
	se[1].insert(0,1.0)
	se[1].grid(row=r, column=1, sticky=W)
	se[5]=IntVar()
	nvary=Checkbutton(mainW, text="fix", variable=se[5]).grid(row=r, column=2, sticky=W)
	
	LS22 = Label(mainW, text="mu_e = ")
	LS22.grid(row=r, column=3,sticky=E)
	se[2] = Entry(mainW, bd =3, width=11)
	se[2].insert(0,20.0)
	se[2].grid(row=r, column=4, sticky=W)
	se[6]=IntVar()
	muevary=Checkbutton(mainW, text="fix", variable=se[6]).grid(row=r, column=5, sticky=W)

	LS32 = Label(mainW, text="R_e [arcsec] = ")
	LS32.grid(row=r,column=6,sticky=E)
	se[3] = Entry(mainW, bd =3, width=11)
	se[3].insert(0,5.0)
	se[3].grid(row=r,column=7, sticky=W)
	se[7]=IntVar()
	revary=Checkbutton(mainW, text="fix", variable=se[7]).grid(row=r, column=8, sticky=W)

	se[4] = Entry(mainW, bd = 3, w=7)
	se[4].insert(0,'r')
	se[4].grid(row=r, column=9, sticky=W)

	if (ns == 0): sersic_entry=se
	else: sersic_entry=vstack([sersic_entry,se])
	ns=ns+1
	r=r+1
	rowv=r
	
def add_corsic():
	global rowv, nc
	global corsic_entry
	ce=[0.0]*14
	r=rowv
	ce[0]=IntVar()
	ce[0].set(1)
	CS=Checkbutton(mainW, text="core-Sersic", variable=ce[0]).grid(row=r, sticky=W)
	r=r+1
	LC12 = Label(mainW, text="mu_pr = ")
	LC12.grid(row=r, column=0,sticky=E)
	ce[1] = Entry(mainW, bd = 3, width=11)
	ce[1].insert(0,15.0)
	ce[1].grid(row=r, column=1, sticky=W)
	ce[8]=IntVar()
	muvary=Checkbutton(mainW, text="fix", variable=ce[8]).grid(row=r, column=2, sticky=W)
	
	LC22 = Label(mainW, text="alpha_cs = ")
	LC22.grid(row=r, column=3,sticky=E)
	ce[2] = Entry(mainW, bd =3, width=11)
	ce[2].insert(0,2.0)
	ce[2].grid(row=r, column=4, sticky=W)
	ce[9]=IntVar()
	aevary=Checkbutton(mainW, text="fix", variable=ce[9]).grid(row=r, column=5, sticky=W)

	LC32 = Label(mainW, text="R_b [arcsec] = ")
	LC32.grid(row=r,column=6,sticky=E)
	ce[4] = Entry(mainW, bd =3, width=11)
	ce[4].insert(0,0.5)
	ce[4].grid(row=r,column=7, sticky=W)
	ce[11]=IntVar()
	revary=Checkbutton(mainW, text="fix", variable=ce[11]).grid(row=r, column=8, sticky=W)
	
	r=r+1
	
	LC32 = Label(mainW, text="Sersic n = ")
	LC32.grid(row=r,column=0,sticky=E)
	ce[6] = Entry(mainW, bd =3, width=11)
	ce[6].insert(0,3.0)
	ce[6].grid(row=r,column=1, sticky=W)
	ce[13]=IntVar()
	revary=Checkbutton(mainW, text="fix", variable=ce[13]).grid(row=r, column=2, sticky=W)
	
	LC22 = Label(mainW, text="gamma_cs = ")
	LC22.grid(row=r, column=3,sticky=E)
	ce[3] = Entry(mainW, bd =3, width=11)
	ce[3].insert(0,0.2)
	ce[3].grid(row=r, column=4, sticky=W)
	ce[10]=IntVar()
	muevary=Checkbutton(mainW, text="fix", variable=ce[10]).grid(row=r, column=5, sticky=W)
	
	LC32 = Label(mainW, text="R_e [arcsec] = ")
	LC32.grid(row=r,column=6,sticky=E)
	ce[5] = Entry(mainW, bd =3, width=11)
	ce[5].insert(0,5.0)
	ce[5].grid(row=r,column=7, sticky=W)
	ce[12]=IntVar()
	revary=Checkbutton(mainW, text="fix", variable=ce[12]).grid(row=r, column=8, sticky=W)

	ce[7] = Entry(mainW, bd = 3, w=7)
	ce[7].insert(0,'r')
	ce[7].grid(row=r, column=9, sticky=W)

	if (nc == 0): corsic_entry=ce
	else: corsic_entry=vstack([corsic_entry,ce])
	nc=nc+1
	r=r+1
	rowv=r

def add_exponential():
	global rowv, ne, exp_entry
	r=rowv
	ee=[0.0]*6
	ee[0]=IntVar()
	ee[0].set(1)
	CD=Checkbutton(mainW, text="Exponential", variable=ee[0]).grid(row=r, sticky=W)
	r=r+1
	LD1 = Label(mainW, text="h [arcsec] = ")
	LD1.grid(row=r, column=0,sticky=E)
	ee[1] = Entry(mainW, bd = 3, width=10)
	ee[1].insert(0,10.0)
	ee[1].grid(row=r, column=1, sticky=W)
	ee[4]=IntVar()
	hvary=Checkbutton(mainW, text="fix", variable=ee[4]).grid(row=r, column=2, sticky=W)

	LD2 = Label(mainW, text="mu_0 = ")
	LD2.grid(row=r, column=3,sticky=E)
	ee[2] = Entry(mainW, bd =3, width=10)
	ee[2].insert(0,20.0)
	ee[2].grid(row=r, column=4, sticky=W)
	ee[5]=IntVar()
	muvary=Checkbutton(mainW, text="fix", variable=ee[5]).grid(row=r, column=5, sticky=W)

	ee[3] = Entry(mainW, bd = 3, w=7)
	ee[3].insert(0,'b')
	ee[3].grid(row=r, column=9, sticky=W)

	if (ne == 0): exp_entry=ee
	else: exp_entry=vstack([exp_entry,ee])
	ne=ne+1
	r=r+1
	rowv=r
	
def add_id():
	global rowv, nid, id_entry
	r=rowv
	es=[0.0]*7
	es[0]=IntVar()
	es[0].set(1)
	CD=Checkbutton(mainW, text=r"Incl. disc", variable=es[0]).grid(row=r, sticky=W)
	r=r+1
	LD1 = Label(mainW, text="z0 [arcsec] = ")
	LD1.grid(row=r, column=0,sticky=E)
	es[1] = Entry(mainW, bd = 3, width=10)
	es[1].insert(0,5.0)
	es[1].grid(row=r, column=1, sticky=W)
	es[4]=IntVar()
	hvary=Checkbutton(mainW, text="fix", variable=es[4]).grid(row=r, column=2, sticky=W)

	LD2 = Label(mainW, text="mu_0z = ")
	LD2.grid(row=r, column=3,sticky=E)
	es[2] = Entry(mainW, bd =3, width=10)
	es[2].insert(0,20.0)
	es[2].grid(row=r, column=4, sticky=W)
	es[5]=IntVar()
	muvary=Checkbutton(mainW, text="fix", variable=es[5]).grid(row=r, column=5, sticky=W)
	
	es[6]=StringVar()
	CaseB=Radiobutton(mainW, text="r=0", variable=es[6], value='r0').grid(row=r, column=6, sticky=E)
	CaseB=Radiobutton(mainW, text="z=0", variable=es[6], value='z0').grid(row=r, column=7, sticky=W)
	es[6].set('r0')
	
	Be = Button(mainW, text ="???", font=("Palatino", 16),command = incl_disc).grid(row=r, column=8, sticky=W)

	es[3] = Entry(mainW, bd = 3, w=7)
	es[3].insert(0,'b')
	es[3].grid(row=r, column=9, sticky=W)

	if (nid == 0): id_entry=es
	else: id_entry=vstack([id_entry,es])
	nid=nid+1
	r=r+1
	rowv=r

def add_gaussian():
	global rowv, ng, gauss_entry
	r=rowv
	ge=[0.0]*8
	ge[0]=IntVar()
	ge[0].set(1)
	CR=Checkbutton(mainW, text="Gaussian", variable=ge[0]).grid(row=r, sticky=W)
	r=r+1
	LR1 = Label(mainW, text="R0r [arcsec] = ")
	LR1.grid(row=r, column=0,sticky=E)
	ge[1] = Entry(mainW, bd = 3, width=10)
	ge[1].insert(0,2.0)
	ge[1].grid(row=r, column=1, sticky=W)
	ge[5]=IntVar()
	rvary=Checkbutton(mainW, text="fix", variable=ge[5]).grid(row=r, column=2, sticky=W)

	LR2 = Label(mainW, text="mu_0r = ")
	LR2.grid(row=r, column=3,sticky=E)
	ge[2] = Entry(mainW, bd =3, width=10)
	ge[2].insert(0,20.0)
	ge[2].grid(row=r, column=4, sticky=W)
	ge[6]=IntVar()
	mvary=Checkbutton(mainW, text="fix", variable=ge[6]).grid(row=r, column=5, sticky=W)

	LR3 = Label(mainW, text="FWHMr [arcsec] = ")
	LR3.grid(row=r, column=6,sticky=E)
	ge[3] = Entry(mainW, bd =3, width=10)
	ge[3].insert(0,1.0)
	ge[3].grid(row=r, column=7, sticky=W)
	ge[7]=IntVar()
	fvary=Checkbutton(mainW, text="fix", variable=ge[7]).grid(row=r, column=8, sticky=W)

	ge[4] = Entry(mainW, bd = 3, w=7)
	ge[4].insert(0,'c')
	ge[4].grid(row=r, column=9, sticky=W)

	if (ng == 0): gauss_entry=ge
	else: gauss_entry=vstack([gauss_entry,ge])
	ng=ng+1

	r=r+1
	rowv=r

def add_psf():
	global rowv, np, psf_entry
	r=rowv
	pe=[0.0]*4
	pe[0]=IntVar()
	pe[0].set(1)
	CP=Checkbutton(mainW, text="PSF", variable=pe[0]).grid(row=r, sticky=W)
	r=r+1
	LP = Label(mainW, text="mu0p = ")
	LP.grid(row=r, column=0,sticky=E)
	pe[1] = Entry(mainW, bd = 3, width=10)
	pe[1].insert(0,15.)
	pe[1].grid(row=r, column=1, sticky=W)
	pe[3]=IntVar()
	mvary=Checkbutton(mainW, text="fix", variable=pe[3]).grid(row=r, column=2, sticky=W)

	pe[2] = Entry(mainW, bd = 3, w=7)
	pe[2].insert(0,'limegreen')
	pe[2].grid(row=r, column=9, sticky=W)
	
	if (np == 0): psf_entry=pe
	else: psf_entry=vstack([psf_entry,pe])
	np=np+1
	
	r=r+1
	rowv=r

def add_ferrer():
	global rowv, nf, ferrer_entry
	r=rowv
	fe=[0.0]*10
	fe[0]=IntVar()
	fe[0].set(1)
	CB=Checkbutton(mainW, text="Ferrer", variable=fe[0]).grid(row=r, sticky=W)
	r=r+1
	LB1 = Label(mainW, text="R_out [arcsec] = ")
	LB1.grid(row=r, column=0, sticky=E)
	fe[1] = Entry(mainW, bd = 3, width=10)
	fe[1].insert(0,30.0)
	fe[1].grid(row=r, column=1, sticky=W)
	fe[6]=IntVar()
	rvary=Checkbutton(mainW, text="fix", variable=fe[6]).grid(row=r, column=2, sticky=W)

	LB2 = Label(mainW, text="mu_0f = ")
	LB2.grid(row=r, column=3, sticky=E)
	fe[2] = Entry(mainW, bd =3, width=10)
	fe[2].insert(0,20.0)
	fe[2].grid(row=r, column=4, sticky=W)
	fe[7]=IntVar()
	mvary=Checkbutton(mainW, text="fix", variable=fe[7]).grid(row=r, column=5, sticky=W)

	LB3 = Label(mainW, text="alpha_F = ")
	LB3.grid(row=r, column=6, sticky=E)
	fe[3] = Entry(mainW, bd =3, width=10)
	fe[3].insert(0,0.5)
	fe[3].grid(row=r, column=7, sticky=W)
	fe[8]=IntVar()
	avary=Checkbutton(mainW, text="fix", variable=fe[8]).grid(row=r, column=8, sticky=W)
	r=r+1
	LB4 = Label(mainW, text="beta_F = ")
	LB4.grid(row=r, column=6, sticky=E)
	fe[4] = Entry(mainW, bd =3, width=10)
	fe[4].insert(0,0.1)
	fe[4].grid(row=r, column=7, sticky=W)
	fe[9]=IntVar()
	bvary=Checkbutton(mainW, text="fix", variable=fe[9]).grid(row=r, column=8, sticky=W)

	fe[5] = Entry(mainW, bd = 3, w=7)
	fe[5].insert(0,'orange')
	fe[5].grid(row=r, column=9, sticky=W)
	
	if (nf == 0): ferrer_entry=fe
	else: ferrer_entry=vstack([ferrer_entry,fe])
	nf=nf+1
	
	r=r+1
	rowv=r
	
def add_td():
	global rowv, ntd, td_entry
	r=rowv
	td=[0.0]*10
	td[0]=IntVar()
	td[0].set(1)
	CB=Checkbutton(mainW, text="Truncated disc", variable=td[0]).grid(row=r, sticky=W)
	r=r+1
	LB1 = Label(mainW, text="mu_0 = ")
	LB1.grid(row=r, column=0, sticky=E)
	td[1] = Entry(mainW, bd = 3, width=10)
	td[1].insert(0,20.0)
	td[1].grid(row=r, column=1, sticky=W)
	td[6]=IntVar()
	rvary=Checkbutton(mainW, text="fix", variable=td[6]).grid(row=r, column=2, sticky=W)

	LB2 = Label(mainW, text="r_b [arcsec] = ")
	LB2.grid(row=r, column=3, sticky=E)
	td[2] = Entry(mainW, bd =3, width=10)
	td[2].insert(0,20.0)
	td[2].grid(row=r, column=4, sticky=W)
	td[7]=IntVar()
	mvary=Checkbutton(mainW, text="fix", variable=td[7]).grid(row=r, column=5, sticky=W)

	LB3 = Label(mainW, text="h1 = ")
	LB3.grid(row=r, column=6, sticky=E)
	td[3] = Entry(mainW, bd =3, width=10)
	td[3].insert(0,5.0)
	td[3].grid(row=r, column=7, sticky=W)
	td[8]=IntVar()
	avary=Checkbutton(mainW, text="fix", variable=td[8]).grid(row=r, column=8, sticky=W)
	r=r+1
	LB4 = Label(mainW, text="h2 = ")
	LB4.grid(row=r, column=6, sticky=E)
	td[4] = Entry(mainW, bd =3, width=10)
	td[4].insert(0,5.0)
	td[4].grid(row=r, column=7, sticky=W)
	td[9]=IntVar()
	bvary=Checkbutton(mainW, text="fix", variable=td[9]).grid(row=r, column=8, sticky=W)

	td[5] = Entry(mainW, bd = 3, w=7)
	td[5].insert(0,'b')
	td[5].grid(row=r, column=9, sticky=W)
	
	if (ntd == 0): td_entry=td
	else: td_entry=vstack([td_entry,td])
	ntd=ntd+1
	
	r=r+1
	rowv=r
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------main GUI window construction--------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
mainW = Tk()
mainW.title("profiler 1.0 : 1-D Surface Brightness Profile Fitting")
from pylab import rcParams
rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

rowv=0

Profiler = Label(mainW, text="profiler", font=("Palatino", 24, "bold italic")).grid(row=rowv, column=4)
B3=Button(mainW, text="Quit", command=quit)
B3.grid(row=rowv,column=9)
rowv=rowv+1
Profiler = Label(mainW, text="v1.0", font=("Palatino", 20)).grid(row=rowv, column=4)
Lspace = Label(mainW, text="------------").grid(row=rowv, column=3)
Lspace = Label(mainW, text="------------").grid(row=rowv, column=5)
rowv=rowv+1
Lspace = Label(mainW, text=" ").grid(row=rowv, column=3)
rowv=rowv+1

Image=Label(mainW, text='Input information:', font=("Palatino", 16, "bold italic")).grid(row=rowv, column=0, sticky=W)

B1 = Button(mainW, text ="Fit & plot", font=("Palatino", 16, "bold italic"),command = fitcomm)
B1.grid(row=rowv,column=4)
rowv=rowv+1
NP=Label(mainW, text='Data file: ').grid(row=rowv, column=0)
IP=Entry(mainW, bd=3, width=20)
plottype=StringVar()
RB2=Radiobutton(mainW, text="Linear plot", variable=plottype, value='lin').grid(row=rowv, column=7, sticky=W)

NPty=Label(mainW, text='Panels: ').grid(row=rowv, column=8, sticky=W)

input_table_type=StringVar()
B0 = Button(mainW, text ="Plot data", font=("Palatino", 16),command = plotcomm)
B0.grid(row=rowv,column=4)
rowv=rowv+1
RBitt2=Radiobutton(mainW, text="isofit table", variable=input_table_type, value='iso').grid(row=rowv, column=2, sticky=W)
RB2=Radiobutton(mainW, text="Log plot", variable=plottype, value='log').grid(row=rowv, column=7, sticky=W)
pan_el = IntVar()
pan_el.set(1)
pan_e=Checkbutton(mainW, text="ellip.", variable=pan_el).grid(row=rowv, column=8, sticky=W)
plottype.set('lin')

IP.insert(0,".dat")
IP.grid(row=rowv,column=0,sticky=E)

rowv=rowv+1
RBitt=Radiobutton(mainW, text="ellipse table", variable=input_table_type, value='ell').grid(row=rowv, column=2, sticky=W)
input_table_type.set('iso')

axistype=StringVar()
RB3=Radiobutton(mainW, text="Major", variable=axistype, value='maj').grid(row=rowv, column=7, sticky=W)
pan_pa = IntVar()
pan_pa.set(1)
pan_pang=Checkbutton(mainW, text="P.A.", variable=pan_pa).grid(row=rowv, column=8, sticky=W)
TitlLab=Label(mainW, text='Plot title: ').grid(row=rowv, column=0)
TitlLabEn=Entry(mainW, bd=3, width=20)

rowv=rowv+1

RBitt3=Radiobutton(mainW, text="I(R) table", variable=input_table_type, value='ir').grid(row=rowv, column=2, sticky=W)
input_table_type.set('ell')
TitlLabEn.insert(0,"")
TitlLabEn.grid(row=rowv,column=0,sticky=E)

RB3=Radiobutton(mainW, text="Equiv.", variable=axistype, value='equ').grid(row=rowv, column=7, sticky=W)
pan_b4 = IntVar()
pan_b4.set(1)
pan_b4b=Checkbutton(mainW, text="B4", variable=pan_b4).grid(row=rowv, column=8, sticky=W)
axistype.set('maj')

rowv=rowv+1
Lspace = Label(mainW, text=" ").grid(row=rowv, column=3)
rowv=rowv+1
PScale=Label(mainW, text='Scale ["/px] = ').grid(row=rowv,column=0,sticky=E)
PSC=Entry(mainW, bd=3, width=10)
PSC.insert(0,1.0)
#PSC.insert(0,"0.12825")
PSC.grid(row=rowv, column=1,sticky=W)
mag0p=Label(mainW, text='zero-pt. mag = ').grid(row=rowv,column=2,sticky=E)
mag0p=Entry(mainW, bd=3, width=11)
mag0p.insert(0,"21.5")
mag0p.grid(row=rowv, column=3,sticky=W)

rlowlab=Label(mainW, text='Data range ["]:').grid(row=rowv, column=4,sticky=E)
rlo=Entry(mainW, bd=3, width=7)
rlo.insert(0,"0.")
rlo.grid(row=rowv, column=5, sticky=E)

rhi=Entry(mainW, bd=3, width=7)
rhi.insert(0,"50.")
rhi.grid(row=rowv, column=6,sticky=W)

rowv=rowv+1
Skytxt=Label(mainW, text=r'Threshold [mag] = ').grid(row=rowv,column=0,sticky=E)
Sky=Entry(mainW, bd=3, width=10)
Sky.insert(0,"99.")
Sky.grid(row=rowv, column=1,sticky=W)
#PSamp=Label(mainW, text='Sampling [px] = ').grid(row=rowv,column=0,sticky=E)
#PSAM=Entry(mainW, bd=3, width=10)
#PSAM.insert(0,"0.3")
#PSAM.grid(row=rowv, column=1,sticky=W)
#mag0p=Label(mainW, text='zero-pt. mag = ').grid(row=rowv,column=2,sticky=E)
#mag0p=Entry(mainW, bd=3, width=11)
#mag0p.insert(0,"24.6949")
#mag0p.grid(row=rowv, column=3,sticky=W)

rflowlab=Label(mainW, text='Fit range ["]:').grid(row=rowv, column=4,sticky=E)
rflo=Entry(mainW, bd=3, width=7)
rflo.insert(0,"0.")
rflo.grid(row=rowv, column=5, sticky=E)

rfhi=Entry(mainW, bd=3, width=7)
rfhi.insert(0,"45.")
rfhi.grid(row=rowv, column=6,sticky=W)

rowv=rowv+1


Lspace = Label(mainW, text=" ").grid(row=rowv, column=3)
rowv=rowv+1
ps=StringVar()
PSFChoice=Label(mainW, text='Point Spread Function : ', font=("Palatino", 16, "bold italic")).grid(row=rowv,column=0, sticky=W)

ellip_label=Label(mainW, text='Global Ellipticity:').grid(row=rowv,column=7, sticky=E)
ellipentry=Entry(mainW, bd=3, width=10)
ellipentry.insert(0,0.0)
ellipentry.grid(row=rowv, column=8, sticky=E)

rowv=rowv+1
RB=Radiobutton(mainW, text="Gaussian", variable=ps, value='gaussian').grid(row=rowv, column=1, sticky=W)
PGF=Label(mainW, text='FWHM [pix] = ').grid(row=rowv,column=2, sticky=E)
PGFV=Entry(mainW, bd=3, width=10)
PGFV.insert(0,"2.0")
PGFV.grid(row=rowv, column=3, sticky=W)
rowv=rowv+1
RB=Radiobutton(mainW, text="Moffat", variable=ps, value='moffat').grid(row=rowv, column=1, sticky=W)
PMF=Label(mainW, text='FWHM [pix] = ').grid(row=rowv,column=2, sticky=E)
PMFV=Entry(mainW, bd=3, width=10)
PMFV.insert(0,"2.0")
PMFV.grid(row=rowv, column=3,sticky=W)
PMB=Label(mainW, text='beta = ').grid(row=rowv,column=4, sticky=E)
PMBV=Entry(mainW, bd=3, width=10)
PMBV.insert(0,"4.765")
PMBV.grid(row=rowv, column=5, sticky=W)
rowv=rowv+1
RB=Radiobutton(mainW, text="numerical", variable=ps, value='numerical').grid(row=rowv, column=1, sticky=W)
PN=Label(mainW, text='PSF file ').grid(row=rowv,column=2, sticky=E)
PNV=Entry(mainW, bd=3, width=10)
PNV.insert(0,".dat")
PNV.grid(row=rowv, column=3, sticky=W)
ps.set('moffat')

rowv=rowv+1
Lspace = Label(mainW, text="  ").grid(row=rowv, column=3)
rowv=rowv+1
text=Label(mainW, text='Choose your weapons:', font=("Palatino", 16, "bold italic")).grid(row=rowv,column=0, sticky=W)

AS = Button(mainW, text = '+ Sersic', command=add_sersic, width=11)
AS.grid(row=rowv, column=1)


AE = Button(mainW, text = '+ Exponential', command=add_exponential, width=11)
AE.grid(row=rowv, column=3)

AG = Button(mainW, text = '+ Gaussian', command=add_gaussian, width=11)
AG.grid(row=rowv, column=4)


AP = Button(mainW, text = '+ PSF', command=add_psf, width=11)
AP.grid(row=rowv, column=5)

AF = Button(mainW, text = '+ Ferrer', command=add_ferrer, width=11)
AF.grid(row=rowv, column=6)

ASi = Button(mainW, text = r'+ Inclined disc', command=add_id, width=11)
ASi.grid(row=rowv, column=7)

ACS = Button(mainW, text = '+ core-Sersic', command=add_corsic, width=11)
ACS.grid(row=rowv, column=2)

ATD = Button(mainW, text = '+ Truncd. disc', command=add_td, width=11)
ATD.grid(row=rowv, column=8)

colour = Label(mainW, text="colour *").grid(row=rowv, column=9)


rowv=rowv+1
Lspace = Label(mainW, text="-----------------------").grid(row=rowv, column=0)
Lspace = Label(mainW, text="------------").grid(row=rowv, column=1)
Lspace = Label(mainW, text="------------").grid(row=rowv, column=2)
Lspace = Label(mainW, text="------------").grid(row=rowv, column=3)
Lspace = Label(mainW, text="------------").grid(row=rowv, column=4)
Lspace = Label(mainW, text="------------").grid(row=rowv, column=5)
Lspace = Label(mainW, text="------------").grid(row=rowv, column=6)
Lspace = Label(mainW, text="------------").grid(row=rowv, column=7)
Lspace = Label(mainW, text="------------").grid(row=rowv, column=8)
#Lspace = Label(mainW, text="------------").grid(row=rowv, column=9)
#Lspace = Label(mainW, text="------------").grid(row=rowv, column=10)
#Lspace = Label(mainW, text="", width=11).grid(row=rowv, column=7)
#Lspace = Label(mainW, text="", width=11).grid(row=rowv, column=8)
#Lspace = Label(mainW, text="", width=11).grid(row=rowv, column=9)
rowv=rowv+1

def incl_disc():

	explW = Tk()
	explW.title("Inclined disc options")
	text = Label(explW, text='This type of component is adapted from the inclined disk surface brightness equation in', font=("Palatino", 16)).grid(row=0, column=0, sticky=W)
	text = Label(explW, text='van der Kruit & Searle (1981):', font=("Palatino", 16)).grid(row=1, column=0, sticky=W)
	Lspace = Label(explW, text="").grid(row=2, column=0, sticky=W)
	txt = Label(explW, text=r'I(r,z) = I(0,0) (r/h) K_1(r/h) [sech(z/z_s)]^2', font=("Palatino", 16, "bold")).grid(row=3, column=0)
	Lspace = Label(explW, text="").grid(row=4, column=0)
	txt2 = Label(explW, text=r'where h is the scale length and z_s the scale height. The above is a 2D description, however. ', font=("Palatino", 16)).grid(row=5, column=0, sticky=W)
	txt3 = Label(explW, text=r'Here it is reduced to 1D by taking the special cases r=0 (minor axis) or z=0 (major axis).', font=("Palatino", 16)).grid(row=6, column=0, sticky=W)
	explW.mainloop()

mainW.mainloop()


