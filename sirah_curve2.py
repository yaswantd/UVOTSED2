import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from os import path
from scipy.interpolate import UnivariateSpline
import scipy.stats as stats
import sys
import glob





#Setup the bolometric spline function of SN2011fe
bolo = pd.read_csv('SN2011fe_bolo.txt',delim_whitespace=True)


plt.figure(1,figsize=(7,5),dpi=250)
plt.subplot(211)
phase = np.asarray(bolo.DaysSinceMax)+3.0
p2 = np.where((-10.0 <= phase) & (phase <= 15.0))
phase[p2]
nux = np.linspace(phase[p2][0],phase[p2][-1],100)

ubol = np.asarray(bolo.U)
#ubol = 18.34-(np.log10(ubol)/0.4)
uspl = UnivariateSpline(phase[p2],ubol[p2],w=1/ubol[p2],s=0.01)
plt.plot(nux,uspl(nux),color='red')
plt.scatter(phase[p2],ubol[p2],label='U',color='red')

phase = np.asarray(bolo.DaysSinceMax)
p2 = np.where((-10.0 <= phase) & (phase <= 15.0))
phase[p2]
nux = np.linspace(phase[p2][0],phase[p2][-1],100)

bbol = np.asarray(bolo.B)
#bbol = 19.11-(np.log10(bbol)/0.4)
#bspl = UnivariateSpline(phase[p2],bbol[p2])
#plt.plot(nux,bspl(nux),color='purple')
plt.scatter(phase[p2],bbol[p2],label='B',color='purple')

bspl = UnivariateSpline(phase[p2],bbol[p2],w=1/bbol[p2],s=0.01)
plt.plot(nux,bspl(nux),color='purple')

phase = np.asarray(bolo.DaysSinceMax)-2.0
p2 = np.where((-10.0 <= phase) & (phase <= 15.0))
phase[p2]
nux = np.linspace(phase[p2][0],phase[p2][-1],100)

vbol = np.asarray(bolo.V)
#vbol = 17.89-(np.log10(vbol)/0.4)
vspl = UnivariateSpline(phase[p2],vbol[p2],w=1/vbol[p2],s=0.01)
plt.plot(nux,vspl(nux),color='brown')
plt.scatter(phase[p2],vbol[p2],label='V',color='brown')
plt.legend(loc='upper right')
plt.xlabel('Phase Relative to Peak')
plt.ylabel('Photon Counts')


plt.subplot(212)
phase = np.asarray(bolo.DaysSinceMax)+3.0
p2 = np.where((-10.0 <= phase) & (phase <= 15.0))
phase[p2]
nux = np.linspace(phase[p2][0],phase[p2][-1],100)

w1bol = np.asarray(bolo.W1)
#w1bol = 17.44-(np.log10(w1bol)/0.4)
w1spl = UnivariateSpline(phase[p2],w1bol[p2],w=1/w1bol[p2],s=0.01)
plt.plot(nux,w1spl(nux),color='green')
plt.scatter(phase[p2],w1bol[p2],label='UVW1',color='green')

phase = np.asarray(bolo.DaysSinceMax)+3.0
p2 = np.where((-10.0 <= phase) & (phase <= 15.0))
phase[p2]
nux = np.linspace(phase[p2][0],phase[p2][-1],100)

m2bol = np.asarray(bolo.M2)
#m2bol = 16.85-(np.log10(m2bol)/0.4)
m2spl = UnivariateSpline(phase[p2],m2bol[p2],w=1/m2bol[p2],s=0.01)
plt.plot(nux,m2spl(nux),color='orange')
plt.scatter(phase[p2],m2bol[p2],label='UVM2',color='orange')

phase = np.asarray(bolo.DaysSinceMax)+3.0
p2 = np.where((-10.0 <= phase) & (phase <= 15.0))
phase[p2]
nux = np.linspace(phase[p2][0],phase[p2][-1],100)

w2bol = np.asarray(bolo.W2countrate)
#w2bol = 17.38-(np.log10(w2bol)/0.4)
w2spl = UnivariateSpline(phase[p2],w2bol[p2],w=1/w2bol[p2],s=0.01)
plt.plot(nux,w2spl(nux),color='blue')
plt.scatter(phase[p2],w2bol[p2],label='UVW2',color='blue')
plt.xlabel('Phase Relative to Peak')
plt.ylabel('Photon Counts')
plt.legend(loc='upper right')
plt.suptitle('SN2011fe Light Curve')
plt.tight_layout()
plt.savefig('sirah/bolometric_curves2.png',facecolor='white')

plt.close()


sys.exit()

# #read in the swift data for the SN of choice
# df = pd.read_csv('sirah/sirahdata/SN2020jgl_uvotB22.1.dat',delim_whitespace=True,comment='#')
# df.columns=['Filter','MJD','Mag','Magerr','3SigErr','0.98Lim','Rate','RateErr','Ap','Frametime','Exp','Telapse']

# #Guess peak data, define filter subsets
# w2data = df.loc[(df.Filter == 'UVW2') & (np.isnan(df.Mag) == False )]
# m2data = df.loc[(df.Filter == 'UVM2') & (np.isnan(df.Mag) == False )]
# w1data = df.loc[(df.Filter == 'UVW1') & (np.isnan(df.Mag) == False )]
# udata = df.loc[(df.Filter == 'U') & (np.isnan(df.Mag) == False )]
# bdata = df.loc[(df.Filter == 'B') & (np.isnan(df.Mag) == False )]
# vdata = df.loc[(df.Filter == 'V') & (np.isnan(df.Mag) == False )]

# pkl = bdata.loc[(bdata.Mag == min(bdata.Mag))]
# pkl = np.asarray(pkl.MJD)[0]
# print('best guess for start date is '+str(pkl))


#Do the fitting 

def sirfit(band, sfile, pk):

	sname = sfile.split('/')[-1]
	sname = sname.split('_')[0]

	if band == 'V':
		pk_x = pk + 2.0
		sdat = vdata
		def sfit(x,a,c):
			return vspl(x/a) *c

	elif band == 'B':
		pk_x = pk 
		sdat = bdata
		def sfit(x,a,c):
			return bspl(x/a) *c

	elif band == 'U':
		pk_x = pk-2.0
		sdat = udata
		def sfit(x,a,c):
			return uspl(x/a) *c

	elif band == 'UVW1':
		pk_x = pk-2.0
		sdat = w1data
		def sfit(x,a,c):
			return w1spl(x/a) *c

	elif band ==  'UVM2':
		pk_x = pk-2.0
		sdat = m2data
		def sfit(x,a,c):
			return m2spl(x/a) *c

	elif band == 'UVW2':
		pk_x = pk-2.0
		sdat = w2data
		def sfit(x,a,c):
			return w2spl(x/a) *c
	else:
		print('Error: Band not defined')
		return

	plt.figure(2,figsize=(7,5),dpi=250)
	plt.title(band)
	plt.subplot(211)
#if uf == 1:
	draw = np.random.normal(pk_x,1.5,1000)
	def sfit(x,a,c):
	    return vspl(x/a) *c
	c1 = 100
	c2 = 0
	c3 = 0
	up1, up2, up3 = np.asarray(sdat.MJD), np.asarray(sdat.Rate), np.asarray(sdat.RateErr)
	chid = []
	lxd = []
	allrs = []
	print(up1)
	for i in range(len(draw)):
	    uphase1 = np.asarray(sdat.MJD) - draw[i]
	    umag1 = np.asarray(sdat.Rate)
	    uerr1 = np.asarray(sdat.RateErr)
	    uphase, umag, uerr = [],[],[]
	    for j in range(len(uphase1)):
	        if uphase1[j] >= -5.0 and uphase1[j] <= 5.0:
	            uphase.append(uphase1[j])
	            umag.append(umag1[j])
	            uerr.append(uerr1[j])

	    uphase = np.asarray(uphase)
	    umag = np.asarray(umag)
	    uerr = np.asarray(uerr)
	    if len(uphase) < 3:
	        continue

	    try:
	    	popt,pcov = curve_fit(sfit,uphase,umag,sigma=uerr,p0=[1,1])
	    except RuntimeError:
	    	print('Not Fit')
	    	continue


	    ynew = sfit(uphase,popt[0],popt[1])
	    chi2 = (sum((umag - ynew)**2 / (uerr)**2))
	    chid.append(chi2)
	    lxd.append(draw[i])
	    #allrs.append(sfit(0.0,popt[0],popt[1]))
	    if chi2 < c1:
	        c1 = chi2
	        c2 = popt
	        c3 = draw[i]
	        c4 = pcov
	        up1, up2, up3= uphase, umag, uerr
	if c1 == 100:
	    plt.text(0.5,0.5,'Model Unable to Converge')
	    to_u, mo_u = 0,0
	    nu_y = 0.0
	else:
	    xnew = np.linspace(up1[0],up1[-1],100)
	    #print(up1,up2,up3)
	    plt.errorbar(up1,up2,yerr=up3,ls='none',marker='o',label='Observations')
	    plt.plot(xnew,sfit(xnew,c2[0],c2[1]),label='Best Fit')
	    plt.plot(xnew,sfit(xnew,1.0,c2[1]),ls='--',label='SN2011fe Template')
	    #plt.plot(xnew,spl(xnew))
	    nu_y = sfit(0.0,c2[0],c2[1])
	    to_u, mo_u = c3, nu_y
	    #plt.gca().invert_yaxis()
	    plt.minorticks_on()
	    plt.tick_params(which='both',direction='in')
	    plt.legend(loc='best')
	    plt.xlabel('Phase Relative to Peak')
	    plt.ylabel('Flux (Photon Counts)')
	    print(c1,c3)
	    print(nu_y,c4[1][1])
	    
	plt.subplot(212)
	jux = np.asarray(lxd)-c3
	juy2 = 1.0/np.asarray(chid)

	jux = np.asarray(jux)
	d = {'col1': jux, 'col2': juy2}
	df = pd.DataFrame(data=d)
	ndf = df.sort_values(by=['col1'])
	plt.scatter(ndf.col1 +c3,ndf.col2)
	plt.xlabel('MJD')
	plt.ylabel('$1/\chi^2$')
	def Gauss(x, mu, sigma,A):
	    #y = A*np.exp(-1*B*x**2)

	    y = A*(1/(sigma * np.sqrt(2 * np.pi))) * np.exp( - 0.5*((x - mu)/sigma)**2)
	    return y
	                                                
	try:
		poptg, pcovg = curve_fit(Gauss, ndf.col1, ndf.col2)
		plt.plot(ndf.col1 +c3,Gauss(ndf.col1,poptg[0],poptg[1],poptg[2]))
		print(poptg)
		t0best = round(c3,2)
		t0errbest = round(poptg[1],2)
		print(t0best,t0errbest)
		plt.text(c3-4.0,0.5*max(ndf.col2),'T0: '+str(t0best)+' err: '+str(t0errbest))

		#plt.text(c3-4.5,0.25*max(ndf.col2),'R0: '+str(round(nu_y,2))+' err: '+str(round(c4[1][1],2)))

	except RuntimeError:
		plt.text(c3-4.5,0.5*max(ndf.col2),'No Optimal Fit')
		t0best, t0errbest = np.nan, np.nan
		poptg = [0.0,0.0,0.0]
	except ValueError:
		plt.ylabel('No Optimal Fit (2)')
		t0best, t0errbest = np.nan, np.nan
		poptg = [0.0,0.0,0.0]
		return [sname, band, np.nan, np.nan, np.nan, np.nan]


	uphase1 = np.asarray(sdat.MJD) - (c3-poptg[1])
	umag1 = np.asarray(sdat.Rate)
	uerr1 = np.asarray(sdat.RateErr)
	uphase, umag, uerr = [],[],[]
	for j in range(len(uphase1)):
	    if uphase1[j] >= -10.0 and uphase1[j] <= 15.0:
	        uphase.append(uphase1[j])
	        umag.append(umag1[j])
	        uerr.append(uerr1[j])

	uphase = np.asarray(uphase)
	umag = np.asarray(umag)
	uerr = np.asarray(uerr)
	if len(uphase) >= 3 and poptg[1] != 0.0:
		try:
			popt,pcov = curve_fit(sfit,uphase,umag,sigma=uerr,p0=[1,1])
			rerr = sfit(0.0,popt[0],popt[1])
			plt.text(c3-4.0,0.25*max(ndf.col2),'F0: '+str(round(nu_y,2))+' err: '+str(round(abs(nu_y-rerr),2)))
		except RuntimeError:
			plt.text(c3-4.5,0.25*max(ndf.col2),'No F0 err Fit (RuntimeError)')
			rerr = np.nan
	else:
		plt.text(c3-4.5,0.25*max(ndf.col2),'No F0 err Fit '+str(len(uphase)))
		rerr = np.nan

	plt.suptitle(str(sname)+' '+str(band)+' Filter')
	plt.tight_layout()
	plt.savefig('sirah/thesis_'+sname+'_'+band+'.png',facecolor='white')
	plt.close()

	# print('T0: '+str(c3)+' err: '+str(popt[1]))
	# print('M0: '+str(nu_y)+' err: '+str(c4[1][1]))

	return [sname, band, t0best, t0errbest, round(nu_y,2), round(abs(nu_y-rerr),2)]



outdf = pd.DataFrame(columns=['SNname','Filter','T0','T0_err','R0','R0_err'])


indf = pd.read_csv('SNPY_Sample_Decline.csv')


mysn = indf.sname

for i in range(len(mysn)):
	if mysn[i] == 'SN2011by':
		df = pd.read_csv('sndata_files/'+mysn[i]+'_uvotB15.1.dat',delim_whitespace=True,comment='#')
		df.columns=['Filter','MJD','Mag','Magerr','3SigErr','0.98Lim','Rate','RateErr','Ap','Frametime','Exp','Telapse']

		#Guess peak data, define filter subsets
		w2data = df.loc[(df.Filter == 'UVW2') & (np.isnan(df.Mag) == False )]
		m2data = df.loc[(df.Filter == 'UVM2') & (np.isnan(df.Mag) == False )]
		w1data = df.loc[(df.Filter == 'UVW1') & (np.isnan(df.Mag) == False )]
		udata = df.loc[(df.Filter == 'U') & (np.isnan(df.Mag) == False )]
		bdata = df.loc[(df.Filter == 'B') & (np.isnan(df.Mag) == False )]
		vdata = df.loc[(df.Filter == 'V') & (np.isnan(df.Mag) == False )]

		pkl = bdata.loc[(bdata.Mag == min(bdata.Mag))]
		pkl = np.asarray(pkl.MJD)[0]
		print('best guess for start date is '+str(pkl))
		fils = ['V','B','U','UVW1','UVM2','UVW2']
		for k in range(len(fils)):

			vres = sirfit(fils[k],mysn[i],pkl)
			outdf = outdf.append({'SNname':vres[0], 'Filter':vres[1], 'T0': vres[2], 'T0_err': vres[3], 'R0': vres[4], 'R0_err': vres[5]}, ignore_index=True)


# outdf.to_csv('sirah/compare/curvefit_small.csv',index=True)



