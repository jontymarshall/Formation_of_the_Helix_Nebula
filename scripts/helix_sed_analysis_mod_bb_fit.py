 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 11:25:55 2021

@author: jonty
"""

import numpy as np
from astropy.io import ascii,fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy import units as u
from skimage.transform import rescale
from scipy.ndimage.interpolation import shift
import matplotlib.pyplot as plt
import copy
#constants
h = 6.626e-34
c = 299792458.0 # m/s
k = 1.38e-23
sb = 5.67e-8 #
au     = 1.495978707e11 # m 
pc     = 3.0857e16 # m
lsol   = 3.828e26 # W
rsol   = 6.96342e8 # m
MEarth = 5.97237e24 # kg

um = 1e-6 #for wavelengths in microns

direc = '../'


#fit blackbody to SED
from astropy.modeling.physical_models import BlackBody as blackbody
from scipy.optimize import curve_fit
from scipy import interpolate

#Spectral energy distribution
#Read in data
#Flux in mJy
converters = {'Wavelength': [ascii.convert_numpy(np.float32)],
              'Flux': [ascii.convert_numpy(np.float32)],
              'Uncertainty': [ascii.convert_numpy(np.float32)]}

data = ascii.read(direc+'../data/photometry/helix_compact_photometry.csv',converters=converters,format='csv',comment='#')

obs_l = data['Wavelength'].data
obs_f = data['Flux'].data
obs_u = data['Uncertainty'].data

dstar = 201.0 #pc

#Read in white dwarf spectrum
converters = {'Lambda':[ascii.convert_numpy(np.float32)],
              'Fstar':[ascii.convert_numpy(np.float32)]}

star = ascii.read(direc+'../data/models/photosphere/helix_photosphere_v2.txt',converters=converters)

star_l = star['Lambda'].data*1e-4
star_f = 3.33564095E+04*star['Fstar'].data*(star['Lambda'].data**2)#*6e-21


#extrapolate WD spectrum
extra_l = np.logspace(np.log10(300), np.log10(10000),num=1000,endpoint=True)
extra_f = star_f[-1]*(star_l[-1]/extra_l)**2

star_l = np.append(star_l,extra_l)
star_f = np.append(star_f,extra_f)

#Read in Spitzer IRS spectrum
converters = {'wavelength':[ascii.convert_numpy(np.float32)],
              'flux':[ascii.convert_numpy(np.float32)],
              'error':[ascii.convert_numpy(np.float32)]}

#Flux in mJy
irs = ascii.read(direc+'../data/spectroscopy/spitzer/cassis/cassis_tbl_spcf_20677888t.tbl',converters=converters)

irs_l = irs['wavelength'].data
irs_f = irs['flux'].data*1e3
irs_u = irs['error'].data*1e3

#Scale WD SED to photometry
f = interpolate.interp1d(star_l,star_f,kind='linear')
a = np.where((obs_l < 8.) &(obs_f/obs_u > 3.))
star_f_int = f(obs_l)

def scale(a,y):
    return a*y

popt, pcov = curve_fit(scale, star_f_int[a], obs_f[a], sigma=obs_u[a],absolute_sigma=True)
photosphere_scaling_factor = popt[0]
star_f = photosphere_scaling_factor*star_f #in mJy

f = interpolate.interp1d(star_l,star_f,kind='linear')
photospheric_contribution = f(obs_l)
wavelengths = [24.,54.,70.,160.,250.,1300.]
phot = f(wavelengths)

#Add modified blackbody to SED
Tdust= 102.0 #Kelvin
beta =   0.0 
lam0 = 210.0 #um 

bb = blackbody(temperature=Tdust*u.K)
wav = star_l*1e-6*u.m
bb_f = bb(wav).value

a = np.where(star_l >= lam0)
bb_f[a] = bb_f[a]*(lam0/star_l[a])**beta
bb_f = (bb_f / np.max(bb_f))
bb_f_save = bb_f


a = np.where((obs_l > 8.0)&(obs_f/obs_u > 3))
f = interpolate.interp1d(star_l,bb_f,kind='linear')
blackbody_model = f(obs_l)


popt, pcov = curve_fit(scale, blackbody_model[a], obs_f[a] - photospheric_contribution[a], sigma=obs_u[a],absolute_sigma=True)

blackbody_scaling_factor = popt[0]
bb_f = blackbody_scaling_factor*bb_f

frac_lum = float(np.sum(bb_f)/np.sum(star_f))

#Compact component of the SED
lam_extend = np.asarray([70.])
flx_extend = np.asarray([222.])
unc_extend = np.asarray([12.]) 

lam_ul_e = np.asarray([160.,250.])
sig_ul_e = np.asarray([120.,60.])

lam_compact = np.asarray([8.0,24.0,70.0])
flx_compact = np.asarray([0.174,48.4,36.0])
unc_compact = np.asarray([0.017,7.3,4.5])

lam_ul_c = np.asarray([54.,160.,1300.])
sig_ul_c = np.asarray([96.,15.,8e-3])

#Plot up SED
fig, ax = plt.subplots(figsize=(8,6))
ax.set_yscale('log')
ax.set_xlim((3e-3,3e3))
ax.set_ylim((1e-4,3e3))
ax.set_xscale('log')
ax.set_xlabel(r'Wavelength ($\mu$m)',size=20)
ax.set_ylabel(r'Flux density (mJy)',size=20)
ax.plot(star_l,star_f,marker=None,linestyle='-',color='black')
ax.errorbar(irs_l, irs_f, xerr=None,yerr=irs_u,marker='.',linestyle=None,color='darkgray',alpha=0.1)


#set up grid
#radii     = np.linspace(50,400,num=11,endpoint=True)#[50.,100.,200.,400.]
#widths    = np.linspace(0.25,0.9,num=11,endpoint=True)#[0.25,0.5,0.75,0.9]
#lam_break = [25.,50.,75.,100.]
#betas     = [1.5,2.0,2.5,3.0,3.5]
#temps = np.linspace(50.,150.,num=101,endpoint=True)
csq_bf = 1e30

accept_lam0_c = []
accept_beta_c = []

accept_lam0_e = []
accept_beta_e = []
accept_temp_e = []

niter=0

posterior_flag_extended = 0
posterior_flag_compact  = 0

while niter <= 100000:
    temp = np.random.uniform(20.,60.)
    beta_e = np.random.uniform(0.0,2.0)
    lam0_e = np.random.uniform(70.0,200.0)
    bb = blackbody(temperature=temp*u.K)
    wav = star_l*1e-6*u.m
    bb_e = bb(wav).value
    lam0 = np.random.uniform(24.0,50.0)
    beta = np.random.uniform(0.0,4.0)
    
    #Compact component
    bb = blackbody(temperature=102*u.K)
    wav = star_l*1e-6*u.m
    bb_f = bb(wav).value
    disc_f = copy.deepcopy(bb_f)
    long = np.where(star_l > lam0)
    disc_f[long] = bb_f[long]*(lam0/star_l[long])**beta
    
    #print(disc_f)
    f = interpolate.interp1d(star_l,disc_f,kind='linear')
    fdisc = f(wavelengths)
    fscale = flx_compact[1] + unc_compact[1]*(2.0*np.random.uniform(0.0,1.0)-1.0)
    disc_f = (fscale/fdisc[1])*disc_f
    
    val_compact = np.asarray([fdisc[0],fdisc[1],fdisc[2]])

    f = interpolate.interp1d(star_l,disc_f,kind='linear')
    fdisc = f(wavelengths)
    
    alpha = 0.01
    #Extended component
    disc_e = copy.deepcopy(bb_e)
    long = np.where(star_l > lam0_e)
    disc_e[long] = bb_e[long]*(lam0_e/star_l[long])**beta_e
    f = interpolate.interp1d(star_l,disc_e,kind='linear')
    fdise = f(wavelengths)
    
    #print(disc_f)
    fscale = flx_extend[0]+unc_extend[0]*(2.0*np.random.uniform(0.0,1.0)-1.0)
    disc_e = (fscale/fdise[2])*disc_e
    
    f = interpolate.interp1d(star_l,disc_e,kind='linear')
    fdise = f(wavelengths)
    
    val_extend = np.asarray([fdise[2]])
    
    inbounds_c = False
    inbounds_e = False
    
    if flx_compact[2]-unc_compact[2] < fdisc[2] < flx_compact[2]+unc_compact[2] \
        and fdisc[1] < 3*sig_ul_c[0] \
        and fdisc[3] < 3*sig_ul_c[1] :
        #and flx_compact[1]-unc_compact[1] < fdisc[0] < flx_compact[1]+unc_compact[1] \
        
        
        inbounds_c = True
                
    if flx_extend[0]-unc_extend[0] < fdise[2] < flx_extend[0]+unc_extend[0] \
        and fdise[1] < 3*sig_ul_c[0] \
        and fdise[3] < 3*sig_ul_e[0] \
        and fdise[4] < 3*sig_ul_e[1] \
        and fdise[5] < 21*sig_ul_c[2] :

        inbounds_e = True
    
    #print(inbounds_c,inbounds_e)
    
    if inbounds_c == True and inbounds_e == True:
        csq = 0.0
        csq = np.sum((flx_compact[1:] - val_compact[1:])**2/unc_compact[1:]**2) + np.sum((flx_extend - val_extend)**2/unc_extend**2)
        print("Here",csq)
        if csq < csq_bf:
            csq_bf = csq
            disc_f_bf = disc_f
            disc_e_bf = disc_e
            best_fit_e = (temp,lam0_e,beta_e)
            best_fit_c = (lam0,beta)
        accept_lam0_e.append(lam0_e)
        accept_beta_e.append(beta_e)
        accept_temp_e.append(temp)    
        accept_beta_c.append(beta)
        accept_lam0_c.append(lam0)
        
    #plot X % of all models
    roulette = np.random.uniform()
    if roulette >= 0.95 and inbounds_c == True:
        if posterior_flag_compact == 0:
            ax.plot(star_l,1e-20*disc_f,marker=None,linestyle='-',color='purple',alpha=1,label='Compact component')
            posterior_flag_compact = 1
        else:
            ax.plot(star_l,disc_f,marker=None,linestyle='-',color='purple',alpha=alpha)
    
    if roulette <= 0.05 and inbounds_e == True:
        if posterior_flag_extended == 0:
            ax.plot(star_l,1e-20*disc_e,marker=None,linestyle='-',color='dodgerblue',alpha=1,label='Extended component')
            posterior_flag_extended = 1
        else: 
            ax.plot(star_l,disc_e,marker=None,linestyle='-',color='dodgerblue',alpha=alpha)
    niter+=1

ax.plot(star_l,disc_f_bf,marker=None,linestyle='--',color='black',alpha=1.0,label='Compact component (max. prob.)')
ax.plot(star_l,disc_e_bf,marker=None,linestyle=':',color='black',alpha=1.0,label='Extended component (max. prob.)')
disc_t_bf = np.asarray(star_f + disc_f_bf + disc_e_bf)
ax.plot(star_l,star_f+disc_f_bf+disc_e_bf,marker=None,linestyle='-',color='black',alpha=1.0,label='Total emission')


a = np.where(obs_l < 24.)
ax.errorbar(obs_l[a],obs_f[a],yerr=obs_u[a],marker='o',mfc='darkgray',mec='black',linestyle='',label='Ancilllary photometry')

ax.errorbar(lam_compact,flx_compact,xerr=None,yerr=unc_compact,linestyle='',marker='o',mec='black',mfc='white',label='Compact photometry')
ax.plot(lam_ul_c,3*sig_ul_c,marker='v',linestyle='',mfc='white',mec='black')

ax.errorbar(lam_extend,flx_extend,xerr=None,yerr=unc_extend,linestyle='',marker='o',color='black',label='Extended photometry')
ax.plot(lam_ul_e,3*sig_ul_e,marker='v',linestyle='',color='black')
#ax.plot(lam_ul_c[2],21*sig_ul_c[2],marker='v',linestyle='',color='black')
ax.plot(lam_compact[0],3*unc_compact[0],marker='v',linestyle='',color='black')



ax.legend(loc='lower left')

fig.savefig('sed_helix_two_component_modbb.pdf',dpi=200,bbox_inches="tight")

print(np.min(accept_lam0_c),np.max(accept_lam0_c))
print(np.min(accept_beta_c),np.max(accept_beta_c))
print(np.min(accept_lam0_e),np.max(accept_lam0_e))
print(np.min(accept_beta_e),np.max(accept_beta_e))
print(np.min(accept_temp_e),np.max(accept_temp_e))

print(best_fit_c)
print(best_fit_e)

from astropy.table import Table

t = Table([star_l, star_f, disc_f_bf, disc_e_bf, disc_t_bf], names=('Wavelength', 'Star', 'Compact', 'Extended', 'Total'))
t.meta['comments'] = ['Wavelength in microns.','Flux density in milliJanskys.']
filename = direc+'helix_best_fit_sed_table.csv'
t.write(filename, format='csv',overwrite=True,comment='#')

#Interpolate best fit model to IRS spectrum wavelengths, then subtract
sort = np.argsort(irs_l)

irs_l = irs_l[sort]
irs_f = irs_f[sort]
irs_u = irs_u[sort]

f = interpolate.interp1d(star_l,disc_t_bf)
model_fluxes = f(irs_l)
difference = irs_f - model_fluxes

fig, ax = plt.subplots(figsize=(8,6))
ax.set_yscale('linear')
ax.set_xlim((5.0,35.0))
ax.set_ylim((-10,10))
ax.set_xscale('linear')
ax.set_xlabel(r'Wavelength ($\mu$m)',size=20)
ax.set_ylabel(r'Flux density (mJy)',size=20)
ax.plot(irs_l,difference/irs_u,marker=None,linestyle='-',color='black')
#ax.errorbar(irs_l, irs_f, xerr=None,yerr=irs_u,marker='.',linestyle=None,color='darkgray',alpha=0.5)

fig.savefig(direc+'irs_helix_model_subtract.pdf',dpi=200,bbox_inches="tight")

radius = 100.
width  =  25.

#Plot up the observations used in the analysis
theta = (phot[0],fdisc[0],radius,width,0.0,0.0)
spitzer_model = make_disc(dstar, theta, pixscale=2.45, npix= 61)
spitzer_synth = convolve(spitzer_model,kernel_24um,normalize_kernel=True)

theta = (phot[1],fdisc[1],radius,width,0.0,0.0)
sofia_model = make_disc(dstar, theta, pixscale=1.0, npix= 61)
sofia_synth = convolve(sofia_model,kernel_54um,normalize_kernel=True)*1e-3              

theta = (phot[2],fdisc[2],radius,width,0.0,0.0)
herschel_model = make_disc(dstar, theta, pixscale=1.6, npix= 61)
herschel_synth = convolve(herschel_model,kernel_70um,normalize_kernel=True)*1e-3              

theta = (phot[3],fdisc[3],radius,width,0.0,0.0)
alma_model = make_disc(dstar, theta, pixscale=0.04, npix= 361)*1e-3*alma_beam/(0.04**2)
alma_synth = convolve(alma_model,kernel_13mm,normalize_kernel=True)
                                
fig, axes = plt.subplots(nrows=4, ncols=4)
ax = axes.ravel()

#Spitzer
dn = 30
ax[0].imshow(spitzer_obsvn, cmap='rainbow',vmin=np.min(spitzer_obsvn),vmax=np.max(spitzer_obsvn))
ax[0].set_title("Observation")

spitzer_model[dn,dn] = 0.0
ax[1].imshow(spitzer_model, cmap='rainbow',vmin=np.min(spitzer_model),vmax=np.max(spitzer_model))
ax[1].set_title("Raw Model")
                
ax[2].imshow(spitzer_synth, cmap='rainbow',vmin=np.min(spitzer_obsvn),vmax=np.max(spitzer_obsvn))
ax[2].set_title("Convolved Model")    

ax[3].imshow(spitzer_obsvn - spitzer_synth, cmap='rainbow',vmin=np.min(spitzer_obsvn),vmax=np.max(spitzer_obsvn))
ax[3].set_title("Residual")

#SOFIA
ax[4].imshow(sofia_obsvn, cmap='rainbow',vmin=np.min(sofia_obsvn),vmax=np.max(sofia_obsvn))

sofia_model[dn,dn] = 0.0
ax[5].imshow(sofia_model, cmap='rainbow',vmin=np.min(sofia_model),vmax=np.max(sofia_model))
                
ax[6].imshow(sofia_synth, cmap='rainbow',vmin=np.min(sofia_obsvn),vmax=np.max(sofia_obsvn))

ax[7].imshow(sofia_obsvn - sofia_synth, cmap='rainbow',vmin=np.min(sofia_obsvn),vmax=np.max(sofia_obsvn))

#Herschel
ax[8].imshow(herschel_obsvn, cmap='rainbow',vmin=np.min(herschel_obsvn),vmax=np.max(herschel_obsvn))

herschel_model[dn,dn] = 0.0
ax[9].imshow(herschel_model, cmap='rainbow',vmin=np.min(herschel_model),vmax=np.max(herschel_model))
                
ax[10].imshow(spitzer_synth, cmap='rainbow',vmin=np.min(herschel_obsvn),vmax=np.max(herschel_obsvn))

ax[11].imshow(herschel_obsvn - herschel_synth, cmap='rainbow',vmin=np.min(herschel_obsvn),vmax=np.max(herschel_obsvn))

#ALMA
dn = 120
ax[12].imshow(alma_obsvn, cmap='rainbow',vmin=np.min(alma_obsvn),vmax=np.max(alma_obsvn))

alma_model[dn,dn] = 0.0
ax[13].imshow(alma_model, cmap='rainbow',vmin=np.min(alma_model),vmax=np.max(alma_model))
                
ax[14].imshow(alma_synth, cmap='rainbow',vmin=np.min(alma_obsvn),vmax=np.max(alma_obsvn))

ax[15].imshow(alma_obsvn - alma_synth, cmap='rainbow',vmin=np.min(alma_obsvn),vmax=np.max(alma_obsvn))

fig.savefig(direc+'Helix_ring_model_'+str(lam0)+'_'+str(beta)+'.png',dpi=200)
plt.close()