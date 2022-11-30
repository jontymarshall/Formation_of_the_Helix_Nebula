#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 18:45:04 2022

@author: jonty
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from RT_Code import RTModel

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

#Read in compact component photometry from csv file
converters = {'Wavelength': [ascii.convert_numpy(np.float32)],
              'Flux': [ascii.convert_numpy(np.float32)],
              'Uncertainty': [ascii.convert_numpy(np.float32)]}

direc = '../'

data = ascii.read(direc+'data/photometry/helix_compact_photometry.csv',converters=converters,format='csv')

obs_l = data['Wavelength'].data
obs_f = data['Flux'].data
obs_u = data['Uncertainty'].data

ul_wave = [54,160,250,1300]
ul_flux = [111,45,24,0.180]

#Small grains, broad size distribution, steady state, compact disc
#Create disk model
disk = RTModel()
RTModel.get_parameters(disk,'Helix_RTModel_Input_File.txt')
disk.obs_wave = obs_l
disk.obs_flux = obs_f
disk.obs_uncs = obs_u

disk.parameters['directory'] = direc
disk.parameters['prefix']    = 'Helix_SED_Model_A_'

#Plot parameters
disk.parameters['lmin'] = 0.01

#Dust
disk.parameters['dtype'] = 'onepl'
disk.parameters['mdust'] = 0.008
disk.parameters['rin']  = 60.
disk.parameters['rout'] = 80.
disk.parameters['alpha_out'] = 2.0
disk.parameters['amin'] = 2.0
disk.parameters['amax'] = 20.
disk.parameters['q'] = -3.5

RTModel.make_star(disk)
#RTModel.scale_star(disk)
RTModel.make_dust(disk)
RTModel.make_disc(disk)
RTModel.calculate_surface_density(disk)
RTModel.read_optical_constants(disk)
RTModel.calculate_qabs(disk)
RTModel.calculate_dust_emission(disk,mode='full',tolerance=1e-3)
RTModel.calculate_dust_scatter(disk)
RTModel.flam_to_fnu(disk)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.loglog(disk.sed_wave, disk.sed_emit, color='black',linestyle='--')
ax.loglog(disk.sed_wave, disk.sed_star, color='black',linestyle='-')
ax.loglog(disk.sed_wave, disk.sed_star + disk.sed_emit, color='black',linestyle='-')

ax.errorbar(disk.obs_wave,disk.obs_flux,yerr=disk.obs_uncs,marker='o',mec='blue',mfc='dodgerblue',linestyle='')
ax.plot(ul_wave,ul_flux,marker='v',mec='black',mfc='white',linestyle='')
ax.set_xlabel(r'Wavelength ($\mu$m)')
ax.set_ylabel(r'Flux density (mJy)')

ax.set_xlim(disk.parameters["lmin"],disk.parameters["lmax"])
ax.set_ylim(3e-2,1e3)

plt.annotate('Model A',(200,300),fontsize=14)
plt.show()
fig.savefig(disk.parameters['directory']+disk.parameters['prefix']+'_sed.pdf',dpi=200)
plt.close(fig)

#Blowout grains, broad size distribution, steady state slope, compact disc
#Create disk model
disk = RTModel()
RTModel.get_parameters(disk,'Helix_RTModel_Input_File.txt')
disk.obs_wave = obs_l
disk.obs_flux = obs_f
disk.obs_uncs = obs_u

disk.parameters['directory'] = direc
disk.parameters['prefix'] = 'Helix_SED_Model_B_'
#Plot parameters
disk.parameters['lmin'] = 0.01

#Dust
disk.parameters['dtype'] = 'onepl'
disk.parameters['mdust'] = 0.01
disk.parameters['rin']  = 30.
disk.parameters['rout'] = 50.
disk.parameters['alpha_out'] = -1.0
disk.parameters['amin']  = 60.0
disk.parameters['amax'] =  1000.
disk.parameters['q'] = -3.5

RTModel.make_star(disk)
#RTModel.scale_star(disk)
RTModel.make_dust(disk)
RTModel.make_disc(disk)
RTModel.calculate_surface_density(disk)
RTModel.read_optical_constants(disk)
RTModel.calculate_qabs(disk)
RTModel.calculate_dust_emission(disk,mode='full',tolerance=1e-3)
RTModel.calculate_dust_scatter(disk)
RTModel.flam_to_fnu(disk)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.loglog(disk.sed_wave, disk.sed_emit, color='black',linestyle='--')
ax.loglog(disk.sed_wave, disk.sed_star, color='black',linestyle='-')
ax.loglog(disk.sed_wave, disk.sed_star + disk.sed_emit, color='black',linestyle='-')

ax.errorbar(disk.obs_wave,disk.obs_flux,yerr=disk.obs_uncs,marker='o',mec='blue',mfc='dodgerblue',linestyle='')
ax.plot(ul_wave,ul_flux,marker='v',mec='black',mfc='white',linestyle='')
ax.set_xlabel(r'Wavelength ($\mu$m)')
ax.set_ylabel(r'Flux density (mJy)')

ax.set_xlim(disk.parameters["lmin"],disk.parameters["lmax"])
ax.set_ylim(3e-2,1e3)

plt.annotate('Model B',(200,300),fontsize=14)
plt.show()
fig.savefig(disk.parameters['directory']+disk.parameters['prefix']+'_sed.pdf',dpi=200)
plt.close(fig)

#Large grains, broad size distribution, steady state slope, large disc 
#Create disk model
disk = RTModel()
RTModel.get_parameters(disk,'Helix_RTModel_Input_File.txt')
disk.obs_wave = obs_l
disk.obs_flux = obs_f
disk.obs_uncs = obs_u

disk.parameters['directory'] = direc
disk.parameters['prefix'] = 'Helix_SED_Model_C_'
#Plot parameters
disk.parameters['lmin'] = 0.01

#Dust
disk.parameters['dtype'] = 'onepl'
disk.parameters['mdust'] = 0.01
disk.parameters['rin']  = 30.
disk.parameters['rout'] = 90.
disk.parameters['alpha_out'] = -1.0
disk.parameters['amin']  = 60.0
disk.parameters['amax'] =  1000.
disk.parameters['q'] = -3.5

RTModel.make_star(disk)
#RTModel.scale_star(disk)
RTModel.make_dust(disk)
RTModel.make_disc(disk)
RTModel.calculate_surface_density(disk)
RTModel.read_optical_constants(disk)
RTModel.calculate_qabs(disk)
RTModel.calculate_dust_emission(disk,mode='full',tolerance=1e-3)
RTModel.calculate_dust_scatter(disk)
RTModel.flam_to_fnu(disk)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.loglog(disk.sed_wave, disk.sed_emit, color='black',linestyle='--')
ax.loglog(disk.sed_wave, disk.sed_star, color='black',linestyle='-')
ax.loglog(disk.sed_wave, disk.sed_star + disk.sed_emit, color='black',linestyle='-')

ax.errorbar(disk.obs_wave,disk.obs_flux,yerr=disk.obs_uncs,marker='o',mec='blue',mfc='dodgerblue',linestyle='')
ax.plot(ul_wave,ul_flux,marker='v',mec='black',mfc='white',linestyle='')
ax.set_xlabel(r'Wavelength ($\mu$m)')
ax.set_ylabel(r'Flux density (mJy)')

ax.set_xlim(disk.parameters["lmin"],disk.parameters["lmax"])
ax.set_ylim(3e-2,1e3)

plt.annotate('Model C',(200,300),fontsize=14)
plt.show()
fig.savefig(disk.parameters['directory']+disk.parameters['prefix']+'_sed.pdf',dpi=200)
plt.close(fig)

#Middle grains, broad size distribution, steep slope, broad disc 
#Create disk model
disk = RTModel()
RTModel.get_parameters(disk,'Helix_RTModel_Input_File.txt')
disk.obs_wave = obs_l
disk.obs_flux = obs_f
disk.obs_uncs = obs_u

disk.parameters['directory'] = direc
disk.parameters['prefix'] = 'Helix_SED_Model_D_'
#Plot parameters
disk.parameters['lmin'] = 0.01

#Dust
disk.parameters['dtype'] = 'onepl'
disk.parameters['mdust'] = 0.01
disk.parameters['rin']  = 30
disk.parameters['rout'] = 90.
disk.parameters['alpha_out'] = -1.0
disk.parameters['amin']  = 10.0
disk.parameters['amax'] =  1000.
disk.parameters['q'] = -3.5

RTModel.make_star(disk)
#RTModel.scale_star(disk)
RTModel.make_dust(disk)
RTModel.make_disc(disk)
RTModel.calculate_surface_density(disk)
RTModel.read_optical_constants(disk)
RTModel.calculate_qabs(disk)
RTModel.calculate_dust_emission(disk,mode='full',tolerance=1e-3)
RTModel.calculate_dust_scatter(disk)
RTModel.flam_to_fnu(disk)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.loglog(disk.sed_wave, disk.sed_emit, color='black',linestyle='--')
ax.loglog(disk.sed_wave, disk.sed_star, color='black',linestyle='-')
ax.loglog(disk.sed_wave, disk.sed_star + disk.sed_emit, color='black',linestyle='-')

ax.errorbar(disk.obs_wave,disk.obs_flux,yerr=disk.obs_uncs,marker='o',mec='blue',mfc='dodgerblue',linestyle='')
ax.plot(ul_wave,ul_flux,marker='v',mec='black',mfc='white',linestyle='')
ax.set_xlabel(r'Wavelength ($\mu$m)')
ax.set_ylabel(r'Flux density (mJy)')

ax.set_xlim(disk.parameters["lmin"],disk.parameters["lmax"])
ax.set_ylim(3e-2,1e3)

plt.annotate('Model D',(200,300),fontsize=14)
plt.show()
fig.savefig(disk.parameters['directory']+disk.parameters['prefix']+'_sed.pdf',dpi=200)
plt.close(fig)
