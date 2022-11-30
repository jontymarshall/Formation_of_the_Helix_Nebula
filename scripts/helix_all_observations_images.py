#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 14:22:20 2022

@author: jonty
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from astropy.io import fits

direc = '/Users/jonty/mydata/fearless/helix/image_grid/'

#Read in big Herschel PACS image
plt.rcParams["figure.figsize"] = [8,10]
plt.rcParams["figure.autolayout"] = True


fig,ax = plt.subplots(2,2)

#SPITZER
fwhm_maj = 6.0
fwhm_min = 6.0
pixscale = 2.45
beam_ang = 0.0

hdul = fits.open(direc+'../spitzer/mips/r4621056/ch1/pbcd/SPITZER_M1_4621056_0000_10_E6427898_maic.fits')
spitzer_obsvn  = hdul[0].data
cy = 509
cx = 669
dn = 18
spitzer_obsvn = spitzer_obsvn[cx-dn:cx+dn+1,cy-dn:cy+dn+1] #in Jy/pixel
spitzer_obsvn *= 23.5045 * 1e-3 #in mJy/arcsec**2

xsu = (2*dn)/6
im = ax[0,0].imshow(spitzer_obsvn,vmin=1.,vmax=1.5, cmap='inferno',origin='lower')
ax[0,0].set_xticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[0,0].set_xticklabels(['+45','+30','+15','0','-15','-30','-45'],fontsize=12)
ax[0,0].set_yticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[0,0].set_yticklabels(['-45','-30','-15','0','+15','+30','+45'],fontsize=12)
#ax[0,0].set_xlabel(r'East Offset ($^{\prime}$)',fontsize=16)
ax[0,0].set_ylabel(r'North Offset ($^{\prime\prime}$)',fontsize=16)
axins1 = inset_axes(ax[0,0],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    borderpad=-1.5,
                    loc='upper center')

ellipse = Ellipse((15/pixscale,15/pixscale),fwhm_maj/pixscale,fwhm_min/pixscale,angle=beam_ang,color='white',alpha=0.8)
ax[0,0].add_patch(ellipse)

cb = fig.colorbar(mappable=im,cax=axins1,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux density (mJy/arcsec$^{2}$)', size=10)

#SOFIA
fwhm_maj = 3.0
fwhm_min = 3.0
pixscale = 1.0
beam_ang = 0.0

hdul = fits.open(direc+'../steve_data/sofia_2017-10-17_HA_F440/p5245/F0440_HA_IMA_0500543_HAWAHWPOpen_CRH_035-053.fits')
sofia_obsvn  = hdul[0].data
cy = 109
cx = 100
dn = 45
sofia_obsvn = sofia_obsvn[cx-dn:cx+dn+1,cy-dn:cy+dn+1] #in Jy/pixel
sofia_obsvn *= 1e3  #in mJy/arcsec**2

xsu = (2*dn)/6
im = ax[0,1].imshow(sofia_obsvn,vmin=-2,vmax=2, cmap='inferno',origin='lower')
ax[0,1].set_xticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[0,1].set_xticklabels(['+45','+30','+15','0','-15','-30','-45'],fontsize=12)
ax[0,1].set_yticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[0,1].set_yticklabels(['-45','-30','-15','0','+15','+30','+45'],fontsize=12)
#ax[0,0].set_xlabel(r'East Offset ($^{\prime}$)',fontsize=16)
#ax[0,0].set_ylabel(r'North Offset ($^{\prime}$)',fontsize=16)
axins1 = inset_axes(ax[0,1],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    borderpad=-1.5,
                    loc='upper center')

ellipse = Ellipse((15/pixscale,15/pixscale),fwhm_maj/pixscale,fwhm_min/pixscale,angle=beam_ang,color='white',alpha=0.8)
ax[0,1].add_patch(ellipse)

cb = fig.colorbar(mappable=im,cax=axins1,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux density (mJy/arcsec$^{2}$)', size=10)


#HERSCHEL
fwhm_maj = 5.4
fwhm_min = 5.4
pixscale = 1.6
beam_ang = 0.0

hdul = fits.open(direc+'../herschel/Archive/helix_ksu_hppjsmapb.fits')
herschel_obsvn  = hdul[1].data
cy = 128
cx = 163
dn = 30
herschel_obsvn = herschel_obsvn[cx-dn:cx+dn+1,cy-dn:cy+dn+1] #in Jy/pixel
herschel_obsvn *= 1e3 / 1.6**2 #in mJy/arcsec**2

xsu = (2*dn)/6
im = ax[1,0].imshow(herschel_obsvn,vmin=-0.1,vmax=0.5, cmap='inferno',origin='lower')
ax[1,0].set_xticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[1,0].set_xticklabels(['+45','+30','+15','0','-15','-30','-45'],fontsize=12)
ax[1,0].set_yticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[1,0].set_yticklabels(['-45','-30','-15','0','+15','+30','+45'],fontsize=12)
ax[1,0].set_xlabel(r'East Offset ($^{\prime\prime}$)',fontsize=16)
ax[1,0].set_ylabel(r'North Offset ($^{\prime\prime}$)',fontsize=16)
axins1 = inset_axes(ax[1,0],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    borderpad=-1.5,
                    loc='upper center')

ellipse = Ellipse((15/pixscale,15/pixscale),fwhm_maj/pixscale,fwhm_min/pixscale,angle=beam_ang,color='white',alpha=0.8)
ax[1,0].add_patch(ellipse)

cb = fig.colorbar(mappable=im,cax=axins1,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux density (mJy/arcsec$^{2}$)', size=10)

#ALMA
fwhm_maj = 0.257
fwhm_min = 0.229
pixscale = 0.040
beam_ang = 8.021694946289E+01
beam_area = 4.*(fwhm_maj*fwhm_min)/(2*np.log(2))
hdul = fits.open(direc+'../alma/WD2226-2_a_06_TE_cont.image.pbcor.fits')
alma_obsvn  = hdul[0].data
cy = 654
cx = 656
dn = 112
alma_obsvn = alma_obsvn[0,0,cx-dn:cx+dn+1,cy-dn:cy+dn+1] #in mJy/beam
alma_obsvn *= 1e3 * beam_area / pixscale**2 #in mJy/arcsec**2

xsu = (2*dn)/6
im = ax[1,1].imshow(alma_obsvn,vmin=-2,vmax=2, cmap='inferno',origin='lower')
ax[1,1].set_xticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[1,1].set_xticklabels(['+4.5','+3.0','+1.5','0','-1.5','-3.0','-4.5'],fontsize=12)
ax[1,1].set_yticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[1,1].set_yticklabels(['-4.5','-3.0','-1.5','0','+1.5','+3.0','+4.5'],fontsize=12)
ax[1,1].set_xlabel(r'East Offset ($^{\prime\prime}$)',fontsize=16)
#ax[1,0].set_ylabel(r'North Offset ($^{\prime}$)',fontsize=16)
axins1 = inset_axes(ax[1,1],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    borderpad=-1.5,
                    loc='upper center')

ellipse = Ellipse((1.5/pixscale,1.5/pixscale),fwhm_maj/pixscale,fwhm_min/pixscale,angle=beam_ang,color='white',alpha=0.8)
ax[1,1].add_patch(ellipse)

cb = fig.colorbar(mappable=im,cax=axins1,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux density (mJy/arcsec$^{2}$)', size=10)


plt.show()
plt.close()
fig.savefig(direc+'Helix_Nebula_All_Observations.pdf',dpi=200)
