#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 14:22:20 2022

@author: jonty
"""

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as patches
from matplotlib.patches import Ellipse
import numpy as np
from astropy.io import fits
import scipy as sp
from scipy.ndimage import filters

direc = '/Users/jonty/mydata/fearless/helix/image_grid/'

#Read in big Herschel PACS image
plt.rcParams["figure.figsize"] = [8,12]
plt.rcParams["figure.autolayout"] = False

pixscale = 1.6
hdul = fits.open(direc+'../herschel/Archive/helix_kpgt_hppjsmapb.fits')
herschel_obsvn2  = hdul[1].data
cy = 909
cx = 1213
dn = 187
herschel_obsvn2 = herschel_obsvn2[cx-dn:cx+dn+1,cy-dn:cy+dn+1] #in Jy/pixel
herschel_obsvn2 *= 1e3 / 1.6**2 #in mJy/arcsec**2


herschel_obsvn2 = filters.gaussian_filter(herschel_obsvn2, 1.6875, mode='constant')

fig,ax = plt.subplots()

xsu = 93.6
im = ax.imshow(herschel_obsvn2,vmin=-2,vmax=1, cmap='inferno',origin='lower')
ax.set_xticks([0,xsu,2*xsu,3*xsu,4*xsu])
ax.set_xticklabels(['+10','+5','0','-5','-10'],fontsize=12)
ax.set_yticks([0,xsu,2*xsu,3*xsu,4*xsu])
ax.set_yticklabels(['-10','-5','0','+5','+10'],fontsize=12)
ax.set_xlabel(r'East Offset ($^{\prime}$)',fontsize=16)
ax.set_ylabel(r'North Offset ($^{\prime}$)',fontsize=16)

#Herschel PACS 70um
fwhm_maj = 5.4
fwhm_min = 5.4
pixscale = 1.6
beam_ang = 0.0

ellipse = Ellipse((75/pixscale,75/pixscale),fwhm_maj/pixscale,fwhm_min/pixscale,angle=beam_ang,color='white',alpha=0.8)
ax.add_patch(ellipse)


rect = patches.Rectangle((159.375,159.375), 56.25, 56.25, linewidth=2, edgecolor='white', facecolor='none')
ax.add_patch(rect)

rect = patches.Rectangle((184.1875,184.1875), 5.625, 5.625, linewidth=2, edgecolor='white', facecolor='none')
ax.add_patch(rect)

#ax.plot(187,187,'*',mfc='white',mec='white',markersize=12,alpha=0.5)

axins2 = inset_axes(ax,
                    width="100%",  # width = 50% of parent_bbox width
                    height="3%",  # height : 5%
                    borderpad=-2.5,
                    loc='upper center')

cb = fig.colorbar(mappable=im,cax=axins2,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux density (mJy/arcsec$^{2}$)', size=10)
plt.show()
plt.close()

fig.savefig(direc+'Helix_Nebula_Big_Observation.pdf',dpi=200)
