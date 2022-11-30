#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 14:22:20 2022

@author: jonty
"""

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from astropy.io import fits
from matplotlib.patches import Ellipse

plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams["figure.autolayout"] = True


direc = '/Users/jonty/mydata/fearless/helix/image_grid/'

#Read in Herschel PACS image
hdul = fits.open(direc+'../herschel/Archive/helix_ksu_hppjsmapb.fits')
herschel_obsvn  = hdul[1].data
cy = 128
cx = 164
dn = 30
herschel_obsvn = 1e3*herschel_obsvn[cx-dn:cx+dn+1,cy-dn:cy+dn+1] #in Jy/pixel

#2D Gaussian model
def makeGaussian2(x_center=0, y_center=0, theta=0, sigma=[5.4,5.4], size=[61,61]):
    # x_center and y_center will be the center of the gaussian, theta will be the rotation angle
    # sigma_x and sigma_y will be the stdevs in the x and y axis before rotation
    # x_size and y_size give the size of the frame 
    
    if len(sigma) == 2:
        sigma_x = sigma[0]
        sigma_y = sigma[1]
    else:
        sigma_x = sigma
        sigma_y = sigma
    
    if len(size) == 2:
        x_size = size[0]
        y_size = size[1]
    else:
        x_size = size
        y_size = size
    
    theta = 2*np.pi*theta/360
    x = np.arange(0,x_size, 1, float)
    y = np.arange(0,y_size, 1, float)
    y = y[:,np.newaxis]
    sx = sigma_x
    sy = sigma_y
    x0 = x_center
    y0 = y_center

    # rotation
    a=np.cos(theta)*x -np.sin(theta)*y
    b=np.sin(theta)*x +np.cos(theta)*y
    a0=np.cos(theta)*x0 -np.sin(theta)*y0
    b0=np.sin(theta)*x0 +np.cos(theta)*y0

    return np.exp(-(((a-a0)**2)/(2*(sx**2)) + ((b-b0)**2) /(2*(sy**2))))

#Make 2D Gaussian model for extended component 
herschel_extmd = makeGaussian2(x_center=26, y_center=34, theta=25, sigma=[28/2.355,24/2.355], size=[61,61])
herschel_extsb = herschel_obsvn - 222*(herschel_extmd/np.sum(herschel_extmd))

#Make 2D Gaussian model for compasct component
herschel_commd = makeGaussian2(x_center=30, y_center=30, theta=0, sigma=[5.4/2.355,5.4/2.355], size=[61,61])
herschel_comsb = herschel_extsb - (36/1.6**2)*(herschel_commd/np.sum(herschel_commd))

#Plot 3x2 figure
fig,ax = plt.subplots(2,3)#,sharex=True,sharey=True)

im = ax[0,0].imshow(herschel_obsvn,vmin=np.min(herschel_obsvn),vmax=np.max(herschel_obsvn),cmap='inferno',origin='lower')
xsu = 10.
ax[0,0].set_xticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[0,0].set_xticklabels(['+45','+30','+15','0','-15','-30','-45'],fontsize=12)
ax[0,0].set_yticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[0,0].set_yticklabels(['-45','-30','-15','0','+15','+30','+45'],fontsize=12)
ax[0,0].set_ylabel(r'North Offset (")',fontsize=16)

axins1 = inset_axes(ax[0,0],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    borderpad=-1.5,
                    loc='upper center')

cb = fig.colorbar(im,cax=axins1,orientation="horizontal",ticklocation='top')
cb.set_label(r'Intensity (mJy/pixel)', size=10)

im = ax[0,1].imshow(herschel_extmd, cmap='inferno',vmin=0.0,vmax=1.0,origin='lower')
ax[0,1].set_xticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[0,1].set_xticklabels(['+45','+30','+15','0','-15','-30','-45'],fontsize=12)
ax[0,1].set_yticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[0,1].set_yticklabels(['-45','-30','-15','0','+15','+30','+45'],fontsize=12)

axins1 = inset_axes(ax[0,1],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    borderpad=-1.5,
                    loc='upper center')

cb = fig.colorbar(mappable=im,cax=axins1,orientation="horizontal",ticklocation='top')
cb.set_label(r'Norm. intensity (arb. units)', size=10)

im = ax[0,2].imshow(herschel_extsb, cmap='inferno',vmin=np.min(herschel_obsvn),vmax=np.max(herschel_obsvn),origin='lower')
ax[0,2].set_xticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[0,2].set_xticklabels(['+45','+30','+15','0','-15','-30','-45'],fontsize=12)
ax[0,2].set_yticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[0,2].set_yticklabels(['-45','-30','-15','0','+15','+30','+45'],fontsize=12)

axins1 = inset_axes(ax[0,2],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    borderpad=-1.5,
                    loc='upper center')

cb = fig.colorbar(mappable=im,cax=axins1,orientation="horizontal",ticklocation='top')
cb.set_label(r'Intensity (mJy/pixel)', size=10)

im = ax[1,0].imshow(herschel_extsb, cmap='inferno',vmin=np.min(herschel_obsvn),vmax=np.max(herschel_obsvn),origin='lower')
ax[1,0].set_xticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[1,0].set_xticklabels(['+45','+30','+15','0','-15','-30','-45'],fontsize=12)
ax[1,0].set_yticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[1,0].set_yticklabels(['-45','-30','-15','0','+15','+30','+45'],fontsize=12)
ax[1,0].set_xlabel(r'East Offset (")',fontsize=16)
ax[1,0].set_ylabel(r'North Offset (")',fontsize=16)

axins1 = inset_axes(ax[1,0],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    borderpad=-1.5,
                    loc='upper center')

cb = fig.colorbar(mappable=im,cax=axins1,orientation="horizontal",ticklocation='top')
cb.set_label(r'Intensity (mJy/pixel)', size=10)

im = ax[1,1].imshow(herschel_commd, cmap='inferno',vmin=0.0,vmax=1.0,origin='lower')
ax[1,1].set_xticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[1,1].set_xticklabels(['+45','+30','+15','0','-15','-30','-45'],fontsize=12)
ax[1,1].set_yticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[1,1].set_yticklabels(['-45','-30','-15','0','+15','+30','+45'],fontsize=12)
ax[1,1].set_xlabel(r'East Offset (")',fontsize=16)

axins1 = inset_axes(ax[1,1],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    borderpad=-1.5,
                    loc='upper center')

cb = fig.colorbar(mappable=im,cax=axins1,orientation="horizontal",ticklocation='top')
cb.set_label(r'Norm. intensity (arb. units)', size=10)


im = ax[1,2].imshow(herschel_comsb, cmap='inferno',vmin=np.min(herschel_obsvn),vmax=np.max(herschel_obsvn),origin='lower')
ax[1,2].set_xticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[1,2].set_xticklabels(['+45','+30','+15','0','-15','-30','-45'],fontsize=12)
ax[1,2].set_yticks([0,xsu,2*xsu,3*xsu,4*xsu,5*xsu,6*xsu])
ax[1,2].set_yticklabels(['-45','-30','-15','0','+15','+30','+45'],fontsize=12)
ax[1,2].set_xlabel(r'East Offset (")',fontsize=16)

axins1 = inset_axes(ax[1,2],
                    width="100%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    borderpad=-1.5,
                    loc='upper center')

ellipse = Ellipse((26,34),28,24,angle=-25,color='white',alpha=0.25)
ax[1,2].add_patch(ellipse)
ellipse = Ellipse((30,30),5.4/1.6,5.4/1.6,angle=0,color='white',alpha=0.5)
ax[1,2].add_patch(ellipse)

cb = fig.colorbar(mappable=im,cax=axins1,orientation="horizontal",ticklocation='top')
cb.set_label(r'Intensity (mJy/pixel)', size=10)
#plt.tight_layout()
fig.savefig(direc+'Helix_Nebula_Herschel_PACS70_Two_Component_Model.pdf',dpi=200)
plt.show()
plt.close()
