pro helix_analysis

;procedure to examine extent and structure of central source in Helix Nebula

dir = '/Users/jonty/ownCloud/helix/'
;read in fits images:

;Kate Su image - 1.6"/pixel, cx = 128, cy = 164
image_1 = mrdfits(dir+'helix_ksu_hppjsmapb.fits','image',img1_head)
image_1_unc = mrdfits(dir+'helix_ksu_hppjsmapb.fits','stDev',unc1_head)
;KPGT image - 3.2"/pixel, cx = 909, cy = 1213
image_2 = mrdfits(dir+'helix_kpgt_hppjsmapb.fits','image',img2_head)
image_2_unc = mrdfits(dir+'helix_kpgt_hppjsmapb.fits','stDev',unc2_head)

;convolve images with SPIRE beam cf GV van de Steene et al., 2015, A&A 574, A134
psf_1 = psf_gaussian(npix=128.,fwhm=11.25d0,/normalize,ndimen=2)
image_1_conv = convolve(image_1,psf_1)
psf_2 = psf_gaussian(npix=128.,fwhm=5.55556d0,/normalize,ndimen=2)
image_2_conv = convolve(image_2,psf_2)

;centre image on WD2226-210 position: 22 29 38.541 -20 50 13.64
image_1 = image_1(128-30:128+30,164-30:164+30)
image_1_unc = image_1_unc(128-30:128+30,164-30:164+30)
image_2 = image_2(909-15:909+15,1213-15:1213+15)
image_2_unc = image_2_unc(909-15:909+15,1213-15:1213+15)

image_1_conv = image_1_conv(128-30:128+30,164-30:164+30)
image_2_conv = image_2_conv(909-15:909+15,1213-15:1213+15)

;print, mean(image_1(where(finite(image_1) eq 1))),median(image_1(where(finite(image_1) eq 1)))
;print, stdev(image_1(where(finite(image_1) eq 1)))
;print, mean(image_2_unc(where(finite(image_2_unc) eq 1))),median(image_2_unc(where(finite(image_2_unc) eq 1)))
;print, stdev(image_2_unc(where(finite(image_2_unc) eq 1)))


loadct,1,ncolors=256,/silent

set_plot,'ps'
!P.MULTI=[0,2,2]
!P.FONT=0
!X.STYLE=1
!Y.STYLE=1
!P.THICK=4
!X.THICK=4
!Y.THICK=4

device,filename=dir+'helix_images.eps',/helvetica,font_size=30,xsize=60,ysize=60,$
       /color,/encapsulated,bits_per_pixel=8


xnames = ['+45','+30','+15','0','-15','-30','-45']
ynames = ['-45','-30','-15','0','+15','+30','+45']

imdisp,image_1,range=[min(image_1),max(image_1)],/axis,$
       xticks=6,yticks=6,xtickname=xnames,ytickname=ynames,$
       xtitle='x-axis offset [arcsec]',ytitle='y-axis offset [arcsec]',$
       /isotropic,positon=[0.20,0.55,0.55,0.95]
xyouts,45,50,'K_SU org',color=256

imdisp,image_2,range=[min(image_2),max(image_2)],/axis,$
       xticks=6,yticks=6,xtickname=xnames,ytickname=ynames,$
       xtitle='x-axis offset [arcsec]',ytitle='y-axis offset [arcsec]',$
       /isotropic,positon=[0.55,0.55,0.95,0.95]
xyouts,22.5,25,'KPGT org',color=256

imdisp,image_1_conv,range=[min(image_1_conv),max(image_1_conv)],/axis,$
       xticks=6,yticks=6,xtickname=xnames,ytickname=ynames,$
       xtitle='x-axis offset [arcsec]',ytitle='y-axis offset [arcsec]',$
       /isotropic,positon=[0.20,0.20,0.55,0.55]
xyouts,45,50,'K_SU cnv',color=256
FITS_WRITE,dir+'helix_ksu_jsmapb_conv.fits',image_1_conv

imdisp,image_2_conv,range=[min(image_2_conv),max(image_2_conv)],/axis,$
       xticks=6,yticks=6,xtickname=xnames,ytickname=ynames,$
       xtitle='x-axis offset [arcsec]',ytitle='y-axis offset [arcsec]',$
       /isotropic,positon=[0.55,0.20,0.95,0.55]
xyouts,22.5,25,'KPGT cnv',color=256
FITS_WRITE,dir+'helix_kpgt_jsmapb_conv.fits',image_2_conv

device,/close

;fit 2-D Gauassian to the convolved and original K Su images
x = ( findgen(61) - (30.) ) / 1.6 & y = x
xx = x # (y*0 + 1)
yy = (x*0 + 1) # y
estimates = [0.0,max(image_1),5.4/1.6,5.4/1.6,0.0,0.0,(!DPI/180.)*30.0]
p = MPFIT2DPEAK(image_1,Ai,xx,yy,/TILT,ERROR=image_1_unc,ESTIMATES=estimates,PERROR=Aiunc,/GAUSSIAN,CHISQ=csq)
print,'---'
print,'Chi-squared:',csq
print,'Background: ',Ai[0],' +- ',Aiunc[0]
print,'Peak      : ',Ai[1],' +- ',Aiunc[1]
print,'FWHMx     : ',2.35*Ai[2],' +- ',2.35*Aiunc[2]
print,'FWHMy     : ',2.35*Ai[3],' +- ',2.35*Aiunc[3]
print,'Peak x    : ',Ai[4],' +- ',Aiunc[4]
print,'Peak y    : ',Ai[5],' +- ',Aiunc[5]
print,'Rotation  : ',90 - (180.*Ai[6]/!DPI),' +- ',(180.*Aiunc[6]/!DPI)
print,'---'
blob_vals = Ai

;subtract model from the K Su image
xp = xx - Ai[4]
widx = abs(Ai[2])
yp = yy - Ai[5]
widy = abs(Ai[3])
c = cos(Ai[6])
s = sin(Ai[6])

u = ( (xp * (c/widx) - yp * (s/widx))^2 + $
                (xp * (s/widy) + yp * (c/widy))^2 )
model = Ai[0] + Ai[1]*exp(-0.5*u)
model1 = model

residual = image_1 - model

beam = !dpi*(5.4*5.4) / (4.*alog(2.))
print,'Total Emission from Model:',(total(model)*(61.*1.6)^2)/beam

FITS_WRITE,dir+'helix_ksu_image.fits',image_1
FITS_WRITE,dir+'helix_ksu_model.fits',model
FITS_WRITE,dir+'helix_ksu_resid.fits',residual

loadct,5,ncolors=256,/silent

set_plot,'PS'
!P.MULTI=[0,3,1]
!P.FONT=0
!X.STYLE=1
!Y.STYLE=1
!P.THICK=4
!X.THICK=4
!Y.THICK=4

device,filename=dir+'helix_residual_plot.eps',/HELVETICA,/ENCAPSULATED,FONT_SIZE=16,$
       XSIZE=30,YSIZE=20,/COLOR

xtnames = ['+48','+32','+16','0','-16','-32','-48']
ytnames = ['-48','-32','-16','0','+16','+32','+48']

IMDISP,image_1,RANGE=[min(image_1),max(image_1)],/axis,TITLE='Image',$
         	     XTICKS=6,YTICKS=6,XTICKNAME=xnames,YTICKNAME=ynames,$
       		     XTITLE='Offset R.A. [arcsec]',YTITLE='Offset Dec. [arcsec]',$
       		     /ISOTROPIC,COLOR=0
IMDISP,model,RANGE=[min(image_1),max(image_1)],/axis,TITLE='Model',$
         	     XTICKS=6,YTICKS=6,XTICKNAME=xnames,YTICKNAME=ynames,$
       		     XTITLE='Offset R.A. [arcsec]',YTITLE='Offset Dec. [arcsec]',$
       		     /ISOTROPIC,COLOR=0
IMDISP,residual,RANGE=[min(image_1),max(image_1)],/axis,TITLE='Residual',$
         	     XTICKS=6,YTICKS=6,XTICKNAME=xnames,YTICKNAME=ynames,$
       		     XTITLE='Offset R.A. [arcsec]',YTITLE='Offset Dec. [arcsec]',$
       		     /ISOTROPIC,COLOR=0	
device,/close

;fit profile to residual image
x = ( findgen(61) - (30.) ) / 1.6 & y = x
xx = x # (y*0 + 1)
yy = (x*0 + 1) # y
estimates = [0.0,max(residual),5.4/1.6,5.4/1.6,0.0,0.0,(!DPI/180.)*30.0]
p = MPFIT2DPEAK(residual,Ar,xx,yy,/TILT,ERROR=image_1_unc,ESTIMATES=estimates,PERROR=Arunc,/GAUSSIAN,CHISQ=csq)

Ar[2] = 5.4/2.35/1.6
Ar[3] = 5.4/2.35/1.6
print,'---'
print,'Chi-squared:',csq
print,'Background: ',Ar[0],' +- ',Arunc[0]
print,'Peak      : ',Ar[1],' +- ',Arunc[1]
print,'FWHMx     : ',2.35*Ar[2],' +- ',2.35*Arunc[2]
print,'FWHMy     : ',2.35*Ar[3],' +- ',2.35*Arunc[3]
print,'Peak x    : ',Ar[4],' +- ',Arunc[4]
print,'Peak y    : ',Ar[5],' +- ',Arunc[5]
print,'Rotation  : ',0.0;90 - (180.*Ar[6]/!DPI),' +- ',(180.*Arunc[6]/!DPI)
print,'---'

;subtract model from the K Su image
xp = xx - Ar[4]
widx = abs(Ar[2])
yp = yy - Ar[5]
widy = abs(Ar[3])
c = cos(Ar[6])
s = sin(Ar[6])

u = ( (xp * (c/widx) - yp * (s/widx))^2 + $
                (xp * (s/widy) + yp * (c/widy))^2 )
model = Ar[0] + Ar[1]*exp(-0.5*u)
model2 = model


resid = residual - model

beam = !dpi*(5.4*5.4) / (4.*alog(2.))
print,'Total emission from Model2:',max(resid)*1000.*beam

FITS_WRITE,dir+'helix_ksu_model_2.fits',model
FITS_WRITE,dir+'helix_ksu_resid_2.fits',resid

device,filename=dir+'helix_residual_plot2.eps',/HELVETICA,/ENCAPSULATED,FONT_SIZE=16,$
       XSIZE=30,YSIZE=20,/COLOR

xtnames = ['+48','+32','+16','0','-16','-32','-48']
ytnames = ['-48','-32','-16','0','+16','+32','+48']

IMDISP,residual,RANGE=[min(residual),max(residual)],/axis,TITLE='Image',$
         	     XTICKS=6,YTICKS=6,XTICKNAME=xnames,YTICKNAME=ynames,$
       		     XTITLE='Offset R.A. [arcsec]',YTITLE='Offset Dec. [arcsec]',$
       		     /ISOTROPIC,COLOR=0
IMDISP,model,RANGE=[min(residual),max(residual)],/axis,TITLE='Model',$
         	     XTICKS=6,YTICKS=6,XTICKNAME=xnames,YTICKNAME=ynames,$
       		     XTITLE='Offset R.A. [arcsec]',YTITLE='Offset Dec. [arcsec]',$
       		     /ISOTROPIC,COLOR=0
IMDISP,resid,RANGE=[min(residual),max(residual)],/axis,TITLE='Residual',$
         	     XTICKS=6,YTICKS=6,XTICKNAME=xnames,YTICKNAME=ynames,$
       		     XTITLE='Offset R.A. [arcsec]',YTITLE='Offset Dec. [arcsec]',$
       		     /ISOTROPIC,COLOR=0	
device,/close

;calcualate radial profiles
rmaj = image_1(30,0:60)
rmin = image_1(0:60,30)
rval = x

mmaj1 = model1(30,0:60)
mmin1 = model1(0:60,30)

mmaj2 = model2(30,0:60)
mmin2 = model2(0:60,30)

set_plot,'PS'
!P.MULTI=[0,2,2]
!P.FONT=0
!X.STYLE=1
!Y.STYLE=1
!P.THICK=4
!X.THICK=4
!Y.THICK=4

device,filename=dir+'helix_radial_profile.eps',/HELVETICA,/ENCAPSULATED,FONT_SIZE=24,$
       XSIZE=40,YSIZE=40,/COLOR

plot,rval,rmaj,psym=4,xtitle='Radial distance [arcsec]',ytitle='Flux [Jy/pixel]'
oplot,rval,mmaj1,linestyle=2
oplot,rval,mmaj2,linestyle=2
oplot,rval,mmaj1+mmaj2,linestyle=0

plot,rval,rmin,psym=5,xtitle='Radial distance [arcsec]',ytitle='Flux [Jy/pixel]'
oplot,rval,mmin1,linestyle=2
oplot,rval,mmin2,linestyle=2
oplot,rval,mmin1+mmin2,linestyle=0

scatter1 = rmaj - (mmaj1+mmaj2) 
print, stdev(scatter1)*beam*1000.
plot,rval,scatter1,yrange=[min(scatter1),max(scatter1)],psym=4,$
     xtitle='Radial distance [arcsec]',ytitle='Flux [Jy/pixel]'
oplot,[min(rval),max(rval)],[0.0,0.0],linestyle=1

scatter2 = rmin - (mmin1+mmin2) 
print, stdev(scatter2)*beam*1000.
plot,rval,scatter2,yrange=[min(scatter2),max(scatter2)],psym=5,$
     xtitle='Radial distance [arcsec]',ytitle='Flux [Jy/pixel]'
oplot,[min(rval),max(rval)],[0.0,0.0],linestyle=1

device,/close

;
; 160 um
;
;Kate Su image - 3.2"/pixel, cx = 64 cy = 84
image_3 = mrdfits(dir+'helix_ksu_hppjsmapr.fits','image',img3_head)
image_3_unc = mrdfits(dir+'helix_ksu_hppjsmapr.fits','stDev',unc3_head)
image_3 = image_3(64-15:64+15,83-15:83+15)
image_3_unc = image_3_unc(64-15:64+15,83-15:83+15)

;fit 2-D Gauassian to the convolved and original K Su images
x = ( findgen(31) - (15.) ) / 3.2 & y = x
xx = x # (y*0 + 1)
yy = (x*0 + 1) # y
estimates = [0.0,max(image_3),12.3/3.2,12.3/3.2,0.0,0.0,(!DPI/180.)*30.0]
p = MPFIT2DPEAK(image_2,Ai,xx,yy,/TILT,ERROR=image_2_unc,ESTIMATES=estimates,PERROR=Aiunc,/GAUSSIAN,CHISQ=csq)
print,'---'
print,'Chi-squared:',csq
print,'Background: ',Ai[0],' +- ',Aiunc[0]
print,'Peak      : ',Ai[1],' +- ',Aiunc[1]
print,'FWHMx     : ',2.35*Ai[2],' +- ',2.35*Aiunc[2]
print,'FWHMy     : ',2.35*Ai[3],' +- ',2.35*Aiunc[3]
print,'Peak x    : ',Ai[4],' +- ',Aiunc[4]
print,'Peak y    : ',Ai[5],' +- ',Aiunc[5]
print,'Rotation  : ',90 - (180.*Ai[6]/!DPI),' +- ',(180.*Aiunc[6]/!DPI)
print,'---'

;subtract model from the K Su image
xp = xx - Ai[4]
widx = abs(Ai[2])
yp = yy - Ai[5]
widy = abs(Ai[3])
c = cos(Ai[6])
s = sin(Ai[6])

u = ( (xp * (c/widx) - yp * (s/widx))^2 + $
                (xp * (s/widy) + yp * (c/widy))^2 )
model = Ai[0] + Ai[1]*exp(-0.5*u)
model3 = model

residual = image_2 - model

beam = !dpi*(12.3*11.7) / (4.*alog(2.))
print,'Total Emission from Model:',(total(model3)*(61.*1.6)^2)/beam

FITS_WRITE,dir+'helix_ksu_image_160.fits',image_3
FITS_WRITE,dir+'helix_ksu_model_160.fits',model3
FITS_WRITE,dir+'helix_ksu_resid_160.fits',residual

;estimate upper limits by inserting extended source into image
base_img = image_3
;rescale extended blob gaussian for 160um
blob_vals[0] = 0.0 ;take background from 160um model fits rather than 70um
blob_vals[1] = blob_vals[1]*4. ;scale peak value due to new pixel area
blob_vals[2] = blob_vals[2]/2. ;extent of blob scaled to new pixel size
blob_vals[3] = blob_vals[3]/2.
blob_vals[4] = blob_vals[4]/2. ;offset from centre scaled to new pixel size
blob_vals[5] = blob_vals[5]/2.

print, blob_vals

x = ( findgen(31) - (15.) )  & y = x
xx = x # (y*0 + 1)
yy = (x*0 + 1) # y

xp = xx - blob_vals[4]
widx = abs(blob_vals[2])
yp = yy - blob_vals[5]
widy = abs(blob_vals[3])
c = cos(blob_vals[6])
s = sin(blob_vals[6])

u = ( (xp * (c/widx) - yp * (s/widx))^2 + $
                (xp * (s/widy) + yp * (c/widy))^2 )
blob_model = blob_vals[1]*exp(-0.5*u)

;fits_write,dir+'helix_160_blob_base.fits',base_img
;fits_write,dir+'helix_160_blob_mdl.fits',blob_model
blob_est1 = base_img - 0.25*blob_model
fits_write,dir+'helix_160_blob_est_m0p25.fits',blob_est1
blob_est2 = base_img - 0.50*blob_model
fits_write,dir+'helix_160_blob_est_m0p50.fits',blob_est2
blob_est3 = base_img - 0.75*blob_model
fits_write,dir+'helix_160_blob_est_m0p75.fits',blob_est3

;
; 250 um
;
;KPGT image - 6.0"/pixel, cx = 705 cy = 590
image_4 = mrdfits(dir+'helix_kpgt_250_psrc.fits','image',img3_head)
image_4_unc = mrdfits(dir+'helix_kpgt_250_psrc.fits','stDev',unc3_head)

print, median(image_4(where(finite(image_4) eq 1))),stdev(image_4(where(finite(image_4) eq 1)))

end