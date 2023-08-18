#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 10:06:00 2021
@author: cvaneck

Breaks large cube into 4 parts.

Call sequence:
    
script.py incube source_list outfile_basename

"""

import astropy.io.fits as pf
import sys
import os
import astropy.wcs as wcs
import numpy as np
import astropy.table as at

padding=35 #How many pixels to increase the size by, as appropriate


basefile=sys.argv[-3]
basefile=basefile.replace('.i.','._.').replace('.q.','._.').replace('.u.','._.')

outfile_base=sys.argv[-1]

#Open I cube, get size:
hdu=pf.open(basefile.replace('._.','.i.'),memmap=True)

header=hdu[0].header
xsize=header['NAXIS1']
ysize=header['NAXIS2']


print('Writing I cubes.')
#Bottom left:
subimage=hdu[0].data[:,:,0:ysize//2+padding,0:xsize//2+padding]
pf.writeto(outfile_base+'.i.bottomleft.fits',subimage,header)

#Bottom right:
subimage=hdu[0].data[:,:,0:ysize//2,xsize//2-padding:]
outheader=header.copy()
outheader['CRPIX1']=header['CRPIX1']-(xsize//2-padding)
pf.writeto(outfile_base+'.i.bottomright.fits',subimage,outheader)

#Top left:
subimage=hdu[0].data[:,:,ysize//2-padding:,0:xsize//2+padding]
outheader=header.copy()
outheader['CRPIX2']=header['CRPIX2']-(ysize//2-padding)
pf.writeto(outfile_base+'.i.topleft.fits',subimage,outheader)

#Top right:
subimage=hdu[0].data[:,:,ysize//2-padding:,xsize//2-padding:]
outheader=header.copy()
outheader['CRPIX1']=header['CRPIX1']-(xsize//2-padding)
outheader['CRPIX2']=header['CRPIX2']-(ysize//2-padding)
pf.writeto(outfile_base+'.i.topright.fits',subimage,outheader)

hdu.close()

print('Writing Q cubes.')
hdu=pf.open(basefile.replace('._.','.q.'),memmap=True)


#Bottom left:
subimage=hdu[0].data[:,:,0:ysize//2+padding,0:xsize//2+padding]
pf.writeto(outfile_base+'.q.bottomleft.fits',subimage,header)

#Bottom right:
subimage=hdu[0].data[:,:,0:ysize//2,xsize//2-padding:]
outheader=header.copy()
outheader['CRPIX1']=header['CRPIX1']-(xsize//2-padding)
pf.writeto(outfile_base+'.q.bottomright.fits',subimage,outheader)

#Top left:
subimage=hdu[0].data[:,:,ysize//2-padding:,0:xsize//2+padding]
outheader=header.copy()
outheader['CRPIX2']=header['CRPIX2']-(ysize//2-padding)
pf.writeto(outfile_base+'.q.topleft.fits',subimage,outheader)

#Top right:
subimage=hdu[0].data[:,:,ysize//2-padding:,xsize//2-padding:]
outheader=header.copy()
outheader['CRPIX1']=header['CRPIX1']-(xsize//2-padding)
outheader['CRPIX2']=header['CRPIX2']-(ysize//2-padding)
pf.writeto(outfile_base+'.q.topright.fits',subimage,outheader)

hdu.close()

print('Writing U cubes.')
hdu=pf.open(basefile.replace('._.','.u.'),memmap=True)


#Bottom left:
subimage=hdu[0].data[:,:,0:ysize//2+padding,0:xsize//2+padding]
pf.writeto(outfile_base+'.u.bottomleft.fits',subimage,header)

#Bottom right:
subimage=hdu[0].data[:,:,0:ysize//2,xsize//2-padding:]
outheader=header.copy()
outheader['CRPIX1']=header['CRPIX1']-(xsize//2-padding)
pf.writeto(outfile_base+'.u.bottomright.fits',subimage,outheader)

#Top left:
subimage=hdu[0].data[:,:,ysize//2-padding:,0:xsize//2+padding]
outheader=header.copy()
outheader['CRPIX2']=header['CRPIX2']-(ysize//2-padding)
pf.writeto(outfile_base+'.u.topleft.fits',subimage,outheader)

#Top right:
subimage=hdu[0].data[:,:,ysize//2-padding:,xsize//2-padding:]
outheader=header.copy()
outheader['CRPIX1']=header['CRPIX1']-(xsize//2-padding)
outheader['CRPIX2']=header['CRPIX2']-(ysize//2-padding)
pf.writeto(outfile_base+'.u.topright.fits',subimage,outheader)

hdu.close()




print('Breaking up source list')

cat=at.Table.read(sys.argv[-2])
extension=sys.argv[-2][-4:]
if extension == '.csv':
    fileformat='csv'
elif extension == '.xml':
    fileformat='votable'
else:
    raise Exception('This script currently only supports CSV or VOTable catalogs.')

cs = wcs.WCS(header,fix=False,naxis=(1,2))
dx, dy = cs.all_world2pix(cat['ra_deg_cont'],cat['dec_deg_cont'],0)
ix = np.rint(dx)
iy = np.rint(dy)

bottomleft=(dx <= xsize//2) & (dy <= ysize//2)
at.Table(cat[bottomleft]).write(outfile_base+'.bottomleft'+extension,format=fileformat)

bottomlright=(dx > xsize//2) & (dy <= ysize//2)
at.Table(cat[bottomlright]).write(outfile_base+'.bottomright'+extension,format=fileformat)
    
topleft=(dx <= xsize//2) & (dy > ysize//2)
at.Table(cat[topleft]).write(outfile_base+'.topleft'+extension,format=fileformat)

topright=(dx > xsize//2) & (dy > ysize//2)
at.Table(cat[topright]).write(outfile_base+'.topright'+extension,format=fileformat)








