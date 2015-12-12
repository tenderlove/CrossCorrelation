import pdb, os
import numpy as np
from astropy import wcs
import astropy.io.fits as pyfits

# Write the SDSS check catalog
#p=pyfits.open('/home/ketron/cross/sdss_DR12.fits')
#d = d[d['RA'] < 20]
#d = d[d['dec'] > 1.3]
#ra = d['RA']
#dec = d['dec']
#tb = (pyfits.BinTableHDU.from_columns([
#    pyfits.Column(name='ra',format='E',array=ra),
#    pyfits.Column(name='dec',format='E',array=dec)]))
#tb.writeto('/home/ketron/cross/check/sdss_DR12_partial_check.fits',clobber=True)

# Read sdss check catalog.
cat = pyfits.open('/home/ketron/cross/check/sdss_DR12_partial_check.fits')
ra = cat[1].data['RA']
dec = cat[1].data['dec']
radec = np.dstack((ra,dec))[0,0:]


# Build the header programatically.
w = wcs.WCS(naxis=2)
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

crval0 = np.median(ra)
crval1 = np.median(dec)
w.wcs.crval = [crval0, crval1]

# 30"/pix
cd = 0.008333333333333333 #degree
w.wcs.cdelt = np.array([-cd, cd])
rangex = int(round((max(ra)-min(ra)) / cd)) + 100.
rangey = int(round((max(dec)-min(dec)) / cd)) + 50.
#if rangex % 2 != 0: rangex += 1
#if rangey % 2 != 0: rangey += 1
crpix0 = rangex / 2
crpix1 = rangey / 2
w.wcs.crpix = [crpix0, crpix1]

#w.wcs.naxis1 = rangex
#w.wcs.naxis2 = rangey


x,y = w.wcs_world2pix(ra,dec,1)
#x,y = hd.wcs_world2pix(radec,1)

x = x.round().astype(int)
y = y.round().astype(int)
#pdb.set_trace()
#im = np.random.normal(0,0.01,(rangex,rangey))
im = np.zeros((rangex,rangey))
im[(x,y)] = 1.0

header = w.to_header()
header.insert(0,('SIMPLE','T'))
header.insert(1,('BITPIX','-32'))
#header.insert(2,('NAXIS','2'))
#header.insert(3,('NAXIS1',rangex))
#header.insert(4,('NAXIS1',rangey))

#pdb.set_trace()
p=pyfits.PrimaryHDU(im, header=header) # DO I NEED TO TRANSPOSE???
p.writeto('/home/ketron/cross/check/sdss_check_map.fits',clobber=True)
print('Wrote $cross/check/sdss_check_map.fits')
