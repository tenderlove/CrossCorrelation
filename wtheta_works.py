from __future__ import division
import pyfits, pdb, pickle, os, time
import numpy as np
import math as m
from astropy import wcs
import matplotlib.pyplot as plt
from matplotlib.pyplot import ion



def get_info(filt):

    ''' 
    Pixel scale, filename, dimensions, etc 
    from filter.
    '''

    cr = os.getenv('cr')
    filt = str(filt)
    filt = filt.lower()
    
    if filt == 'plw':
        pscale = 12.
        fname = cr + 'herschel/HELMS_image_500_SANEPIC_v0.2.fits'
    elif filt == 'pmw':
        pscale = 8.3333
        fname = cr + 'herschel/HELMS_image_350_SANEPIC_v0.2.fits'
    elif filt == 'psw':
        pscale = 6.
        fname = cr + 'herschel/HELMS_image_250_SANEPIC_v0.2.fits'
    else:
        print "You didn't choose an appropriate Herschel filter."
        raise Exception(filt)

    hdu = pyfits.open(fname)
    im = hdu[1].data
    # replace nans with zeros?
    inds = np.where(np.isnan(im))
    im[inds] = 0.

    # Need to make array square.
    d1 = hdu[1].data.shape[0]
    d2 = hdu[1].data.shape[1]
    if d1 > d2:
        increase = d1 - d2
        im = np.pad(im, ((0,0),(0,increase)),'edge')
    if d2 > d1:
        increase = d2 - d1
        im = np.pad(im, ((0,increase),(0,0)),'edge')

    # Astropy wcs doesn't accept comments for some dumb reason.
    while 'COMMENT' in hdu[1].header: 
        hdu[1].header.remove('COMMENT')

    hd = wcs.WCS(hdu[1].header)
    d1 = hdu[1].data.shape[0]
    d2 = hdu[1].data.shape[1]
    
    info = ({'pscale':pscale, 'd1':d1, 'd2':d2, 'fname':fname, 
             'im':im, 'hd':hd, 'nbins': 20})
    return info



def wsum(ra, dec, filt='plw', **kwargs):
    
    ''' 
    Given a catalog of RA and dec, go to those positions and add
    up all the values in the array within some annulus, defined by the r
    binning in linbins.
    '''

    class WTheta:
        """
        Empty object container for each ra, dec pair. These will go into
        wthetas list below.
        """
        def __init__(self): 
            self.mean = None
            self.std = None
            self.median = None
            self.numel = None
            self.numel_nonz = None
            self.sum = None
            self.max = None
            self.min = None
            self.r = None

    wthetas = [] #there will be len(radec) WTheta classes in here.
            
    info = get_info(filt)
    im = info['im']    
    #im = np.ones((10,20)) ### TEST
    d1 = im.shape[0]
    d2 = im.shape[1]
    annulus_width = 5000.


    # ra dec -> x y
    radec = np.dstack((ra,dec))[0,0:]
    xy = info['hd'].wcs_world2pix(radec,1)
    x,y = xy[:,0], xy[:,1]
    
    # Loop through each x, y pair
    for i in range(len(x)):                

        # Meshgrid
        xvec = np.arange(d1) - x[i] #- inc1 #these shifts are paramount to this exercise. 
        yvec = np.arange(d2) - y[i] #- inc2

        xgrid, ygrid = np.meshgrid(xvec, yvec)

        # Radial coords
        r = abs(xgrid + 1j*ygrid)
        working_mask = np.ones(im.shape, bool)
        rmax = r[working_mask].max()
        #rmax = min(im.shape)
        dr = np.abs([xgrid[0,0] - xgrid[0,1]]) * annulus_width
        radial = np.arange(rmax/dr)*dr + dr/2.
        nrad = len(radial)

        wtheta = WTheta()
        wtheta.mean = np.zeros(nrad)
        wtheta.sum = np.zeros(nrad)
        #wtheta.std = np.zeros(nrad)
        #wtheta.median = np.zeros(nrad)
        wtheta.numel = np.zeros(nrad, dtype=int)
        wtheta.numel_nonz = np.zeros(nrad)
        #wtheta.max = np.zeros(nrad)
        #wtheta.min = np.zeros(nrad)
        wtheta.r = radial
        
        # Time
        if i == 0: time0 = time.time()
        for irad in range(nrad): #= 1:numel(radial)
            minrad = irad*dr
            maxrad = minrad + dr
            thisindex = (r>=minrad) * (r<maxrad) * working_mask
            
            if not thisindex.ravel().any():
                wtheta.mean[irad] = np.nan
                wtheta.sum[irad] = np.nan
                #wtheta.std[irad]  = np.nan
                #wtheta.median[irad] = np.nan
                wtheta.numel[irad] = 0
                wtheta.numel_nonz[irad] = 0
                #wtheta.max[irad] = np.nan
                #wtheta.min[irad] = np.nan
                #print 'notelse'
                #pdb.set_trace()
            else:
                wtheta.mean[irad] = im[thisindex].mean()
                wtheta.sum[irad] = im[thisindex].sum()#[r<maxrad].sum()
                #wtheta.std[irad]  = im[thisindex].std()
                #wtheta.median[irad] = np.median(im[thisindex])
                wtheta.numel[irad] = im[thisindex].size
                nonzindex = np.where(im[thisindex] != 0)
                wtheta.numel_nonz[irad] = nonzindex[0].size
                #wtheta.max[irad] = im[thisindex].max()
                #wtheta.min[irad] = im[thisindex].min()
                #plt.imshow(thisindex*1.)
                #ion()
                #plt.show()#;plt.colorbar()
                #print 'else'
                
            #pdb.set_trace()

        wthetas.append(wtheta)
        print 'finished '+ str(i) +' of '+str(len(x))
        if i == 0:
            oneloop = time.time() - time0
            print("%1.5f seconds for one radial loop." % (oneloop))
            ttot = len(x) * oneloop / 3600.
            print("This is going to take ~ %4.2f hours!" % (ttot))
    return wthetas

hdu = pyfits.open(os.getenv('cr') + 'cats/nvss.fits')
ra = hdu[1].data['ALPHA_J2000']
dec = hdu[1].data['DELTA_J2000']
radec = np.dstack((ra,dec))[0,0:]
s = wsum(ra, dec)

savename='nvss_wthetas_plw.pickle'
pickle.dump(s, open(savename, 'wb'))

print 'Saved wthetas as '+savename

pdb.set_trace()
 
