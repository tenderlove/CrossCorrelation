from __future__ import division
import pyfits, pdb, pickle, os
import numpy as np
import math as m
from astropy import wcs
import matplotlib.pyplot as plt
from matplotlib.pyplot import ion
import radial_data


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

    # Astropy wcs doesn't accept comments for some dumb reason.
    while 'COMMENT' in hdu[1].header: 
        hdu[1].header.remove('COMMENT')

    hd = wcs.WCS(hdu[1].header)
    dim1 = hdu[1].data.shape[0]
    dim2 = hdu[1].data.shape[1]
    
    info = ({'pscale':pscale, 'dim1':dim1, 'dim2':dim2, 'fname':fname, 
             'im':im, 'hd':hd, 'nbins': 20})
    return info


def logbins(**kwargs):

    '''
    Logarithmic binning
    '''

    dim1, dim2, nbins = kwargs['dim1'], kwargs['dim2'], kwargs['nbins']
    min, max = 10, m.floor(m.sqrt(dim1**2 + dim2**2))
    step = pow(max / min, 1. / (nbins - 1))
    lbins = min * pow(step, np.arange(nbins))
    np.insert(lbins, 0, 1.) #first element
    return lbins


def linbins(**kwargs):

    '''
    linear binning
    '''

    dim1, dim2, nbins = kwargs['dim1'], kwargs['dim2'], kwargs['nbins']
    len = m.floor(m.sqrt(dim1**2 + dim2**2))
    lbins = np.arange(1, len, len/nbins)
    return lbins



# def wsum(ra, dec, filt='plw', **kwargs):
    
#     ''' 
#     Given a catalog of RA and dec, go to those positions and add
#     up all the values in the array within some annulus, defined by the r
#     binning in linbins.
#     '''

    

#     info = get_info(filt)
#     bins = linbins(**info)
#     im = info['im']
#     d1 = info['dim1']
#     d2 = info['dim2']

#     # ra dec -> x y
#     radec = np.dstack((ra,dec))[0,0:]
#     xy = info['hd'].wcs_world2pix(radec,1)
#     #x, y = xy[:,0], xy[:,1]
    
#     # The array needs to be padded with zeros, so if I'm in 
#     # a corner, the radial measurement doesn't go off the array.
#     if d1 % 2 != 0: 
#         d1 +=1
#     if d2 % 2 != 0:
#         d2 += 1
    
        


#     #    im = np.pad(im, ((d1/2, d1/2), (d2/2, d2/2)), 'edge')
#     #    hdu = pyfits.PrimaryHDU(im)
#     #    hdu.writeto('t.fits',clobber=True)
#     #    pdb.set_trace()

#     # Meshgrid
#     ### THIS SHOULD BE SQUARE, SIZE = MAX(IM.SHAPE)
#     working_mask = np.ones(im.shape,bool)
#     xvec = np.arange(d1) - d1/2. 
#     yvec = np.arange(d2) - d2/2.
#     xgrid, ygrid = np.meshgrid(xvec, yvec)

#     # Radial
#     annulus_width = 500.
#     r = abs(xgrid + 1j*ygrid)
#     #####################################################3
#     ### MAKE THE ARRAY SQUARE AND EVEN!!!
#     ####################################################
#     pdb.set_trace()
#     rmax = r[working_mask].max()
#     #rmax = min(im.shape)
#     dr = np.abs([x[0,0] - x[0,1]]) * annulus_width
#     radial = np.arange(rmax/dr)*dr + dr/2.
    
#     pdb.set_trace()

#     # Loop through each x, y pair
#     for i in range(len(xy[0,:])):
#         xpos, ypos = xy[i,:][0], xy[i,:][1]

#         # make a grid and compute sum(annulus) at each bin.
# #        for j in range(len(bins)-1):
# #            if j == 0: 
# #                rmin, rmax = 0., bins[j]
                
# #            else:
# #                rmax, rmin = bins[j], bins[j+1]
            

#             # rmax = min(im.shape)

#             # dr = ny.abs([x[0,0] - x[0,1]]) * annulus_width

#             # # shift to the ra, dec pos we're at.
#             # xgrid += xpos
#             # ygrid += ypos
            
#             # r = abs(xgrid + 1j*ygrid)
#             # dr = np.abs([xvec[0,0] - xvec[0,1]]) * annulus_width


#             # r = radial_data.radial_data(im, x=xx,y=yy,)
#             # pdb.set_trace()   
            

#hdu = pyfits.open(os.getenv('cr') + 'cats/nvss.fits')
#ra = hdu[1].data['ALPHA_J2000']
#dec = hdu[1].data['DELTA_J2000']
#radec = np.dstack((ra,dec))[0,0:]
#s = wsum(ra, dec)
    
########

#def get_r(**kwargs):

#    """ Make a grid from the center of the image, then offset it?
#    """
    
    
def wsum(*args, **kwargs):
    
    '''
    
    '''

    im = kwargs['im']
    im = np.zeros((100,50))
    
    # Compute theta modes via fourier modes. 
    # Dimensions don't matter at the moment.
    lx = np.fft.fftfreq(im.shape[0]*2)
    #lx = np.roll(lx, double(im.shape[0] / 2) # shift to center.

    ly = np.fft.fftfreq(im.shape[1]*2)
    #ly = np.roll(ly, im.shape[1] / 2)

    gx, gy = np.meshgrid(ly, lx)
    #gx *= 2. * m.pi / (pscale / 3600. * m.pi / 180.) / 2.
    #gy *= 2. * m.pi / (pscale / 3600. * m.pi / 180.) / 2.

    # why can't I do this in one line?
    gxsq, gysq = (gx)**2, (gy)**2
    ladd = np.add(gxsq,gysq)
    ell = np.sqrt(ladd)

    pdb.set_trace()
    thetax = 1 / lx * pscale * 3600.
    thetay = 1 / ly * pscale * 3600.


    #zero = np.where(~np.isnan(refim))
    #sx = np.linspace(-array1.shape[0]/2., 
    #sy = array1.shape[1]


info = get_info('plw')
w=wsum(**info)

ddir = '/data-2/cross/helms/correlations/'
savename = ddir + 'plw_x_nvss.pickle'


# ctheta = c_theta(im1, im2, pscale, savename=savename)

# pdb.set_trace()


# pickle.dump(corr, open(savename, 'wb'))
# print 'Saved correlation as '+savename

# #y, x = np.unravel_index(np.argmax(corr), corr.shape) # find the match
# ion()
# plt.imshow(corr)
# cbar = plt.colorbar()
# plt.show()
# pdb.set_trace()


# does numpy.fftconvolve shift the arrays?
#g = np.random.normal([100,100])
#check = signal.fftconvolve(g, g, mode='same')



################################3
## PADDING:
################################3
# Pad the array with zeros.
    #nans = np.where(np.isnan(im))
    #if nans: im[nans] = 0.
    # increase = round(pow(d1**2+d2**2,0.5))
    # increase = 0.
    # if increase % 2 != 0: increase += 1
    # if d1 == d2:
    #     inc1 = increase/2
    #     inc2 = increase/2
    # elif d2 > d1:
    #     inc1 = np.floor((d2 - d1)/2 + increase/2)
    #     inc2 = increase/2
    # else:
    #     inc1 = increase/2
    #     inc2 = np.floor((d2 - d1)/2 + increase/2)
        
    # if d1 % 2 != 0 and d2 % 2 != 0: 
    #     padim = np.pad(im, ((inc1,inc1+1), (inc2,inc2+1)), 'edge')
    # elif d2 % 2 != 0 and d1 % 2 == 0: 
    #     padim = np.pad(im, ((inc1,inc1+1), (inc2,inc2)), 'edge')
    # elif d1 %2 !=0 and d2 % 2 == 0:
    #     padim = np.pad(im, ((inc1,inc1), (inc2,inc2+1)), 'edge')
    # else:
    #     padim = np.pad(im, ((inc1,inc1), (inc2,inc2)), 'edge')
    # print padim.shape
