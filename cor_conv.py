from __future__ import division
import pyfits, pdb, pickle
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.pyplot import ion


#def get_modes(dim1, dim2, pscale):
    


def c_theta(array1, array2, pscale, **kwargs):
    
    """ Correlation is same as convolution but 
    with one input reversed."""

    
    # Arrays must be same size
    if array1.shape != array2.shape:
        print "Input arrays must be of the same size."
        raise 


    # Normalize and compute discrete cross-correlation.
    array1 = (array1 - np.mean(array1)) / (np.std(array1) * len(array1))
    array2 = (array2 - np.mean(array2)) / np.std(array2)
    corr = signal.fftconvolve(array1, array2[::-1, ::-1], mode='same')

    pdb.set_trace()
    # Compute theta modes.
    lx = np.fft.fftfreq(array1.shape[0])
    lx = np.roll(lx, array1.shape[0] / 2) # shift to center.

    ly = np.fft.fftfreq(array1.shape[1])
    ly = np.roll(ly, array1.shape[1] / 2)

    gx, gy = np.meshgrid(lx, ly)
    gx *= 2. * m.pi / (pscale / 3600. * m.pi / 180.) / 2.
    gy *= 2. * m.pi / (pscale / 3600. * m.pi / 180.) / 2.

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



pscale = 12. # "/pix
im1 = pyfits.open('helms_plw.fits')
im2 = pyfits.open('nvss.fits')

im1 = im1[0].data
im2 = im2[0].data

ddir = '/data-2/cross/helms/correlations/'
savename = ddir + 'plw_x_nvss.pickle'

ctheta = c_theta(im1, im2, pscale, savename=savename)

pdb.set_trace()


pickle.dump(corr, open(savename, 'wb'))
print 'Saved correlation as '+savename

#y, x = np.unravel_index(np.argmax(corr), corr.shape) # find the match
ion()
plt.imshow(corr)
cbar = plt.colorbar()
plt.show()
pdb.set_trace()


# does numpy.fftconvolve shift the arrays?
#g = np.random.normal([100,100])
#check = signal.fftconvolve(g, g, mode='same')
