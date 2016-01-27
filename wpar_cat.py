import pdb, os, time, sharedmem, glob, socket
import pickle, sys, tempfile, shutil, argparse
import numpy as np
import math as m
import astropy.io.fits as pyfits
from astropy import wcs
import randMapDiscrete, randPosDiscrete, filterCatalog, checkXY, loadCat
from astroML.correlation import bootstrap_two_point_angular
import cross_correlation
import treecorr
import matplotlib.pyplot as plt


if __name__ == '__main__':
    parser = (argparse.ArgumentParser(
        description='Compute w(theta) for 2 catalogs. This just sets up '+
        'the catalogs and calls the astroML bootstrap_two_point_angular '+
        'routine. For an auto-correlation, set cat1=cat2.'))
    (parser.add_argument('-cat1', '--cat1', type=str,
                         help='Name of first catalog. (See loadCat.py input.)'))
    (parser.add_argument('-cat2','--cat2', type=str, 
                         help='Name of second catalog. (See loadCat.py input.)'))
    (parser.add_argument('-nbins', '--nbins', type=int, default=24, 
                         help='Number of bins...'))
    (parser.add_argument('-filt','--filt',type=str,default='plw',
                         help='I load this to get the area, and to '+
                         'convert RA, dec -> x, y for stuff.'))
    
    options = parser.parse_args()
    cat1    = options.cat1
    cat2    = options.cat2
    nbins   = options.nbins
    filt    = options.filt
    
    #def compute(cat1, cat2, nbins=24, filt='plw'):
    
    # '''
    # Computes the angular cross-correlation of cat1 against cat2.
    # For an auto-correlation, set cat1=cat2. This just sets up the
    # catalogs (given the Herschel map area & [dust] masks), and 
    # calls the astroML bootstrap_two_point_angular routine.
    # '''
    
    cr = os.getenv('cr')
    crmaps = cr + 'herschel/'
    crcats = cr + 'cats/'
    crcheck = '/home/ketron/cross/check/'


    ###
    #Get Herschel FOV
    fname = crmaps + 'helms_' + filt.lower() + '.fits'
    print 'Getting Herschel FOV from:' + fname
    hdu = pyfits.open(fname)
    exten = 0 #if filt == 'check' else 1
    dataIm = hdu[exten].data
    Nx, Ny = dataIm.shape[0], dataIm.shape[1]
    # Astropy wcs doesn't like comments for some reason.
    while 'COMMENT' in hdu[exten].header: hdu[exten].header.remove('COMMENT')
    hd = wcs.WCS(hdu[exten].header)

    mask = dataIm * 0.
    mask += 1.
    mask[np.where(dataIm == 0)] = 0.
    # Remove NaNs
    #dataIm[np.where(np.isnan(dataIm))] = 0.
    # INCORPORATE DUST MASK HERE
    print('You should incorporate the dust mask?')


    ###
    #Load catalogs
    ra1, dec1 = loadCat.loadCat(cat1, crcats)
    if cat1 == cat2:
        ra2, dec2 = ra1, dec1
    else:
        ra2, dec2 = loadCat.loadCat(cat2, crcats)

    # Generate random catalogs. randMapDiscrete considers FOV.
    Ngal = len(ra1)#*5
    randIm1 = randMapDiscrete.randMapDiscrete(Nx, Ny, Ngal, nanmask=mask)
    ind1 = np.where(randIm1 == 1)
    randRA1,randDec1 = hd.wcs_pix2world(ind1[1],ind1[0],1)

    Ngal = len(ra2)#*5
    randIm2 = randMapDiscrete.randMapDiscrete(Nx, Ny, Ngal, nanmask=mask)
    ind2 = np.where(randIm2 == 1)
    randRA2, randDec2 = hd.wcs_pix2world(ind2[1],ind2[0],1)

    #hdu1 = pyfits.PrimaryHDU(randIm1)
    #hdu1.writeto('rand1.fits')



    ###
    #Convert to ra, dec to (x,y) on Herschel map. 
    x1, y1 = hd.wcs_world2pix(ra1, dec1, 1)
    x1, y1 = x1.round().astype(int), y1.round().astype(int)
    
    x2, y2 = hd.wcs_world2pix(ra2, dec2, 1)
    x2, y2 = x2.round().astype(int), y2.round().astype(int)
    
    # Why does wcs swap x and y? I have some strange header maybe.
    if max(x1) > Nx or max(y1) > Ny:
        print('x1 and y1 are swapped from wcs_world2pix for some reason. '+
                    "I'm swapping them.")
        x1, y1 = y1, x1
    if max(x2) > Nx or max(y2) > Ny:
        print('x2 and y2 are swapped from wcs_world2pix for some reason. '+
                    "I'm swapping them.")
        x2, y2 = y2, x2 

    # If I'm running on my rotated cropped image the catalog will be
    # slightly larger than it should, so those sources need to be
    # removed via checkXY().
    x1, y1   = checkXY.checkXY(x1,   y1, Nx, Ny)
    x2, y2   = checkXY.checkXY(x2,   y2, Nx, Ny)

    # Remove sources in catalog that should be behind the mask.
    x1, y1 = filterCatalog.filterCatalog(x1, y1, Nx, Ny, nanmask=mask)
    x2, y2 = filterCatalog.filterCatalog(x2, y2, Nx, Ny, nanmask=mask)


    ###
    # Go back to RA, dec and feed to bootstrap.
    #ra1, dec1 = hd.wcs_pix2world(x1, y1, 1)
    #ra2, dec2 = hd.wcs_pix2world(x2, y2, 1)
    ra1 = ra1[0:10000]
    dec1 = dec1[0:10000]

    print('Computing W(theta) for '+cat1+' x '+cat2+'...')
    c1 = treecorr.Catalog(ra=ra1, dec=dec1, ra_units='deg', dec_units='deg')
    c2 = treecorr.Catalog(ra=ra2, dec=dec2, ra_units='deg', dec_units='deg')    

    r1 = treecorr.Catalog(ra=randRA1, dec=randDec1, ra_units='deg', dec_units='deg')
    r2 = treecorr.Catalog(ra=randRA2, dec=randDec2, ra_units='deg', dec_units='deg')

    bsz, min_sep, max_sep, sep_units, bin_slop = 0.3, 0.1, 20, 'deg', 0.5
    nn = treecorr.NNCorrelation(bin_size=bsz, min_sep=min_sep, \
                                max_sep=max_sep,sep_units=sep_units)#, bin_slop=bin_slop)
    
    rr = treecorr.NNCorrelation(bin_size=bsz, min_sep=min_sep, \
                                max_sep=max_sep,sep_units=sep_units)#, bin_slop=bin_slop)
    
    dr = treecorr.NNCorrelation(bin_size=bsz, min_sep=min_sep, \
                                max_sep=max_sep,sep_units=sep_units)#, bin_slop=bin_slop)
    
    rd = treecorr.NNCorrelation(bin_size=bsz, min_sep=min_sep, \
                                max_sep=max_sep,sep_units=sep_units)#, bin_slop=bin_slop)


    
    nn.process(c1, c2)
    rr.process(r1, r2)
    dr.process(c1, r2)
    rd.process(c2, r1)

    xi,varxi = nn.calculateXi(rr,dr)#,rd)

    plt.errorbar(np.exp(nn.logr), xi, varxi, fmt='o', ecolor='gray',lw=10)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$w(\theta)$',fontsize=15)
    plt.xlabel(r'$\theta\ (deg)$',fontsize=15)
    plt.ion()
    plt.show()

    pdb.set_trace()

    


    bins = 10 ** np.linspace(np.log10(1. / 60.), np.log10(60), nbins)
    wtheta = [bins]
    Nbootstraps=10
    #if cat1 == cat2:
    #print('Auto-correlation...')
    wtheta += bootstrap_two_point_angular(ra1, dec1,
                                          bins=bins,
                                          method='landy-szalay',
                                          Nbootstraps=Nbootstraps)
    
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    bins, b_corr, b_corr_err, b_bootstraps = wtheta[0], wtheta[1], wtheta[2],wtheta[3]
    pdb.set_trace()
    plt.figure()
    plt.xscale('log'); plt.yscale('log')
    plt.errorbar(bin_centers, b_corr, b_corr_err,fmt='o', ecolor='gray', lw=10)
    plt.show()
