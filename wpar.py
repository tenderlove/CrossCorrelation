from __future__ import division
import pdb, os, time, sharedmem, glob, socket
import pickle, sys, tempfile, shutil, argparse
import numpy as np
import math as m
import astropy.io.fits as pyfits
from astropy import wcs
from joblib import Parallel, delayed, load, dump
import randMapDiscrete, randPosDiscrete, filterCatalog, checkXY, loadCat


def randPosDiscrete(Nx, Ny, Ngal):
    ''' 
    Return Ngal random integer points in 2D array 
    of size [Nx, Ny].
    '''

    randx = np.random.random_integers(0, Nx-1, Ngal)
    randy = np.random.random_integers(0, Ny-1, Ngal)
    return (randx, randy)
 

def randMapDiscrete(Nx, Ny, Ngal, nanmask=1.):
    '''
    Apply the ROI and source mask (if any) to the random
    map output from randPosDiscrete.
    '''

    im = np.zeros((Nx,Ny))
    randx, randy = randPosDiscrete(Nx, Ny, Ngal)
    im[randx, randy] = 1.0
    im *= nanmask
    return im


#def randMapDiffuse(Nx,Ny,Nal): for background=False. Not there yet.


def filterCatalog(xx, yy, Nx, Ny, nanmask=1.0):

    '''
    I'll be using catalogs with sources outside the Herschel area (mask),
    so inject them into a map, apply a mask, and return
    the sources that remain!
    '''

    xymap = np.zeros((Nx, Ny))
    xymap[xx,yy] = 1.0
    xymap *= nanmask
    remainders = np.where(xymap == 1)
    xcat, ycat = remainders[0], remainders[1]
    return (xcat, ycat)

def checkXY(xxx, yyy, Nx, Ny):

    '''
    x, y arrays should be between 0 and N?-1. This is broken when specifically
    when I use catalogs generated from separate masks for Herschel
    auto-correlation.
    '''

    ptsInside = ((xxx >= 0) & (yyy >= 0) & (xxx < Nx) & (yyy < Ny))
    xxx = xxx[ptsInside]
    yyy = yyy[ptsInside]
    return (xxx, yyy)


def loadCat(cname, crcats):

    '''
    Loads catalogs, which are presumably unique.
    This is not general and applies to only my catalogs.
    '''

    #NVSS radio catalog
    if cname == 'nvss':
        print('Reading NVSS catlog:')
        print(crcats + 'nvss_cat.fits')
        hdu = pyfits.open(crcats+'nvss_cat.fits')
        ra = hdu[1].data['ALPHA_J2000']
        dec = hdu[1].data['DELTA_J2000']

    # SDSS optical catalog -- 4 million+ entries ==> magnitude clip!
    elif cname == 'u' or cname == 'g' or cname == 'r' or cname == 'z':
        print('Reading optical catalog:')
        print(crcats + 'sdss_master.fits')
        hdu = pyfits.open(crcats + 'sdss_master.fits')
        data = hdu[1].data
        if cname == 'u':
            mag = data['psfMag_u']
        elif cname == 'g':
            mag = data['psfMag_g']
        elif cname == 'r':
            mag = data['psfMag_r']
        elif cname == 'i':
            mag = data['psfMag_i']
        elif cname == 'z':
            mag = data['psfMag_z']
        data = data[(mag > 17) & (mag <21)]
        ra = data['ra']
        dec = data['dec']
        print('Clipped sources for only 17 < ' + cname.upper()  + ' < 21')

    # Check against a cropped SDSS region and compare against astroML code.
    # It works.
    elif cname == 'check':
        print('Reading check catlog:')
        print('$cross/check/sdss_DR12_partial_check.fits')
        hdu = pyfits.open(crcheck + 'sdss_DR12_partial_check.fits')
        ra = hdu[1].data['ra']
        dec = hdu[1].data['dec']

    # Herschel auto-correlation if catname = filt.
    #I set ra, dec = hra, hdec below.
    elif cname == 'plw' or cname == 'pmw' or cname == 'psw':
        ra, dec = 0,0
    else:
        raise Exception('Which catalog?')

    return (ra, dec)


def getCorrelationData(filt, catname, wtype, nbins, nomask, background):

    '''
    Loads Herschel arrays and exterior catalog from kwds. Returns them with Pixel scale,
    filename, dimensions, etc. If wtype == 'RR' then random maps and catalogs are generated
    instead. If wtype='DR' a random map is made and x,y positions are from catalog. If
    background is set, the code correlates the input catalog name positions with
    the CFIRB in Herschel maps. Otherwise it generates maps from Herschel catalogs.

    '''

    class correlationData:

        '''
        Empty data structure.
        '''

        def __init__(self):
            # Input arguments.
            self.filt       = None
            self.catname    = None
            self.wtype      = None
            self.nomask     = None
            self.background = None
            self.im         = None
            self.xy         = None
            self.nbins      = None
            self.pscale     = None
            self.crop       = False

    filt = str(filt)
    filt = filt.lower()

    hostname = socket.gethostname()
    if hostname != 'thegord':
        cr = '/data/users/ketron/'
        crcats = cr
        crmaps = cr
        crcheck = cr
    else:
        cr = os.getenv('cr')
        crmaps = cr + 'herschel/'
        crcats = cr + 'cats/'
        crcheck = '/home/ketron/cross/check/'


    # Herschel catalog should be loaded no matter what.
    cname = crmaps + 'helms_v0.2_' + filt.upper() + '_SXT.fits'
    print 'Reading Herschel catalog:' + cname
    hcat = pyfits.open(cname)
    hdata = hcat[1].data
    hdata = hdata[hdata['flux'] > 30]
    hra, hdec = hdata['ra'], hdata['dec']
    print('Doing herschel flux cut > 30 mJy.')


    # Load catalog.
    ra, dec = loadCat(catname, crcats)
    if filt == 'check':
        print('Setting hra, hdec = checkRA, checkDEC...')
        hra, hdec = ra, dec
    if catname == 'plw' or catname == 'pmw' or catname == 'psw':
        print('Setting ra, dec = hra, hdec')
        ra, dec = hra, hdec


    # Read [Herschel] array.
    if filt == 'plw' or filt == 'pmw' or filt == 'psw':
        fname = crmaps + 'helms_' + filt + '.fits'
        pscale = 12.
    elif filt == 'check':
        pscale = 30.
        fname = crcheck + 'sdss_check_map.fits'
    else:
        print "You didn't choose an appropriate Herschel filter."
        raise Exception(filt)

    print('Reading '+fname)
    hdu = pyfits.open(fname)
    exten = 0 #if filt == 'check' else 1


    dataIm = hdu[exten].data
    Nx, Ny = dataIm.shape[0], dataIm.shape[1]
    # Astropy wcs doesn't like comments for some reason.
    while 'COMMENT' in hdu[exten].header: hdu[exten].header.remove('COMMENT')
    hd = wcs.WCS(hdu[exten].header)


    # Mask the simulation maps the same as the real image via mask.
    nanmask = dataIm * 0.
    nanmask += 1.
    nanmask[np.where(dataIm == 0)] = 0.
    # Remove NaNs
    dataIm[np.where(np.isnan(dataIm))] = 0.


    # Incorporate the source mask into the nanmask?
    if not nomask:
        hdumask = pyfits.open(crmaps + 'mask_'+filt+'.fits')
        mask = hdumask[0].data
        nanmask *= mask
        print('Applying source mask to the nanmap map.')


    # I want 3x more random sources and the nanmap will
    # decrease this, so find the percentage I need to increase 3 by.
    ones = np.where(nanmask == 1)
    zeros = np.where(nanmask == 0)
    nz, no = len(zeros[0]), len(ones[0])
    increase = float(nz)/float(no)
    Rfact = (3 * (.91 + increase))

    data = correlationData()
    #CROP! (Needs to go after Rfact calculation.)
    #nanmask *= 0.
    #nanmask[2500:3500, 2000:3500] = 1.0
    #p = pyfits.PrimaryHDU(nanmask)
    #p.writeto('nanmask.fits',clobber=True)
    #print('Cropping!')
    #data.crop=True


    # Get the array-positional info from RA dec
    y,x = hd.wcs_world2pix(ra, dec, 1)
    x,y = x.round().astype(int), y.round().astype(int)
    hy, hx = hd.wcs_world2pix(hra, hdec, 1)
    hx, hy = hx.round().astype(int), hy.round().astype(int)

    # In bounds?
    x, y   = checkXY(x,   y, Nx, Ny)
    hx, hy = checkXY(hx, hy, Nx, Ny)


    # Get the image to correlate on. If I'm not doing background correlation
    # then I just need to read in the Herschel catalog and inject ones into
    # the boolean map (depending on wtype).
    if not background:
        # I should have 2 catalogs. Inject the sources into the map from the
        # largest. Keep track which is DR1 and DR2
        if len(x) > len(hx):
            x1, y1 = x, y
            x2, y2 = hx, hy
        else:
            x1, y1 = hx, hy
            x2, y2 = x, y

        if wtype == 'DD':
            # I should have 2 catalogs. Inject the sources into the
            # map from the largest catalog.

            im = dataIm * 0.
            im[x1, y1] = 1.0
            im *= nanmask
            x, y = x2, y2
            print('Injecting x1, y1 into map and using x, y = x2, y2 for wtype=DD')

        elif wtype == 'DR1':
            Ngal = round(len(x1) * Rfact)
            im = randMapDiscrete(Nx, Ny, Ngal, nanmask=nanmask)
            print('Replacing map with len(x1)*3 random positions with x,y=x1,y1 (wtype=DR1).')
            x, y = x1, y1

        elif wtype == 'DR2':
            Ngal = round(len(x2) * Rfact)
            im = randMapDiscrete(Nx, Ny, Ngal, nanmask=nanmask)
            print('Replacing map with len(x2)*3 random positions with x,y=x2,y2(wtype=DR2).')
            x, y = x2, y2

        elif wtype == 'RR':
            Ngal = round(len(x1) * Rfact) #x1 is largest
            im = randMapDiscrete(Nx, Ny, Ngal, nanmask=nanmask)
            Ngal2 = round(len(x2) * Rfact)
            #x, y = randPosDiscrete(Nx, Ny, Ngal2) need to incorporate mask shape,
            # so this ^ is insufficient. Just take the x, y positions from a masked
            # random map!
            im2 = randMapDiscrete(Nx, Ny, Ngal, nanmask=nanmask)
            inds = np.where(im2 == 1)
            x, y = inds[0], inds[1]
            print(('Replacing map with len(x1)*3 random positions and catalog '+
                'with ~ len(x2)*3 random positions.'))

        else:
            print(('You need to set wtype = RR|DD|DR1|DR2. DR is not accepted '+
                'with background=False.'))
            raise Exception(wtype)

    ##########################################
    else:    # Background = True
        # Just apply the source mask as default?
        pdb.set_trace()
        if wtype == 'DD':
            # I have one catalog and one map. Catalog should be x,y, not hx, hy.
            im = dataIm # not times zero...
            im *= nanmask
            print('Injecting x1, y1 into map and using x, y = x2, y2 for wtype=DD')

        elif wtype == 'DR1':
            Ngal = round(len(x1) * Rfact)
            im = randMapDiscrete(Nx, Ny, Ngal, nanmask=nanmask)
            print('Replacing map with len(x1)*3 random positions with x,y=x1,y1 (wtype=DR1).')
            x, y = x1, y1

        elif wtype == 'DR2':
            Ngal = round(len(x2) * Rfact)
            im = randMapDiscrete(Nx, Ny, Ngal, nanmask=nanmask)
            print('Replacing map with len(x2)*3 random positions with x,y=x2,y2(wtype=DR2).')
            x, y = x2, y2

        elif wtype == 'RR':
            Ngal = round(len(x1) * Rfact) #x1 is largest
            im = randMapDiscrete(Nx, Ny, Ngal, nanmask=nanmask)
            Ngal2 = round(len(x2) * Rfact)
            #x, y = randPosDiscrete(Nx, Ny, Ngal2) need to incorporate mask shape,
            # so this ^ is insufficient. Just take the x, y positions from a masked
            # random map!
            im2 = randMapDiscrete(Nx, Ny, Ngal, nanmask=nanmask)
            inds = np.where(im2 == 1)
            x, y = inds[0], inds[1]
            print(('Replacing map with len(x1)*3 random positions and catalog '+
                'with ~ len(x2)*3 random positions.'))

        else:
            print(('You need to set wtype = RR|DD|DR1|DR2. DR is not accepted '+
                'with background=False.'))
            raise Exception(wtype)



#        else:
#            print("**** If I'm not running check, the background should just be gaussian noise? *****")
#            pdb.set_trace()
        #p=pyfits.PrimaryHDU(im, hdu[exten].header)
        #p.writeto('check/RR_check.fits',clobber=True)
        #pdb.set_trace()

    # replace nans with zeros.
    im[np.where(np.isnan(im))] = 0.


    # Apply the mask?
    if not nomask:
        hdumask = pyfits.open(crmaps + 'mask_'+filt+'.fits')
        mask = hdumask[0].data
        im *= mask
        print('Applying source mask to Herschel map.')


    # Remove sources in catalog that should be behind the mask.
    x, y = filterCatalog(x, y, Nx, Ny, nanmask=nanmask)
    xy = np.dstack((x,y))[0,0:]

    # ra, dec = hd.wcs_pix2world(x, y, 1)
    # pickle.dump(ra, open('ra_plw_crop.pickle', 'wb'), -1)
    # pickle.dump(dec, open('dec_plw_crop.pickle', 'wb'), -1)
    # pdb.set_trace()


    #Fill the data class and return it.
    data.filt       = filt
    data.catname    = catname
    data.wtype      = wtype
    data.nomask     = nomask
    data.background = background
    data.im         = im
    data.xy         = xy
    data.nbins      = nbins
    data.pscale     = pscale

    return data


def wsum(xy, shared_im, rbins, rbins_edges, savename, index, tic):

    '''
    Given a list of x and y positions, go to those positions and add
    up all the values in the array within some annulus, defined by the r
    binning in rbins*.
    '''

    # If this loop has already been done, skip this iteration.
    # Useful for stopping and starting.
    thissave = savename+str(index)+'.pickle'
#    if os.path.isfile(thissave):
#        sys.stdout.write('Skipping %d \n' % (index))
#        return


    nrad = len(rbins)
    wtheta = ({'sum':np.zeros(nrad), 'r':np.zeros(nrad)})

    xc, yc = xy[index,0].round(), xy[index,1].round()
    x,y = np.ogrid[:shared_im.shape[0], :shared_im.shape[1]]
    r2 = (x - xc)**2 + (y - yc)**2
    #print shared_im[xc,yc]

    j = 0
    for irad in range(0, nrad*2, 2):
        minrad = rbins_edges[irad]
        maxrad = rbins_edges[irad+1]
        mask = (r2 <= (maxrad-0.5)**2) & (r2 >= (minrad-0.5)**2)
        #mask = (r2 <= maxrad**2) & (r2 >= minrad**2)
        region = np.asarray(shared_im[np.where(mask)])
        if region.size > 0:
            wtheta['sum'][j] = region.sum()

        # Is the code doing what I think it is?
        #checkim = np.array(shared_im)
        # checkim[np.where(mask)] -= 1.
        # checkp = pyfits.PrimaryHDU(checkim)
        # checkp.writeto('check/checkim_'+str(j)+'.fits',clobber=True)
        # print 'x,y = ',xc,yc
        # print 'rminmax = ',minrad, maxrad
        # print 'sum = ', region.sum()
        #pdb.set_trace()

        wtheta['r'][j] = rbins[j]
        j += 1


    #Save wtheta.
    pickle.dump(wtheta, open(thissave, 'wb'), -1)
    #print("[Worker %d] Saved %s" % (os.getpid(), thissave))

    # Check progress
    #if index > 0 and index < 30: sys.stdout.write(".")
    if index == 30:
        oneloop = time.time() - tic
        num = len(glob.glob(thissave[0:thissave.rfind('/')] + '/*.pickle'))
        print(' ')
        print("~ %1.5f seconds for one set." % (oneloop/num))
        ttot = len(xy) * oneloop / 60. / 30#num
        print("This is going to take less than ~ %4.2f minutes!" % (ttot))

    if index > 10:
        if index % 500 == 0:
            pct = (float(index)/float(len(xy)))*100.
            print('%1.2f%%' % (pct))



## Invoke Parallel processing for each x,y pair.
if __name__ == '__main__':
    # Parse input parameters
    parser = (argparse.ArgumentParser(description=
                                      'Compute a w(theta) for each source position.'))
    (parser.add_argument('-wtype',
                         '--wtype', default='DD', type=str,
                         help = "We're using the standard estimator, so this is either the "+
                         "data-data ('DD'), or random-random ('RR') for w(theta)=2DD/RR-1."))
    (parser.add_argument('-filt',
                         '--Hfilter', default='plw', type=str,
                         help='Herschel filter. Default is PLW. Set -filt=check ' +
                         'to run on the $cross/check/sdss_check_map.fits to check'+
                         ' against the Landy-Szalay estimator (sdss_check_landy.png output of check.py)'))
    (parser.add_argument('-cat',
                         '--catname', default='nvss', type=str,
                         help='Input positional catalog. Default is NVSS'))
    (parser.add_argument('-nbins',
                         '--nbins', default=15, type=int,
                         help='Number of bins...'))
    (parser.add_argument('-ncores',
                         '--ncores', default=8, type=int,
                         help='Number of cores to parallelize over.'))
    (parser.add_argument('-nomask',
                         '--nomask', default=True, type=bool,
                         help='Set to False to apply the herschel source mask to the Herschel map.'+
                         ' To wit - probably want to run with background=True.'))
    (parser.add_argument('-background',
                         '--background', default=False, type=bool,
                         help='Set True to run cross-correlation of catname against Herschel background'))



    options    = parser.parse_args()
    wtype      = options.wtype.upper()
    Hfilter    = options.Hfilter.lower()
    catname    = options.catname.lower()
    ncores     = options.ncores
    nbins      = options.nbins
    nomask     = options.nomask
    background = options.background


    if nomask == False and background == False: # e.g. masking without doing background? nah.
        print 'Why would you want a masked image while not doing background computations?'

    checkrun = False
    if Hfilter == 'check':
        checkrun = True
        catname = 'check'


    if wtype != 'DD' and wtype != 'RR' and wtype != 'DR1' and wtype != 'DR2':
        print('You need to choose either "RR", "DD", "DR1" or "DR2" for wtype!')
        raise Exception(wtype)


    # First of all, if I'm runnin on the cluster everything should be saved to /scratch. dumb. this is.
    hostname = socket.gethostname()
    if hostname != 'thegord':
        scratch = os.path.exists('/scratch/ketron/')
        if not scratch:
            os.mkdir('/scratch/ketron/')
            print('Made directory /scratch/ketron/')


    # Naming conventions for save variables.
    cr = os.getenv('cr')
    pdir = cr + 'pickle/' + catname + '_' + Hfilter  # pickle save directory/name.
    if nomask:
        pdir += '_nomask'
    if checkrun: pdir = cr + 'pickle/check/'
    if pdir[-1] != '/': pdir += '/'
    check_pdir = os.path.exists(pdir)
    print pdir


    # Check for pickle save directory. Make one if not.
    if not check_pdir:
        os.system('mkdir -p '+pdir)
        print('Made directory "'+pdir+'"')

    savename = pdir + wtype + '_'
    ###pdb.set_trace()


    # Load the data image (Herschel) and exterior catalog.
    data = getCorrelationData(Hfilter, catname, wtype, nbins, nomask, background)
    im = data.im
    xy = data.xy
    d1 = im.shape[1]
    d2 = im.shape[0]

    data.nbins = nbins # set to keyword value!

    # Check stuff
    if data.crop:
        xx, yy = xy[:,0], xy[:,1]
        im = im[min(xx):max(xx)+1, min(yy):max(yy)+1]
        x = xx - min(xx)
        y = yy - min(yy)
        xy = np.dstack((x,y))[0,0:]
        data.im = im
        data.xy = xy
    o = np.where(im == 1)
    no = float(len(o[0]))
    ncat = no / float(len(xy))
    print 'N(im) / N(cat) = %1.4f' %ncat
    p=pyfits.PrimaryHDU(im)
    p.writeto('checkim.fits',clobber=True)
    print('Writing im to checkim.fits')
    newnp = np.zeros((im.shape[0], im.shape[1]))
    newnp[xy[:,0], xy[:,1]] = 1.0
    ppp=pyfits.PrimaryHDU(newnp)
    ppp.writeto('checkcat.fits', clobber=True)
    print 'wrote checkim, checkcat'


    # Binning
    minb = 2 #if checkrun else 10
    maxb = 2.4*max(im.shape) #m.floor(m.sqrt(dim1**2 + dim2**2))
    radial = 10 ** np.linspace(np.log10(minb), np.log10(maxb), data.nbins+2)
    nrad = len(radial)
    rbins_edges = np.zeros((nrad-2)*2); rbins = np.zeros(nrad-2); j=0
    for i in range(nrad-2):
        rbins[i] = np.mean([radial[i], radial[i+1]])
        rbins_edges[j]= radial[i]
        rbins_edges[j+1] = radial[i+1]
        j += 2

    nrad = len(rbins)
    outto = radial[-1] * data.pscale / 3600.
    into  = radial[0]  * data.pscale / 3600.
    print(('%d radial integrations at log widths between %1.1f and %1.1f degrees'
           %(nrad, into, outto)))


    # Naming again.
    if hostname != 'thegord':
        folder = tempfile.mkdtemp(prefix=cr)
    else:
        folder = tempfile.mkdtemp()

    print('Saving shared temp files to '+folder)
    im_name = os.path.join(folder, 'im')
    #mask_name = os.path.join(folder, 'mask')


    try:

        # Dump the input data to the disk to free up some memory.
        dump(im, im_name)

        # Release the reference on the original in memory array and replace it
        # by a reference to the memmap array so that the garbage collector can
        # release the memory before forking. gc.collect() is internally called
        # in Parallel just before forking.
        im = load(im_name, mmap_mode='r')

        tic = time.time()

        ## Parallel command.
        print('Starting parallelization with '+str(ncores)+' cores for '+str(len(xy)) + ' iterations.')

        Parallel(n_jobs=ncores)(delayed(wsum)(xy, im, rbins, rbins_edges, savename, index, tic)
                                for index in range(len(xy)))

        print('Total time was %1.1f minutes' % ((time.time() - tic)/60.))

    finally:

        try:
            print('Removing '+folder) # sharedmem folder.
            shutil.rmtree(folder)
        except:
            print('Failed to delete: ' + folder)

        if hostname != 'thegord':
            try:
                #print('Copying '+pdir+'*.pickle to my home directory.')
                #thomedir=tempfile.mkdtemp(prefix=wtype+'_tmp',dir='/data/users/ketron/')
                if not background:
                    pfix = '/data/users/ketron/' + catname + '_' + Hfilter + '_'
                else:
                    pfix = '/data/users/ketron/' + catname + '_' + Hfilter + '_bkg_'
                tempdir = tempfile.mkdtemp(prefix=pfix)
                print('Copying ' + pdir +'/*.pickle to '+tempdir)
                os.system('cp '+pdir + '/*.pickle ' + tempdir)#*.pickle ' + thomedir)
                print('Removing '+pdir)
                os.system('rm -r '+pdir)
                #shutil.rmtree(pdir)
                #os.system('scp -r '+ pdir
            except:
                print('Failed to copy files home from '+pdir)
                print("That directory still exists and it shouldn't")
