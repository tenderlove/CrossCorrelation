from __future__ import division
import pdb, os, time, sharedmem, glob, socket
import pickle, sys, tempfile, shutil, argparse
import numpy as np
import math as m
import astropy.io.fits as pyfits
from astropy import wcs
from joblib import Parallel, delayed , load, dump


class correlationData:
    def __init__(self, **kwargs)#filt, catname, wtype, im, xy, nomask=True, background=False, nbins=20):
        # Input arguments.
        self.filt       = filt
        self.catname    = catname
        self.wtype      = wtype
        self.nomask     = nomask
        self.background = background
        self.im         = im
        self.xy         = xy
        self.nbins      = nbins
    #pdb.set_trace()

    def get_data(self):#filt, catname, wtype, nomask=False, background=False):


        '''
        Loads Herschel arrays and exterior catalog from kwds. Returns them with Pixel scale,
        filename, dimensions, etc. If wtype == 'RR' then random maps and catalogs are generated
        instead. If wtype='DR' a random map is made and x,y positions are from catalog. If
        background is set, the code correlates the input catalog name positions with
        the CFIRB in Herschel maps. Otherwise it generates maps from Herschel catalogs.

        '''

        hostname = socket.gethostname()
        if hostname != 'elgordo':
            cr = '/data/users/ketron/'
            crcats = cr
            crcheck = cr
        else:
            cr = os.getenv('cr')
            #cr += '/herschel/'
            crcats = cr + 'cats/'
            crcheck = '/home/ketron/cross/check/'

        # First load catalogs.
        if catname == 'nvss':
            print('Reading NVSS catlog:')
            print(crcats+'nvss_cat.fits')
            hdu = pyfits.open(crcats+'nvss_cat.fits')
            ra = hdu[1].data['ALPHA_J2000']
            dec = hdu[1].data['DELTA_J2000']
            radec = np.dstack((ra,dec))[0,0:]
        elif catname == 'check':
            print('Reading check catlog:')
            print('$cross/check/sdss_DR12_partial_check.fits')
            hdu = pyfits.open(crcheck+'sdss_DR12_partial_check.fits')
            ra = hdu[1].data['ra']
            dec = hdu[1].data['dec']
            radec = np.dstack((ra,dec))[0,0:]
        else:
            raise Exception('Which catalog?')


        # Now read [Herschel] array.
        filt = str(filt)
        filt = filt.lower()

        if filt == 'plw':
            pscale = 12.
            fname = cr + 'HELMS_image_500_SANEPIC_v0.2.fits'
        elif filt == 'pmw':
            pscale = 8.3333
            fname = cr + 'HELMS_image_350_SANEPIC_v0.2.fits'
        elif filt == 'psw':
            pscale = 6.
            fname = cr + 'HELMS_image_250_SANEPIC_v0.2.fits'
        elif filt == 'check':
            pscale = 30.
            fname = crcheck + 'sdss_check_map.fits'
        else:
            print "You didn't choose an appropriate Herschel filter."
            raise Exception(filt)

        print('Reading '+fname) 
        hdu = pyfits.open(fname)
        exten = 0 if filt == 'check' else 1

        # Astropy wcs doesn't like comments for some dumb reason.
        while 'COMMENT' in hdu[exten].header:
            hdu[exten].header.remove('COMMENT')

            im = hdu[exten].data
            hd = wcs.WCS(hdu[exten].header)

            # Get the array-positional info from RA dec
            x, y = hd.wcs_world2pix(ra, dec, 1)


            # If wtype = 'RR' or 'DR' then we need to replace the signal in the map with random
            # sources, and the clustered galaxy positions with random positions.
            # **** If I'm not running check, the background should just be gaussian noise? *****
            # DEFAULT is to use 3x the number of galaxies. This means RR calculation is 3x as long :/
            if wtype == 'RR' or wtype == 'DR':
                #pdb.set_trace()
                if filt == 'check':
                    # This is a catalog-catalog simulation, so RR Ngal = 3x DD Ngal. But for background
                    # I'm just replacing with a normal distribution.
                    Ngal = len(radec)
                    Nx, Ny = im.shape[0], im.shape[1]
                    randx = np.random.random_integers(0, Nx-1, Ngal*3)
                    randy = np.random.random_integers(0, Ny-1, Ngal*3)
                    im *= 0.
                    im[randx, randy] = 1.0
                    print('Replacing map with random positions (wtype=RR|DR).')

                if wtype == 'RR':
                    x, y = randx, randy
                    print('Setting x, y = randx, randy for wtype=RR')

                else:
                    print("**** If I'm not running check, the background should just be gaussian noise? *****")
                    pdb.set_trace()
            #p=pyfits.PrimaryHDU(im, hdu[exten].header)
            #p.writeto('check/RR_check.fits',clobber=True)
            #pdb.set_trace()

        # replace nans with zeros.
        im[np.where(np.isnan(im))] = 0.


        # Apply the mask?
        if ~nomask:
            hdumask = pyfits.open(cr + 'mask_'+filt+'.fits')
            mask = hdumask[0].data
            im *= mask
            print('Applying source mask to Herschel map.')


        #if filt == 'check':
    #    # catalog is slightly shifted...
    #    print('Using where(im != 0) as check catalog.')
    #    inds = np.where(im == 1)
    #    xy = np.dstack(inds)[0,0:]
    #else:
    #    xy = np.dstack((x,y))[0,0:]

    xy = np.dstack((x,y))[0,0:]
    data = ({'im':im, 'xy':xy, 'nbins': 0, 'fname':fname})

    #return data




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
                         '--nbins', default=20, type=int,
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



    options = parser.parse_args()
    wtype   = options.wtype.upper()
    Hfilter = options.Hfilter.lower()
    catname = options.catname.lower()
    ncores  = options.ncores
    nbins   = options.nbins
    nomask  = options.nomask


    if nomask == False and background == False: # e.g. masking without doing background? nah.
        print 'Why would you want a masked image while not doing background computations?'

    checkrun = False
    if Hfilter == 'check':
        checkrun = True
        catname = 'check'


    if wtype != 'DD' and wtype != 'RR' and wtype != 'DR':
        print('You need to choose either "RR", "DD" or "DR" for wtype, biotch!')
        raise Exception(wtype)


    # First of all, if I'm runnin on the cluster everything should be saved to /scratch. dumb. this is.
    hostname = socket.gethostname()
    if hostname != 'elgordo':
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

    pdb.set_trace()
    # Check for pickle save directory. Make one if not.
    if not check_pdir:
        os.system('mkdir -p '+pdir)
        print('Made directory "'+pdir+'"')

    savename = pdir + wtype + '_'
    ###pdb.set_trace()

    # Load the data image (Herschel) and exterior catalog.
    data = get_data(Hfilter, wtype)
    im = data['im']
    xy = data['xy']

    d1 = im.shape[1]
    d2 = im.shape[0]

    data['nbins'] = nbins # set to keyword value!


    # Binning
    min, max = 2., 3.*max(im.shape) #m.floor(m.sqrt(dim1**2 + dim2**2))
    radial = 10 ** np.linspace(np.log10(min), np.log10(max), data['nbins']+2)
    nrad = len(radial)
    rbins_edges = np.zeros((nrad-2)*2); rbins = np.zeros(nrad-2); j=0
    for i in range(nrad-2):
        rbins[i] = np.mean([radial[i], radial[i+1]])
        rbins_edges[j]= radial[i]
        rbins_edges[j+1] = radial[i+1]
        j += 2

    nrad = len(rbins)
    outto = radial[len(radial)-1] * data['pscale'] / 3600.
    print(('%d radial integrations at log widths out to %1.1f degrees'
           %(nrad, outto)))


    # Naming again.
    if hostname != 'elgordo':
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

        if hostname != 'elgordo':
            try:
                #print('Copying '+pdir+'*.pickle to my home directory.')
                #thomedir=tempfile.mkdtemp(prefix=wtype+'_tmp',dir='/data/users/ketron/')
                homepdir = os.path.exists('/data/users/ketron/'+pdir)
                print 'Moving ' + pdir +' to home'
                os.system('mv '+pdir + ' /data/users/ketron/ ')#*.pickle ' + thomedir)
                print('Removing '+pdir)
                shutil.rmtree(pdir)
                #os.system('scp -r '+ pdir
            except:
                print('Failed to copy files home from '+pdir)
                print("That directory still exists and it shouldn't")
