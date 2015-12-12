from __future__ import division
import pdb, os, time, sharedmem, glob, socket
import pickle, sys, tempfile, shutil, argparse
import numpy as np
import math as m
import astropy.io.fits as pyfits
from astropy import wcs
from joblib import Parallel, delayed , load, dump 
#from joblib.pool import has_shareable_memory
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import ion


def get_info(filt):

    ''' 
    Pixel scale, filename, dimensions, etc 
    from filter.
    '''

    #cr = os.getenv('cr')
    cr = '/data/users/ketron/'
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
    else:
        print "You didn't choose an appropriate Herschel filter."
        raise Exception(filt)
    print('Reading '+fname)
    hdu = pyfits.open(fname)
    im = hdu[1].data

    # Apply the mask! (kwd for this?)
    hdumask = pyfits.open(cr + 'mask_'+filt+'.fits')
    mask = hdumask[0].data
    im *= mask

    # replace nans with zeros?
    inds = np.where(np.isnan(im))
    im[inds] = 0.
    p = pyfits.PrimaryHDU(im)
    p.writeto('im.fits')
    print('Wrote im.fits')

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
    
    rmax = min([d1,d2])

    info = ({'pscale':pscale, 'd1':d1, 'd2':d2, 'fname':fname, 
             'im':im, 'hd':hd, 'nbins': 20, 'Linear':False, 'rmax':rmax})
    return info




def wsum(xy, shared_xgrid, shared_ygrid, shared_im, shared_mask, rbins, rbins_edges, savename, index, tic): 

    ''' 
    Given a list of x and y positions, go to those positions and add
    up all the values in the array within some annulus, defined by the r
    binning in rbins*.
    '''

    # If this loop has already been done, skip this iteration.
    # (Useful if the fucker crashes, which it inevitably will).
    thissave = savename+str(index)+'.pickle'
    if os.path.isfile(thissave): 
        sys.stdout.write('Skipping %d ' % (index))
        return

    x, y = xy[index,0].round(), xy[index,1].round()
    #x, y = xy[0].round(), xy[1].round()
    

    # Radial coords
    r = abs( (shared_xgrid - x) + 1j*(shared_ygrid - y) ) 
    nrad = len(rbins)    

    wtheta = ({'sum':np.zeros(nrad), 'numel':np.zeros(nrad), 'r':np.zeros(nrad),
               'min':np.zeros(nrad), 'max':np.zeros(nrad), 'std':np.zeros(nrad)}) 
                #'mean':np.zeros(nrad)
    j = 0
    for irad in range(0, nrad*2, 2):
        minrad = rbins_edges[irad]
        maxrad = rbins_edges[irad+1]#minrad + dr
        thisindex = (r>=minrad) * (r<maxrad) * shared_mask
        
        if not thisindex.ravel().any():
            wtheta['numel'][j] = 0
            wtheta['sum'][j] = np.nan
            wtheta['min'][j] = np.nan
            wtheta['max'][j] = np.nan            
            wtheta['std'][j] = np.nan
            #wtheta['mean'][j] = np.nan            
        else:
            wtheta['sum'][j] = shared_im[thisindex].sum()#[r<maxrad].sum()
            nonzindex = np.where(shared_im[thisindex] != 0)
            wtheta['numel'][j] = nonzindex[0].size
            wtheta['min'][j] = shared_im[thisindex][nonzindex].min()
            wtheta['max'][j] = shared_im[thisindex][nonzindex].min()            
            wtheta['std'][j] = shared_im[thisindex][nonzindex].std()

            #wtheta['mean'][j] = shared_im[thisindex].min()
            #ion() #CHECK
            #plt.imshow(thisindex*1)
            #pdb.set_trace()
        wtheta['r'][j] = rbins[j]
        j += 1
            

    #Save wtheta. (Should I pass it back, or somehow do this differently?)
    pickle.dump(wtheta, open(thissave, 'wb'), -1)
    #print("[Worker %d] Saved %s" % (os.getpid(), thissave))


    # Check progress
    if index > 0 and index < 30: sys.stdout.write(".")
    if index == 30:
        oneloop = time.time() - tic
        num = len(glob.glob(thissave[0:thissave.rfind('/')] + '/*.pickle'))
        print(' ')
        print("~ %1.5f seconds for one set." % (oneloop/num))
        ttot = len(xy) * oneloop / 3600. / num
        print("This is going to take less than ~ %4.2f hours!" % (ttot))
        
    if index > 10:
        if index % 200 == 0:
            pct = (float(index)/float(len(xy)))*100.
            print('%1.2f%%' % (pct))
                
                

## Invoke Parallel processing for each x,y pair.
if __name__ == '__main__':
    # Parse input parameters
    parser = (argparse.ArgumentParser(description=
                                      'Compute a w(theta) for each source position.'))
    (parser.add_argument('-filt', 
                         '--Hfilter', default='plw', type=str, 
                         help='Herschel filter. Default is PLW'))
    (parser.add_argument('-cat', 
                         '--catname', default='nvss', type=str, 
                         help='Input positional catalog. Default is NVSS'))
    (parser.add_argument('-ncpu',
                         '--ncores', default=8, type=int, 
                         help='Number of cores to parallelize over.'))
    (parser.add_argument('-test',
                         '--testrun', default=False, type=bool, 
                         help='Set to run on smaller map with 2 CPUs.'))
    (parser.add_argument('-check',
                         '--checkrun', default=False, type=bool, 
                         help='Set to run on the $cross/check/sdss_check_map.fits to check 
                         against the Landy-Szalay estimator (sdss_check_landy.png output of check.py)'))

    options = parser.parse_args()
    Hfilter = options.Hfilter.lower()
    catname = options.catname.lower()
    ncores  = options.ncores
    testrun = options.testrun
    checkrun = options.checkrun()
    if checkrun: catname = 'check'

    # First of all, if I'm runnin on the cluster everything should be saved to /scratch.
    # If /scratch/ketron doesn't exist, I can't make /scratch/ketron/whateverdir.
    hostname = socket.gethostname()
    if hostname != 'elgordo.ps.uci.edu':
        scratch = os.path.exists('/scratch/ketron/')
        if not scratch:
            os.mkdir('/scratch/ketron/')
            print('Made directory /scratch/ketron/')


    # Naming conventions for save variables.
    cr = os.getenv('cr')
    pdir = cr + 'pickle_' + catname+'_'+Hfilter  # pickle save directory.
    if pdir[-1] != '/': pdir += '/'
    check_pdir = os.path.exists(pdir)


    # Read catalogs. Presumably column-specific for each filter. 
    if not testrun: 
        info = get_info(Hfilter)
        im = info['im']
        Linear = info['Linear']

        # NVSS
        if catname == 'nvss':
            # Check for pickle save directory. Make one if not.
            if not check_pdir:
                os.system('mkdir '+pdir)
                print('Made directory "'+pdir+'"')
            savename = pdir+'wthetas_'

            #hdu = pyfits.open(os.getenv('cr') + 'nvss_cat.fits')
            hdu = pyfits.open('/data/users/ketron/nvss_cat.fits')
            ra = hdu[1].data['ALPHA_J2000']
            dec = hdu[1].data['DELTA_J2000']
            radec = np.dstack((ra,dec))[0,0:]
            radec = np.dstack((ra,dec))[0,0:]
            xy = info['hd'].wcs_world2pix(radec,1)
        else:
            raise Exception('Which catalog?')


    else: # Test run
        pdir = cr + 'pickleT/'
        check_pdir = os.path.exists(pdir)
        if not check_pdir:
            os.mkdir(pdir)
            print('Made directory "'+pdir+'"')
        ncores = 1
        print('Doing a test run on a small image with 1 core')
        im = np.ones((100,100)) ### TEST
        x = np.arange(5)
        y = np.arange(5)
        xy = np.dstack((x,y))[0,0:]
        #x = [25,30]
        #y = [25,30]
        info = {'pscale':12.,'nbins':20}
        Linear = False
        savename = pdir+'t_'


    d1 = im.shape[0]
    d2 = im.shape[1]


    # Meshgrid
    xvec = np.arange(d1) #- x[0]  #these shifts are paramount to this exercise. (now done in r)
    yvec = np.arange(d2) #- y[0] 
    xgrid, ygrid = np.meshgrid(xvec, yvec)
    mask = np.ones(im.shape, bool)


    # Binning
    if Linear: 
        rmax = min(im.shape)
        dr = np.abs([xgrid[0,0] - xgrid[0,1]]) * annulus_width
        pdb.set_trace()
        radial = np.arange(rmax/dr)*dr + dr/2.
        rbins = []
        for irad in range(nrad):
            minrad = irad*dr
            maxrad = minrad + dr
            rbins.append(minrad)
            rbins.append(maxrad)
        outto = radial[len(radial)-1] * info['pscale'] / 3600.
        nrad = len(rbins)
        print(('%d radial integrations at log widths of %1.1f pixels out to %1.1f degrees' 
               %(nrad, annulus_width, outto)))
        print('Need to update rbins and rbins_edges for linear case.')
        pdb.set_trace()
    else:
        min, max = 3, min(im.shape) #m.floor(m.sqrt(dim1**2 + dim2**2))
        steps = pow(max / min, 1. / (info['nbins']))
        radial = min * pow(steps, np.arange(info['nbins']+1))
        radial = radial.round().astype('int')
        radial = np.insert(radial, 0, 0)
        nrad = len(radial)
        rbins_edges = np.zeros((nrad-2)*2); rbins = np.zeros(nrad-2); j=0
        #pdb.set_trace()
        for i in range(nrad-2):
            #pdb.set_trace()
            rbins[i] = np.mean([radial[i], radial[i+1]])
            rbins_edges[j]= radial[i]
            rbins_edges[j+1] = radial[i+1]
            j += 2
        #pdb.set_trace()
        nrad = len(rbins)
        outto = radial[len(radial)-1] * info['pscale'] / 3600.
        print(('%d radial integrations at log widths out to %1.1f degrees' 
               %(nrad, outto)))
        #radial = np.insert(lbins, 0, 1.) #first element
    
        
    rng = np.random.RandomState(42)

    if hostname != 'elgordo.ps.uci.edu':
        folder = tempfile.mkdtemp(prefix=cr)
    else:
        folder = tempfile.mkdtemp()

    print('Saving shared temp files to '+folder)
    mask_name = os.path.join(folder, 'mask')
    im_name = os.path.join(folder, 'im')
    xgrid_name = os.path.join(folder, 'xgrid')
    ygrid_name = os.path.join(folder, 'ygrid')


    try:

        # Dump the input data to the disk to free up some memory.        
        dump(im, im_name)
        dump(mask, mask_name)
        dump(xgrid, xgrid_name)
        dump(ygrid, ygrid_name)


        # Release the reference on the original in memory array and replace it
        # by a reference to the memmap array so that the garbage collector can
        # release the memory before forking. gc.collect() is internally called
        # in Parallel just before forking.
        im    = load(im_name, mmap_mode='r')
        mask  = load(mask_name, mmap_mode='r')
        xgrid = load(xgrid_name, mmap_mode='r')
        ygrid = load(ygrid_name, mmap_mode='r')

        tic = time.time()

        
        ## Parallel command.
        if not testrun:
            ncores = 2
            
            print('Starting parallelization with '+str(ncores)+' cores for '+str(range(len(xy))) + ' iterations.')
            
            Parallel(n_jobs=ncores)(delayed(wsum)(xy, xgrid, ygrid, im, mask, rbins, rbins_edges, savename, index, tic)
                                    for index in range(50))
                                    #for index in range(len(xy)))
            
            print('Total time was %1.1f hours' % ((time.time() - tic)/3600.))
            
            #Ntest = 100
            #print('Starting Test with %1s iterations' %Ntest)
            #for jdogs in range(100):
            #    wsum(xy, xgrid, ygrid, im, mask, rbins, rbins_edges, savename, jdogs, tic)
            
        else:
            Ntest = 100
            print('Starting Test with %1s iterations' %Ntest)
            for jdogs in range(100):
                wsum(xy, xgrid, ygrid, im, mask, rbins, rbins_edges, savename, jdogs, tic)

    finally:
        try:
            print('Removing '+folder)
            shutil.rmtree(folder)
        except:
            print('Failed to delete: ' + folder)

        if hostname != 'elgordo.ps.uci.edu':
            try:
                print('Copying '+pdir+' to my home directory.')
                os.system('mv '+pdir + ' ~')
            except:
                print('Failed to copy files home from '+pdir)
                print("That directory still exists and it shouldn't")
                
            #try:
            #    print('Removing '+pdir)+'.'
            #    os.system('rm -r '+pdir)
            #except:
            #    print('Failed to remove directory '+pdir)
