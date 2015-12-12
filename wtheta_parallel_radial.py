from __future__ import division
import pyfits, pdb, os, time, sharedmem, multiprocessing, glob, pickle, sys
import numpy as np
import math as m
from astropy import wcs
import matplotlib.pyplot as plt
from matplotlib.pyplot import ion
#from joblib import Parallel, delayed  
#from joblib.pool import has_shareable_memory


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


def calc_radial(irad, r, dr, working_mask, im):

    #print 'here0'

    #print 'here1'
    
    radial_i = {'sum':0., 'numel_nonz':0.}#,'r':0., 'min':0., 'max':0., 'mean':0.,  'numel':0.}

    minrad = irad*dr
    maxrad = minrad + dr
    thisindex = (r>=minrad) * (r<maxrad) * working_mask[0]
    thisindex = thisindex.astype(bool)
    #print type(thisindex)
    #print type(im[0])
                
    if not thisindex.ravel().any():
        #radial_i['mean'] = np.nan
        radial_i['sum'] = np.nan
        #radial_i['numel'] = 0
        radial_i['numel_nonz'] = 0
        #radial_i['min'] = np.nan
        #radial_i['max'] = np.nan
        #radial_i['r'] = np.nan
        
    else:
        good = im[thisindex]
        #radial_i['mean'] = good.mean()
        radial_i['sum'] = good.sum()#[r<maxrad].sum()
        #radial_i['min'] = good.min()
        #radial_i['max'] = good.max()
        #radial_i['numel'] = good.size
        nonzindex = np.where(good != 0)
        radial_i['numel_nonz'] = nonzindex[0].size
        #radial_i['r'] = irad
        
    #print radial_i

    #return radial_i
    # ^that doesn't pass the class back to processes. All I can quickly think of to do is pickle :(
    pickle.dump(radial_i, open('radial_'+str(irad)+'.pickle', 'wb'), -1)


''' 
Given a list of x and y positions, go to those positions and add
up all the values in the array within some annulus, defined by the r
binning in linbins.
'''
def split_work(): 

    class Wtheta_i:
        '''
        Empty object container for one single radial computation.
        '''
        def __init__(self):
            #self.mean = None
            #self.std = None
            #self.median = None
            #self.numel = None
            self.numel_nonz = None
            self.sum = None
            #self.max = None
            #self.min = None
            self.r = None
    
    
    wthetas = []          # List of wtheta_i classes

    info = get_info('plw')
    im = info['im']
    #im = np.ones((10,20)) ### TEST
    d1 = im.shape[0]
    d2 = im.shape[1]
    annulus_width = 500. # BINNING!

    #NVSS
    hdu = pyfits.open(os.getenv('cr') + 'cats/nvss.fits')
    ra = hdu[1].data['ALPHA_J2000']
    dec = hdu[1].data['DELTA_J2000']
    radec = np.dstack((ra,dec))[0,0:]
    
    # ra dec -> x y
    radec = np.dstack((ra,dec))[0,0:]
    xy = info['hd'].wcs_world2pix(radec,1)
    x, y = xy[:,0].round(), xy[:,1].round()
    kinds = np.where(x <= 1500)
    x, y = x[kinds], y[kinds]
    
    # Meshgrid
    xvec = np.arange(d1) #- x[i]  #these shifts are paramount to this exercise. (now done in r)
    yvec = np.arange(d2) #- y[i] 
    xgrid, ygrid = np.meshgrid(xvec, yvec)

    # Define r for position[0], just so we can keep the same working_mask throughout.
    # It's in the corner, so this will be max rmax anyway. Probably too large (26 degrees).
    r = abs( (xgrid + x[0]) + 1j*(ygrid + y[0]) )
    working_mask = np.ones(im.shape, bool)
    #rmax = rad[working_mask].max()
    rmax = min(im.shape)
    dr = np.abs([xgrid[0,0] - xgrid[0,1]]) * annulus_width
    radial = np.arange(rmax/dr)*dr + dr/2.
    nrad = len(radial)
    
    # Multiprocessing
    shared1 = sharedmem.empty(working_mask.shape)
    shared1[:] = working_mask
    shared2 = sharedmem.empty(working_mask.shape)
    shared2[:] = im
    del im, working_mask
    
    # Now we have a list of positions and the square image array.
    # Loop through each x, y pair
    for i in range(len(x)):                
        
        # Radial coords
        if i == 0:
            if i == 0: tic = time.time()
            # (r already defined above)
        else:
            r = abs( (xgrid + x[i] - x[i-1]) + 1j*(ygrid + y[i] - y[i-1]) ) 
            
        shared3 = sharedmem.empty(r.shape)
        shared3[:] = r
        del r
        

        # multiprocessing     
        processes = ([multiprocessing.Process(
            target=calc_radial, args=(irad, shared3, dr, shared1, shared2)) for irad in range(nrad)])
        for pp in processes:
            pp.start()
        for pp in processes:
            pp.join()
                
                
        ## Keep track of this set of pickled radial measurements. There will be numel(RA) of these...

        # --- Initialize the data container ----
        wtheta_i = Wtheta_i() # each class has nrad elements.
        wtheta_i.r = np.zeros(nrad)
        wtheta_i.sum = np.zeros(nrad)
        wtheta_i.numel_nonz = np.zeros(nrad)
        
        picklefiles = glob.glob('radial*.pickle')    
        for picklefile in picklefiles:
            this_pickle = pickle.load(open(picklefile,'r'))
            
            # Outer radius in pixels. I will eventually be confused by how I did this,.
            irad = int(picklefile[str.find(picklefile,'_')+1:str.find(picklefile,'.')])
            
            wtheta_i.r[irad] = radial[irad] 
            wtheta_i.sum[irad] = this_pickle['sum']
            wtheta_i.numel_nonz[irad] = this_pickle['numel_nonz']
            #wtheta_i
        
        wthetas.append(wtheta_i)
        
        #print 'finished '+ str(i) +' of '+str(len(x))
        if i == 0: print('Checking how long this is going to take...')
        if i > 0 and i < 10: sys.stdout.write(".")
        if i == 10:
            oneloop = time.time() - tic
            print(' ')
            print("%1.5f seconds for one radial loop." % (oneloop/10.))
            ttot = len(x) * oneloop / 3600. / 10.
            print("This is going to take ~ %4.2f hours!" % (ttot))


        if i > 10:
            if i % 20 == 0:
                pct = (float(i)/float(len(x)))*100.
                print('%1.2f%%' % (pct))
    print('Total time was %1.2f' % ((time.time() - tic)/3600.))
    return wthetas
        #print i
    ## Just save the wthetas periodically in case it fucks up.
#    if i % 2000 == 0:
#        if os.path.isfile('wthetas_nvss.fits'):
            #append existing
#        else:
            #open new
#        tb = pyfits.PrimaryHDU(
#        pickle.dump(wthetas, open('wthetas_nvss.pickle','w'), -1)
#        print('Wrote wthetas_nvss.fits at %d' %i)

if __name__ == '__main__':
    split_work()
    pdb.set_trace()

# joblib command
#results = Parallel(n_jobs=12,verbose=5)(delayed(calc_radial)(irange, dr, working_mask) for irange in range(nrad))

#for irad in range(nrad): #= 1:numel(radial)
    

    #     wtheta = WTheta()
    #     wtheta.mean = np.zeros(nrad)
    #     wtheta.sum = np.zeros(nrad)
    #     #wtheta.std = np.zeros(nrad)
    #     #wtheta.median = np.zeros(nrad)
    #     wtheta.numel = np.zeros(nrad, dtype=int)
    #     wtheta.numel_nonz = np.zeros(nrad)
    #     #wtheta.max = np.zeros(nrad)
    #     #wtheta.min = np.zeros(nrad)
    #     wtheta.r = radial
        
    #     #call parallel shit here.

    # wthetas = [] #there will be len(radec) WTheta classes in here.



    #     print 'finished '+ str(i) +' of '+str(len(x))
    #     if i == 0:
            
    #         oneloop = time.time() - time0
    #         print("%1.5f seconds for one radial loop." % (oneloop))
    #         ttot = len(x) * oneloop / 3600.
    #         print("This is going to take ~ %4.2f hours!" % (ttot))
    #         return wthetas
            


