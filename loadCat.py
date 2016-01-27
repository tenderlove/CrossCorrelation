import astropy.io.fits as pyfits
import os

def loadCat(cname, crcats):

    '''
    Loads catalogs, which are presumably unique.
    '''

    cname = cname.lower()
    cr = os.getenv('cr')

    #NVSS radio catalog
    if cname == 'nvss':
        print('Reading NVSS catlog:')
        print(crcats + 'nvss_cat.fits')
        hdu = pyfits.open(crcats+'nvss_cat.fits')
        ra = hdu[1].data['ALPHA_J2000']
        dec = hdu[1].data['DELTA_J2000']

    # SDSS optical catalog -- 4+ million entries ==> magnitude clip!
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
        data = data[data['z'] > 0.08]
        data = data[data['z'] < 0.12]
        ra = data['ra']
        dec = data['dec']
        print('Clipped sources for only 17 < ' + cname.upper()  + ' < 21')
        print('and  0.08 > z > 0.12...')

    # HELMS Herschel data.
    elif cname == 'plw' or cname == 'pmw' or cname == 'psw':
        print('Reading Herschel HELMS '+cname+' map.')
        hcat = pyfits.open(cr+'herschel/helms_v0.2_'+cname.upper() + '_SXT.fits')
        hdata = hcat[1].data
        #hdata = hdata[hdata['flux'] > 30]
        ra, dec = hdata['ra'], hdata['dec']

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
