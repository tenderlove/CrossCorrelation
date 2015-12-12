import montage_wrapper as montage
from astropy.io import fits

#field = 'helms'
#cr = os.getenv('cr')
#plw = cr + field + 'herschel/HELMS
ref = '/data-2/cross/helms/herschel/HELMS_image_500_SANEPIC_v0.2.fits'
hdu = fits.open(ref)
hd = hdu[1].header
montage.mosaic('/data-2/cross/helms/wise/w1', '/home/ketron/cross/mont')
