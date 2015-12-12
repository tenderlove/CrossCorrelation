from astropy import wcs
from astropy.io import fits
import numpy as np
import treecorr, pdb, os
import matplotlib.pyplot as pl
from astroML.datasets import fetch_sdss_specgals


## Set up SDSS data.
data = fetch_sdss_specgals()
m_max = 17.7
# redshift and magnitude cuts
data = data[data['z'] > 0.08]
data = data[data['z'] < 0.12]
data = data[data['petroMag_r'] < m_max]

# RA/DEC cuts
RAmin, RAmax = 200,220 #140, 220
DECmin, DECmax = 5, 15 #5, 45
data = data[data['ra'] < RAmax]
data = data[data['ra'] > RAmin]
data = data[data['dec'] < DECmax]
data = data[data['dec'] > DECmin]


ur = data['modelMag_u'] - data['modelMag_r']
flag_red = (ur > 2.22)
flag_blue = ~flag_red

data_red = data[flag_red]
data_blue = data[flag_blue]


print "data size:"
print "  red gals: ", len(data_red)
print "  blue gals:", len(data_blue)


ra1 = data_red['ra']; dec1 = data_red['dec']
ra2 = data_blue['ra']; dec2 = data_blue['dec']



## Generate treecorr catalog from map pixels
cat1 = treecorr.Catalog(ra=ra1, dec=dec1,ra_units='deg', dec_units='deg')

cat2 = treecorr.Catalog(ra=ra2, dec=dec2, ra_units='deg', dec_units='deg')

nn = treecorr.NNCorrelation(config)#bin_size=0.1, min_sep=.05, max_sep=10, sep_units='arcmin', bin_slop=0.5)

nn.process(cat1,cat2)

pdb.set_trace()
pl.errorbar(np.exp(nn.logr), nn.xip, nn.varxi, c='b',marker='o')
pl.xscale('log')
pl.ylabel(r'$\xi$',fontsize=20)
pl.xlabel(r'$\theta$ (arcmin)',fontsize=15)
#pl.size('large')
pl.xlim((.1,20))
pl.title('E x B')
#matplotlib.rcParams.update({'font.size': 22})
#pl.axis([
#pl.errorbar(np.exp(kk.logr), kk.xi, kk.xi, c='k')
# pl.xscale('log')
# pl.ylabel(r'$\xi$')
# pl.xlabel('R (arcmin)')

## ____________________________________
## Cross-correlation with E-mode maps.


pdb.set_trace()

