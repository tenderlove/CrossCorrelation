# Author: Jake VanderPlas
# License: BSD
#   The figure produced by this code is published in the textbook
#   "Statistics, Data Mining, and Machine Learning in Astronomy" (2013)
#   For more information, see http://astroML.github.com
#   To report a bug or issue, use the following forum:
#    https://groups.google.com/forum/#!forum/astroml-general
import numpy as np
import pdb
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt

from astroML.decorators import pickle_results
from astroML.datasets import fetch_sdss_specgals
from astroML.correlation import bootstrap_two_point_angular

#----------------------------------------------------------------------
# This function adjusts matplotlib settings for a uniform feel in the textbook.
# Note that with usetex=True, fonts are rendered with LaTeX.  This may
# result in an error if LaTeX is not installed on your system.  In that case,
# you can set usetex to False.
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=12, usetex=False)

#------------------------------------------------------------
# Get data and do some quality cuts
data = pyfits.open('../sdss_DR12.fits')
data = data[1].data
#data = fetch_sdss_specgals()
m_max = 17.7

# redshift and magnitude cuts
data = data[data['redshift'] > 0.08]
data = data[data['redshift'] < 0.12]
data = data[data['r'] < 21]#m_max]
data = data[data['r'] > 17]

# RA/DEC cuts
# RAmin, RAmax = 0,19
# DECmin, DECmax = -8, 8
# data1 = data[data['ra'] < RAmax]
# data1 = data1[data1['ra'] > RAmin]
# data1 = data1[data1['dec'] < DECmax]
# data1 = data1[data1['dec'] > DECmin]

# RAmin, RAmax = 349.5, 360
# data2 = data[data['ra'] < RAmax]
# data2 = data2[data2['ra'] > RAmin]
# data2 = data2[data2['dec'] < DECmax]
# data2 = data2[data2['dec'] > DECmin]

# data = np.concatenate((data1,data2),axis=0)

#pdb.set_trace()

ur = data['u'] - data['r']
flag_red = (ur > 2.22)
flag_blue = ~flag_red

data_red = data#[flag_red]
data_blue = data#[flag_blue]


# Save positions to make map for my sim.
rr = data_red['ra']
dr = data_red['dec']
rb = data_blue['ra']
db = data_blue['dec']
col_ra_red = pyfits.Column(name='RA', format='E', array=rr)
col_dec_red = pyfits.Column(name='dec', format='E', array=dr)
col_ra_blue = pyfits.Column(name='RA', format='E', array=rb)
col_dec_blue = pyfits.Column(name='dec', format='E', array=db)

cols_red = pyfits.ColDefs([col_ra_red, col_dec_red])
cols_blue = pyfits.ColDefs([col_ra_blue, col_dec_blue])

tbhdu_r = pyfits.BinTableHDU.from_columns(cols_red)
tbhdu_b = pyfits.BinTableHDU.from_columns(cols_blue)

tbhdu_r.writeto('sdss_red.fits', clobber=True)
tbhdu_b.writeto('sdss_blue.fits',clobber=True)
print('Wrote sdss_red/blue.fits\n')

print "data size:"
print "  red gals: ", len(data_red)
print "  blue gals:", len(data_blue)

#data_blue = data_red
#print('***** Data blue = data red')


#------------------------------------------------------------
# Set up correlation function computation
#  This calculation takes a long time with the bootstrap resampling,
#  so we'll save the results.
@pickle_results("correlation_functions.pkl")
def compute_results(Nbins=16, Nbootstraps=10,  method='landy-szalay', rseed=0):
    np.random.seed(rseed)
    bins = 10 ** np.linspace(np.log10(1. / 60.), np.log10(6), 16)

    results = [bins]

    results += bootstrap_two_point_angular(data_blue['ra'],
                                           data_blue['dec'],
                                           bins=bins,
                                           method=method,
                                           Nbootstraps=Nbootstraps)

    return results

bins, b_corr, b_corr_err, b_bootstraps = compute_results()

bin_centers = 0.5 * (bins[1:] + bins[:-1])

#------------------------------------------------------------
# Plot the results
corr = [ b_corr]
corr_err = [ b_corr_err]
bootstraps = [ b_bootstraps]
labels = ['$u-r < 2.22$\n$N=%i$' % len(data_blue)]

fig = plt.figure()#figsize=(5, 5))
#fig.subplots_adjust(bottom=0.2, top=0.9,
#                    left=0.13, right=0.95)
i=0
#ax = fig(121 + i, xscale='log', yscale='log')
plt.xscale('log'); plt.yscale('log')
plt.errorbar(bin_centers, corr[i], corr_err[i],
             fmt='.k', ecolor='gray', lw=3)

t = np.array([0.01, 10])
plt.plot(t, 10 * (t / 0.01) ** -0.8, ':k', linewidth=3)

plt.text(9, 9., labels[i],
        ha='right', va='top')#, transform=ax.transAxes)
plt.xlabel(r'$\theta\ (deg)$',fontsize=15)
if i == 0:
    plt.ylabel(r'$w(\theta)$',fontsize=15)

plt.ion()
plt.show()
