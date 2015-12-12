import numpy as np
import pdb, os, pickle
import astropy.io.fits as pyfits
from astropy import wcs
from matplotlib import pyplot as plt
from astroML.decorators import pickle_results
from astroML.datasets import fetch_sdss_specgals
from astroML.correlation import bootstrap_two_point_angular
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=12, usetex=False)
# Author: Jake VanderPlas
# License: BSD
#   The figure produced by this code is published in the textbook
#   "Statistics, Data Mining, and Machine Learning in Astronomy" (2013)
#   For more information, see http://astroML.github.com
#   To report a bug or issue, use the following forum:
#    https://groups.google.com/forum/#!forum/astroml-general
#------------------------------------------------------------
# Now with ketron edits.


#os.system('rm correlation_functions.pkl')

do_RR=False

# Get data and do some quality cuts
#data = pyfits.open('../sdss_DR12.fits')
if not do_RR:
    ## small area: ##
    p = pyfits.open('/home/ketron/cross/sdss_DR12.fits')
    data = p[1].data
    data = data[data['RA'] < 20]
    data = data[data['dec'] > 1.3]
    data = data[data['class'] == 'GALAXY']
    data = data[data['r'] < 21]#m_max]
    data = data[data['r'] > 17]

    ## large area ##
    #p=pyfits.open('/data-2/cross/helms/cats/sdss_master.fits')
    #data = p[1].data
    #data = data[data['psfMag_r'] < 21]
    #data = data[data['psfMag_r'] > 17]

    # redshift and magnitude cuts
    #data = data[data['redshift'] > 0.08]
    #data = data[data['redshift'] < 0.12]
    
    ra = data['RA']
    dec = data['dec']
    tb = (pyfits.BinTableHDU.from_columns([
        pyfits.Column(name='ra',format='E',array=ra),
        pyfits.Column(name='dec',format='E',array=dec)]))
    tb.writeto('/home/ketron/cross/check/sdss_DR12_partial_check.fits',clobber=True)
    print('Wrote $cross/check/sdss_DR12_partial_check.fits')
    pklname = "correlation_functions.pkl"
else:
    # Check that my random map/catalog is indeed random.
    p = pyfits.open('RR_check.fits')
    array = p[0].data
    hd = wcs.WCS(p[0].header)
    inds = np.where(array == 1)
    ra,dec = hd.all_pix2world(inds[0],inds[1], 1)
    data = {'ra':ra,'dec':dec}
    pklname = 'correlation_RR.pkl'

#ur = data['u'] - data['r']
#flag_red = (ur > 2.22)
#flag_blue = ~flag_red


data_blue = data#[flag_blue]
#data_red = data#[flag_red]

#pdb.set_trace()

print "data size:", len(data_blue['ra'])
#print "  red gals: ", len(data_red)
#print "  blue gals:", len(data_blue)

#data_blue = data_red
#print('***** Data blue = data red')


#------------------------------------------------------------
# Set up correlation function computation
#  This calculation takes a long time with the bootstrap resampling,
#  so we'll save the results.
@pickle_results(pklname)
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
labels = ['$17 < r < 21$\n$N=%i$' % len(data_blue)]
#labels = ['N=%i' % len(data_blue)]

fig = plt.figure()#figsize=(5, 5))
#fig.subplots_adjust(bottom=0.2, top=0.9,
#                    left=0.13, right=0.95)
i=0
#ax = fig(121 + i, xscale='log', yscale='log')
plt.xscale('log'); plt.yscale('log')
plt.errorbar(bin_centers, corr[i], corr_err[i],lw=3,fmt='o',color='gray')

xx = bin_centers
yy = corr[i]
er = corr_err[i]
pickle.dump(xx, open('xx.pkl', 'wb'), -1)
pickle.dump(yy, open('yy.pkl', 'wb'), -1)
pickle.dump(er, open('er.pkl', 'wb'), -1)

t = np.array([0.01, 10])
#plt.plot(t, 10 * (t / 0.01) ** -0.8, ':k', linewidth=3)

plt.text(7, .9, labels[i],
         ha='right', va='top',fontsize=15)#, transform=ax.transAxes)
plt.xlabel(r'$\theta\ (deg)$',fontsize=15)
if i == 0:
    plt.ylabel(r'$w(\theta)$',fontsize=15)

plt.ion()
plt.show()
