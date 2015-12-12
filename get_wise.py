from astropy.io import fits
import numpy as np
import glob, os, pdb, math

""" With $cr/field/PLW maps as a reference, makes a point source file
for WISE IRSA input to get 2 deg^2 tiles. Need subsequent mosaicing."""

#field = 'helms'
cr = os.getenv('cr')
#plw1 = glob.glob(cr + field + '/herschel/*PLW*.fits')
#plw = glob.glob(cr + field + '/herschel/*500*.fits')

#hdu = fits.open(plw[0])
#refim = hdu['IMAGE'].data

#find min, max of non-nan values
#inds = np.where(~np.isnan(refim))
#x, y = inds[0], inds[1]

# fuck this. Do it manually.
#x0, y0, x1, y1 = 19.0, -8.85, 349.6, 8.67

#for x, concatenate from 19:0, 360:350. Incriment is 2 deg
inc = 1.95
ra = np.concatenate((np.linspace(20,0,41), np.linspace(358,348,21)))
dec = np.linspace(-8,8,33)
#grid = np.meshgrid(ra, dec)
#grid = grid[0]
#ra = grid[0,:]
f = open('/home/ketron/cross/wise_search_table.txt', 'w')
f.write('|              ra|              dec|\n')
f.write('|          double|           double|\n')
f.write('|             deg|              deg|\n')
for i in ra:
    for j in dec:
        f.write("%2f %2f\n" % (i, j))

f.close()

