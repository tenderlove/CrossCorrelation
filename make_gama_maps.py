# This makes maps from L1 *filtered* data. They're written to
# /home/cthacke1/Atlas/G*/maps/


##################################
# Written by: Ketron
# This version modified by: Cameron Thacker, UC Irvine 2012
# Last modified April 2012
###################################
#from herschel.spire.ia.pipeline.phot.scanmap import MadScanMapperTask


G15 = True
G12 = False
G09 = False

halfMaps = True # set to make 2 tile maps (for overlaps) with 1 timestream from
		# each --so up for one and over for the other
SingleTiles = False
CoaddedTiles = False 

#Possible NaN removing way of making maps

#scans = baselineRemovalPolynomial(obs.level1, wholeTimeline=True,\
#        polyDegree=10)
#chanNoise = obs.calibration.phot.chanNoiseList.getProduct(\
#        obs.level1.getProduct(0).meta["biasMode"].value,\
#        obs.level1.getProduct(0).startDate)
#image = madScanMapper(input=scans, array="PSW", chanNoise=chanNoise)


corrlen = 10000
###################
PLW_size = 12.0 #arcsec/pix
PMW_size = 8.333
PSW_size = 6.0
###################

if G15:
	mult1_plw = 5480
	mult2_plw = 2168
	mult1_psw = 10960
	mult2_psw = 4336
if G12:
	mult1_plw = 3056
	mult2_plw = 2088
	mult1_psw = 6112
	mult2_psw = 4176
if G09:
	mult1_plw = 5392
	mult2_plw = 2164
	mult1_psw = 10784
	mult2_psw = 4328	

Lnaxis1 = long(round(0.002777777777777778 * mult1_plw / (PLW_size/3600)))
Lnaxis2 = long(round(0.002777777777777778 * mult2_plw / (PLW_size/3600)))
Lcrpix1 = long(round(Lnaxis1 / 2.))
Lcrpix2 = long(round(Lnaxis2 / 2.))
Lcdelt1 = PLW_size / 3600. * -1
Lcdelt2 = PLW_size / 3600. 

Mnaxis1 = long(round(0.002777777777777778 * mult1_plw / (PMW_size/3600)))
Mnaxis2 = long(round(0.002777777777777778 * mult2_plw / (PMW_size/3600)))
Mcrpix1 = long(round(Mnaxis1 / 2.))
Mcrpix2 = long(round(Mnaxis2 / 2.))
Mcdelt1 = PMW_size / 3600. * -1
Mcdelt2 = PMW_size / 3600. 

Snaxis1 = long(round(0.002777777777777778 * mult1_psw / (PSW_size/3600)))
Snaxis2 = long(round(0.002777777777777778 * mult2_psw / (PSW_size/3600)))
Scrpix1 = long(round(Snaxis1 / 2.))
Scrpix2 = long(round(Snaxis2 / 2.))
Scdelt1 = PSW_size / 3600. * -1
Scdelt2 = PSW_size / 3600. 

#G09
if G09:
	tile_names = ['26A2','4C95','4E22','4E3B','26A3','4C96','4E23','4E3C']
	#writedir = '/data-1/ATLAS/G09/maps/'
	writedir = '/data-1/ATLAS/SingleTile/G09/'

	wcs_plw = Wcs(crval1 = 134.68098, crval2 = 0.51089, \
		equinox = 2000.00, epoch = 2000.00, \
		cunit1 = "DEGREES",cunit2 = "DEGREES", ctype1 = "RA---TAN", \
		ctype2 = "DEC--TAN")
	wcs_pmw = Wcs(crval1 = 134.68098, crval2 = 0.51089, \
		equinox = 2000.00, epoch = 2000.00, \
		cunit1 = "DEGREES",cunit2 = "DEGREES", ctype1 = "RA---TAN", \
		ctype2 = "DEC--TAN")
	wcs_psw = Wcs(crval1 = 134.68098, crval2 = 0.51089, \
		equinox = 2000.00, epoch = 2000.00, \
		cunit1 = "DEGREES",cunit2 = "DEGREES", ctype1 = "RA---TAN", \
		ctype2 = "DEC--TAN")
	n1 = "/data-1/ATLAS/G09/Filtered/GAMA9-"
	n2 = "-L1-Tcorr-hpf-20101206-noturn-bolcorr-pdt-AstCorr.fits"

#G12
if G12:
	tile_names = ['5930','5939','5931','593A']
	#writedir = '/data-1/ATLAS/G12/maps/'
	writedir = '/data-1/ATLAS/SingleTile/G12/'

	wcs_plw = Wcs(crval1 = 176.52198, crval2 = -0.402658, \
		equinox = 2000.00, epoch = 2000.00, \
		cunit1 = "DEGREES",cunit2 = "DEGREES", ctype1 = "RA---TAN", \
		ctype2 = "DEC--TAN")
	wcs_pmw = Wcs(crval1 = 176.52198, crval2 = -0.402658, \
		equinox = 2000.00, epoch = 2000.00, \
		cunit1 = "DEGREES",cunit2 = "DEGREES", ctype1 = "RA---TAN", \
		ctype2 = "DEC--TAN")
	wcs_psw = Wcs(crval1 = 176.52198, crval2 = -0.402658 , \
		equinox = 2000.00, epoch = 2000.00, \
		cunit1 = "DEGREES",cunit2 = "DEGREES", ctype1 = "RA---TAN", \
		ctype2 = "DEC--TAN")
	n1 = "/data-1/ATLAS/G12/Filtered/GAMA12-"
	n2 = '-L1-Tcorr-hpf-20101130-noturn-pdt-AstCorr.fits'

#G15	
if G15:
	tile_names = ['305D','307B','3115','32AE','305E','307C','3116','32AF']
	writedir = '/data-1/ATLAS/SingleTile/G15/'
	
	wcs_plw = Wcs(crval1 = 217.5 , crval2 = 0.5, \
		equinox = 2000.00, epoch = 2000.00, \
		cunit1 = "DEGREES",cunit2 = "DEGREES", ctype1 = "RA---TAN", \
		ctype2 = "DEC--TAN")
	wcs_pmw = Wcs(crval1 = 217.5 , crval2 = 0.5, \
		equinox = 2000.00, epoch = 2000.00, \
		cunit1 = "DEGREES",cunit2 = "DEGREES", ctype1 = "RA---TAN", \
		ctype2 = "DEC--TAN")
	wcs_psw = Wcs(crval1 = 217.5 , crval2 = 0.5, \
		equinox = 2000.00, epoch = 2000.00, \
		cunit1 = "DEGREES",cunit2 = "DEGREES", ctype1 = "RA---TAN", \
		ctype2 = "DEC--TAN")
	n1 = '/data-1/ATLAS/G15/Filtered/GAMA15-'
	n2 = '-L1-Tcorr-hpf-20101110-noturn-bolcorr-pdt-astcorr.fits'


wcs_plw.naxis1 = Lnaxis1
wcs_plw.naxis2 = Lnaxis2
wcs_plw.cdelt1 = Lcdelt1
wcs_plw.cdelt2 = Lcdelt2
wcs_plw.crpix1 = Lcrpix1
wcs_plw.crpix2 = Lcrpix2

wcs_pmw.naxis1 = Mnaxis1
wcs_pmw.naxis2 = Mnaxis2
wcs_pmw.cdelt1 = Mcdelt1
wcs_pmw.cdelt2 = Mcdelt2
wcs_pmw.crpix1 = Mcrpix1
wcs_pmw.crpix2 = Mcrpix2

wcs_psw.naxis1 = Snaxis1
wcs_psw.naxis2 = Snaxis2
wcs_psw.cdelt1 = Scdelt1
wcs_psw.cdelt2 = Scdelt2
wcs_psw.crpix1 = Scrpix1
wcs_psw.crpix2 = Scrpix2



if CoaddedTiles:
	scans = SpireListContext()
	name = '_'
	print 'Making coadded tile with scans '+ str(tile_names)
	for i in range(len(tile_names)):
		data = fitsReader(file = n1 + tile_names[i] + n2)
		data=baselineRemovalPolynomial(data, wholeTimeline=True,\
			polyDegree=5)
		scans.addProduct(data)
		name += tile_names[i] + '_'
	
	chanNoise=fitsReader("/data-1/ATLAS/Auxilary/SCalPhotChanNoise_Nominal_20110811_v4.0.fits")

	mapPsw=madScanMapper(scans, array="PSW", wcs = wcs_psw,chanNoise=chanNoise, \
		correlationLength = corrlen, resolution = PSW_size)
	simpleFitsWriter(product=mapPsw,file=writedir + 'mad_hpf_noNoise'+ name +'PSW.fits')
	del(mapPsw)

	mapPmw=madScanMapper(scans, array="PMW", wcs = wcs_pmw,chanNoise=chanNoise, \
		correlationLength = corrlen, resolution = PMW_size)
	simpleFitsWriter(product=mapPmw,file=writedir + 'mad_hpf_noNoise'+ name +'PMW.fits')
	del(mapPmw)


	mapPlw=madScanMapper(scans, array="PLW", wcs = wcs_plw,chanNoise=chanNoise, \
		correlationLength = corrlen, resolution = PLW_size)
	simpleFitsWriter(product=mapPlw,file=writedir + 'mad_hpf_noNoise'+ name +'PLW.fits')
	del(mapPlw)
	del(scans)


if halfMaps:
	if G09:
		tile_horiz = ['4C95','26A3','26A3','4E3B','4E22','4C95']
		tile_vert  = ['26A2','4C96','4E3C','26A2','4C96','4E23']
	if G12:
		tile_horiz = ['5930','5939']
		tile_vert = ['593A','5931']
	if G15:
		tile_horiz = ['32AE','3115','307B','3115','305D','307B']
		tile_vert  = ['3116','32AF','3116','307C','307C','305E']
	writedir = '/data-1/ATLAS/Overlap/G15/PSW/'
	if G12:
		tile_names_h1 = tile_names[0:2] #we only have part of data set...
		tile_names_h2 = tile_names[2:4]
	else:
		tile_names_h1 = tile_names[0:4]
		tile_names_h2 = tile_names[4:8]
	scans_h1 = SpireListContext()
	scans_h2 = SpireListContext()
	name_h1 = ''
	name_h2 = ''

	for i in range(len(tile_horiz)):
		scans_h1 = SpireListContext()
		name_h1 = ''

		data_h1 = fitsReader(file = n1 + tile_horiz[i] + n2)
		data_h1=baselineRemovalPolynomial(data_h1, wholeTimeline=True,\
			polyDegree=5)
		scans_h1.addProduct(data_h1)
		name_h1 += tile_horiz[i] + '_'

		data_h1 = fitsReader(file = n1 + tile_vert[i] + n2)
		data_h1=baselineRemovalPolynomial(data_h1, wholeTimeline=True,\
			polyDegree=5)
		scans_h1.addProduct(data_h1)
		name_h1 += tile_vert[i] + '_'


	#half 1
		mapPsw=madScanMapper(scans_h1, array="PSW", wcs = wcs_psw, \
			correlationLength = corrlen, resolution = PSW_size)
		simpleFitsWriter(product=mapPsw,file=writedir + name_h1 +'PSW.fits')
		del(mapPsw)
#		mapPmw=madScanMapper(scans_h1, array="PMW", wcs = wcs_pmw, \
#			correlationLength = corrlen, resolution = PMW_size)
#		simpleFitsWriter(product=mapPmw,file=writedir + name_h1 +'PMW.fits')
#		del(mapPmw)
#		mapPlw=madScanMapper(scans_h1, array="PLW", wcs = wcs_plw, \
#			correlationLength = corrlen, resolution = PLW_size)
#		simpleFitsWriter(product=mapPlw,file=writedir + name_h1 +'PLW.fits')
#		del(mapPlw)
		del(scans_h1)


#	for i in range(len(tile_names_h1)):
#		scans_h1 = SpireListContext()
#		scans_h2 = SpireListContext()
#		name_h1 = ''
#		name_h2 = ''
#		data_h1 = fitsReader(file = n1 + tile_names_h1[i] + n2)
#		scans_h1.addProduct(data_h1)
#		name_h1 += tile_names_h1[i] + '_'
#
#		data_h2 = fitsReader(file = n1 + tile_names_h2[i] + n2)
#		scans_h2.addProduct(data_h2)
#		name_h2 += tile_names_h2[i] + '_'

#	#half 1
#		mapPsw=madScanMapper(scans_h1, array="PSW", wcs = wcs_psw, \
#			correlationLength = corrlen, resolution = PSW_size)
#		mapPmw=madScanMapper(scans_h1, array="PMW", wcs = wcs_pmw, \
#			correlationLength = corrlen, resolution = PMW_size)
#		mapPlw=madScanMapper(scans_h1, array="PLW", wcs = wcs_plw, \
#			correlationLength = corrlen, resolution = PLW_size)
#		simpleFitsWriter(product=mapPlw,file=writedir + name_h1 +'PLW.fits')
#		simpleFitsWriter(product=mapPmw,file=writedir + name_h1 +'PMW.fits')
#		simpleFitsWriter(product=mapPsw,file=writedir + name_h1 +'PSW.fits')
#        
#		del(mapPlw,mapPmw,mapPsw)
#
#	#half 2
#		mapPsw=madScanMapper(scans_h2, array="PSW", wcs = wcs_psw, \
#			correlationLength = corrlen, resolution = PSW_size)
#		mapPmw=madScanMapper(scans_h2, array="PMW", wcs = wcs_pmw, \
#			correlationLength = corrlen, resolution = PMW_size)
#		mapPlw=madScanMapper(scans_h2, array="PLW", wcs = wcs_plw, \
#			correlationLength = corrlen, resolution = PLW_size)
#		simpleFitsWriter(product=mapPlw,file=writedir + name_h2 +'PLW.fits')
#		simpleFitsWriter(product=mapPmw,file=writedir + name_h2 +'PMW.fits')
#		simpleFitsWriter(product=mapPsw,file=writedir + name_h2 +'PSW.fits')
#        
#		del(scans_h1, scans_h2)


if SingleTiles:
	for i in range(len(tile_names)/2):
		print 'making tile for '+ tile_names[i]
		scans = SpireListContext()
		data = fitsReader(file = n1 + tile_names[i] + n2)
		#data=baselineRemovalPolynomial(data, wholeTimeline=True,\
			#polyDegree=5)
		scans.addProduct(data)
		data = fitsReader(file = n1 + tile_names[i+len(tile_names)/2] + n2)
#		data=baselineRemovalPolynomial(data, wholeTimeline=True,\
#			polyDegree=5)
		scans.addProduct(data)
	# Set Channel Noise for madmap making
		chanNoise=fitsReader("/data-1/ATLAS/Auxilary/SCalPhotChanNoise_Nominal_20110811_v4.0.fits")

	### madMapper
		#Note I took away wcs = wcs_psw, wcs = wcs_pmw, wcs = wcs_plw
		# this was done to make single tile maps smaller
		# Also changed filename to PXW_big.fits instead of PXW_astr.fits as temporary change


		mapPsw=madScanMapper(input=scans, array="PSW",wcs = wcs_psw, \
			chanNoise=chanNoise,correlationLength=corrlen, resolution = PSW_size)
		simpleFitsWriter(product=mapPsw,file=writedir + 'mad_'+ tile_names[i] +'_PSW_big.fits')
		del(mapPsw)

		mapPmw=madScanMapper(input=scans, array="PMW", wcs = wcs_pmw, \
			chanNoise=chanNoise,correlationLength=corrlen, resolution = PMW_size)
		simpleFitsWriter(product=mapPmw,file=writedir + 'mad_'+ tile_names[i] +'_PMW_big.fits')
		del(mapPmw)

		mapPlw=madScanMapper(input=scans, array="PLW", wcs = wcs_plw, \
			chanNoise=chanNoise,correlationLength=corrlen, resolution = PLW_size)
		simpleFitsWriter(product=mapPlw,file=writedir + 'mad_'+ tile_names[i] +'_PLW_big.fits')
		del(mapPlw)
		del(scans)


	### Destriper - maybe use in the future. Best results, but new and uses naivemapper
		#scanss, mapPsw,diagprod,tod,signalminusmap = \
			#destriper(level1=scans, array='PSW',pixelSize=6.0,polyDegree=0,\
			#withMedianCorrected=False,iterThresh=1.0E-8)
		#del(mapPsw)
		#mapPsw=naiveScanMapper(input=scanss, array="PSW", wcs = wcs_psw, \
			#resolution = PSW_size)
		#simpleFitsWriter(product=mapPsw,file=writedir + 'mad_'+ tile_names[i] +'_PSW_astr.fits')
		#del(mapPsw,scansp)


		#scansm, mapPmw,diagprod,tod,signalminusmap = \
			#destriper(level1=scans, array='PMW',pixelSize=10.0,polyDegree=0,\
			#withMedianCorrected=False,iterThresh=1.0E-8)
		#del(mapPmw)
		#mapPmw=naiveScanMapper(input=scansm, array="PMW", wcs = wcs_pmw, \
			#resolution = PMW_size)
		#simpleFitsWriter(product=mapPmw,file=writedir + 'mad_'+ tile_names[i] +'_PMW_astr.fits')
		#del(mapPmw,scansm)


		#scansl, mapPlw,diagprod,tod,signalminusmap = \
			#destriper(level1=scans, array='PLW',pixelSize=14.0,polyDegree=0,\
			#withMedianCorrected=False,iterThresh=1.0E-8)
		#del(mapPlw)
		#mapPlw=naiveScanMapper(input=scansl, array="PLW", wcs = wcs_plw, \
			#resolution = PLW_size)
		#simpleFitsWriter(product=mapPlw,file=writedir + 'mad_'+ tile_names[i] +'_PLW_astr.fits')
		#del(mapPlw,scansl)


 




