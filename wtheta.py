import glob, pickle, pdb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import ion

hfilt = 'psw'
pscale = 30.

#pickledir = '/data-2/cross/helms/pickle/nomask/pickle_nvss_' + hfilt + '/'
pickledir = '/data-2/cross/helms/pickle/check/'
dd_pickle = glob.glob(pickledir + 'DD*.pickle')
rr_pickle = glob.glob(pickledir + 'RR*.pickle')
dr_pickle = glob.glob(pickledir + 'DR*.pickle')

rlen, dlen, drlen = len(rr_pickle), len(dd_pickle), len(dr_pickle)
if dlen != drlen:
    print('What the fuuuuck?')
    raise Exception(rlen)

goto = dlen
goto_RR = rlen


#goto = rlen if rlen < dlen else dlen
print('Using %i data files.' %goto)
print('Using %i random files.' %goto_RR)


# Check the binning and params.
for i, dpick in enumerate([dd_pickle[0],rr_pickle[0], dr_pickle[0]]):
    thisp = pickle.load(open(dpick,'r'))
    nbins = len(thisp['r'])
    if i == 0: 
        print 'DD object has these keys:' 
    elif i == 1: 
        print 'RR object has these keys:'
    else:
        print 'DR object has these keys:'

    print thisp.keys()



# Compute w(theta) from DD and RR pickles.
dd_sum = []#np.zeros(goto)
rr_sum = []#np.zeros(goto)
dr_sum = []
dd_n = []
rr_n = []
dr_n = []
for i in range(goto):
    dname = pickledir + 'DD_'+str(i)+'.pickle'
    drname= pickledir + 'DR_'+str(i)+'.pickle'

    this_dd = pickle.load(open(dname,'r'))
    this_dr = pickle.load(open(drname,'r'))

    dd_sum.append(this_dd['sum'])
    dr_sum.append(this_dr['sum'])

    #dd_n.append(this_dd['n'])
    #dr_n.append(this_dr['n'])

for i in range(goto_RR):
    rname = pickledir + 'RR_'+str(i)+'.pickle'
    this_rr = pickle.load(open(rname,'r'))
    rr_sum.append(this_rr['sum'])
    #rr_n.append(this_rr['n'])


DD = zip(*dd_sum)
RR = zip(*rr_sum)
DR = zip(*dr_sum)
#DD_n = zip(*dd_n)
#RR_n = zip(*rr_n)
#DR_n = zip(*dr_n)

DDi, RRi, DRi = [], [], []
varDD, varRR, varDR = [], [], []
#fgal, fran = [], [] # Ngal, Nrandom factors for each
#for tRR, tDD, tDR, tRR_n, tDD_n, tDR_n in zip(RR, DD, DR, RR_n, DD_n, DR_n):
for tDD, tDR in zip(DD, DR):
    DD_avg = np.asarray(tDD)
    DR_avg = np.asarray(tDR)

    DDi_ = np.sum(DD_avg)
    DRi_ = np.sum(DR_avg)

    DDi.append(DDi_)
    DRi.append(DRi_)

    v1 = np.std(DD_avg)
    v2 = np.std(DR_avg)

    varDD.append(v1)
    varDR.append(v2)

for tRR in zip(RR):
    RR_avg = np.asarray(tRR)

    RRi_ = np.sum(RR_avg)

    RRi.append(RRi_)

    v3 = np.std(RR_avg)
    varRR.append(v3)


DD = np.asarray(DDi)
RR = np.asarray(RRi)
DR = np.asarray(DRi)

varDD = np.asarray(varDD)
varDR = np.asarray(varDR)
varRR = np.asarray(varRR)



fuckit = np.where(~np.isnan(DD) & ~np.isnan(RR) & ~np.isnan(DR))
if len(fuckit) <= 0:
    print 'fucking NaNs'
    pdb.set_trace()



factor = 3. # 3x as many random galaxies as data.
corr = (factor ** 2 * DD - 2 * factor * DR + RR) / RR 
corr *= 0.5 # wtf

# Error propagation for (DD - 2 * DR + RR) / RR
top = corr * RR
bot = RR

dtop = (varDD**2 + 4*varDR**2 + varRR**2) ** 0.5
dbot = varRR

dcorr = (corr * ((dtop/top)**2 + (dbot/bot)**2)) ** 0.5

ugh = np.where(corr > 0)[0]
if len(ugh) == 0:
    print 'Everything is zero'
    pdb.set_trace()


#Nr, Nd = 1., 1.
#corr = Nr/Nd * DD/RR - 1   # DD/DR - 1 wtf
theta = thisp['r'] * 30./3600
print factor

plt.ion()
plt.figure()

plt.xscale('log'); plt.yscale('log')
#plt.errorbar(r, wthetas, dw,lw=3,fmt='o',color='gray')

#plt.plot(theta, corr, 'o',color='gray',lw=3)
plt.errorbar(theta,corr, dcorr, color='gray',lw=3,fmt='o', label='Mine (map-catalog)')
plt.xlabel(r'$\theta\ (degrees)$', fontsize=18)
plt.ylabel(r'$w(\theta)$', fontsize=18)
plt.xlim([.01,10]); plt.ylim([1e-4,1])

#plt.ticklabel_format(style='sci', axis='y', scilimits=(1,4))
plt.rcParams.update({'font.size': 11})
plt.gcf().subplots_adjust(bottom=0.15)

xx=pickle.load(open("check/xx.pkl",'r'))
yy=pickle.load(open("check/yy.pkl",'r'))
er=pickle.load(open("check/er.pkl",'r'))
plt.errorbar(xx,yy,er,lw=3,fmt='o',color='red', label='Landy-Szalay (catalog-catalog)')

plt.text(7, .2, '17 < $r$ < 21',ha='right', va='top',fontsize=15)
plt.legend(frameon=False)

#newx, newy = np.loadtxt('/home/ketron/cross/check/sdss_digitized.txt',unpack=True)
#plt.plot(newx, newy, 'o', color = 'blue')
pdb.set_trace()

#plt.errorbar(r,wtheta,yerr=dw)

